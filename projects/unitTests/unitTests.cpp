/*
This file is part of HOBAK.

Permission is hereby granted to use this software solely for non-commercial applications
and purposes including academic or industrial research, evaluation and not-for-profit media
production.
 
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES, YALE, OR
THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <cmath>
#include <iostream>
#include <float.h>

#define CATCH_CONFIG_MAIN

#include "util/MATRIX_UTIL.h"
#include "Hyperelastic/Volume/ARAP.h"
#include "Hyperelastic/Volume/NEO_HOOKEAN_BW.h"
#include "Hyperelastic/Volume/SNH.h"
#include "Hyperelastic/Volume/SNH_WITH_BARRIER.h"
#include "Hyperelastic/Volume/STVK.h"
#include "Hyperelastic/Volume/LINEAR.h"
#include "Hyperelastic/Volume/ANISOTROPIC_ARAP.h"
#include "Hyperelastic/Volume/ANISOTROPIC_STVK.h"
#include "Hyperelastic/Volume/ANISOTROPIC_FUNG.h"
#include "Hyperelastic/Volume/ANISOTROPIC_DIRICHLET.h"
#include "Hyperelastic/Volume/VERTEX_FACE_COLLISION.h"
#include "Hyperelastic/Volume/VERTEX_FACE_SQRT_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_SQRT_COLLISION.h"
#include "Hyperelastic/Volume/MCADAMS_COLLISION.h"
#include "Damping/Volume/GREEN_DAMPING.h"
#include "Geometry/TET_MESH.h"
#include "Geometry/TET_MESH_FASTER.h"
#include "Geometry/CAPSULE.h"
#include "Geometry/CYLINDER.h"
#include "Geometry/CUBE.h"
#include "Geometry/SPHERE.h"
#include "Geometry/FIELD_3D.h"
#include "Geometry/LINE_INTERSECT.h"
#include "Timestepper/QUASISTATIC.h"
#include "Timestepper/NEWMARK.h"
#include "Timestepper/BACKWARD_EULER_VELOCITY.h"
#include "util/COLLISION_UTIL.h"
#include "catch_amalgamated.hpp"

using namespace HOBAK;
using namespace std;

#include "COLLISION_FORCE_TESTS.h"
#include "COLLISION_DETECTION_TESTS.h"
#include "HYPERELASTIC_TESTS.h"
#include "DAMPING_TESTS.h"
#include "SCENE_TESTS.h"
#include "FIELD_TESTS.h"

//////////////////////////////////////////////////////////////////////////////
// see if the optimizations in stiffness assembly broke anything
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Stiffness assembly optimization", "[Stiffness optimization]" )
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING TET_MESH_FASTER Stiffness assembly is still consistent" << endl;
  cout << "=============================================================== " << endl;
  // read in the tet mesh file
  string testFile("../data/cube_10.tobj");
  vector<VECTOR3> inputVertices;
  vector<VECTOR4I> inputTets;
  bool success = TET_MESH::readTobjFile(testFile, inputVertices, inputTets);
  inputVertices = TET_MESH::normalizeVertices(inputVertices);
  REQUIRE(success);

  // build the tet mesh objects
  TET_MESH        tetMesh(inputVertices, inputTets);
  TET_MESH_FASTER tetMeshFast(inputVertices, inputTets);
  VOLUME::ARAP    hyperelastic(1.0, 1.0);

  // scramble the vertices
  vector<VECTOR3>& vertices = tetMesh.vertices();
  std::mt19937 gen(456);
  std::uniform_real_distribution<REAL> dist(0.0, 1.0);

  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    vertices[x][0] = dist(gen);
    vertices[x][1] = dist(gen);
    vertices[x][2] = dist(gen);
  }
  tetMeshFast.vertices() = tetMesh.vertices();

  tetMesh.computeFs();
  tetMesh.computeSVDs();
  tetMeshFast.computeFs();
  tetMeshFast.computeSVDs();

  SPARSE_MATRIX slowA = tetMesh.computeHyperelasticClampedHessian(hyperelastic);
  SPARSE_MATRIX fastA = tetMeshFast.computeHyperelasticClampedHessian(hyperelastic);

  SPARSE_MATRIX diff = fastA - slowA;
  REAL diffNorm = diff.norm();
  cout << " Diff against TET_MESH is: " << diffNorm << " (should be zero) " << endl;

  REQUIRE(diffNorm < 1e-8);
}

//////////////////////////////////////////////////////////////////////////////
// Testing dihedral bending, but for surface faces on a tet mesh,
// this corresponds to the SIMULATION_SCENE test DIHEDRAL_TEST
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Dihedral Bending on Tet Meshes", "[Bending Energies]" )
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING TET_MESH dihedral bending computation" << endl;
  cout << "=============================================================== " << endl;

  string testFile("../data/single_tet.tobj");
  vector<VECTOR3> vertices;
  vector<VECTOR4I> tets;
  bool success = TET_MESH::readTobjFile(testFile, vertices, tets);
  REQUIRE(success);
  TET_MESH* tetMesh = new TET_MESH(vertices, tets);

  REAL angle;
  VECTOR3 original = tetMesh->vertices()[3];

  // smoosh the two vertices together so the angle starts at M_PI
  tetMesh->vertices()[3] = tetMesh->vertices()[0];
  angle = tetMesh->surfaceFaceDihedralAngle(1,3);
  angle *= 360.0 / (2.0 * M_PI);
  cout << "Angle " << angle << " should be 180"<< endl;
  REQUIRE(fabs(angle - 180.0) < 1e-8);

  tetMesh->vertices()[3][0] = 0.0;
  angle = tetMesh->surfaceFaceDihedralAngle(1,3);
  angle *= 360.0 / (2.0 * M_PI);
  cout << "Angle " << angle << " should be 135"<< endl;
  REQUIRE(fabs(angle - 135.0) < 1e-8);

  tetMesh->vertices()[3][0] = original[0];
  angle = tetMesh->surfaceFaceDihedralAngle(1,3);
  angle *= 360.0 / (2.0 * M_PI);
  cout << "Angle " << angle << " should be 90"<< endl;
  REQUIRE(fabs(angle - 90.0) < 1e-8);
  
  tetMesh->vertices()[3][1] = 0.0;
  angle = tetMesh->surfaceFaceDihedralAngle(1,3);
  angle *= 360.0 / (2.0 * M_PI);
  cout << "Angle " << angle << " should be 45"<< endl;
  REQUIRE(fabs(angle - 45.0) < 1e-8);
  
  tetMesh->vertices()[3] = -1.0 * tetMesh->vertices()[0];
  angle = tetMesh->surfaceFaceDihedralAngle(1,3);
  angle *= 360.0 / (2.0 * M_PI);
  cout << "Angle " << angle << " should be 0"<< endl;
  REQUIRE(fabs(angle - 0.0) < 1e-8);
  
  tetMesh->vertices()[3][0] = 0.0;
  angle = tetMesh->surfaceFaceDihedralAngle(1,3);
  angle *= 360.0 / (2.0 * M_PI);
  cout << "Angle " << angle << " should be -45"<< endl;
  REQUIRE(fabs(angle + 45.0) < 1e-8);
  
  tetMesh->vertices()[3] = tetMesh->vertices()[0];
  tetMesh->vertices()[3][1] *= -1;
  angle = tetMesh->surfaceFaceDihedralAngle(1,3);
  angle *= 360.0 / (2.0 * M_PI);
  cout << "Angle " << angle << " should be -90"<< endl;
  REQUIRE(fabs(angle + 90.0) < 1e-8);

  tetMesh->vertices()[3] = tetMesh->vertices()[0];
  tetMesh->vertices()[3][1] = 0;
  angle = tetMesh->surfaceFaceDihedralAngle(1,3);
  angle *= 360.0 / (2.0 * M_PI);
  cout << "Angle " << angle << " should be -135"<< endl;
  REQUIRE(fabs(angle + 135.0) < 1e-8);

  delete tetMesh;
}
