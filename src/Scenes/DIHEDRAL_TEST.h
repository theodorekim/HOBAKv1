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
#ifndef DIHEDRAL_TEST_H
#define DIHEDRAL_TEST_H

#include "SIMULATION_SCENE.h"

namespace HOBAK {

class DIHEDRAL_TEST : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Single tet, and we wind one vertex around to test if the dihedral " << endl;
    cout << " angle being computed is correct. " << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "dihedral_test";

    // read in the tet mesh file
    //string testFile("../data/two_tets_edge_edge.tobj");
    _tetMeshFilename = string("../data/single_tet.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }
    _gravity = VECTOR3(0, -0.3, 0);

    REAL E = 1.0;
    REAL nu = 0.45; // lambda \approx 10

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::SNH(mu, lambda);

    // build the time integrator
    _solver = new TIMESTEPPER::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic);
    _solver->setDt(1.0 / 30.0);

    VECTOR3 center0(0.0,-5.90, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _solver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _solver->addKinematicCollisionObject(_kinematicShapes[0]);

    // collision constants
    const REAL collisionMu = 10.0;
    _solver->collisionStiffness() = collisionMu;
    _solver->collisionDampingBeta() = 0.01;

    _solver->vertexFaceSelfCollisionsOn() = false;
    _solver->edgeEdgeSelfCollisionsOn() = true;

    _eye    = VECTOR3(-4.05987, 0.675697, 3.46079);
    _lookAt = VECTOR3(-3.35421, 0.395048, 2.81018);
    _up     = VECTOR3(0.237082, 0.958807, -0.156453);

    // smoosh the two vertices together so the angle starts at zero
    _tetMesh->vertices()[3] = _tetMesh->vertices()[0];

    return true;
  }

  virtual void stepSimulation(const bool verbose = true) override {
    // rotate the vertex about z axis, see what angle we get
    VECTOR3& v3 = _tetMesh->vertices()[3];
    MATRIX3 R = Eigen::AngleAxisd(0.01, VECTOR3::UnitZ()).toRotationMatrix();
    v3 = R * v3;

    REAL angle = _tetMesh->surfaceFaceDihedralAngle(1,3);
    angle *= 360.0 / (2.0 * M_PI);

    if (verbose)
      cout << " Dihedral: " << angle << endl;

  };

  // drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    SIMULATION_SCENE::drawScene();

    glDisable(GL_DEPTH_TEST);
    //drawSurfaceFace(*_tetMesh, _arrowCounter);
    glPointSize(10.0);
    glColor4f(1.0, 0.0, 0.0, 10.0);
    drawVertex(*_tetMesh, _arrowCounter);
    glEnable(GL_DEPTH_TEST);
  };
#endif
};

}

#endif
