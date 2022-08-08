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
#ifndef DROPPED_C_CCD_TEST_H
#define DROPPED_C_CCD_TEST_H

#include "SIMULATION_SCENE.h"

namespace HOBAK {

class DROPPED_C_CCD_TEST : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " A C-shaped mesh is dropped on the ground, and the collisions " << endl;
    cout << " are sufficiently severe that the top lip goes past the collision " << endl;
    cout << " epsilon in a single timestep" << endl;
    cout << endl;
    cout << " CCD is (presumably) needed for this scenario, but NOT IMPLEMENTED." << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "dropped_c_ccd_test";

    // read in the tet mesh file
    //string testFile("../data/c_8.tobj");
    //string testFile("../data/c_10.tobj");
    _tetMeshFilename = string("../data/c_20.tobj");
    //string testFile("../data/c_40.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    _gravity = VECTOR3(0, -1.0, 0);

    REAL E = 10.0;
    REAL nu = 0.45; // lambda \approx 10

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::ARAP(mu, 1.0);
    _tetMesh->setCollisionEps(0.01);

    const vector<REAL>& areas = _tetMesh->surfaceTriangleAreas();
    REAL smallest = areas[0];
    REAL largest = areas[0];
    for (unsigned int x = 1; x < areas.size(); x++)
    {
      if (areas[x] > largest) largest  = areas[x];
      if (areas[x] < largest) smallest = areas[x];
    }
    cout << "Largest triangle area:  "  << largest << endl;
    cout << "Smallest triangle area: "  << smallest << endl;

    // build the time integrator
    _solver = new TIMESTEPPER::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic);
   
    VECTOR3 center0(0.0,-4.77, 0.0);
    center0[1] = -5.0;
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _solver->addKinematicCollisionObject(_kinematicShapes[0]);

    // collision constants
    const REAL collisionMu = 1000.0;
    _solver->collisionStiffness() = collisionMu;
    _solver->collisionDampingBeta() = 0.01;

    _solver->vertexFaceSelfCollisionsOn() = true;
    _solver->edgeEdgeSelfCollisionsOn() = true;

    _eye    = VECTOR3(-0.755866, 0.174748, 0.402933);
    _lookAt = VECTOR3(0.24092, 0.108045, 0.447299);
    _up     = VECTOR3(0.0, 1.0, 0.0);  
    return true;
  }
  
};

}

#endif
