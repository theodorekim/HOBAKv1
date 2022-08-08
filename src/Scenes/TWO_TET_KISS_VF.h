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
#ifndef TWO_TET_KISS_VF_H
#define TWO_TET_KISS_VF_H

#include "SIMULATION_SCENE.h"

namespace HOBAK {

class TWO_TET_KISS_VF : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " One tet is just barely kissing another in an vertex-edge collision." << endl;
    cout << " The integrator should just barely push them apart according to the" << endl;
    cout << " collision epsilon" << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "two_tet_kiss_vf";

    // read in the tet mesh file
    _tetMeshFilename = string("../data/two_tets_vertex_face_kiss.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    _gravity = VECTOR3(0, 0, 0);

    REAL E = 10.0;
    REAL nu = 0.45; // lambda \approx 10

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH(vertices, tets);
    _hyperelastic = new VOLUME::ARAP(mu, lambda);
    _tetMesh->setCollisionEps(0.03);

    // build the time integrator
    _solver = new TIMESTEPPER::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic);
    _solver->setDt(1.0 / 30.0);

    VECTOR3 center0(0.0,-4.90, 0.0);
    VECTOR3 center1(0.0,6.9, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _kinematicShapes.push_back(new CUBE(center1, 10.0));
    _solver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _solver->attachKinematicSurfaceConstraints(_kinematicShapes[1]);
    _solver->addKinematicCollisionObject(_kinematicShapes[0]);

    _solver->setRayeligh(0,0);

    // testing settings - why doesn't this work?
    const REAL collisionMu = mu * 10.0;
    const REAL dampingFactor = 0.01;

    // collision constants
    _solver->collisionStiffness() = collisionMu;
    _solver->collisionDampingBeta() = dampingFactor;

    _solver->vertexFaceSelfCollisionsOn() = true;
    _solver->edgeEdgeSelfCollisionsOn() = false;

    _eye    = VECTOR3(-2.13557, 0.947023, 2.39579);
    _lookAt = VECTOR3(-1.35707, 0.964614, 1.76839);
    _up     = VECTOR3(-0.00368317, 0.999718, 0.0234604);

    _pauseFrame = 100;

    return true;
  }
  
};

}

#endif
