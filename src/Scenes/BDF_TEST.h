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
#ifndef BDF_TEST_H
#define BDF_TEST_H

#include "SIMULATION_SCENE.h"
#include "Timestepper/BDF_1.h"
#include "Timestepper/BDF_2.h"

namespace HOBAK {

class BDF_TEST : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Dropping a cylinder on the kinematic ground, making sure that it" << endl;
    cout << " reboundes correctly " << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "bdf_test";

    // read in the tet mesh file
    _tetMeshFilename = string("../data/cylinder_10.tobj");

    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    _gravity = VECTOR3(0, -1.0, 0);

    const REAL E = 3.0;  // default
    const REAL nu = 0.48; // lambda \approx 25  // default
    const REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    const REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::SNH(mu, lambda);

    // build the time integrator
    //_solver = new TIMESTEPPER::BDF_1(*_tetMesh, *_hyperelastic);
    _solver = new TIMESTEPPER::BDF_2(*_tetMesh, *_hyperelastic);

    _solver->vertexFaceSelfCollisionsOn() = false;
    _solver->edgeEdgeSelfCollisionsOn() = false;

    const VECTOR3 center0(0.0,-5.0, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _solver->addKinematicCollisionObject(_kinematicShapes[0]);

    _eye    = VECTOR3(-2.5996, 0.52361, 0.286395);
    _lookAt = VECTOR3(-1.60313, 0.44686, 0.32046);
    _up     = VECTOR3(0.0765102, 0.997036, 0.00830762);

    // give it a big ballistic velocity
    VECTOR velocity(vertices.size() * 3);
    velocity.setZero();
    for (int x = 0; x < velocity.size() / 3; x++)
      velocity[3 * x + 1] = -3.0;

    _solver->velocity() = velocity;

    // if it's BDF-1, need to set the old velocity too
    TIMESTEPPER::BDF_1* setOld = dynamic_cast<TIMESTEPPER::BDF_1*>(_solver);
    if (setOld != NULL)
      setOld->velocityOld() = velocity;

    // if it's BDF-2, need to set the old velocity too
    TIMESTEPPER::BDF_2* setOlder = dynamic_cast<TIMESTEPPER::BDF_2*>(_solver);
    if (setOlder != NULL)
    {
      setOlder->velocityOld() = velocity;
      setOlder->velocityOlder() = velocity;
    }

    _pauseFrame = 250;

    return true;
  }
 
// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene() override
  {
    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    glDisable(GL_DEPTH_TEST);
    glPointSize(10.0);
    drawVertices(*_tetMesh, _solver->constrainedNodes());
    drawPlaneConstraints(_tetMesh, _solver);
  };
#endif

};

}

#endif
