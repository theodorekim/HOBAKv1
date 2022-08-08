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
#ifndef CRUSH_TEST_H
#define CRUSH_TEST_H

#include "SIMULATION_SCENE.h"
#include "Damping/Volume/GREEN_DAMPING.h"
#include "Hyperelastic/Volume/SNH_WITH_BARRIER.h"

namespace HOBAK {

class CRUSH_TEST : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Smashing a cube down to see what weird shapes we can get out of " << endl;
    cout << " the Newton solver." << endl;
    cout << endl;
    cout << " This one's just for fun -- there's no right or wrong setting. " << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "crush_test";
    // read in the tet mesh file
    _tetMeshFilename = string("../data/cube_40.tobj");

    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    _gravity = VECTOR3(0, 0, 0);

    REAL E = 3.0;
    //REAL E = 10.0;
    //REAL nu = 0.495; // lambda \approx 100
    //REAL nu = 0.49; // lambda \approx 50
    //REAL nu = 0.4875; // lambda \approx 40
    //REAL nu = 0.48; // lambda \approx 25
    //REAL nu = 0.475; // lambda \approx 20
    //REAL nu = 0.45; // lambda \approx 10
    REAL nu = 0.4; // lambda \approx 4.2

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    //_hyperelastic = new VOLUME::LINEAR(10.0, 1.0);
    //_hyperelastic = new VOLUME::ARAP(10.0, 1.0);
    //_hyperelastic = new VOLUME::STVK(10.0, 100.0);
    //_hyperelastic = new VOLUME::SNH(mu, lambda);
    _hyperelastic = new VOLUME::SNH_WITH_BARRIER(mu, lambda);

    //damping = new VOLUME::GREEN_DAMPING(0.001);
    //damping = new VOLUME::GREEN_DAMPING(0.005);
    _damping = new VOLUME::GREEN_DAMPING(0.02);
    // build the time integrator
    //_solver = new TIMESTEPPER::NEWMARK(*_tetMesh, *_hyperelastic, *_damping);
    _solver = new TIMESTEPPER::NEWMARK(*_tetMesh, *_hyperelastic);
    //solver = new TIMESTEPPER::BACKWARD_EULER(*tetMesh, *hyperelastic);

    //solver->dt() *= 0.001;
    //solver->dt() *= 0.01;
    _solver->setDt(_solver->dt() * 0.01);
    //((TIMESTEPPER::NEWMARK*)solver)->maxNewtonIterations() = 20;
    ((TIMESTEPPER::NEWMARK*)_solver)->maxNewtonIterations() = 3;

    /*
    // crush it
    const REAL scalingAmount = 0.01;
    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      tetMesh->vertices()[x] *= scalingAmount;
      tetMesh->vertices()[x] += VECTOR3(1,1,1);
    }
    */

    /*
    // reflect it
    const REAL scalingAmount = -1;
    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      tetMesh->vertices()[x][1] *= scalingAmount;
      tetMesh->vertices()[x] += VECTOR3(1,1,1);
    }
    */

    // inflate it
    const REAL scalingAmount = 2.0;
    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      _tetMesh->vertices()[x] *= scalingAmount;
      //_tetMesh->vertices()[x] += VECTOR3(1,1,1);
    }

    _solver->position() = _tetMesh->getDisplacement();

    _eye    = VECTOR3(-2.10574, 0.885076, -0.74685);
    _lookAt = VECTOR3(-1.203, 0.830111, -0.32019);
    _up     = VECTOR3(0.0512588, 0.998481, 0.020177);

    _eye    = VECTOR3(-0.779069, 1.08268, 0.326908);
    _lookAt = VECTOR3(0.150359, 1.0276, 0.691778);
    _up     = VECTOR3(0.0512582, 0.998482, 0.0201768);

    _eye    = VECTOR3(-0.715487, 0.682568, 0.663179);
    _lookAt = VECTOR3(0.218776, 0.627493, 1.01548);
    _up     = VECTOR3(0.0512521, 0.998481, 0.0201745);

    _eye    = VECTOR3(-1.49481, 1.88272, -0.164984);
    _lookAt = VECTOR3(-0.623691, 1.58973, 0.229106);
    _up     = VECTOR3(0.269791, 0.956093, 0.114442);

    /*
    // point
    _eye    = VECTOR3(-0.146893, 1.10591, 0.619131);
    _lookAt = VECTOR3(0.782535, 1.05083, 0.984001);
    _up     = VECTOR3(0.0512522, 0.998482, 0.0201745);

    // line
    _eye    = VECTOR3(-0.20828, 1.56163, 0.573036);
    _lookAt = VECTOR3(0.721148, 1.50655, 0.937906);
    _up     = VECTOR3(0.0512523, 0.998482, 0.0201746);
    */

    // crush it

    return true;
  }
#if 0
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "crush_test";
    // read in the tet mesh file
    _tetMeshFilename = string("../data/cube_6.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    // build the tet mesh object
    _tetMesh = new TET_MESH(vertices, tets);
    _hyperelastic = new VOLUME::ARAP(1.0, 1.0);

    // build the time integrator
    _solver = new TIMESTEPPER::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic);

    // add the kinematics
    const VECTOR3 center0(0.56,0.48, -0.25);
    const VECTOR3 center1(0.56,0.48,1.23);
    _kinematicShapes.push_back(new CUBE(center0, 1.0));
    _kinematicShapes.push_back(new CUBE(center1, 1.0));

    // attach it to the rest of the scene
    _solver->attachKinematicSurfaceConstraints(_kinematicShapes[1]);
    _solver->addKinematicCollisionObject(_kinematicShapes[0]);

    _eye    = VECTOR3(-1.61233, 0.617816, -0.169647);
    _lookAt = VECTOR3(-0.630224, 0.573968, 0.0135059);
    _up     = VECTOR3(0.0418256, 0.999014, 0.0148915);
    return true;
  }
#endif
 
  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    _solver->externalForces().setZero();
    _solver->addGravity(_gravity);
    _solver->solve(verbose);

    _frameNumber++;
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

private:
  VOLUME::DAMPING* _damping;
};

}

#endif
