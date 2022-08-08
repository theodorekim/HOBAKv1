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
#ifndef ARM_BEND_H
#define ARM_BEND_H

#include "SIMULATION_SCENE.h"

namespace HOBAK {

class ARM_BEND : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Bending a cylinder like an arm, to see what the self-collisions" << endl;
    cout << " look like. " << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "arm_bend";

    // read in the tet mesh file
    //_tetMeshFilename = string("../data/arm_bend/arm_bend_0.01.tobj");
    _tetMeshFilename = string("../data/arm_bend/arm_bend_0.02.tobj"); // small enough to fit in github
    //_tetMeshFilename = string("../data/arm_bend/arm_bend_0.03.tobj");
    //_tetMeshFilename = string("../data/arm_bend/arm_bend_0.035.tobj");
    //_tetMeshFilename = string("../data/cylinder_0.03.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);
    vertices = TET_MESH::normalizeVertices(vertices);

    _normalizedVertices = true;

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::SNH(1.0, 10.0);
    _tetMesh->setCollisionEps(0.005);

    // build the time integrator
    _solver = new TIMESTEPPER::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic);
    //_solver = new TIMESTEPPER::BACKWARD_EULER_VELOCITY(*_tetMesh, *_hyperelastic);
    //_solver = new TIMESTEPPER::NEWMARK(*_tetMesh, *_hyperelastic);

    // radius of the bones in the arm
    const REAL boneRadius = 0.05;
    VECTOR3 shoulderCenter(0.64, 0.71, 0.64);
    VECTOR3 forearmCenter(0.64, 1.3, 0.64);
    //const REAL boneLength = 0.25; // seems to work well
    const REAL boneLength = 0.35; // does not create a GIA-necessary kinematic situation
    //const REAL boneLength = 0.4;// cylinders still overlap at the end
    //const REAL boneLength = 0.5; // bugs become ambiguous -- is is because GIA is needed?
    _kinematicShapes.push_back(new CYLINDER(shoulderCenter, boneRadius, boneLength));
    _kinematicShapes.push_back(new CYLINDER(forearmCenter, boneRadius, boneLength));

    // attach it to the rest of the scene
    _solver->attachKinematicConstraints(_kinematicShapes[0]);
    _solver->attachKinematicConstraints(_kinematicShapes[1]);

    const REAL collisionMu = 1000.0;  // default
    _solver->collisionStiffness() = collisionMu;
    _solver->collisionDampingBeta() = 0.01; // default

    _solver->edgeEdgeSelfCollisionsOn() = true;
    _solver->edgeEdgeSelfCollisionsOn() = true;

    // look at the overall arm
    _eye    = VECTOR3(2.15236, 1.03626, 0.813244);
    _lookAt = VECTOR3(1.15282, 1.02777, 0.784225);
    _up     = VECTOR3(-0.029005, -0.00184828, 0.999576);

    _pauseFrame = 300;

    // for this one, I'd like to see the bounding box
    VECTOR3 mins, maxs;
    _tetMesh->getBoundingBox(mins, maxs);
    cout << " Bounding box mins: " << mins.transpose() << endl;
    cout << " Bounding box maxs: " << maxs.transpose() << endl;

    cout << " center: " << ((mins + maxs) * 0.5).transpose() << endl;

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    static int timestep = 0;
    if (timestep < 250)
    {
      KINEMATIC_SHAPE& forearm = *_kinematicShapes[1];
      forearm.rotation() = Eigen::AngleAxisd(0.01, VECTOR3::UnitX()) * forearm.rotation();
      VECTOR3 world = forearm.localVertexToWorld(VECTOR3(0,-0.25, 0));
      VECTOR3 pin(0.64, 1.05, 0.64);
      forearm.translation() -= world - pin;
    }

    _solver->externalForces().setZero();
    _solver->addGravity(_gravity);
    _solver->solve(verbose);

    timestep++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene() override
  {
    setToPreviousTimestep();

    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
   
    glDisable(GL_DEPTH_TEST);
    glColor4f(0.2, 0.0, 1.0, 0.5);
    drawCylinder(*(CYLINDER*)(_kinematicShapes[0]));
    glColor4f(1.0, 0.0, 0.2, 0.5);
    drawCylinder(*(CYLINDER*)(_kinematicShapes[1]));
    
    glPointSize(10.0);
    glColor4f(0.0, 1.0, 0.0, 1.0);
    drawVertices(*_tetMesh, _solver->constrainedNodes());

    glEnable(GL_DEPTH_TEST);
    if (_drawFeature)
      drawVertexFacePairs(*_tetMesh, _arrowCounter);

    //drawSurfaceFaceCollisionCell(*_tetMesh, 718);
    //drawSurfaceFaceCollisionCell(*_tetMesh, 727);
    
    restoreToCurrentTimestep();
  };
#endif

};

}

#endif
