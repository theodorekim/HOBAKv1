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
#ifndef SIMULATION_SCENE_H
#define SIMULATION_SCENE_H

#ifndef GL_DISABLED
#include "util/DRAW_GL.h"
#endif

#include "Geometry/CAPSULE.h"
#include "Geometry/CUBE.h"
#include "Geometry/SPHERE.h"
#include "Geometry/TET_MESH.h"
#include "Geometry/TET_MESH_FASTER.h"
#include "Hyperelastic/Volume/SNH.h"
#include "Hyperelastic/Volume/STVK.h"
#include "Hyperelastic/Volume/ARAP.h"
#include "Hyperelastic/Volume/LINEAR.h"
#include "Timestepper/BACKWARD_EULER_VELOCITY.h"
#include "Timestepper/BACKWARD_EULER_POSITION.h"
#include "Timestepper/NEWMARK.h"
#include "Timestepper/QUASISTATIC.h"
#include "util/TIMER.h"

namespace HOBAK {

class SIMULATION_SCENE {
public:

  // initialize the scene
  SIMULATION_SCENE() { 
    _pauseFrame = -2; 
    _arrowCounter = -1;
    _leftArrow = false; 
    _rightArrow = false; 
    _drawFeature = false; 
    _frameNumber = 0;  
    _normalizedVertices = false;
    _sceneName = std::string("default");
    _initialA = MATRIX3::Identity();
    _initialTranslation = VECTOR3::Zero();

    _solver = NULL;
    _tetMesh = NULL;
    _hyperelastic = NULL;
    _gravity.setZero();
  };

  virtual ~SIMULATION_SCENE() 
  {
    delete _tetMesh;
    delete _solver;
    delete _hyperelastic;

    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      delete _kinematicShapes[x];
  };

  // TODO: Build the actual scene. You have to implement this!
  virtual bool buildScene() = 0;

  // TODO: Describe the scene build built. You have to do this!
  virtual void printSceneDescription() = 0;

  const VECTOR3& eye() const         { return _eye; };
  const VECTOR3& lookAt() const      { return _lookAt; };
  const VECTOR3& up() const          { return _up; };
  const VECTOR3& worldCenter() const { return _worldCenter; };
  const VECTOR3& gravity() const  { return _gravity; };
  const int& pauseFrame() const   { return _pauseFrame; };
  const int& frameNumber() const  { return _frameNumber; };
  const bool& drawFeature() const { return _drawFeature; };
  const int& arrowCounter() const { return _arrowCounter; };
  const string& tetMeshFilename() const { return _tetMeshFilename; };
  const vector<KINEMATIC_SHAPE*>& kinematicShapes() const { return _kinematicShapes; };
  const bool& normalizedVertices() const { return _normalizedVertices; };

  const MATRIX3& initialA() const           { return _initialA; };
  const VECTOR3& initialTranslation() const { return _initialTranslation; };
  
  VECTOR3& eye()         { return _eye; };
  VECTOR3& lookAt()      { return _lookAt; };
  VECTOR3& up()          { return _up; };
  VECTOR3& worldCenter() { return _worldCenter; };
  VECTOR3& gravity()     { return _gravity; };
  int& pauseFrame()      { return _pauseFrame; };
  int& frameNumber()     { return _frameNumber; };
  bool& drawFeature()    { return _drawFeature; };
  int& arrowCounter()    { return _arrowCounter; };
  bool& leftArrow()      { return _leftArrow; };
  bool& rightArrow()     { return _rightArrow; };
  string& tetMeshFilename()  { return _tetMeshFilename; };
  bool& normalizedVertices() { return _normalizedVertices; };
  vector<KINEMATIC_SHAPE*>& kinematicShapes() { return _kinematicShapes; };

  const string sceneName() const { return _sceneName; };
  const string movieName() const { return _sceneName + std::string(".mov"); };
  const string jsonName() const  { return _sceneName + std::string(".json"); };

  const TIMESTEPPER::TIMESTEPPER* solver() const { return _solver; };
  TIMESTEPPER::TIMESTEPPER* solver()             { return _solver; };
  const TET_MESH* tetMesh() const { return _tetMesh; };
  TET_MESH* tetMesh()             { return _tetMesh; };

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) {
    _solver->externalForces().setZero();
    _solver->addGravity(_gravity);
    _solver->solve(verbose);

    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
  };
#endif

protected:
  // set the positions to previous timestep, in case the user wants to 
  // look at that instead of the current step
  void setToPreviousTimestep()
  {
    const VECTOR& old = _solver->positionOld();
    _tetMesh->setDisplacement(old);
  }
  
  // restore positions from previous timestep, in case the user just drew
  // the previous timestep, but now we want the state to be consistent
  // when drawing the next frame
  void restoreToCurrentTimestep()
  {
    const VECTOR& current = _solver->position();
    _tetMesh->setDisplacement(current);
  }

  // scene geometry
  TET_MESH* _tetMesh;
  vector<KINEMATIC_SHAPE*> _kinematicShapes;

  // solver and materials
  TIMESTEPPER::TIMESTEPPER* _solver;
  VOLUME::HYPERELASTIC* _hyperelastic;

  // simulation parameters
  VECTOR3 _gravity;

  // initial rotation-scale and translation of tet mesh
  MATRIX3 _initialA;
  VECTOR3 _initialTranslation;

  // drawing parameters
  int _pauseFrame;

  // counter that can be incremented and decremented by the arrow keys
  int _arrowCounter;

  // bools that can be toggled by the arrow keys
  bool _leftArrow;
  bool _rightArrow;

  // flag for whether or not to draw some user-specific feature
  bool _drawFeature;

  // what frame are we on?
  int _frameNumber;

  // tet mesh filename
  string _tetMeshFilename;

  // did we normalize the vertices when we read them in?
  bool _normalizedVertices;

  // camera parameters
  VECTOR3 _eye;
  VECTOR3 _lookAt;
  VECTOR3 _up;
  VECTOR3 _worldCenter;

  // scene name, used to determine the JSON and MOV filenames
  std::string _sceneName;
};

}

#endif
