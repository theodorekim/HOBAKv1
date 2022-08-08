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
#include <iostream>
#include <random>
#include "util/FFMPEG_MOVIE.h"
#include "util/FILE_IO.h"

#include "Scenes/BUNNY_DROP.h"
#include "Scenes/BUNNY_DROP_DEBUG.h"
#include "Scenes/PINNED_C.h"
#include "Scenes/DROPPED_E.h"
#include "Scenes/TWO_TET_DROP_VF.h"
#include "Scenes/TWO_TET_DROP_EE.h"
#include "Scenes/TWO_TET_KISS_VF.h"
#include "Scenes/TWO_TET_KISS_EE.h"
#include "Scenes/DIHEDRAL_TEST.h"
#include "Scenes/FORK_CCD_TEST.h"
#include "Scenes/DROPPED_C_CCD_TEST.h"
#include "Scenes/ARM_BEND.h"
#include "Scenes/KINEMATICS_TEST.h"
#include "Scenes/CRUSH_TEST.h"
#include "Scenes/DROP_TEST.h"
#include "Scenes/BDF_TEST.h"
#include "Scenes/NEWMARK_STRETCH.h"
#include "Scenes/QUASISTATIC_COLLISIONS.h"
#include "Scenes/QUASISTATIC_STRETCH.h"
#include "util/DRAW_GL.h"

#include <snapshot.h>

using namespace HOBAK;
using namespace std;

GLVU glvu;
FFMPEG_MOVIE movie;

bool animate = true;
bool singleStep = false;
int timestep = -1;

// animation should pause on this frame;
int pauseFrame = -1;

// scene being simulated
SIMULATION_SCENE* scene = NULL;

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glvu.BeginFrame();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    scene->drawScene();
    if (animate)
      movie.addFrameGL();
  glvu.EndFrame();
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate) 
  {
    scene->stepSimulation();
    recordFrameToJSON(*scene);
    if (singleStep) 
    {
      animate = false;
      singleStep = false;
    }
    timestep++;
  }
  
  if (timestep == scene->pauseFrame() && animate)
  {
    cout << " Hit pause frame specific in scene file " << endl;
    animate = false;
  }

  glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 27:  // escape key
  case 'Q':
    TIMER::printTimings();
    exit(0);
    break;
  case 'q':
    TIMER::printTimings();
    movie.streamWriteMovie(scene->movieName().c_str());
    writeSceneJSON(scene->jsonName().c_str());
    exit(0);
    break;
  case 'a':
    animate = !animate;
    break;
  case 'd':
    scene->drawFeature() = !scene->drawFeature();
    break;
  case ' ':
    animate = true;
    singleStep = true;
    TIMER::printTimings();
    break;
  case 'v': {
      Camera *camera = glvu.GetCurrentCam();
      glvuVec3f eye;
      glvuVec3f ref;
      glvuVec3f up;
      camera->GetLookAtParams(&eye, &ref, &up);
      cout << "    _eye    = VECTOR3(" << eye[0] << ", " << eye[1] << ", " << eye[2] << ");" << endl;
      cout << "    _lookAt = VECTOR3(" << ref[0] << ", " << ref[1] << ", " << ref[2] << ");" << endl;
      cout << "    _up     = VECTOR3(" << up[0] << ", " << up[1] << ", " << up[2] << ");" << endl;
    } 
    break;
  default:
    break;
  }

  glvu.Keyboard(key, x, y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key) {
  case GLUT_KEY_LEFT:
    scene->leftArrow() = true;
    break;
  case GLUT_KEY_RIGHT:
    scene->rightArrow() = true;
    break;
  case GLUT_KEY_UP:
    scene->arrowCounter()++;
    cout << "Arrow counter: " << scene->arrowCounter() << endl;
    break;
  case GLUT_KEY_DOWN:
    scene->arrowCounter()--;
    cout << "Arrow counter: " << scene->arrowCounter() << endl;
    break;
  case GLUT_KEY_HOME:
    break;
  case GLUT_KEY_END:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y) 
{ 
  glvu.Mouse(button, state, x, y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y) 
{ 
  glvu.Motion(x, y);
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  const string title = string("Simulating scene: ") + scene->sceneName();
  char buffer[title.length() + 1];
  strcpy(buffer, title.c_str());
  glvu.Init(buffer, GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE, 0, 0, 800, 800);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = { 10.0, 4.0, 10.0, 1.0 };
  GLfloat lightZeroColor[] = { 1.0, 1.0, 0.5, 1.0 };
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);

  GLfloat lightPosition1[] = { -4.0, 1.0, 1.0, 1.0 };
  GLfloat lightColor1[] = { 1.0, 0.0, 0.0, 1.0 };
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT1, GL_POSITION, lightPosition1);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);
  glClearColor(1, 1, 1, 0);

  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);

  glEnable(GL_MULTISAMPLE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glvuVec3f ModelMin(-10, -10, -10), ModelMax(10, 10, 10);
  glvuVec3f Eye, LookAtCntr, Up, WorldCenter;

  const VECTOR3& eye    = scene->eye();
  const VECTOR3& lookAt = scene->lookAt();
  const VECTOR3& up     = scene->up();
  const VECTOR3& worldCenter = scene->worldCenter();

  for (unsigned int x = 0; x < 3; x++)
  {
    Eye[x] = eye[x];
    LookAtCntr[x] = lookAt[x];
    Up[x] = up[x];
    WorldCenter[x] = worldCenter[x];
  }

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov, Aspect, Near, Far);
  glvu.SetWorldCenter(WorldCenter);

  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  // self-collision tests
  //scene = new ARM_BEND();
  scene = new BUNNY_DROP();
  //scene = new BUNNY_DROP_DEBUG();
  //scene = new PINNED_C();
  //scene = new DROPPED_E();

  // unit tests
  //scene = new TWO_TET_DROP_VF();
  //scene = new TWO_TET_DROP_EE();
  //scene = new TWO_TET_KISS_VF();
  //scene = new TWO_TET_KISS_EE();
  
  // debugging scenes
  //scene = new DIHEDRAL_TEST();

  // CCD tests that we still don't pass
  //scene = new FORK_CCD_TEST();
  //scene = new DROPPED_C_CCD_TEST();

  // kinematics tests
  //scene = new KINEMATICS_TEST();
  //scene = new CRUSH_TEST();
  //scene = new DROP_TEST();
  //scene = new BDF_TEST();
  //scene = new NEWMARK_STRETCH();
  //scene = new QUASISTATIC_COLLISIONS();
  //scene = new QUASISTATIC_STRETCH();

  // make sure to print out what we should expect from this scene
  scene->printSceneDescription();

  bool success = scene->buildScene();

  // prepare scene for JSON output
  initializeSceneJSON(*scene);

  if (!success)
  {
    cout << " Failed to build scene-> Exiting ..." << endl;
    return 1;
  }

  glutInit(&argc, argv);
  glvuWindow();

  return 0;
}
