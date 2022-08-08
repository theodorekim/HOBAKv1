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

#include "Scenes/JSON_SCENE.h"
#include "util/DRAW_GL.h"

#include <snapshot.h>

using namespace HOBAK;
using namespace std;

GLVU glvu;
FFMPEG_MOVIE movie;

bool animate = false;
bool singleStep = false;
int timestep = -1;

//bool bodyCentered = true;
bool bodyCentered = false;

bool drawingFirstMesh = true;
bool drawingSecondMesh = true;

// animation should pause on this frame;
int pauseFrame = -1;

// scene being replayed
JSON_SCENE* scene0 = NULL;
JSON_SCENE* scene1 = NULL;

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glvu.BeginFrame();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if (!bodyCentered)
    {
      if (drawingFirstMesh)
        scene0->drawScene();

      if (drawingSecondMesh && drawingFirstMesh)
        scene1->drawTetMesh();

      // make sure to draw the kinematics if the first scene isn't
      if (drawingSecondMesh && !drawingFirstMesh)
        scene1->drawScene();

      // draw differences explicitly
      drawDiff(*scene0->tetMesh(), *scene1->tetMesh());
    }
    else
    {
      if (drawingFirstMesh)
        scene0->drawBodyCenteredScene();
      if (drawingSecondMesh)
        scene1->drawBodyCenteredScene(VECTOR3(0,0,0.5), VECTOR3(1,1,1));
    }
  glvu.EndFrame();
}

///////////////////////////////////////////////////////////////////////
// print out the diff norms
///////////////////////////////////////////////////////////////////////
void printDiffs()
{
  const int frame = scene0->drawFrame();
  const VECTOR positionDiff = scene0->positions()[frame] - scene1->positions()[frame];
  const VECTOR velocityDiff = scene0->velocities()[frame] - scene1->velocities()[frame];
  cout << " position diff: " << positionDiff.norm() / scene0->positions()[frame].norm() << endl;
  cout << " velocity diff: " << velocityDiff.norm() / scene0->velocities()[frame].norm() << endl;
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate) 
  {
    scene0->stepSimulation();
    scene1->stepSimulation();

    printDiffs();

    if (singleStep) 
    {
      animate = false;
      singleStep = false;
    }
    timestep++;
  }
  
  if (timestep == scene0->pauseFrame() && animate)
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
  case 27:
  case 'Q':
  case 'q':
    exit(0);
    break;
  case 'a':
    animate = !animate;
    break;
  case 'd':
    scene0->drawFeature() = !scene0->drawFeature();
    scene1->drawFeature() = !scene1->drawFeature();
    break;
  case '1':
    drawingFirstMesh = !drawingFirstMesh;
    return;
    break;
  case '2':
    drawingSecondMesh = !drawingSecondMesh;
    return;
    break;
  case 'f':
    if ((!drawingFirstMesh && !drawingSecondMesh) ||
        (drawingFirstMesh && drawingSecondMesh)) 
    {
      drawingFirstMesh = false;
      drawingSecondMesh = true;
    }
    drawingFirstMesh = !drawingFirstMesh;
    drawingSecondMesh = !drawingSecondMesh;

    if (drawingFirstMesh)
      cout << " Drawing scene 0 " << endl;
    else
      cout << " Drawing scene 1 " << endl;
    return;
    break;
  case ' ':
    animate = true;
    singleStep = true;
    //TIMER::printTimings();
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
    if (animate) animate = false;
    scene0->drawFrame()--;
    scene1->drawFrame()--;

    if (scene0->drawFrame() < 0)
    {
      scene0->drawFrame() = 0;
      scene1->drawFrame() = 0;
    }
    break;
  case GLUT_KEY_RIGHT:
    if (animate) animate = false;
    scene0->drawFrame()++;
    scene1->drawFrame()++;
    break;
  case GLUT_KEY_UP:
    scene0->arrowCounter()++;
    scene1->arrowCounter()++;
    break;
  case GLUT_KEY_DOWN:
    scene0->arrowCounter()--;
    scene1->arrowCounter()--;
    break;
  case GLUT_KEY_HOME:
    break;
  case GLUT_KEY_END:
    break;
  }
  cout << " Arrow counter: " << scene0->arrowCounter() << endl;
  cout << " Draw frame: " << scene0->drawFrame() << endl;
  scene0->updatePositions();
  scene1->updatePositions();
  printDiffs();
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
  char title[] = "Replay Scene";
  glvu.Init(title, GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE, 0, 0, 800, 800);

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

  //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glPolygonMode(GL_FRONT, GL_FILL);
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

  const VECTOR3& eye    = scene0->eye();
  const VECTOR3& lookAt = scene0->lookAt();
  const VECTOR3& up     = scene0->up();
  const VECTOR3& worldCenter = scene0->worldCenter();

  for (unsigned int x = 0; x < 3; x++)
  {
    Eye[x] = eye[x];
    LookAtCntr[x] = lookAt[x];
    Up[x] = up[x];
    WorldCenter[x] = worldCenter[x];
  }

  // if we're just body-centered, snap to origin
  if (bodyCentered)
  {
    //Eye = glvuVec3f(1, 0, 0);
    //LookAtCntr = glvuVec3f(0, 0, 0);
    //Up = glvuVec3f(0, 1, 0);

    Eye = glvuVec3f(-0.971139, 0.715275, 1.80284);
    LookAtCntr = glvuVec3f(-0.531403, 0.385029, 0.967634);
    Up = glvuVec3f(0.107162, 0.942592, -0.316286);
    WorldCenter = glvuVec3f(0, 0, 0);
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
bool validateScenes()
{
  // do the tet meshes have the same number of vertices?
  const int vertices0 = scene0->tetMesh()->vertices().size();
  const int vertices1 = scene1->tetMesh()->vertices().size();

  if (vertices0 != vertices1)
  {
    cout << " Vertex counts don't match! Can't do a diff. " << endl;
    cout << " Scene 0 vertices: " << vertices0 << endl;
    cout << " Scene 1 vertices: " << vertices1 << endl;
    return false;
  }

  // do the scenes have the same number of frames?
  const int frames0 = scene0->positions().size();
  const int frames1 = scene1->positions().size();

  if (frames0 != frames1)
  {
    cout << " Number of frames don't match! Can't do a diff. " << endl;
    cout << " Scene 0 frames: " << frames0 << endl;
    cout << " Scene 1 frames: " << frames1 << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    cout << " USAGE: " << argv[0] << " <JSON output from simulateScene 0> <JSON output from simulateScene 1>" << endl;
    return 0;
  }

  // load up the first scene
  scene0 = new JSON_SCENE();
  scene0->printSceneDescription();
  cout << "=====================================================================" << endl;
  cout << " Loading scene 0: " << argv[1] << endl;
  readSceneJSON(argv[1], *scene0);

  bool success = scene0->buildScene();
  if (!success)
  {
    cout << " Failed to build scene 0 -> Exiting ..." << endl;
    return 1;
  }

  // load up the second scene
  scene1 = new JSON_SCENE();
  cout << "=====================================================================" << endl;
  cout << " Loading scene 1: " << argv[2] << endl;
  readSceneJSON(argv[2], *scene1);
  success = scene1->buildScene();
  if (!success)
  {
    cout << " Failed to build scene 1-> Exiting ..." << endl;
    return 1;
  }

  // make sure the scenes are even comparable
  success = validateScenes();

  glutInit(&argc, argv);
  glvuWindow();

  return 0;
}
