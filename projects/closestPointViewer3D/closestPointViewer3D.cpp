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

#include "FIELD_2D.h"
#include "FIELD_3D.h"
#include "VECTOR3_FIELD_3D.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h> // OpenGL itself.
#include <GL/glu.h> // GLU support library.
#include <GL/glut.h> // GLUT support library.
#endif

#include <iostream>
using namespace std;
using namespace HOBAK;

// the original field data
VECTOR3_FIELD_3D vectorField;
VECTOR3_FIELD_3D displayVectorField;

// the distance field data
FIELD_3D distanceField;

// z slice currently viewing viewed
int zSlice = -1;

string windowLabel("Closest Point Viewer");

// scaled, biased and manipulated field
FIELD_2D viewingField;

// FFT of the original field
FIELD_2D fftReal;
FIELD_2D fftIm;

int xScreenRes = 800;
int yScreenRes = 800;
int xMouse, yMouse;
int mouseButton;
int mouseState;
int xField = -1;
int yField = -1;
float zoom = 1.0;

float scale = 1.0;
float bias = 0.0;

bool drawingGrid = false;
bool useAbsolute = false;
bool useLog = false;
bool normalized = false;
float oldScale = 1.0;
float oldBias = 0.0;

VECTOR3 eyeCenter(0.5, 0.5, 1);
VECTOR3_FIELD_3D centers;
VECTOR3_FIELD_3D displayCenters;

//////////////////////////////////////////////////////////////////////
// draw to GL, but assume the field represents closest points,
// so don't add the origin to the endpoint
//////////////////////////////////////////////////////////////////////
void drawClosestPointZSlice(const VECTOR3_FIELD_3D& display, 
                            const VECTOR3_FIELD_3D& origins, 
                            const int zSlice, const REAL scale, const int stride)
{
  float radius = 0.00125;
  const int xRes = display.xRes();
  const int yRes = display.yRes();
  const int zRes = display.zRes();
  
  glBegin(GL_LINES);
  for (int z = 0; z < zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < yRes; y+=stride)
      for (int x = 0; x < xRes; x+=stride)
      {
        const VECTOR3& origin = origins(x,y,z);
        VECTOR3 endpoint = display(x,y,z);

        //if (norm2(endpoint) < 1e-5)
        if (endpoint.norm() < 1e-5)
          continue;

        //glColor4f(1,1,1,1);
        glColor4f(0.5,0.5,0.5,1);
        glVertex3f(origin[0], origin[1], origin[2]);
        glVertex3f(endpoint[0], endpoint[1], 0);
      }
  }
  glEnd();

  for (int z = 0; z < zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < yRes; y+=stride)
      for (int x = 0; x < xRes; x+=stride)
      {
        const VECTOR3& origin = origins(x,y,z);
        const VECTOR3 endpoint = display(x,y,z);

        //if (norm2(endpoint) < 1e-5)
        if (endpoint.norm() < 1e-5)
          continue;

        glColor4f(1,0,0,1);
        glPushMatrix();
        glTranslatef(origin[0], origin[1], origin[2]);
        glutSolidSphere(radius, 10,10);
        glPopMatrix();
        
        glColor4f(0,0,1,1);
        glPushMatrix();
        glTranslatef(endpoint[0], endpoint[1], endpoint[2]);
        glutSolidSphere(radius, 10,10);
        glPopMatrix();

      }
  }
}

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
  //glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  for (unsigned int x = 0; x < output.size(); x++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
}

///////////////////////////////////////////////////////////////////////
// initialize the GL texture
///////////////////////////////////////////////////////////////////////
void initTexture(FIELD_2D& texture)
{
  float* rgb = new float[3 * texture.totalCells()];

  for (int x = 0; x < texture.totalCells(); x++)
  {
    rgb[3 * x] = rgb[3 * x + 1] = rgb[3 * x + 2] = 0;

    if (distanceField[zSlice * distanceField.slabSize() + x] > 0)
      rgb[3 * x] = texture[x];
    else
      rgb[3 * x + 2] = texture[x];
  }

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, 
      texture.xRes(), 
      texture.yRes(), 0, 
      GL_RGB, GL_FLOAT, 
      rgb);

  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glEnable(GL_TEXTURE_2D);

  delete[] rgb;
}

///////////////////////////////////////////////////////////////////////
// update the texture being viewed
///////////////////////////////////////////////////////////////////////
void updateViewingTexture()
{
  viewingField = distanceField.zSlice(zSlice);

  if (useAbsolute)
    viewingField.abs();

  viewingField += bias;
  viewingField *= scale;

  if (useLog)
    viewingField.log(10.0);

  initTexture(viewingField);
}

///////////////////////////////////////////////////////////////////////
// draw a grid over everything
///////////////////////////////////////////////////////////////////////
void drawGrid()
{
  glColor4f(0.1, 0.1, 0.1, 1.0);

  int xRes = distanceField.xRes();
  int yRes = distanceField.yRes();

  float dx = 1.0 / xRes;
  float dy = 1.0 / yRes;

  if (xRes < yRes)
    dx *= (REAL)xRes / yRes;
  if (xRes > yRes)
    dy *= (REAL)yRes / xRes;

  glBegin(GL_LINES);
  for (int x = 0; x < distanceField.xRes() + 1; x++)
  {
    glVertex3f(x * dx, 0, 1);
    glVertex3f(x * dx, 1, 1);
  }
  for (int y = 0; y < distanceField.yRes() + 1; y++)
  {
    glVertex3f(0, y * dy, 1);
    glVertex3f(1, y * dy, 1);
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the velocity field
///////////////////////////////////////////////////////////////////////
void drawVectorField()
{
  // cache the cell centers
  REAL xLength = 1.0;
  REAL yLength = 1.0;

  int xRes = distanceField.xRes();
  int yRes = distanceField.yRes();

  if (xRes < yRes)
    xLength = (REAL)xRes / yRes;
  if (yRes < xRes)
    yLength = (REAL)yRes / xRes;

  glPushMatrix();

    glTranslatef(0.5 - (1 - xLength) * 0.5, 0.5 - (1 - yLength) * 0.5, 0);

    drawClosestPointZSlice(displayVectorField, displayCenters, zSlice, 0.1, 1);
    //displayVectorField.drawClosestPointZSlice(displayCenters, zSlice, 0.1, 1);
    //displayVectorField.drawClosestPointZSlice(displayCenters, zSlice, 0.1, 2);
    //displayVectorField.drawClosestPointZSlice(displayCenters, zSlice, 0.1, 4);
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  float halfZoom = zoom * 0.5;

  glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);

  // set the matrix mode back to modelview
  glMatrixMode(GL_MODELVIEW);

  // set the lookat transform
  glLoadIdentity();
  gluLookAt(eyeCenter[0], eyeCenter[1], 1,  // eye
            eyeCenter[0], eyeCenter[1], 0,  // center 
            0, 1, 0);   // up

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  REAL xLength = 1.0;
  REAL yLength = 1.0;

  int xRes = distanceField.xRes();
  int yRes = distanceField.yRes();
  if (xRes < yRes)
    xLength = (REAL)xRes / yRes;
  if (yRes < xRes)
    yLength = (REAL)yRes / xRes;

  glEnable(GL_TEXTURE_2D); 
  glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex3f(0.0, 0.0, 0.0);
    glTexCoord2f(0.0, 1.0); glVertex3f(0.0, yLength, 0.0);
    glTexCoord2f(1.0, 1.0); glVertex3f(xLength, yLength, 0.0);
    glTexCoord2f(1.0, 0.0); glVertex3f(xLength, 0.0, 0.0);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  if (drawingGrid)
    drawGrid();

  drawVectorField();

  // if there's a valid field index, print it
  if (xField >= 0 && yField >= 0 &&
      xField < distanceField.xRes() && yField < distanceField.yRes())
  {
    glLoadIdentity();

    // must set color before setting raster position, otherwise it won't take
    glColor4f(0.0f, 1.0f, 0.0f, 1.0f);

    // normalized screen coordinates (-0.5, 0.5), due to the glLoadIdentity
    float halfZoom = 0.5 * zoom;
    glRasterPos3f(-halfZoom* 0.95, -halfZoom* 0.95, 0);

    // build the field value string
    char buffer[256];

    // add the integer index
    int index = xField + yField * distanceField.xRes() + zSlice * distanceField.slabSize();
    string fieldValue;
    sprintf(buffer, "(%i", xField);
    fieldValue = fieldValue + string(buffer);
    sprintf(buffer, "%i", yField);
    fieldValue = fieldValue + string(", ") + string(buffer);
    sprintf(buffer, "%i", zSlice);
    fieldValue = fieldValue + string(", ") + string(buffer) + string(") = \n");

    VECTOR3 centerValue = centers[index];
    sprintf(buffer, "%.3f", centerValue[0]);
    fieldValue = fieldValue + string("P(") + string(buffer);
    sprintf(buffer, "%.3f", centerValue[1]);
    fieldValue = fieldValue + string(",") + string(buffer);
    sprintf(buffer, "%.3f", centerValue[2]);
    fieldValue = fieldValue + string(",") + string(buffer) + string("), ");

    // add the vector field value
    VECTOR3 vectorValue = vectorField[index];
    sprintf(buffer, "%.3f", vectorValue[0]);
    fieldValue = fieldValue + string("CP(") + string(buffer);
    sprintf(buffer, "%.3f", vectorValue[1]);
    fieldValue = fieldValue + string(",") + string(buffer);
    sprintf(buffer, "%.3f", vectorValue[2]);
    fieldValue = fieldValue + string(",") + string(buffer) + string(")");

    // add the global position
    printGlString(fieldValue);
  }

  glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  glutDisplay();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void printCommands()
{
  cout << " q           - quits" << endl;
  cout << " left mouse  - pan around" << endl;
  cout << " right mouse - zoom in and out " << endl;
  cout << " right arrow - increase bias " << endl;
  cout << " left arrow  - decrease bias " << endl;
  cout << " up arrow    - increase scale " << endl;
  cout << " down arrow  - decrease scale " << endl;
  cout << " n           - normalize (auto-set scale and bias) " << endl;
  cout << " g           - throw a grid over the pixels " << endl;
  cout << " a           - take absolute value of cells " << endl;
  cout << " l           - take log of cells " << endl;
  cout << " r           - look at real component of FFT" << endl;
  cout << " i           - look at imaginary component of FFT" << endl;
  cout << " s           - look at spatial (non-FFT)" << endl;
}

///////////////////////////////////////////////////////////////////////
// normalize the texture
///////////////////////////////////////////////////////////////////////
void normalize()
{
  viewingField = distanceField.zSlice(zSlice);

  if (useAbsolute)
    viewingField.abs();

  float minFound = viewingField.min();
  float maxFound = viewingField.max();

  // cache the values in case we want to undo the normalization
  oldScale = scale;
  oldBias = bias;

  // go ahead and compute the normalized version
  bias = -minFound;
  scale = 1.0 / (maxFound - minFound);

  updateViewingTexture();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
      bias -= 0.01;
      break;
    case GLUT_KEY_RIGHT:
      bias += 0.01;
      break;
    case GLUT_KEY_UP:
      scale += 0.01;
      break;
    case GLUT_KEY_DOWN:
      scale -= 0.01;
      break;
  }
  cout << " scale: " << scale << " bias: " << bias << endl;

  updateViewingTexture();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'a':
      useAbsolute = !useAbsolute;
      updateViewingTexture();
      break;
    case 'l':
      useLog = !useLog;
      updateViewingTexture();
      break;
    case 'g':
      drawingGrid = !drawingGrid;
      break;
    case 'n':
      if (normalized)
      {
        scale = oldScale;
        bias = oldBias;
        updateViewingTexture();
      }
      else
        normalize();

      normalized = !normalized;
      break;
    case '?':
      printCommands();
      break;
    case 'z':
      zSlice++;
      if (zSlice > distanceField.zRes() - 1)
        zSlice = distanceField.zRes() - 1;
      updateViewingTexture();
      break;
    case 'Z':
      zSlice--;
      if (zSlice < 0)
        zSlice = 0;
      updateViewingTexture();
      break;
    case 'q':
      exit(0);
      break;
    default:
      break;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  xMouse = x;  
  yMouse = y;

  mouseButton = button;
  mouseState = state;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  float xDiff = x - xMouse;
  float yDiff = y - yMouse;
  float speed = 0.001;
  
  if (mouseButton == GLUT_LEFT_BUTTON) 
  {
    eyeCenter[0] -= xDiff * speed;
    eyeCenter[1] += yDiff * speed;
  }
  if (mouseButton == GLUT_RIGHT_BUTTON)
  {
    zoom -= yDiff * speed;
  }

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutPassiveMouseMotion(int x, int y)
{
  // make the lower left the origin
  y = yScreenRes - y;

  float xNorm = (float)x / xScreenRes;
  float yNorm = (float)y / yScreenRes;

  float halfZoom = 0.5 * zoom;
  float xWorldMin = eyeCenter[0] - halfZoom;
  float xWorldMax = eyeCenter[0] + halfZoom;

  // get the bounds of the field in screen coordinates
  //
  // if non-square textures are ever supported, change the 0.0 and 1.0 below
  float xMin = (0.0 - xWorldMin) / (xWorldMax - xWorldMin);
  float xMax = (1.0 - xWorldMin) / (xWorldMax - xWorldMin);

  float yWorldMin = eyeCenter[1] - halfZoom;
  float yWorldMax = eyeCenter[1] + halfZoom;

  float yMin = (0.0 - yWorldMin) / (yWorldMax - yWorldMin);
  float yMax = (1.0 - yWorldMin) / (yWorldMax - yWorldMin);

  int xRes = distanceField.xRes();
  int yRes = distanceField.yRes();

  REAL xScale = 1.0;
  REAL yScale = 1.0;

  if (xRes < yRes)
    xScale = (REAL)yRes / xRes;
  if (xRes > yRes)
    yScale = (REAL)xRes / yRes;

  // index into the field after normalizing according to screen
  // coordinates
  xField = xScale * xRes * ((xNorm - xMin) / (xMax - xMin));
  yField = yScale * yRes * ((yNorm - yMin) / (yMax - yMin));
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE| GLUT_RGBA);
  glutInitWindowSize(xScreenRes, yScreenRes); 
  glutInitWindowPosition(10, 10);
  glutCreateWindow(windowLabel.c_str());

  // set the viewport resolution (w x h)
  glViewport(0, 0, (GLsizei) xScreenRes, (GLsizei) yScreenRes);

  glClearColor(0.1, 0.1, 0.1, 0);
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glShadeModel(GL_SMOOTH);

  updateViewingTexture();

  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if (argc < 3)
  {
    cout << endl;
    cout << " USAGE: " << argv[0] << " <field file> <distance field file> <window header (optional)> <flip (optional)>" << endl;
    cout << "        for the optional flip parameter, specify xy or yz" << endl;
    cout << "        to choose the flip axis" << endl << endl;
    cout << "        e.g. " << argv[0] << " field.field3d distance.field3d MyField xy" << endl;
    cout << endl;
    return 0;
  }

  vectorField.read(argv[1]);
  distanceField.read(argv[2]);

  cout << " Vector field dims: " << vectorField.dims() << endl;
  cout << " Distance field dims: " << distanceField.dims() << endl;

  assert(vectorField.xRes() == distanceField.xRes());
  assert(vectorField.yRes() == distanceField.yRes());
  assert(vectorField.zRes() == distanceField.zRes());

  // see if there's a window name
  if (argc >= 4)
    windowLabel = string(argv[3]);

  // see if there's a flip parameter
  string flip;
  if (argc >= 5)
    flip = string(argv[4]);
 
  if (flip.compare("xy") == 0 || flip.compare("XY") == 0 ||
      flip.compare("yx") == 0 || flip.compare("YX") == 0)
  {
    vectorField = vectorField.flipXY();
    distanceField = distanceField.flipXY();
  }
  if (flip.compare("yz") == 0 || flip.compare("YZ") == 0 ||
      flip.compare("zy") == 0 || flip.compare("ZY") == 0)
  {
    vectorField = vectorField.flipZY();
    distanceField = distanceField.flipZY();
  }

  // create a scaled copy for display
  displayVectorField = vectorField;

  zSlice = vectorField.zRes() / 2;
  //viewingField = vectorField.zSlice(zSlice);
  viewingField = distanceField.zSlice(zSlice);

  centers = VECTOR3_FIELD_3D::cellCenters(distanceField);
  displayCenters = centers;

  // cache the cell centers
  REAL xLength = distanceField.lengths()[0];
  REAL yLength = distanceField.lengths()[1];

  /*
  int xRes = distanceField.xRes();
  int yRes = distanceField.yRes();
  if (xRes < yRes)
    xLength = (REAL)xRes / yRes;
  if (yRes < xRes)
    yLength = (REAL)yRes / xRes;

  cout << " distance center: " << distanceField.center() << endl;
  cout << " distance length: " << distanceField.lengths() << endl;
  cout << " vector field center: " << vectorField.center() << endl;
  cout << " vector field length: " << vectorField.lengths() << endl;

  cout << " (104,112,100) center: " << centers(104,112,100) << endl;

  REAL scale = 1;
  */
  scale = (xLength > 1) ? 1.0 / xLength : 1;
  scale = (yLength > 1) ? 1.0 / yLength : scale;

  cout << " scaling by: " << scale << endl;

  for (int x = 0; x < displayCenters.totalCells(); x++)
  {
    displayCenters[x][0] *= scale;
    displayCenters[x][1] *= scale;
    displayVectorField[x][0] *= scale;
    displayVectorField[x][1] *= scale;
  }

  glutInit(&argc, argv);

  glvuWindow();
  return 1;
}
