//------------------------------------------------------------------------------
// File : glvu.cpp
//------------------------------------------------------------------------------
// GLVU : Copyright 1997 - 2002 
//        The University of North Carolina at Chapel Hill
//------------------------------------------------------------------------------
// Permission to use, copy, modify, distribute and sell this software and its 
// documentation for any purpose is hereby granted without fee, provided that 
// the above copyright notice appear in all copies and that both that copyright 
// notice and this permission notice appear in supporting documentation. 
// Binaries may be compiled with this software without any royalties or 
// restrictions. 
//
// The University of North Carolina at Chapel Hill makes no representations 
// about the suitability of this software for any purpose. It is provided 
// "as is" without express or implied warranty.

//============================================================================
// glvu.cpp : OpenGL/GLUT -based viewer
//============================================================================

#include "glvu.h"
#include "trackball.h"

GLVU *GLVU::GLVUs[MAX_GLVUS];
static int GLVU_Initialized = 0;

/// The constructor
/**
 * Create a new GLVU object with default settings.
 *
 * @note To actually do anything useful with it you will need to call
 * either Init() or InitWin(), and probably SetAllCams() too.
 *
 * @see Init, InitWin 
 */
GLVU::GLVU()
{
  if (!GLVU_Initialized) {
    for (int i=0; i<MAX_GLVUS; i++)
      GLVUs[i]=NULL;
    GLVU_Initialized = 1;
  }

  WindowID = -1;

  WorldNavMode=0;
  InsideLookingOutMode=0;
  NumCams=4;
  Cams = new Camera[NumCams];
  CamDisplayOn = new int[NumCams];
  for (int i=0;i<NumCams;i++) CamDisplayOn[i]=0;
  Cam = &(Cams[0]);
  MainMenuID=0;
  LeftButtonDown=0;
  RightButtonDown=0;
  MiddleButtonDown=0;
  OldX=0; OldY=0; NewX=0; NewY=0;
  moveSpeed = 1;
  CtrlPressed=0;
  AltPressed=0;
  ShiftPressed=0;
  InertiaOn=0, 
  InertiaEnabled=1;
  InertiaDelay=10;
  SetInertiaFunc(NULL);
  calcFPS = 0;
  lastFPS = 0.0;
  lastFPSCount = 0;
  RecordingOn = 0;
  PlaybackOn = 0;
  RecordPlaybackFP = 0;
}

/// The destructor
/**
 * Virtual since GLVU has virtual methods and allocates some members on the heap. 
 */
GLVU::~GLVU()
{
  delete[] Cams;
  delete[] CamDisplayOn;
}

/// Initialize a viewer
/**
 * Create a GLVU window and initialize it, via GLUT calls. This must
 * be done \em before calling \c glutMainLoop(). OpenGL calls can be made
 * after this to set state for lights, textures, materials, etc.  This
 * version creates its own GLUT window with a call to
 * \c glutCreateWindow().  If you already have a glut window handy for
 * whatever reason, you can call InitWin() instead, to wrap a GLVU
 * around that window.
 *
 * Also set up default keyboard, mouse, and reshape callback with GLUT.
 *
 * @param WindowTitle What should appear on the window's title bar
 * (and in the task manager on Windows, or on the icon on an XWindows
 * desktop).
 *
 * @param VisualMode a set of GLUT flags combined with logical or'ing.
 * These are passed to glutIinitDisplayMode before creating the
 * window.  Possible values include \c GLUT_DEPTH, \c, GLUT_DOUBLE, \c
 * GLUT_ALPHA, \c GLUT_ACCUM, and on SGIs the all-important \c
 * GLUT_MULTISAMPLE.  See GLUT documentation for more information.
 *
 * @param WindowStartX Desired X position for the window's upper left corner.
 * @param WindowStartX Desired Y position for the window's upper left corner.
 * @param WindowWidth Desired width of window in pixels
 * @param WindowHeight Desired height of window in pixels
 * 
 * @note \em Every GLVU object should have this method called on it
 *       (this or InitWin()).
 *
 * @see InitWin */
int GLVU::Init(char *WindowTitle, 
               unsigned int VisualMode,
               int WindowStartX, int WindowStartY,
               int WindowWidth, int WindowHeight)
{
  // OPENGL/GLUT STUFF STUFF TO CREATE WINDOW
  glutInitDisplayMode(VisualMode);
  glutInitWindowPosition(WindowStartX,WindowStartY);
  glutInitWindowSize(WindowWidth,WindowHeight);
  WindowID = glutCreateWindow(WindowTitle);  // STARTS AT WindowID=1

  // DO GLVU INITIALIZATION PART
  return InitWin(WindowID);
}

/// Initialize a viewer with an existing GLUT window
/**
 * Set a GLVU window to use the existing GLUT window referred to by \p
 * wID.  This is something you way wish to do if you are integrating
 * with a legacy GLUT application, for instance, or if you need to use
 * some special GLUT functionality to create a window with a
 * particular display mode.
 *
 * Also set up default keyboard, mouse, and reshape callback with GLUT.
 * 
 * @param wID The window identifier obtained from GLUT.
 *
 * @see Init
 */
int GLVU::InitWin(int wID)
{
  // MAKE SURE GLUT LIKES THIS ID
  glutSetWindow(wID);

  // IF WE GET HERE THEN GLUT LIKED THE WINDOW ID
  WindowID = wID;

  // ADD WINDOW TO STATIC GLVU ARRAY; FIRST CHECK IF TOO MANY ARE BEING ADDED
  if (WindowID>=MAX_GLVUS) { 
    printf("ERROR: too many windows added!\n");
    glutDestroyWindow(WindowID);
    WindowID=0;
    return(WindowID); 
  }
  GLVUs[WindowID] = this;

  // REGISTER PER-WINDOW DEFAULT CALLBACKS
  glutKeyboardFunc(DefaultKeyboardFunc);
  glutMouseFunc(DefaultMouseFunc);
  glutMotionFunc(DefaultMotionFunc);
  glutReshapeFunc(DefaultReshapeFunc);
  glutDisplayFunc(DefaultDisplayFunc);
  SetInertiaFunc(DefaultInertiaFunc);

  InitMenu();

#if defined(DEBUG) | defined(_DEBUG)
  PrintVisualInfo();
#endif

  return(WindowID);
}

/// The default GLUT display function implementation
/**
 * This just calls GetGLVU()->Display(), thus allowing Object-Oriented
 * people to customize GLVU's behavior by overriding the Display()
 * method instead of by dealing directly with GLUTs callbacks.
 *
 * @see Display, DefaultReshapeFunc, DefaultInertiaFunc, DefaultKeyboardFunc,
 *     DefaultMouseFunc, DefaultMotionFunc         
 */
void GLVU::DefaultDisplayFunc()
{
  GetGLVU()->Display();
}

/// Handler for redrawing application OpenGL content.
/**
 * This method can be overridden to perform the application's actual
 * OpenGL drawing.  Typically one begins by calling BeginFrame() as
 * the first call, and ends with EndFrame() as the last call.  Those
 * two methods handle camera setup (path playback and recording),
 * buffer swapping for double buffered windows, and frame rate timing
 * calculations.
 *
 * The default implementation does nothing.
 *
 * In the default configuration, this method is called indirectly via
 * GLUTs \c glutDisplayFunc callback.  Users can override this method,
 * or register their own glut display callback with GLUT directly.  
 *
 * @see DefaultDisplayFunc, DefaultReshapeFunc, DefaultInertiaFunc, 
 *     DefaultKeyboardFunc, DefaultMouseFunc, DefaultMotionFunc
 */
void GLVU::Display()
{
}

/// The default GLUT display function implementation
/**
 * This just calls GetGLVU()->Reshape(), thus allowing Object-Oriented
 * people to customize GLVU's behavior by overriding the Reshape()
 * method instead of by dealing directly with GLUTs callbacks.
 *
 * @see Reshape, DefaultDisplayFunc, DefaultInertiaFunc, DefaultKeyboardFunc,
 *     DefaultMouseFunc, DefaultMotionFunc 
*/
void GLVU::DefaultReshapeFunc(int WW, int WH)
{
  GetGLVU()->Reshape(WW, WH);
}

/// Handler for changed window size
/**
 * Typically this method is overridden to handle setting up the
 * perspective matrix and viewport matrix.  It corresponds to the \c
 * \c glutReshapeFunc(), and is indirectly called via that very mechanism.  
 *
 * The default implementation is quite serviceable.  It sets the
 * viewport to be the whole visible window area, and adjusts all the
 * Cameras so that the aspect ratio of the perspective transformation
 * is not all whacked out.  One reason you might have to override this
 * is if you are implementing stereo using multiple viewports, or if
 * you are using multiple viewports for any other reason.
 *
 * Rather than overriding this method, you could also call
 * glutReshapeFunc directly, to register your function with GLUT.  If
 * you do so then this method will no longer get called. To invoke the
 * default functionality from your callback you can call
 * GetCurrent()->Reshape() or GLVU::DefaultReshapeFunc().
 *
 * @see DefaultReshapeFunc, DefaultDisplayFunc, DefaultInertiaFunc, 
 *     DefaultKeyboardFunc, DefaultMouseFunc, DefaultMotionFunc
 */
void GLVU::Reshape(int WW, int WH)
{
  glViewport(0, 0, WW, WH);  
  AllCamsPerspectiveAspectChange((GLfloat)WW/(GLfloat)WH);
}

/// Dump info about the selected visuals to standard out
/*
 * Currently, this is called automatically when a glut window is
 * initialized with Init() or InitWin().  Can't you tell that we're
 * researchers?  We use stdout, \em and we spit a bunch of geek
 * information out to it even if you don't explicitly ask.
 */
void GLVU::PrintVisualInfo()
{
  GLint i;
  GLboolean j;
  printf("GRAPHICS VISUAL INFO (# bits of each):\n");
  glGetIntegerv(GL_RED_BITS,&i);    printf("RGBA: %d ", i);
  glGetIntegerv(GL_GREEN_BITS,&i);  printf("%d ", i);
  glGetIntegerv(GL_BLUE_BITS,&i);   printf("%d ", i);
  glGetIntegerv(GL_ALPHA_BITS,&i);  printf("%d\n", i);
  glGetIntegerv(GL_ACCUM_RED_BITS,&i);    printf("Accum RGBA: %d ", i);
  glGetIntegerv(GL_ACCUM_GREEN_BITS,&i);  printf("%d ", i);
  glGetIntegerv(GL_ACCUM_BLUE_BITS,&i);   printf("%d ", i);
  glGetIntegerv(GL_ACCUM_ALPHA_BITS,&i);  printf("%d\n", i);
  glGetIntegerv(GL_INDEX_BITS,&i);  printf("Color Index: %d\n", i);
  glGetIntegerv(GL_DEPTH_BITS,&i);  printf("Depth: %d\n", i);
  glGetIntegerv(GL_STENCIL_BITS,&i);  printf("Stencil: %d\n", i);
  glGetBooleanv(GL_DOUBLEBUFFER,&j); printf("Double Buffer?  %d\n", j);
  glGetBooleanv(GL_STEREO,&j); printf("Stereo Buffer?  %d\n", j);
  glGetIntegerv(GL_AUX_BUFFERS,&i);  printf("# Aux Buffers: %d\n", i);
}

/// Print information about current OpenGL errors. (Debug only)
void GLVU::CheckForGLError( char *msg )
{
#if defined(DEBUG) | defined(_DEBUG)
 GLenum errCode;
 const GLubyte *errStr;
 if ((errCode = glGetError()) != GL_NO_ERROR) 
 {
    errStr = gluErrorString(errCode);
    fprintf(stderr,"OpenGL ERROR: %s: %s\n", errStr, msg);
 }
#endif
}

//-------------------------------------------------------------------------------
// Storage for the current font needed by the Text() routine. Global state for
// the text module. Access functions are provided. The following fonts are available
// and are defined in glut.h (fontids are on the left):
//  0: GLUT_BITMAP_8_BY_13        : 13 pixels high
//  1: GLUT_BITMAP_9_BY_15        : 15 pixels high
//  2: GLUT_BITMAP_TIMES_ROMAN_10 : 10 pixels high
//  3: GLUT_BITMAP_TIMES_ROMAN_24 : 24 pixels high
//  4: GLUT_BITMAP_HELVETICA_10   : 10 pixels high
//  5: GLUT_BITMAP_HELVETICA_12   : 12 pixels high
//  6: GLUT_BITMAP_HELVETICA_18   : 18 pixels high
//-------------------------------------------------------------------------------
//void *CurrentFont = GLUT_BITMAP_9_BY_15;
void *CurrentFont = GLUT_BITMAP_HELVETICA_12;

void *Fonts[] = { GLUT_BITMAP_8_BY_13, 
                  GLUT_BITMAP_9_BY_15,
                  GLUT_BITMAP_TIMES_ROMAN_10, 
                  GLUT_BITMAP_TIMES_ROMAN_24,
                  GLUT_BITMAP_HELVETICA_10, 
                  GLUT_BITMAP_HELVETICA_12,
                  GLUT_BITMAP_HELVETICA_18 };
int FontHeights[] = { 13, 15, 10, 24, 10, 12, 18 };
int NumFonts = 7;
#include "text.h"

//----------------------------------------------------------------------------
/// Set all cameras to the same viewing parameters
/**
// Most every GLVU app will call this routine soon after Init().
// It's the best way to get the camera perspecitve matrix parameters all set
// up.
//
// This method also saves these camera parameters into a separate
// member, so that these camera settings can be restored with a call
// to AllCamsResetToOrig()
//
// @param ModelMin, ModelMax  axis-aligned bounding box of the world
// @param Eye   starting location of the camera origin (eye)
// @param LookAtCntr   starting world point to look towards (usually model center)
// @param ViewUp   starting world up vector (usually [0,1,0])
// @param Yfov   vertical field-of-view in degrees.
// @param Aspect   pixel width/height. 
// @param Nearfactor, Farfactor  determine the distances to the near and 
//        far planes as a factor times the world bounding sphere radius 
//        (sphere that surrounds the given (WorldMin,WorldMax).
// @note \p Eye cannot equal \p LookAtCntr!! Earth-shattering Kaboom!!
//
// @see AllCamsPerspectiveChange, AllCamsPerspectiveAspectChange,
//     AllCamsResetToOrig
*/
//----------------------------------------------------------------------------
void GLVU::SetAllCams(const glvuVec3f& WorldMin, const glvuVec3f& WorldMax, 
                  const glvuVec3f& Eye, const glvuVec3f& LookAtCntr, const glvuVec3f& viewup,
                  float Yfov, float Aspect, float NearFactor, float FarFactor)
{
  // STORE WORLD BOUNDING SPHERE
  WorldCenter = (WorldMax+WorldMin)*0.5;
  glvuVec3f Diagonal = WorldMax-WorldCenter;
  WorldRadius = Diagonal.Length();

  // SET CAMERA PROJECTION
  float NearDist = WorldRadius * NearFactor;
  float FarDist = WorldRadius * FarFactor;
  AllCamsPerspectiveChange(Yfov,Aspect,NearDist,FarDist);

  // STORE THE UP VECTOR
  ViewUp = viewup;

  // SET CAMERA MODELVIEW
  for (int i=0; i<NumCams; i++)
    Cams[i].LookAt(Eye,LookAtCntr,ViewUp);
    
  SetOrigCam(&Cams[0]);
}


//----------------------------------------------------------------------------
/// Change some parameters of all camearas
/**
// 
// This is an easy way to change some key parameters that define the
// perspective transformation being used by every camera, without
// changing the location or orientation of the cameras.
//
// @param Yfov   vertical field-of-view in degrees.
// @param Aspect window width/height in pixels. 
// @param NDist, FDist  Near and far plane distances.
//
// @see SetAllCams, AllCamsPerspectiveAspectChange, AllCamsResetToOrig
*/
//----------------------------------------------------------------------------
void GLVU::AllCamsPerspectiveChange(float Yfov, float Aspect, float Ndist, float Fdist)
{
  for (int i=0; i<NumCams; i++)
    Cams[i].Perspective(Yfov, Aspect, Ndist, Fdist);
}


//----------------------------------------------------------------------------
/// Change the aspect ratio of all cameras
/**
// Change the aspect ratio of the perspective transformation for every
// camera without changing any of the other parameters of the perspective 
// matrix.
//
// This is called by the default implementation of Reshape().
//
// @param NewAspect A new aspect ratio (width over height)
//
// @see SetAllCams, AllCamsPerspectiveChange, AllCamsResetToOrig, Reshape
*/
//----------------------------------------------------------------------------
void GLVU::AllCamsPerspectiveAspectChange(float NewAspect)
{
  // PRESERVE OTHER PROJECTION SETTINGS, ONLY CHANGE ASPECT RATIO FOR
  // NEW WINDOW SIZES
  float Yfov, OldAspect, Ndist, Fdist;
  for (int i=0; i<NumCams; i++)
  {
    Cams[i].GetPerspectiveParams(&Yfov,&OldAspect,&Ndist,&Fdist);
    Cams[i].Perspective(Yfov,NewAspect,Ndist,Fdist);
  }
}

//----------------------------------------------------------------------------
/// Change the near and far distances of all cameras
/**
// Change the near and far distances of the perspective frustum for every
// camera without changing any of the other parameters of the perspective 
// matrix.
//
// @param Ndist,Fdist The new near and far distances
//
// @see SetAllCams, AllCamsPerspectiveChange, AllCamsResetToOrig
*/
//----------------------------------------------------------------------------
void GLVU::AllCamsPerspectiveNearFarChange(float Ndist, float Fdist)
{
  // PRESERVE OTHER PROJECTION SETTINGS, ONLY CHANGE NEAR/FAR
  float Yfov, Aspect, OldNdist, OldFdist;
  for (int i=0; i<NumCams; i++)
  {
    Cams[i].GetPerspectiveParams(&Yfov,&Aspect,&OldNdist,&OldFdist);
    Cams[i].Perspective(Yfov,Aspect,Ndist,Fdist);
  }
}

//----------------------------------------------------------------------------
/// Set all the cameras back to the initial view
/**
// Set all the cameras back to the initial camera that was set with the last
// call to SetAllCams().
//
// @see SetAllCams, AllCamsPerspectiveChange
*/
//----------------------------------------------------------------------------
void GLVU::AllCamsResetToOrig()
{
  for (int i=0; i<NumCams; i++)
    Cams[i].Copy(OrigCam);
  float WW = glutGet((GLenum)GLUT_WINDOW_WIDTH);
  float WH = glutGet((GLenum)GLUT_WINDOW_HEIGHT);
  AllCamsPerspectiveAspectChange(WW/WH);
}

//----------------------------------------------------------------------------
/// Return the next modelview matrix
/**
// This gets the next modelview matrix that should be used for displaying.
// BeginFrame() calls this automatically, but if you aren't calling
// BeginFrame() from your Display() method, then you can call this
// directly to get the next modelview matrix to pass to OpenGL.  The
// matrix is in OpenGL format, so you can pass it right to \c glLoadMatrix().  
// Here is what BeginFrame() does with the matrix stack for example:
//
//  @code
    GLfloat M[16];
    glMatrixMode(GL_PROJECTION);
     glLoadMatrixf( GetProjectionMatrix(M) );
    glMatrixMode(GL_MODELVIEW);
     glLoadMatrixf( GetModelviewMatrix(M) );
    @endcode
//
// This innocuous-sounding getter function also plays a role in camera
// path playback.  If playback mode is engaged, then the matrix
// returned to you will not be based on the user's mouse wigglings,
// but on recorded camera data loaded in from a file on demand.  So
// this isn't your typical do-nothing getter function, it actually
// does some work to figure out what matrix to give you.
//
// @see GetProjectionMatrix, BeginFrame
*/
//----------------------------------------------------------------------------
float* GLVU::GetModelviewMatrix(float* M)
{
  // WRITE ALL CAMERAS TO FILE IF RECORDING ON
  if (RecordingOn) 
  {
    for (int i=0; i<NumCams; i++)
      Cams[i].WriteToFile(RecordPlaybackFP);
  }

  // READ ALL CAMERAS FROM FILE IF PLAYBACK ON
  else if (PlaybackOn)
  {
    int ShouldEndPlayback=0;
    for (int i=0; i<NumCams; i++)
      if ( !(Cams[i].ReadFromFile(RecordPlaybackFP)) )
        ShouldEndPlayback=1;
    if (ShouldEndPlayback)
      EndPlayback();
  }

  return( Cam->GetModelviewMatrix(M) );
}

/// Return the OpenGL-style projection matrix from the current camera
/**
 * The matrix is an array of 16 floats, suitable for passing to
 * \c glLoadMatrixf()
 * 
 * @see GetModelviewMatrix, BeginFrame
 */
float* GLVU::GetProjectionMatrix(float* M)
{
  return( Cam->GetProjectionMatrix(M) );
}

/// Prepare the current frame for drawing in display routine
/**
 * Sets up the modelview and projection matrices based on the current 
 * camera, and updates the frame rate timer.
 * 
 * Unless you need to do some fancy customized manipulation of the
 * OpenGL matrix stack, just calling BeginFrame() is a convenient and effective way 
 * to set things up so you can just do your application-specific OpenGL drawing. 
 *
 * @see EndFrame, GetModelviewMatrix, GetProjectionMatrix, Display
 */
void GLVU::BeginFrame()
{
  GLfloat M[16];
  glMatrixMode(GL_PROJECTION);
   glLoadMatrixf( GetProjectionMatrix(M) );
  glMatrixMode(GL_MODELVIEW);
   glLoadMatrixf( GetModelviewMatrix(M) );
  UpdateFPS();
}

/// Finish up the current OpenGL display frame.  
/**
 * This should be, or at least can be, called at the end of your Display() 
 * method.  
 *
 * Unless you need to swap buffers in the middle of a frame for some
 * reason, this is usually sufficient.  
 *
 * This method is also responsible for drawing the lines in space to
 * represent other cameras, so that you can see from one camera where
 * the others are located in space.  The default keys for switching
 * cameras are 1, 2, 3, and 4.  The default keys for toggling display
 * of camera frusta are Shift-1, Shift-2, Shift-3, and Shift-4.
 * 
 * @see BeginFrame, GetModelviewMatrix, GetProjectionMatrix, Display 
 */
void GLVU::EndFrame()
{
  // DRAW CAMERAS (IF DISPLAY IS ON)
  for (int i=0; i<NumCams; i++)
    if (GetCamDisplayOn(i)) 
      GetCam(i)->Display();

  glutSwapBuffers();
}

//----------------------------------------------------------------------------
/// Return the world space ray corresponding to a particular screen pixel
/** 
// Given a screen pixel location (sx,sy) w/ (0,0) at the lower-left, and the
// screen dimensions, return the ray (start,dir) of the ray in world coords.
//  
// A very handy routine for implementing ray-tracers.
//
// @see  Camera::GetPixelRay()
*/
//----------------------------------------------------------------------------
void GLVU::GetPixelRay(int sx, int sy, int ww, int wh, 
                       glvuVec3f *Start, glvuVec3f *Dir) const
{
  Cam->GetPixelRay((float)sx,(float)sy,ww,wh,Start,Dir);
}


//----------------------------------------------------------------------------
// MOUSE-TO-CAMERA MOVEMENT ROUTINE: TRANSLATE
//----------------------------------------------------------------------------
/// One of the core routines used to translate mouse motion into camera motion
//----------------------------------------------------------------------------
void GLVU::TranslateX(int OldX, int NewX, int WW)
{
  float RelX = (NewX-OldX)/(float)WW;
  float dX;

  if (InsideLookingOutMode)
  {
    dX = RelX * (-WorldRadius*0.25f);
  }
  else // FOR OUTSIDE-LOOKING-IN
  {
    glvuVec3f EyeToWorldCntr = WorldCenter - Cam->Orig;
    float DistToCntr = EyeToWorldCntr.Length();
    if (DistToCntr<WorldRadius) DistToCntr=WorldRadius;

    float vpX = RelX * (Cam->wR-Cam->wL);
    dX = DistToCntr * vpX / -Cam->Near;
  }

  float M[16];
  glvuVec3f Trans = Cam->X * dX * moveSpeed;
  Translate16fv(M,Trans.x,Trans.y,Trans.z);
  Cam->Xform(M);
}

//----------------------------------------------------------------------------
/// One of the core routines used to translate mouse motion into camera motion
//----------------------------------------------------------------------------
void GLVU::TranslateY(int OldY, int NewY, int WH)
{
  float RelY = (NewY-OldY)/(float)WH;
  float dY;

  if (InsideLookingOutMode)
  {
    dY = RelY * (WorldRadius*0.25f);
  }
  else // FOR OUTSIDE-LOOKING-IN
  {
    glvuVec3f EyeToWorldCntr = WorldCenter - Cam->Orig;
    float DistToCntr = EyeToWorldCntr.Length();
    if (DistToCntr<WorldRadius) DistToCntr=WorldRadius;

    float vpY = RelY * (Cam->wT-Cam->wB);
    dY = DistToCntr * vpY / Cam->Near;
  }

  float M[16];
  glvuVec3f Trans = Cam->Y * dY * moveSpeed;
  Translate16fv(M,Trans.x,Trans.y,Trans.z);
  Cam->Xform(M);
}

//----------------------------------------------------------------------------
// MOUSE-TO-CAMERA MOVEMENT ROUTINE: DRIVE MODE
//----------------------------------------------------------------------------
/// One of the core routines used to translate mouse motion into camera motion
//----------------------------------------------------------------------------
void GLVU::DriveY(int OldY, int NewY, int WH)
{
  float RelY = (NewY-OldY)/(float)WH;
 
  glvuVec3f EyeToWorldCntr = WorldCenter - Cam->Orig;
  float DistToCntr = EyeToWorldCntr.Length();
  if (DistToCntr<WorldRadius) DistToCntr=WorldRadius;

  float M[16];
  glvuVec3f Trans = Cam->Z * (DistToCntr*RelY*2.0f) * moveSpeed;
  Translate16fv(M,Trans.x,Trans.y,Trans.z);
  Cam->Xform(M);
}


//----------------------------------------------------------------------------
// MOUSE-TO-CAMERA MOVEMENT ROUTINE: LOOK MODE
//----------------------------------------------------------------------------
/// One of the core routines used to translate mouse motion into camera motion
//----------------------------------------------------------------------------
void GLVU::LookX(int OldX, int NewX, int WW)
{
  float RelX = (NewX-OldX)/(float)WW;
 
  float M[16], N[16];
  Translate16fv(M,Cam->Orig.x,Cam->Orig.y,Cam->Orig.z);
  Rotate16fv(N,-RelX*90,&(ViewUp.x));
  Mult16fv(M,M,N);
  Translate16fv(N,-Cam->Orig.x,-Cam->Orig.y,-Cam->Orig.z);
  Mult16fv(M,M,N);
  Cam->Xform(M);
}

//----------------------------------------------------------------------------
/// One of the core routines used to translate mouse motion into camera motion
//----------------------------------------------------------------------------
void GLVU::LookY(int OldY, int NewY, int WH)
{
  float RelY = (NewY-OldY)/(float)WH;
 
  float M[16], N[16];
  Translate16fv(M,Cam->Orig.x,Cam->Orig.y,Cam->Orig.z);
  Rotate16fv(N,-RelY*90,&(Cam->X.x));
  Mult16fv(M,M,N);
  Translate16fv(N,-Cam->Orig.x,-Cam->Orig.y,-Cam->Orig.z);
  Mult16fv(M,M,N);
  Cam->Xform(M);
}


//----------------------------------------------------------------------------
// MOUSE-TO-CAMERA MOVEMENT ROUTINE: TRACKBALL MODE
//----------------------------------------------------------------------------
/// One of the core routines used to translate mouse motion into camera motion
//----------------------------------------------------------------------------
void GLVU::TrackBallX(int OldX, int NewX, int WW)
{
  float RelX = (NewX-OldX)/(float)WW;
  glvuVec3f Up = InsideLookingOutMode ? ViewUp : Cam->Y;
  float M[16], N[16];
  Translate16fv(M,WorldCenter.x,WorldCenter.y,WorldCenter.z);
  Rotate16fv(N,-RelX*180,&(Up.x));
  Mult16fv(M,M,N);
  Translate16fv(N,-WorldCenter.x,-WorldCenter.y,-WorldCenter.z);
  Mult16fv(M,M,N);
  Cam->Xform(M);
}

//----------------------------------------------------------------------------
/// One of the core routines used to translate mouse motion into camera motion
//----------------------------------------------------------------------------
void GLVU::TrackBallY(int OldY, int NewY, int WH)
{
  float RelY = (NewY-OldY)/(float)WH;
  float M[16], N[16];
  Translate16fv(M,WorldCenter.x,WorldCenter.y,WorldCenter.z);
  Rotate16fv(N,-RelY*180,&(Cam->X.x));
  Mult16fv(M,M,N);
  Translate16fv(N,-WorldCenter.x,-WorldCenter.y,-WorldCenter.z);
  Mult16fv(M,M,N);
  Cam->Xform(M);
}


//----------------------------------------------------------------------------
// MOUSE-TO-CAMERA MOVEMENT ROUTINE: SGI "HYPERBOLIC SHEET" TRACKBALL MODE
//----------------------------------------------------------------------------
/// One of the core routines used to translate mouse motion into camera motion
//----------------------------------------------------------------------------
void GLVU::HyperBall(int OldX, int OldY, int NewX, int NewY, int WW, int WH)
{
  // SCALE SCREEN COORDS TO [-1.0,1.0 range] for trackball() call
  float RelX0 = 2.0*(OldX/(float)WW)-1.0;
  float RelX1 = 2.0*(NewX/(float)WW)-1.0;
  float RelY0 = 1.0-2.0*(OldY/(float)WH);
  float RelY1 = 1.0-2.0*(NewY/(float)WH);
  float q[4];
  trackball(q, RelX0, RelY0, RelX1, RelY1);
  float Q[4][4], M[16], N[16];
  
  Translate16fv(M,WorldCenter.x,WorldCenter.y,WorldCenter.z);
  invViewing16fv(N,&(Cam->X.x),&(Cam->Y.x),&(Cam->Z.x),&(glvuVec3f::ZERO.x));
  Mult16fv(M, M, N);
  build_rotmatrix(Q, q);
  Transpose16fv((float*)Q);
  Mult16fv(M, M, (float*)Q);

  Viewing16fv(N,&(Cam->X.x),&(Cam->Y.x),&(Cam->Z.x),&(glvuVec3f::ZERO.x));
  glvuVec3f tmp(Cam->Orig), nV;
  tmp -= WorldCenter;
  Mult16fv3fv(&(nV.x), N, &(tmp.x));

  Cam->LoadIdentityXform();
  Cam->Orig.Set((const glvuVec3f)nV);
  Cam->Xform(M);
}

//----------------------------------------------------------------------------
// CAMERA RECORDING AND PLAYBACK
//----------------------------------------------------------------------------
/**
 * Begins path recording.  Records the position and orientation of all four 
 * Cameras to the file specified with SetPathFilename(), or to "path0.dat" if 
 * none is specified.
 * @param FileName The name of the file to open and write to.
 * @see EndRecording, StartPlayback, EndPlayback, SelectCam
 */
void GLVU::StartRecording(const char *FileName)
{
  if (RecordingOn) { printf("Already recording!\n"); return; }
  if (PlaybackOn) EndPlayback();
  RecordPlaybackFP = fopen(FileName,"w");
  if (!RecordPlaybackFP) { printf("ERROR (StartRecording) fp==NULL\n"); return; }
  RecordingOn=1;
}

/**
 * Stops path recording and closes the file being recorded to.
 * @see StartRecording, StartPlayback, EndPlayback, SelectCam
 */
void GLVU::EndRecording()
{
  if (!RecordingOn) { printf("Not recording!\n"); return; }
  if (PlaybackOn) {   // NOT SUPPOSED TO HAPPEN EVER
    printf("Playback is on! Not recording!\n"); return; 
  }
  fclose(RecordPlaybackFP);
  RecordingOn=0;
}

/**
 * Begins path playback from the path file specified, or from "path0.dat"
 * if none is specified.
 * Updates the positions and orientations of all four Cameras.  This relies
 * on the use of GetModelviewMatrix() in the display routine.
 * @see EndPlayback, StartRecording, EndRecording
 */
void GLVU::StartPlayback(const char *FileName)
{
  if (RecordingOn) EndRecording();
  if (PlaybackOn) EndPlayback();
  RecordPlaybackFP = fopen(FileName,"r");
  if (!RecordPlaybackFP) { printf("ERROR (StartPlayback) fp==NULL\n"); return; }
  PlaybackOn=1;

  printf("CAM PLAYBACK...\n");
  PlaybackTime = clock();

  PathPlaybackTimerFunc(WindowID);
}

/**
 * Ends camera path playback, and closes the file that was being read from.
 * @see StartPlayback, StartRecording, EndRecording
 */
void GLVU::EndPlayback()
{
  if (!PlaybackOn) { printf("Not playing back!\n"); return; }
  if (RecordingOn) { printf("Recording is on! Not playing back!\n"); return; } // NOT SUPPOSED TO HAPPEN
  fclose(RecordPlaybackFP);
  PlaybackOn=0;

  printf("%.2f seconds\n", (clock()-PlaybackTime)/(float)CLOCKS_PER_SEC);
}

/**
 * Stops both camera path recording and playback.
 * Equivalent to calling both EndPlayback() and EndRecording().  
 * Just a convenience method.
 * @see StartPlayback, StartRecording, EndPlayback, EndRecording
 */
void GLVU::StopRecordingPlayback()
{
  if (PlaybackOn) EndPlayback();
  else if (RecordingOn) EndRecording();
  else printf("Not recording or playing back!\n");
}

void GLVU::PathPlaybackTimerFunc(int winID)
{
  GLVU *g = GLVUs[winID];
  if (g->PlaybackOn) {
    if (glutGetWindow()!=winID) glutSetWindow(winID);
    glutPostRedisplay();
    glutTimerFunc(g->InertiaDelay, PathPlaybackTimerFunc, g->WindowID);
  }
}


//----------------------------------------------------------------------------
// FRAME TIMING
//----------------------------------------------------------------------------
inline float elapsed_ftime(timeb *tstart, timeb *tend)
{
  return (float)(tend->time - tstart->time)
    + ((float)(tend->millitm - tstart->millitm))/1000.0f;
}

/**
 * If FPS calculation is enabled (StartFPSClock()), then this function updates
 * the current Frames Per Second calculation.  Call only if you are not using
 * the default implementation of BeginFrame(), because it calls this for you.
 * @see StartFPSClock, StopFPSClock, BeginFrame, DrawFPS
 */
void GLVU::UpdateFPS()
{
  if (calcFPS) {
    lastFPSCount++;
    struct timeb newClock;
    ftime(&newClock);
    float tdiff = elapsed_ftime(&lastFPSClock,&newClock);
    if (tdiff >= FPS_INTEGRATE_INTERVAL) {
      lastFPS = (float)(lastFPSCount)/tdiff;
      lastFPSClock.time = newClock.time;
      lastFPSClock.millitm = newClock.millitm;
      lastFPSCount = 0;
    }
  }
}

/**
 * Draw the current frame rate in frames per second at position \p x, \p y.
 *
 * This uses the current OpenGL drawing color, whatever that is.   Probably
 * best to set it explicitly yourself before calling this method.
 *
 * Frame rate calculation requires calling StartFPSClock() once, and either 
 * BeginFrame()/EndFrame() ever frame or manual calls to UpdateFPS().
 *
 * @param x,y The position for the lower left corner of the text, relative to the 
 *            lower left corner of the screen.
 *
 * @see StartFPSClock, StopFPSClock, BeginFrame, UpdateFPS
 */
void GLVU::DrawFPS(int x, int y)
{
  // Draw the current fps number in the lower left corner of the screen.
  //glColor3f(1.0, 1.0, 1.0);
  glColor3f(0,0,0);
  Text(x, y, "FPS: %.1f", GetFPS());
}
void GLVU::Print(int x, int y, char* buffer)
{
  Text(x, y, buffer);
}

/// The default GLUT keyboard function implementation
/**
 * This just calls GetGLVU()->Keyboard(), thus allowing Object-Oriented
 * people to customize GLVU's behavior by overriding the Keyboard()
 * method instead of by dealing directly with GLUT's callbacks.
 *
 * @see Keyboard, DefaultDisplayFunc, DefaultInertiaFunc, DefaultReshapeFunc,
 *     DefaultMouseFunc, DefaultMotionFunc 
 */
void GLVU::DefaultKeyboardFunc(unsigned char Key, int x, int y)
{
  GetGLVU()->Keyboard(Key, x, y);
}

/// Handler for keyboard events
/**
 * The default implementation handles all the GLVU default key bindings.
 * Override this method to add your own key bindings, but if you don't handle 
 * the key, be sure to call the superclass (i.e. call GLVU::Keyboard()) 
 *
 * This only handles "normal" key events, i.e. those that correspond
 * to a key with an ascii character. GLUT also has \c glutSpecialFunc
 * that can be called to set a handler for other key events like the
 * F-keys, or arrow keys.  There are also "Up" versions of the key
 * events in GLUT that can be used to detect key release as well as 
 * key press.
 *
 * Users not interested in creating an object-oriented app can simply
 * call GLUT's \c glutKeyboardFunc to set a callback directly.  If you
 * do so, however, you should call GetCurrent()->Keyboard() or 
 * GLVU::DefaultKeyboardFunc(), in your handler when you don't use the 
 * key.
 *
 * @see DefaultKeyboardFunc, DefaultReshapeFunc, DefaultDisplayFunc, 
 *     DefaultInertiaFunc, DefaultMouseFunc, DefaultMotionFunc 
 */
void GLVU::Keyboard(unsigned char Key, int x, int y)
{
  int iModifiers = glutGetModifiers();
  CtrlPressed    = (iModifiers & GLUT_ACTIVE_CTRL) != 0;   // IF CTRL BIT SET
  AltPressed     = (iModifiers & GLUT_ACTIVE_ALT) != 0;    // IF ALT BIT SET
  ShiftPressed   = (iModifiers & GLUT_ACTIVE_SHIFT) != 0;  // IF SHIFT BIT SET

  switch(Key)
  {
    case 27: EscapeHandler(-1); break;

    case 'z': MainMenuHandler(0); break;
    case 'h': MainMenuHandler(1); break;
    case 'x': MainMenuHandler(2); break;
    case 'c': MainMenuHandler(3); break;
    case 'v': MainMenuHandler(4); break;

    case '=': vuOptionsMenuHandler(1); break;
    case 'o': vuOptionsMenuHandler(2); break;
    case 'r': vuOptionsMenuHandler(3); break;
    case 'i': vuOptionsMenuHandler(4); break;

    case 'w': glOptionsMenuHandler(2); break;
    case 'b': glOptionsMenuHandler(3); break;
    case 'l': glOptionsMenuHandler(8); break;
    case 'n': glOptionsMenuHandler(9); break;
    case 'm': glOptionsMenuHandler(10); break;

    case '1': CamViewMenuHandler(0); break;
    case '2': CamViewMenuHandler(1); break;
    case '3': CamViewMenuHandler(2); break;
    case '4': CamViewMenuHandler(3); break;

    case '!': CamViewDisplayMenuHandler(0); break;
    case '@': CamViewDisplayMenuHandler(1); break;
    case '#': CamViewDisplayMenuHandler(2); break;
    case '$': CamViewDisplayMenuHandler(3); break;
  };
}

#include "snapshot.h"  // JUST FOR SCREEN SNAPSHOT MENU OPTION

void GLVU::MainMenuHandler(int value)
{
  switch(value)
  {
  case NAV_MODE_TRACKBALL:
  case NAV_MODE_HYPERBALL:
  case NAV_MODE_DRIVE:
  case NAV_MODE_TRANSLATE:
  case NAV_MODE_LOOK:

/*
  WVB: This messes up anyone who sets their own mouse function
  (pretty common).  Usually they will call DefaultMouseFunc
  themselves in the override, if they want this type of NAV_MODE
  behavior.  If they DON'T set their own mouse func then this isn't
  necessary at all.

  glutMouseFunc(DefaultMouseFunc);   // SET DEFAULT MOUSE FUNCS
  glutMotionFunc(DefaultMotionFunc);
*/
      GLVU *g = GetGLVU();
      g->WorldNavMode=value;
      g->SetInertiaOn(0);
      break;
  };
}

void GLVU::vuOptionsMenuHandler(int value)
{
  Camera *Cam;
  glvuVec3f Eye, ViewRefPt, ViewUp;
  float Yfov, Aspect, Ndist, Fdist;

  switch(value)
  {
    case 1: 
      SnapShot(); 
      break;
    case 2:  // TOGGLE BETWEEN INSIDE-LOOKING-OUT AND OUTSIDE-LOOKING-IN NAV MODES
      GetGLVU()->ToggleInOutMode();
      break;
    case 3:  // RESTORE ALL CAMS TO ORIGINAL (STARTING) CAM POS
      GetGLVU()->AllCamsResetToOrig();
      glutPostRedisplay();
      break;
    case 4: {
      GLVU *g = GetGLVU();
      if (g->GetInertiaEnabled()) { g->SetInertiaEnabled(0); g->SetInertiaOn(0); }
      else g->SetInertiaEnabled(1); }
      break;
    case 5: 
      Cam = GetGLVU()->GetCurrentCam();
      Cam->GetLookAtParams(&Eye,&ViewRefPt,&ViewUp);
      Cam->GetPerspectiveParams(&Yfov,&Aspect,&Ndist,&Fdist);
      printf("--- CURRENT CAM PARAMS ---\n");
      printf("       Eye: "); Eye.Print();
      printf("LookAtCntr: "); ViewRefPt.Print();
      printf("    ViewUp: "); ViewUp.Print();
      printf("     Y FOV: %f\n", Yfov);  
      printf("    Aspect: %f\n", Aspect);
      printf("      Near: %f\n", Ndist);
      printf("       Far: %f\n", Fdist);
      break;
  };
}

void GLVU::glOptionsMenuHandler(int value)
{
  int State;
  int PolygonState[2];

  switch(value)
  {
    case 2: 
      glGetIntegerv(GL_POLYGON_MODE,PolygonState);
      if (PolygonState[0]==GL_POINT) glPolygonMode(GL_FRONT_AND_BACK,GL_FILL); 
      else if (PolygonState[0]==GL_FILL) {
        glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE,&State);
        glPolygonMode(State?GL_FRONT_AND_BACK:GL_FRONT,GL_LINE);
      }
      else glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      break;
    case 3:
      if (glIsEnabled(GL_CULL_FACE)) glDisable(GL_CULL_FACE); else glEnable(GL_CULL_FACE);
      break;
    case 6:
      glGetIntegerv(GL_SHADE_MODEL,&State);
      if (State==GL_SMOOTH) glShadeModel(GL_FLAT); else glShadeModel(GL_SMOOTH);
      break;
    case 8:
      if (glIsEnabled(GL_LIGHTING)) glDisable(GL_LIGHTING); else glEnable(GL_LIGHTING);
      break;
    case 7:
      glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE,&State);
      if (State == GL_TRUE) 
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
      else
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
      break;
    case 1:
      glGetIntegerv(GL_LIGHT_MODEL_LOCAL_VIEWER,&State);
      if (State == GL_TRUE) 
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
      else
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
      break;
    case 9:
      glGetIntegerv(GL_CULL_FACE_MODE,&State);
      if (State==GL_BACK) glCullFace(GL_FRONT); else glCullFace(GL_BACK);
      break;
    case 10:
      if (glIsEnabled(GL_COLOR_MATERIAL)) glDisable(GL_COLOR_MATERIAL); else glEnable(GL_COLOR_MATERIAL);
      break;
  };
  glutPostRedisplay();
}

void GLVU::PathPlaybackMenuHandler(int value)
{
  switch(value)
  {
    case 1: 
      GetGLVU()->StartRecording("cam_record0.dat"); 
      break;
    case 2: 
      GetGLVU()->StartPlayback("cam_record0.dat"); 
      break;
    case 3: 
      GetGLVU()->StopRecordingPlayback(); 
      break;
  };
}

void GLVU::CamViewMenuHandler(int value)
{
  GetGLVU()->SelectCam(value);
  glutPostRedisplay();
}

void GLVU::CamViewDisplayMenuHandler(int value)
{
  GLVU *g = GetGLVU();
  if (g->CamDisplayOn[value]) g->CamDisplayOn[value]=0; 
  else g->CamDisplayOn[value]=1;
  glutPostRedisplay();
}

void GLVU::EscapeHandler(int value)
{
  exit(1);
}

#define GLVU_PATH_PLAY_RECORD    1
#define GLVU_PATH_PLAY_PLAY      2
#define GLVU_PATH_PLAY_STOP      3

#define GLVU_CAMVIEW_0           0
#define GLVU_CAMVIEW_1           1
#define GLVU_CAMVIEW_2           2
#define GLVU_CAMVIEW_3           3 

#define GLVU_VU_OPTIONS_SNAPSHOT    1
#define GLVU_VU_OPTIONS_INOUT_MODE  2
#define GLVU_VU_OPTIONS_RESET_VIEW  3
#define GLVU_VU_OPTIONS_INERTIA     4
#define GLVU_VU_OPTIONS_DUMP_CAM    5

#define GLVU_GL_OPTIONS_POLYMODE    2


/**
 * Sets up all of the default GLVU menu choices.  These are bound to the 
 * GLUT right-click menu.  You can subclass GLVU to override this method
 * if you need to implement a custom menu.
 * You can also get the GLUT ID for the menu by calling GetMainMenuID().
 *
 * This method is  called automatically from Init() and InitWin().
 */
void GLVU::InitMenu()
{
  int i;
  int PathPlayMenu = glutCreateMenu(PathPlaybackMenuHandler);
   glutAddMenuEntry("Record",1);
   glutAddMenuEntry("Play",  2);
   glutAddMenuEntry("Stop",  3);

  char Buffer[10];
  int CamViewMenu = glutCreateMenu(CamViewMenuHandler);
   for (i=0; i<NumCams; i++)
     { sprintf(Buffer,"Cam (%d)",i);
       glutAddMenuEntry(Buffer,i); }

  int CamViewDisplayMenu = glutCreateMenu(CamViewDisplayMenuHandler);
   glutAddMenuEntry("Cam 0 (!)",0);
   glutAddMenuEntry("Cam 1 (@)",1);
   glutAddMenuEntry("Cam 2 (#)",2);
   glutAddMenuEntry("Cam 3 ($)",3);

  int vuOptionsMenu = glutCreateMenu(vuOptionsMenuHandler);
   glutAddMenuEntry("SnapShot (=)",1);
   glutAddMenuEntry("InOut/OutIn (o)",2);
   glutAddMenuEntry("Reset Orig Views (r)",3);
   glutAddMenuEntry("Inertia (i)",4);
   glutAddMenuEntry("Dump Cam Params",5);

  int glOptionsMenu = glutCreateMenu(glOptionsMenuHandler);
   glutAddMenuEntry("Solid/Wire/Pt (w)",2);
   glutAddMenuEntry("Smooth/Flat Shading", 6);
   glutAddMenuEntry("Lighting (l)",8);
   glutAddMenuEntry("One/Two Sided Lighting",7);
   glutAddMenuEntry("Inf/Local Viewer",1);
   glutAddMenuEntry("Materials (m)",10);
   glutAddMenuEntry("Face Culling (b)",     3);
   glutAddMenuEntry("Back/Front Cull Face (n)",9);

  MainMenuID = glutCreateMenu(MainMenuHandler);
   glutAddMenuEntry("TrackBall (z)",    NAV_MODE_TRACKBALL);
   glutAddMenuEntry("HyperBall (h)",    NAV_MODE_HYPERBALL);
   glutAddMenuEntry("Drive (x)",        NAV_MODE_DRIVE);
   glutAddMenuEntry("Translate (c)",    NAV_MODE_TRANSLATE);
   glutAddMenuEntry("Look (v)",         NAV_MODE_LOOK);

   glutAddSubMenu("Viewer Options",vuOptionsMenu);
   glutAddSubMenu("OpenGL Options",glOptionsMenu);
   glutAddSubMenu("PathPlay",PathPlayMenu);
   glutAddSubMenu("Cam View",CamViewMenu);
   glutAddSubMenu("Cam Display",CamViewDisplayMenu);

  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

// SUPPORT FOR INERTIAL SYSTEM

#ifndef ABS
template <class T>
static bool ABS(T x)
{
  return (((x)<0)?(-(x)):x);  
}
#endif
#define MINMOVE 2  // AMOUNT (in pixels) TO MOVE TO INVOKE INERTIAL SYSTEM

/// The default GLUT mouse button callback implementation
/**
 * This just calls GetGLVU()->Mouse(), thus allowing Object-Oriented
 * people to customize GLVU's behavior by overriding the Mouse()
 * method instead of by dealing directly with GLUTs callbacks.
 *
 * @see Mouse, DefaultDisplayFunc, DefaultInertiaFunc, DefaultKeyboardFunc,
 *     DefaultMotionFunc, DefaultReshapeFunc 
 */
void GLVU::DefaultMouseFunc(int button, int state, int x, int y)
{
  GetGLVU()->Mouse(button, state, x, y);
}

/// Handler for mouse clicks (called when a button is pressed or released)
/**
 * This method can be overridden to perform application specific
 * funcionality (e.g. picking).  The default implementation does some
 * important work for handling camera manipulation (world
 * navigation), so if you override this method you should always call
 * GLVU::Mouse() from your override.
 * 
 * Users not interested in creating an object-oriented app can simply
 * call GLUT's \c glutMouseFunc to set a callback directly.  If you
 * do so you can still call GetCurrent()->Mouse() or 
 * GLVU::DefaultMouseFunc(), in your handler to get the default 
 * behavior back.
 *
 *  @param button One of \c GLUT_LEFT_BUTTON, \c GLUT_MIDDLE_BUTTON, 
 *         or \c GLUT_RIGHT_BUTTON
 *  @param state One of \c GLUT_UP or \c GLUT_DOWN
 *  @param x,y The pointer location when the event took place 
 *         in window coordinates.
 * 
 *  @see DefaultMouseFunc, Motion, SetWorldNavMode
 */
void GLVU::Mouse(int button, int state, int x, int y)
{
  int iModifiers = glutGetModifiers();
  CtrlPressed    = (iModifiers & GLUT_ACTIVE_CTRL) != 0;   // IF CTRL BIT SET
  AltPressed     = (iModifiers & GLUT_ACTIVE_ALT) != 0;    // IF ALT BIT SET
  ShiftPressed   = (iModifiers & GLUT_ACTIVE_SHIFT) != 0;  // IF SHIFT BIT SET

  // SET APPROPRIATE FLAGS FOR A LEFT-BUTTON MOUSE EVENT
  if (button==GLUT_LEFT_BUTTON)
  {
    // STORE THE NEW MOUSE POS
    NewX=x; NewY=y;

    if (state==GLUT_DOWN)  // LEFT-BUTTON DOWN
    {
      OldX=x; OldY=y;      // RESET THE OLD TO THE CURRENT (starting over)
      LeftButtonDown=true; // SET LEFT-BUTTON-DOWN FLAG
      SetInertiaOn(0);     // TURN-OFF INERTIA WHEN USER PRESSES LEFT-BUTTON
    }
    else if (LeftButtonDown)    // LEFT-BUTTON UP after LEFT-BUTTON DOWN
    {
      LeftButtonDown=false;    // UNSET LEFT-BUTTON-DOWN FLAG

      // INVOKE THE INERTIAL SYSTEM IF THE LEFT BUTTON IS BEING RELEASED, THE
      // AMOUNT OF MOVEMENT IF "ENOUGH" AS DEFINED BY MINMOVE (IN PIXELS), AND
      // THE GLOBAL InertiaEnabled FLAG IS SET (CAN BE SET BY SetInertiaEnabled).
      if ((ABS(NewX-OldX) >= MINMOVE) || (ABS(NewY-OldY) >= MINMOVE)) 
        if (InertiaEnabled)
          SetInertiaOn(1);
    }
  }
  else if (button==GLUT_RIGHT_BUTTON) {
    if (state==GLUT_DOWN) RightButtonDown = true;
    else                  RightButtonDown = false;
  }
  else if (button==GLUT_MIDDLE_BUTTON) {
    if (state==GLUT_DOWN) MiddleButtonDown = true;
    else                  MiddleButtonDown = false;
  }
  if (state == GLUT_DOWN) {
    OldX = x;
    OldY = y;
  }
}

/// The default GLUT motion function implementation
/**
 * This just calls GetGLVU()->Motion(), thus allowing Object-Oriented
 * people to customize GLVU's behavior by overriding the Motion()
 * method instead of by dealing directly with GLUTs callbacks.
 *
 * @see Motion, DefaultDisplayFunc, DefaultInertiaFunc, DefaultKeyboardFunc,
 *     DefaultMouseFunc, DefaultReshapeFunc 
*/
void GLVU::DefaultMotionFunc(int x, int y)
{
  GetGLVU()->Motion(x, y);
}

/// Handler for 'active' mouse drag events, i.e. dragging with a button pressed.
/**
 *  This method can be overridden to perform application specific
 *  funcionality (e.g. direct manipulation of scene objects).  The
 *  default implementation does some important work for handling
 *  camera manipulation (world navigation), so if you override this
 *  method you should always call GLVU::Mouse() from your override if
 *  you wish to preserve the built-in navigation capabilities.  
 *
 *  The exact effect the default implementation has on the current 
 *  camera depends upon the current world navigation mode.  
 *  See SetWorldNavMode().
 *
 *  Users not interested in creating an object-oriented app can simply
 *  call GLUT's \c glutMotionFunc to set a callback directly.  If you
 *  do so you can still call GetCurrent()->Motion() or 
 *  GLVU::DefaultMotionFunc(), in your handler to get the default 
 *  behavior back.
 *
 *  This is hooked up to the \c glutMotionFunc().  GLUT also provides
 *  a \c glutPassiveMotionFunc() which can be used to handle mouse
 *  motion when there are no mouse buttons pressed.  GLVU does not
 *  have a wrapper for that one, however.
 * 
 *  @param x,y The most recent location of the pointer.
 * 
 *  @see DefaultMouseFunc, Motion, SetWorldNavMode */
void GLVU::Motion(int x, int y)
{
  // STORE PREVIOUS NEW MOUSE POSITION (OLD)
  OldX=NewX; OldY=NewY;
  if (LeftButtonDown)
  {
    // STORE THE NEW MOUSE POSITION
    NewX=x;    NewY=y;

    int WW = glutGet((GLenum)GLUT_WINDOW_WIDTH);  // GET THE WINDOW DIMENSIONS
    int WH = glutGet((GLenum)GLUT_WINDOW_HEIGHT);

    switch(GetWorldNavMode())
    {
      // -------  WORLD NAVIGATION -------
      case NAV_MODE_TRACKBALL:
        if (CtrlPressed) { TrackBallX(OldX,NewX,WW); DriveY(OldY,NewY,WH); }
        else { TrackBallX(OldX,NewX,WW); TrackBallY(OldY,NewY,WH); }
        break;
      case NAV_MODE_HYPERBALL:
        if (CtrlPressed) { 
          HyperBall(OldX,OldY,NewX,OldY,WW,WH); DriveY(OldY,NewY,WH); }
        else { HyperBall(OldX,OldY,NewX,NewY,WW,WH); }
        break;
      case NAV_MODE_DRIVE:
        if (CtrlPressed) { TranslateX(NewX,OldX,WW); TranslateY(NewY,OldY,WH); }
        else if (AltPressed) { LookX(OldX,NewX,WW); LookY(OldY,NewY,WH); }
        else { LookX(OldX,NewX,WW); DriveY(OldY,NewY,WH); }
        break;
      case NAV_MODE_TRANSLATE:
        TranslateX(OldX,NewX,WW); 
        if (CtrlPressed) DriveY(OldY,NewY,WH); else TranslateY(OldY,NewY,WH);
        break;
      case NAV_MODE_LOOK:
        if (CtrlPressed) { TranslateX(OldX,NewX,WW); TranslateY(OldY,NewY,WH); }
        else if (AltPressed) { LookX(OldX,NewX,WW); DriveY(OldY,NewY,WH); }
        else { LookX(OldX,NewX,WW); LookY(OldY,NewY,WH); }
        break;
    };
    glutPostRedisplay();
  }
}


//----------------------------------------------------------------------------
// INERTIAL SYSTEM
//----------------------------------------------------------------------------

/// The default intertia function implementation
/**
 * This just calls GetGLVU()->Inertia(), thus allowing Object-Oriented
 * people to customize GLVU's behavior by overriding the Inertia()
 * method instead of by dealing with callbacks.
 *
 * Unlike the other methods of its ilk, this does not correspond to any GLUT
 * callback, it is purely a GLVU invention.
 *
 * @see Inertia, DefaultDisplayFunc, DefaultReshapeFunc, DefaultKeyboardFunc,
 *     DefaultMouseFunc, DefaultMotionFunc 
 */
void GLVU::DefaultInertiaFunc(int x, int y)
{
  GetGLVU()->Inertia(x, y);
}

/// Handler for inertia events
/**
 * The default implementation calls the Motion() method, causing the
 * camera to move a little more in the direction it was moving when
 * inertia kicked in.
 *
 * This method only gets called when inertia has been triggered by an
 * appropriate mouse drag and release action.
 *
 * @param x, y 
 * @see DefaultInertiaFunc, Motion, Mouse, Keyboard, Reshape, Display 
 */
void GLVU::Inertia(int x, int y)
{
  GLVU::Motion(x, y);
}

/// Call the inertia handler after some setup
/**
 * The way inertia is handled is to trick GLVU into thinking it is
 * getting the exact same mouse motion again and again, i.e, that the
 * mouse was dragged from x1,y1 to x2,y2 repeatedly.  This method munges the 
 * various internal data members as necessary to pull this trick off, then 
 * calls the inertia callback.
 * 
 * Returns TRUE (1) if inertia is enabled, FALSE (0) otherwise.
 * @see SetInertiaOn, Inertia
 */
//----------------------------------------------------------------------------
int GLVU::DoInertia()
{
  if (InertiaOn)
  {
    // SIMPLY REPEATEDLY CALL THE MOUSE MOVE HANDLER (MUST SET APPROPRIATE STATES)
    int tNewX=NewX, tNewY=NewY;          // COPY NEW VALUES TO TEMPS
    NewX=OldX; NewY=OldY;                // COPY OLD TO NEW (MOUSE MOVE COPIES)
    LeftButtonDown=true;                 // "PRETEND" BUTTON IS PRESSED
    if (UserInertiaFunc)
      UserInertiaFunc(tNewX, tNewY);     // CALL MOUSE MOVE HANDLER (INDIRECTLY)
    LeftButtonDown=false;                // RESET BUTTON STATE
    return(1);
  }
  return(0);
}

/// Turn inertia on or off.
/**
 * This is not about whether inertia is \em enabled or not, but whether 
 * it is currently active.  Usually called internally only.
 *
 * @see Inertia, DoInertia, SetInertiaEnabled
 */
void GLVU::SetInertiaOn(int onOrOff)
{
  InertiaOn = onOrOff;
  if (InertiaOn)
    InertialTimerFunc(WindowID);
}

void GLVU::InertialTimerFunc(int winID)
{
  GLVU *g = GLVUs[winID];
  if (glutGetWindow()!=winID) glutSetWindow(winID);
  g->DoInertia();
  if (g->InertiaOn) 
    glutTimerFunc(g->InertiaDelay, InertialTimerFunc, g->WindowID);
}
