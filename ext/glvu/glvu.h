//------------------------------------------------------------------------------
// File : glvu.hpp
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
// glvu.hpp : OpenGL/GLUT -based viewer
//============================================================================

#ifndef _GLVU_
#define _GLVU_

#include <stdlib.h>
#include <time.h>            // FOR TIMED PATH PLAYBACK (glvu_camview)
#include <sys/timeb.h>       // FOR FRAME TIMING (glvu_camview)
#if __APPLE__
#include <GL/glut.h>
#elif __linux__
#include <GL/glut.h>
#endif
#include <glvuVec3f.h>
#include <camera.h>

#ifndef ABS
#define ABS(x) (((x)<0)?(-(x)):x)    // HANDY UNIVERSAL ABSOLUTE VALUE FUNC
#endif
#define INERTIA_MINMOVE 2            // AMOUNT (in pixels) TO MOVE TO INVOKE INERTIAL SYSTEM
#define MAX_GLVUS 33                 // MAX NUMBER OF WINDOWS PLUS 1 SUPPORTED IN ONE APP
#define FPS_INTEGRATE_INTERVAL 1.0f // MINIMUM NUMBER OF SECONDS TO INTEGRATE FRAME-RATE OVER

/**
 * @class GLVU
 * @brief A convenient OpenGL/GLUT based 3D viewer.
 *
 * GLVU is part of a collection of C++ libraries for creating OpenGL
 * applications based on the cross-platform GLUT library.  There are
 * many sub-libraries that comprise GLVU, but the actual GLVU class is
 * the one that encapsulates all of the GLUT user interface
 * functionality.  Most of the other libraries do not depend on GLUT,
 * or depend on it in trivial ways.  Another class that is based on
 * the Qt toolkit, QGLVU, can also be used for more advanced user
 * interface capabilities.  Note however that Qt comes with some
 * significant licensing restrictions, whereas GLUT does not.
 *
 * The viewer library was developed mainly for creating small research
 * applications, and as such implements much of the functionality
 * necessary in a small viewer for prototyping rendering algorithms.
 * Among these features are a flexible camera class [Camera], several
 * camera manipulation, or navigation, modes, frame rate calculation,
 * and the ability to record and play back camera paths.  All of these
 * are features likely to be needed at one time or another by most
 * anyone that wants to look at 3D geometric scene data with OpenGL,
 * and each requires some effor to implement.  In particular,
 * developing the mouse-based camera manipulation routines can be
 * quite time consuming, as how to construct a mapping from 2D to 3D
 * (6D, really) that is intuitive for the user is far from obvious.
 * 
 * Another nice feature supported by GLVU is "inertia", or the ability
 * to start a model spinning and have it keep spinning even after the
 * user releases the mouse.  GLVU can also support multiple windows
 * easily, and the inertia works independently in each via GLUT
 * timers.  GLVU does not use the one global GLUT idle callback at
 * all, leaving that instead for the user to do with as desired.
 *
 * Since GLUT is really more of a windowing toolkit than a user
 * interface toolkit, users of the GLVUT-based version of GLVU may
 * find an external GLUT-based UI toolkit to be useful.  For quick
 * prototyping, GLUI is recommended
 * (http://www.cs.unc.edu/~rademach/glui).  For more sophisticated
 * interfaces, use Qt and QGLVU rather than this class.
 *
 * Part of the philosophy in the design of GLVU was to make it as easy
 * as possible for current GLUT users to start using it right away.
 * For this reason all of the GLUT callbacks are still accessible to
 * the GLVU programmer.  GLVU also supports a more object-oriented
 * approach as well, by allowing you to subclass GLVU and override
 * handler methods.
 *
 * @see QGLVU, Camera 
 */
class GLVU
{
  public:

    /// The world navigation modes.  See SetWorldNavMode()
    /** These are the different modes for manipulating the current
     *  camera with mouse input.
     */
    enum WorldNavMode{
      NAV_MODE_TRACKBALL, 
      NAV_MODE_HYPERBALL,
      NAV_MODE_DRIVE, 
      NAV_MODE_TRANSLATE, 
      NAV_MODE_LOOK
    };

    /// Camera IDs.  There are 4 predefined cameras in 
    /// a GLVU that can be switched between.
    enum CameraID
    {
      CAMERA_1,
      CAMERA_2,
      CAMERA_3,
      CAMERA_4,
      CAMERA_NUM_CAMS
    };
  
    /// A function type used by the 'inertia' handling system
    typedef void (*InertiaFunc)(int x, int y);

    GLVU();
    virtual ~GLVU();
    int Init(char *WindowTitle, 
             unsigned int VisualMode,
             int WindowStartX, int WindowStartY,
             int WindowWidth, int WindowHeight);
    virtual int InitWin(int aPreviouslyCreatedGlutWindow);
    void BeginFrame();
    void EndFrame();

    // CAMVIEW
    int GetPlaybackOn() const;
    int GetRecordingOn() const;
     
    // CAMVIEW
    void SetAllCams(const glvuVec3f& WorldMin, const glvuVec3f& WorldMax, 
                    const glvuVec3f& Eye, const glvuVec3f& LookAtCntr, 
                    const glvuVec3f& viewup,
                    float Yfov, float Aspect, 
                    float NearFactor, float FarFactor);
    void AllCamsPerspectiveChange(float Yfov, float Aspect, 
                                 float Ndist, float Fdist);
    void AllCamsPerspectiveAspectChange(float Aspect);
    void AllCamsPerspectiveNearFarChange(float Ndist, float Fdist);
    void AllCamsResetToOrig();
    float* GetModelviewMatrix(float* M);
    float* GetProjectionMatrix(float* M);
    void GetPixelRay(int sx, int sy, int ww, int wh, glvuVec3f *Start, glvuVec3f *Dir)
      const;
    void SelectCam(int WhichCam);
    int GetCurrentCamId() const;
    Camera* GetCurrentCam();
    void SetCurrentCam(Camera *NewCam);
    Camera* GetCam(int WhichCam);
    int GetInOutMode() const;
    void SetInOutMode(int Bool);
    void ToggleInOutMode();
    void SetOrigCam(Camera *Cam);
    const glvuVec3f& GetViewUp() const ;
    const glvuVec3f& GetWorldCenter() const ;
    void SetWorldCenter(const glvuVec3f& newCenter);
    float GetWorldRadius() const ;
    void SetWorldRadius(float newRadius) ;
    void SetViewUp(glvuVec3f viewup) ;
    void SetMoveSpeed(float speed) ;
    float GetMoveSpeed() const ;
    void TranslateX(int OldX, int NewX, int WW);
    void TranslateY(int OldY, int NewY, int WH);
    void DriveY(int OldY, int NewY, int WH);
    void LookX(int OldX, int NewX, int WW);
    void LookY(int OldY, int NewY, int WH);
    void TrackBallX(int OldX, int NewX, int WW);
    void TrackBallY(int OldY, int NewY, int WH);
    void HyperBall(int OldX, int OldY, int NewX, int NewY,int WW, int WH);
    void StartRecording(const char *FileName = "path0.dat");
    void EndRecording();
    FILE* GetPathFile() const;
    void StartPlayback(const char *FileName = "path0.dat");
    void EndPlayback();
    void StopRecordingPlayback();
    void StartFPSClock() ;
    void StopFPSClock() ;
    void DrawFPS(int xpos = 5, int ypos = 5);
    float GetFPS() const ;
    void UpdateFPS();  // USE ONLY IF YOU DO NOT CALL GLVU::BeginFrame()
    void Print(int x, int y, char* buffer);

    // MENU
    int GetMainMenuID() const ;
    int GetWorldNavMode() const ;
    void SetWorldNavMode(int mode) ;
    int GetCamDisplayOn(int WhichCam) const ;
    int IsWorldNavMode() const ;

    // MOUSE
    bool GetLeftButtonDown() const ;
    bool GetMiddleButtonDown() const ;
    bool GetRightButtonDown() const ;
    bool GetAltDown() const;
    bool GetShiftDown() const;
    bool GetCtrlDown() const;
    int GetMouseDeltaX(int curx) const;
    int GetMouseDeltaY(int cury) const;
    int GetInertiaOn() const ; 
    int GetInertiaEnabled() const ;
    void SetInertiaOn(int Bool);
    void SetInertiaEnabled(int Bool) ;
    int GetInertiaDelay() const ; 
    void SetInertiaDelay(int msec) ;
    void SetInertiaFunc(InertiaFunc f) ; 
    InertiaFunc GetInertiaFunc() const;

    // DEFAULT CALLBACK FUNCTIONS (MOSTLY JUST CALL THE CALLBACK METHODS)
    static void DefaultMouseFunc(int button, int state, int x, int y);
    static void DefaultMotionFunc(int x, int y);
    static void DefaultReshapeFunc(int WW, int WH);
    static void DefaultDisplayFunc();
    static void DefaultKeyboardFunc(unsigned char Key, int x, int y);
    static void DefaultInertiaFunc(int x, int y);
  
    // CALLBACK METHODS
    virtual void Mouse(int button, int state, int x, int y);
    virtual void Motion(int x, int y);
    virtual void Reshape(int WW, int WH);
    virtual void Display() ;
    virtual void Keyboard(unsigned char Key, int x, int y);
    virtual void Inertia(int x, int y);

    // MISC
    int GetWindowID() const;
    void MakeCurrent();
    static GLVU* GetCurrent();
    static GLVU* GetGLVU(int WindowID);
    static GLVU* GetGLVU(); // deprecated
    static void PrintVisualInfo();
    static void CheckForGLError( char *msg );


  protected:

    // GL VIEWER STATE VARIABLES
    Camera *Cams;                // ARRAY OF 4 VIEWER CAMS
    Camera *Cam;                 // PTR TO CURRENT CAM (DEFAULT IS CAM 0)
    Camera OrigCam;              // THE ORIGINAL VIEW FOR RESETTING
    int RecordingOn, PlaybackOn; // CAMERA RECORDING/PLAYBACK FLAGS
    FILE *RecordPlaybackFP;      // RECORDING/PLAYBACK FILE POINTER
    glvuVec3f WorldCenter;           // WORLD BOUNDING-SPHERE CENTER
    float WorldRadius;           // WORLD BOUNDING-SPHERE RADIUS
    glvuVec3f ViewUp;                // THE WORLD UP-VECTOR
    int InsideLookingOutMode;    // NAVIGATION MODE (IN->OUT OR OUT->IN)
    clock_t PlaybackTime;        // FOR PATH PLAYBACK TIMING
    int NumCams;
    struct timeb lastFPSClock;
    int calcFPS;
    float lastFPS;
    int lastFPSCount;

    bool LeftButtonDown, MiddleButtonDown, RightButtonDown;
    mutable int OldX, OldY;
    int NewX, NewY;  // CURRENT OLD AND NEW MOUSE POSITIONS
    float moveSpeed;
    bool CtrlPressed, AltPressed, ShiftPressed;
    int InertiaOn, InertiaEnabled;
    int InertiaDelay;

    int MainMenuID;              // GLUT MAIN MENU ID
    int WorldNavMode;            // WORLD NAVIGATION MODE
    int *CamDisplayOn;           // DISPLAY FLAG FOR EACH CAMERA

    // MENU
    static void MainMenuHandler(int value);
    static void vuOptionsMenuHandler(int value);
    static void glOptionsMenuHandler(int value);
    static void PathPlaybackMenuHandler(int value);
    static void CamViewMenuHandler(int value);
    static void CamViewDisplayMenuHandler(int value);
    static void EscapeHandler(int value);

    // MOUSE
    int DoInertia();
    InertiaFunc UserInertiaFunc;

    static GLVU *GLVUs[MAX_GLVUS];

    // glut WindowID
    int WindowID;
   
    static void InertialTimerFunc(int value);
    static void PathPlaybackTimerFunc(int value);

    virtual void InitMenu();

};

//----------------------------------------------------------------------------

/// Return whether or not playback mode is currently active
/**
 * When playback mode is active, the modelview matrices for viewing
 * are read from a previously recorded file instead of coming from the 
 * user's mousing.
 *
 * @see GetRecordingOn, StartRecording, EndRecording, 
 *     StartPlayback, EndPlayback, GetModelviewMatrix
 */
inline int GLVU::GetPlaybackOn() const
{
  return(PlaybackOn); 
}

/// Return whether or not playback mode is currently active
/**
 * When record mode is active, every camera view change made by the user
 * (through mouse manipulation) is written to a file.  The resulting
 * path can be played back later using StartPlayback().
 *
 * @see GetPlaybackOn, StartRecording, EndRecording, 
 *     StartPlayback, EndPlayback, GetModelviewMatrix
 */
inline int GLVU::GetRecordingOn() const
{
  return(RecordingOn); 
}

/// Return the file pointer being used for playback or recording
/**
 * Return the file pointer being used for playback or recording.
 * This is 0 if neither path playback nor recording is currently active.
 *
 * @see GetPlaybackOn, StartRecording, EndRecording, 
 *     StartPlayback, EndPlayback, GetModelviewMatrix
 */
inline FILE* GLVU::GetPathFile() const
{
  return RecordPlaybackFP;
}


/// Switch to the specified camera
/** 
 * Changes the camera from which the \c GLVU renders the scene.
 *
 * By default this method is bound to the 1,2,3, and 4 keys on the
 * keyboard.
 *
 * @param WhichCam can be one of the #CameraID enum values.
 * @see GetModelviewMatrix, BeginFrame, EndFrame 
*/
inline void GLVU::SelectCam(int WhichCam)
{
  Cam = &Cams[WhichCam]; 
}

/// Return the index of the active camera.
inline int GLVU::GetCurrentCamId() const
{
  return (Cam - &Cams[0]);
}


/// Return a pointer to the active camera.
inline Camera* GLVU::GetCurrentCam() 
{
  return(Cam); 
}

/// Set the active camera.
/**
 *  @param NewCam the new camera to set.
 *
 *  @note The new camera will not be owned by GLVU.  This sets the
 *  current camera \em temporarily to \p NewCam.
 *
 *  @see SelectCam
 */
inline void GLVU::SetCurrentCam(Camera *NewCam) 
{
  Cam=NewCam; 
}

/// Return a pointer to the specified camera
/**
 * @param WhichCam can be one of the \c CameraID identifiers \c
 * CAMERA_1, \c CAMERA_2, \c CAMERA_3, or \c CAMERA_4, or just a
 * number between zero and 3.  
 */
inline Camera* GLVU::GetCam(int WhichCam) 
{
  return(&Cams[WhichCam]);
}

/**
 * Returns a boolean specifying whether the navigation mode is "inside looking 
 * out" (true) or "outside looking in" (false, the default).  
 * @see SetInOutMode, ToggleInOutMode
 */   
inline int GLVU::GetInOutMode() const 
{
  return(InsideLookingOutMode); 
}

/**
 * Sets the navigation mode to either "inside looking out" (true) or
 * "outside looking in" (false, the default).  
 * @see SetInOutMode, ToggleInOutMode.
 */
inline void GLVU::SetInOutMode(int Bool) 
{
  InsideLookingOutMode=Bool; 
}

/**
 * Toggles the display mode between "inside looking out" 
 * and "outside looking in".
 * @see SetInOutMode, GetInOutMode
 */
inline void GLVU::ToggleInOutMode() 
{
  if (InsideLookingOutMode) InsideLookingOutMode=0;
  else InsideLookingOutMode=1; 
}

/**
 * Stores a copy of \p Cam so that Cameras can be reset to its parameters
 * at a later time using AllCamsResetToOrig().  Called from SetAllCams().
 * @see AllCamsResetToOrig, SetAllCams
 */
inline void GLVU::SetOrigCam(Camera *Cam) 
{
  OrigCam.Copy(*Cam); 
}

/** 
 * Returns the global "up" vector, as set by SetViewUp() or SetAllCams().
 * This has some effect on the operation of certain mouse navigation modes.
 * @see SetWorldNavMode
 */
inline const glvuVec3f& GLVU::GetViewUp() const 
{
  return(ViewUp); 
}

/**
 * Returns the global "center" of the world, as set by SetWorldCenter() or
 * (indirectly) by SetAllCams().
 */
inline const glvuVec3f& GLVU::GetWorldCenter() const 
{
  return WorldCenter; 
}

/**
 * Sets the global "center" of the world.  Also set (indirectly) by SetAllCams().
 * @see GetWorldCenter, SetWorldNavMode
 */
inline void GLVU::SetWorldCenter(const glvuVec3f& center)
{
  WorldCenter = center;
}

/**
 * Returns the global "radius" of the world, as set by SetWorldRadius() or
 * (indirectly) by SetAllCams().
 */
inline float GLVU::GetWorldRadius() const 
{
  return WorldRadius; 
}

/**
 * Sets the global "radius" of the world.
 * @see GetWorldRadius
 */
inline void GLVU::SetWorldRadius(float newRadius) 
{
  WorldRadius = newRadius; 
}

/**
 * Sets the global "up" vector.
 * @see GetViewUp
 */
inline void GLVU::SetViewUp(glvuVec3f viewup) 
{
  ViewUp=viewup; 
}

/**
 * Sets the gain factor used in translating mouse motion in pixels
 * into world units.
 * @see GetMoveSpeed, SetWorldNavMode
 */
inline void GLVU::SetMoveSpeed(float speed) 
{
  moveSpeed = speed; 
}

/**
 * Gets the gain factor used in translating mouse motion in pixels
 * into world units.
 * @see SetMoveSpeed, SetWorldNavMode
 */
inline float GLVU::GetMoveSpeed() const
{
  return moveSpeed; 
}


/**
 * Starts the frame timer and resets the frame counter.  You must call
 * this once after creating the GLVU if you would like to take
 * advantage of the built-in frame timing capabilities.
 * @see StopFPSClock, UpdateFPS, BeginFrame 
 */
inline void GLVU::StartFPSClock() 
{
  calcFPS = 1; ftime(&lastFPSClock); lastFPSCount=0;
}

/**
 *   Stops the frame rate timer.
 * @see StartFPSClock
 */
inline void GLVU::StopFPSClock() 
{
  calcFPS = 0; 
}

/**
 * Returns the last calculated Frames Per Second measurement.
 * @see StartFPSClock, StopFPSClock, UpdateFPS
 */
inline float GLVU::GetFPS() const 
{
  return lastFPS; 
}

/**
 * Returns the GLUT ID for the right-click menu automatically installed by GLVU. 
 * This can be used with GLUT calls to modify the contents of the menu.
 * @see InitMenu
 */
inline int GLVU::GetMainMenuID() const 
{
  return(MainMenuID); 
}

/**
 * Returns the current world navigation mode, (i.e. camera control mode).
 * Return value is one of the #WorldNavMode enum constants.
 * @see SetWorldNavMode
 */
inline int GLVU::GetWorldNavMode() const 
{
  return(WorldNavMode); 
}

/**
 * Sets the current world navigation mode, (i.e. camera control mode).
 * @param mode One of the #WorldNavMode enum constants.
 *
 * - \c NAV_MODE_TRACKBALL: a simple trackball.  This is primarily for
 *           rotating the model, but by holding down CTRL you can move
 *           in and out as well.  The rotation generated by the
 *           trackball mode depends only upon the relative motion of
 *           the mouse, not the absolute location.  The drawback is
 *           that this mode has no good way to make a model upright.
 *           The Hyperball is better for that.
 *
 * - \c NAV_MODE_HYPERBALL: an "SGI-style" trackball. The effect of
 *           this trackball is different depending on where on the
 *           screen you drag the mouse.  On the edge of the screen,
 *           rotation happens mostly in the screen plane, but near
 *           the middle of the window, rotation happens perpendicular
 *           to the screen plane.  Hold down CTRL to move in and out.
 *
 * - \c NAV_MODE_DRIVE: Left / Right to steer, Up / Down to move fore and aft.
 * - \c NAV_MODE_TRANSLATE: translate the camera parallel to the view plane.
 * - \c NAV_MODE_LOOK: Rotation about a fixed eye position.
 *
 * The outside-looking-in rotational modes (\c NAV_MODE_TRACKBALL,
 * \c NAV_MODE_HYPERBALL) use the "world center" as the center of rotation.
 * See SetWorldCenter() and SetAllCams() for more information about the 
 * world center and camera settings.
 *
 * GLVU uses the following default keyboard accelerators to switch between modes:
 *
 * - 'z' : \c NAV_MODE_TRACKBALL
 * - 'x' : \c NAV_MODE_DRIVE
 * - 'c' : \c NAV_MODE_TRANSLATE
 * - 'v' : \c NAV_MODE_LOOK
 * - 'h' : \c NAV_MODE_TRACKBALL
 *
 * @see GetWorldNavMode */
inline void GLVU::SetWorldNavMode(int mode) 
{
  WorldNavMode=mode; 
}

/**
 *  Returns the display status of the specified Camera.
 *  Camera display refers to rendering of some lines that show the extents of a 
 *  camera's frustum.   GLVU has some functionality to do this automically.
 *
 *  @see EndFrame
*/
inline int GLVU::GetCamDisplayOn(int WhichCam) const 
{
  return(CamDisplayOn[WhichCam]); 
}
/**
 *  Returns true (nonzero) if the current world nav mode set is valid
 *  false (0) otherwise.
 *
 *  @see SetWorldNavMode
*/
inline int GLVU::IsWorldNavMode() const 
{
  return(WorldNavMode>=0 && WorldNavMode<=3); 
}

/**
 * Returns true if the left button was down
 * last time the Mouse() or Motion() callback got called.
 */
inline bool GLVU::GetLeftButtonDown() const 
{
  return(LeftButtonDown); 
}

/**
 * Returns true if the middle button was down
 * last time the Mouse() or Motion() callback got called.
 */
inline bool GLVU::GetMiddleButtonDown() const 
{
  return(MiddleButtonDown); 
}

/**
 * Returns true if the right button was down
 * last time the Mouse() or Motion() callback got called.
 */
inline bool GLVU::GetRightButtonDown() const 
{
  return(RightButtonDown); 
}

/**
 * Returns true if the Alt key was down
 * last time the Mouse() or Keyboard() callback got called.
 */
inline bool GLVU::GetAltDown() const
{
  return (AltPressed);
}
/**
 * Returns true if the Shift key was down
 * last time the Mouse() or Keyboard() callback got called.
 */
inline bool GLVU::GetShiftDown() const
{
  return (ShiftPressed);
}
/**
 * Returns true if the Control key was down
 * last time the Mouse() or Keyboard() callback got called.
 */
inline bool GLVU::GetCtrlDown() const
{
  return (CtrlPressed);
}

/**
 * Returns the difference between the current x value, curx, passed in
 * and the last mouse X value GLVU got.
 */
inline int GLVU::GetMouseDeltaX(int curx) const
{
  int dx = curx - OldX;
  OldX = curx;
  return dx;
}

/**
 * Returns the difference between the current y value, cury, passed in
 * and the last mouse Y value GLVU got (inverts normal screen coordinate
 * system so that upward mouse motion gives a positive delta).
 */
inline int GLVU::GetMouseDeltaY(int cury) const
{
  int dy = OldY - cury;
  OldY = cury;
  return dy;
}

/**
 * Returns true (nonzero) if inertia is currently active.
 * Note that this is different from whether or not it is \em enabled. 
 * @see SetInertiaEnabled
 */
inline int GLVU::GetInertiaOn() const 
{
  return(InertiaOn); 
} 

/// Return whether inertia is enabled
/**
 * Returns TRUE (1) if inertia is enabled, FALSE (0) otherwise.
 * @see SetInertiaEnabled 
 */
inline int GLVU::GetInertiaEnabled() const 
{
  return(InertiaEnabled);
}

/// Enable or disable the GLVU's inertia feature
/**
 * @see GetInertiaEnabled
 */
inline void GLVU::SetInertiaEnabled(int Bool) 
{
  InertiaEnabled=Bool;
}


/// Return the number of milliseconds between inertia callbacks
/**
 * Inertia callbacks are only made by GLVU when inertia is active and
 * enabled.  But when triggered, inertia callbacks are made repeatedly
 * at regular intervals to animate the camera.
 *
 * @see SetInertiaDelay, GetInertiaOn, GetInertiaEnabled 
 */
inline int GLVU::GetInertiaDelay() const 
{
  return(InertiaDelay);
} 

/// Set the number of milliseconds to wait between inertia callbacks
/**
 * Inertia callbacks are only made by GLVU when inertia is active and
 * enabled.  But when triggered, inertia callbacks are made repeatedly
 * at regular intervals to animate the camera.
 *
 * The implementation relies on \c glutTimerFunc().
 *
 * @see GetInertiaDelay, GetInertiaOn, GetInertiaEnabled 
 */
inline void GLVU::SetInertiaDelay(int msec) 
{
  InertiaDelay=msec;
}

/// Set the function to use as an inertia callback
/**
 * The default implemntation works fine.  There's really no
 * reason to call this.
 */
inline void GLVU::SetInertiaFunc(GLVU::InertiaFunc f) 
{
  UserInertiaFunc=f;
} 

/// Return the function being uses as the inertia callback
/**
 * @see SetInertiaFunc
 */
inline GLVU::InertiaFunc GLVU::GetInertiaFunc() const
{
  return UserInertiaFunc;
}

/// Return the GLUT window ID associated with this GLVU
/**
 * Sometimes it is necessary to get the GLUT window ID for a GLVU
 * window, and when you do, you can call this method.  If all you want
 * to do is make this the active GLUT window (i.e. call \c
 * glutSetWindow() on it) then you can call MakeCurrent() instead.
 *
 * @see MakeCurrent, GetGLVU(int), GetGLVU 
 */
inline int GLVU::GetWindowID() const
{
  return WindowID;
}

/// Make \c this the currently active GL context
/**
 * Equivalent to \c glutSetWindow(GetWindowID()).  
 *
 * @note This bites most people at some point or other in their \c
 * glutIdleFunc.  Many new GLUT users aren't aware that the
 * glutIdleFunc is global -- one per application.  All the other GLUT
 * callbacks are per window, but not the idle func.  When your idle
 * func gets called there is no guarantee which glut window will be
 * current.  If you do a \c glutPostRedisplay() that redisplay event
 * could be sent to a window other than the one you were intending
 * unless you call MakeCurrent() first.
 *
 * None of the other GLUT callbacks have this problem.  They are all
 * per-window, and GLUT guarantees that that window will be current
 * when your callback is called.
 *
 *  @see GetGLVU(void) 
 */
inline void GLVU::MakeCurrent()
{
  glutSetWindow(WindowID);
}

/// Returns the currently active GLVU window
/**
 * See MakeCurrent() for a discussion of what it means to be the active 
 * window, and why you should care.
 *
 * @see MakeCurrent, GetGLVU(int)
 */
inline GLVU* GLVU::GetCurrent()
{
  return GLVUs[ glutGetWindow() ];
}

/// Returns a pointer to the GLVU that corresponds to the specified GLUT window ID.
/**
 * There's a one-to-one mapping between GLVUs and GLUT windows.  This
 * does the reverse lookup to find the GLVU for a given GLUT window
 * ID.  GetWindowID() does the forward lookup.
 *
 * @see MakeCurrent, GetCurrent, GetWindowID
 */
inline GLVU* GLVU::GetGLVU(int WindowID)
{
  return GLVUs[ WindowID ];
}

/// Returns the currently active GLVU window
/**
 * @deprecated
 *
 * This does the same thing as GetCurrent().  You should use that instead.
 *
 * @see MakeCurrent
 */
inline GLVU* GLVU::GetGLVU(void)
{
  return GLVUs[ glutGetWindow() ];
}


#endif






