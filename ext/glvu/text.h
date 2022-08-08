//------------------------------------------------------------------------------
// File : text.hpp
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
// text.hpp : bitmap text drawing routine; works like "printf"
//============================================================================

#ifndef _GLVU_TEXT_H_
#define _GLVU_TEXT_H_

#if __APPLE__
#include <GL/freeglut.h>
#else
#include <GL/freeglut.h>
#include <GL/glu.h>
#endif
#include <stdio.h>
#include <stdarg.h>

inline void* GetCurrentFont() { return(CurrentFont); }
inline void SetCurrentFont(void *font) { CurrentFont=font; }
inline void SetCurrentFont(int fontid) { CurrentFont=Fonts[fontid]; }

//-------------------------------------------------------------------------------
// Get the height in pixels for a specified font. Needed for calculating the
// Yoffset to go to the next line. There are several different ways:
// use the current font, specify font by id or explicitly.
//-------------------------------------------------------------------------------
inline int GetFontHeight(void *font)
{
       if (font==GLUT_BITMAP_8_BY_13) return(13);
  else if (font==GLUT_BITMAP_9_BY_15) return(15);
  else if (font==GLUT_BITMAP_TIMES_ROMAN_10) return(10);
  else if (font==GLUT_BITMAP_TIMES_ROMAN_24) return(24);
  else if (font==GLUT_BITMAP_HELVETICA_10) return(10);
  else if (font==GLUT_BITMAP_HELVETICA_12) return(12);
  else if (font==GLUT_BITMAP_HELVETICA_18) return(18);
  printf("ERROR (Text): Not a valid font!\n"); 
  return(0);
}

inline int GetFontHeight() 
{ 
  return( GetFontHeight(CurrentFont) ); 
}

inline int GetFontHeight(int fontid)
{
  return( FontHeights[fontid] );
}

//-------------------------------------------------------------------------------
// BITMAP TEXT ROUTINE: draws the printf style format string (with args) with
// the bottom-left corner of the string at pixel position (x,y) in the current 
// window with the current color. The bottom-left corner of the window is at 
// (0,0). Many different fonts are available (see CurrentFont above):
// The current raster position is updated to the end of the string; check with:
//  float CurrentRasterPos[4];
//  glGetFloatv(GL_CURRENT_RASTER_POSITION,CurrentRasterPos);
//  x=CurrentRasterPos[0]+0.5;
//  y=CurrentRasterPos[1]+0.5;
// Setting negative values for x or y will avoid the update of the raster pos;
// this is useful for continuing more text on the same line after a previous
// call to Text(). Newline chars ('\n') in the string advance the raster position
// to the next line down (Newline chars are ignored in string continuing on same line).
//-------------------------------------------------------------------------------
inline void Text(int x, int y, char *format, ...)
{
  va_list args;
  char buffer[256];
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);

  int WW = glutGet((GLenum)GLUT_WINDOW_WIDTH);
  int WH = glutGet((GLenum)GLUT_WINDOW_HEIGHT);

  glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
  glPushAttrib(GL_COLOR_BUFFER_BIT);
    glDisable(GL_DITHER);     
  glPushAttrib(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

  glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0, 0, WW, WH);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, WW, 0, WH);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    glLoadIdentity();

  int FontHeight = GetFontHeight(CurrentFont);

  if (x>=0 && y>=0)
    glRasterPos2i(x,y);

  for (char *p=buffer; *p; p++)
    if (*p=='\n')
      { if (x>=0 && y>=0) { y-=FontHeight; glRasterPos2i(x,y); } }
    else
      { glutBitmapCharacter(CurrentFont,*p); }

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glPopAttrib();

  glPopAttrib();
  glPopAttrib();
  glPopAttrib();
}

//-------------------------------------------------------------------------------
// BITMAP TEXT ROUTINE: same as Text, but give 3D coordinates in modelspace.
//-------------------------------------------------------------------------------
inline void Text3D(float x, float y, float z, char *format, ...)
{
  va_list args;
  char buffer[256];
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);

  //int WW = glutGet((GLenum)GLUT_WINDOW_WIDTH);
  //int WH = glutGet((GLenum)GLUT_WINDOW_HEIGHT);

  glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
  glPushAttrib(GL_COLOR_BUFFER_BIT);
    glDisable(GL_DITHER);     
  glPushAttrib(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

  int FontHeight = GetFontHeight(CurrentFont);

  glRasterPos3f(x,y,z);

  for (char *p=buffer; *p; p++) {
    if (*p=='\n') { 
      y-=FontHeight; glRasterPos3f(x,y,z); 
    } else { 
      glutBitmapCharacter(CurrentFont,*p); 
    }
  }

  glPopAttrib();
  glPopAttrib();
  glPopAttrib();
}

#endif
