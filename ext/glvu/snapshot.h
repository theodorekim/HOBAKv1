//------------------------------------------------------------------------------
// File : snapshot.hpp
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
// snapshot.hpp : TAKES A "SNAPSHOT" OF THE WINDOW (WRITES THE FRAMEBUFFER
// TO A .PPM FILE)
//============================================================================

#ifndef _GLVU_SNAPSHOT_H_ 
#define _GLVU_SNAPSHOT_H_ 

/**
 * Writes a snapshot of the current screen in ppm format. 
 * The file written is a file in the current working directory 
 * with a name of the form "snap*.ppm".  Each time SnapShot is 
 * called a search is performed to find a name of that for that 
 * doesn't already exist, and that is what is written to.
 *
 * Implemented using glReadPixels.
 */
inline void SnapShot(void)
{
  FILE *fp;

  // GENERATE A UNIQUE FILENAME
  char FileName[20];
  static int SnapShotNum=100;
  int UniqueFound=0;  
  do { sprintf(FileName,"snap%d.ppm",SnapShotNum);
       if ((fp=fopen(FileName,"r"))) fclose(fp); else UniqueFound=1;
       SnapShotNum++;
     } while(!UniqueFound);

  GLint OldReadBuffer;
  glGetIntegerv(GL_READ_BUFFER,&OldReadBuffer);
  glReadBuffer(GL_FRONT);

  GLint OldPackAlignment;
  glGetIntegerv(GL_PACK_ALIGNMENT,&OldPackAlignment); 
  glPixelStorei(GL_PACK_ALIGNMENT,1);

  int WW = glutGet((GLenum)GLUT_WINDOW_WIDTH);
  int WH = glutGet((GLenum)GLUT_WINDOW_HEIGHT);
  int NumPixels = WW*WH;
  GLubyte* Pixels = new GLubyte[NumPixels*3];
  if (Pixels==NULL) { printf("UNABLE TO ALLOC PIXEL READ ARRAY!\n"); return; }
  glReadPixels(0,0,WW,WH,GL_RGB,GL_UNSIGNED_BYTE,Pixels);

  fp = fopen(FileName, "wb");
  fprintf(fp, "P6\n%d %d\n255\n", WW, WH);
  fwrite(Pixels,1,NumPixels*3,fp);
  fclose(fp);
  delete[] Pixels;

  glPixelStorei(GL_PACK_ALIGNMENT,OldPackAlignment);
  glReadBuffer((GLenum)OldReadBuffer);
}


#endif
