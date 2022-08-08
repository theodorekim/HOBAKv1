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
#ifndef FFMPEG_MOVIE_H
#define FFMPEG_MOVIE_H

#include <vector>

// enables OpenGL screengrabs
#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GL/glut.h>
#elif __linux__
#include <GL/glut.h>
#endif

class FFMPEG_MOVIE {
public:
  FFMPEG_MOVIE() {
    _width = -1;
    _height = -1;
    _totalFrames = 0;
    _fps = 30;
  };

  ~FFMPEG_MOVIE() {
    // clean up all the stored frames
    for (unsigned int x = 0; x < _frames.size(); x++)
      delete[] _frames[x];
  };

  ////////////////////////////////////////////////////////////////////////
  // Grab the current OpenGL frame and add it to the movie
  //
  // The original screengrab code is from GLVU:
  // http://www.cs.unc.edu/~walk/software/glvu/
  ////////////////////////////////////////////////////////////////////////
  void addFrameGL()
  {
    GLint OldReadBuffer;
    glGetIntegerv(GL_READ_BUFFER,&OldReadBuffer);
    glReadBuffer(GL_FRONT);

    GLint OldPackAlignment;
    glGetIntegerv(GL_PACK_ALIGNMENT,&OldPackAlignment); 
    glPixelStorei(GL_PACK_ALIGNMENT,1);

    // get the screen pixels
    int width = glutGet((GLenum)GLUT_WINDOW_WIDTH);
    int height = glutGet((GLenum)GLUT_WINDOW_HEIGHT);
    int NumPixels = width * height;

    // use this to write out movies to PPM
    /*
    GLubyte* Pixels = new GLubyte[NumPixels*3];
    if (Pixels==NULL) { printf("UNABLE TO ALLOC PIXEL READ ARRAY!\n"); return; }
    glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,Pixels);
    */
    // use this if you're just streaming from memory
    GLubyte* Pixels = new GLubyte[NumPixels*4];
    if (Pixels==NULL) { printf("UNABLE TO ALLOC PIXEL READ ARRAY!\n"); return; }
    glReadPixels(0,0,width,height,GL_RGBA,GL_UNSIGNED_BYTE,Pixels);

    // store the screen pixels
    if (_width == -1 && _height == -1)
    {
      _width = width;
      _height = height;
    }
    _frames.push_back(Pixels);

    glPixelStorei(GL_PACK_ALIGNMENT,OldPackAlignment);
    glReadBuffer((GLenum)OldReadBuffer);

    _totalFrames++;
  };

  ////////////////////////////////////////////////////////////////////////
  // dump out the final movie, plus PPMs
  ////////////////////////////////////////////////////////////////////////
  void writeMovie(const char* moviename) 
  {
    using namespace std;

    if (_frames.size() == 0) 
    {
      cout << " No frames to output to movie! " << endl;
      return;
    }

    // we can create a temp directory, right?
    const string mkdir("mkdir temp");
    system(mkdir.c_str());

    // write all the frames out as PPM
    cout << " Writing frame " << flush;
    for (unsigned int x = 0; x < _frames.size(); x++)
    {
      cout << x << " " << flush;
      char filename[1024];
      sprintf(filename, "./temp/frame_%04i.ppm", x);
      FILE* fp = fopen(filename, "wb");
      if (fp == NULL)
      {
        cout << " Couldn't open file " << filename << "!!!" << endl;
        exit(0);
      }
      fprintf(fp, "P6\n%d %d\n255\n", _width, _height);
      fwrite(_frames[x], 1, _width * _height * 3, fp);
      fclose(fp);
    }
    cout << " done." << endl;

    // call FFMPEG
    string ffmpeg("ffmpeg -framerate 60 -i ./temp/frame_%04d.ppm -pix_fmt yuv420p -vf vflip -vcodec libx264 -y -b:v 20000k ");
    ffmpeg = ffmpeg + string(moviename); 
    system(ffmpeg.c_str());

    // delete all the PPM files
    cout << " Cleaning up ... " << flush;
    for (unsigned int x = 0; x < _frames.size(); x++)
    {
      char filename[1024];
      sprintf(filename, "./temp/frame_%04i.ppm", x);
      string rm("rm " + std::string(filename));
      system(rm.c_str());
    }
    cout << " done. " << endl;
  };

  ////////////////////////////////////////////////////////////////////////
  // stream out the final movie, without using files
  ////////////////////////////////////////////////////////////////////////
  void streamWriteMovie(const char* moviename) 
  {
    using namespace std;

    if (_frames.size() == 0) 
    {
      cout << " No frames to output to movie! " << endl;
      return;
    }

    char buffer[1024];
    sprintf(buffer, "ffmpeg -r 60 -f rawvideo -pix_fmt rgba -s %ix%i -i - ", _width, _height);
    string ffmpegCmd = string(buffer) + string("-threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip -vcodec libx264 ");
    ffmpegCmd = ffmpegCmd + string(moviename);

    FILE* ffmpeg = popen(ffmpegCmd.c_str(), "w");

    for (unsigned int x = 0; x < _frames.size(); x++)
      fwrite(_frames[x], sizeof(int) * _width * _height, 1, ffmpeg);

    pclose(ffmpeg);
  }

private:
  std::vector<GLubyte*> _frames;
  int _totalFrames;
  int _width;
  int _height;
  int _fps;
};

#endif

