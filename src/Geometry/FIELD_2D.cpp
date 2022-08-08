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
#include "FIELD_2D.h"
#include "TIMER.h"
#include <iostream>

using namespace std;

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D::FIELD_2D(const int& rows, const int& cols) :
  _xRes(rows), _yRes(cols)
{
  _totalCells = _xRes * _yRes;
  _data = new REAL[_totalCells];

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

FIELD_2D::FIELD_2D(const FIELD_2D& m) :
  _xRes(m.xRes()), _yRes(m.yRes())
{
  _totalCells = _xRes * _yRes;
  _data = new REAL[_totalCells];

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];
}

FIELD_2D::FIELD_2D() :
  _xRes(0), _yRes(0), _totalCells(0), _data(NULL)
{
}


FIELD_2D::FIELD_2D(const MATRIX& m) :
  _xRes(m.cols()), _yRes(m.rows())
{
  _totalCells = _xRes * _yRes;
  _data = new REAL[_totalCells];

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      (*this)(x,y) = m(_yRes - 1 - y, x);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D::~FIELD_2D()
{
  delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::clear()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::write(string filename) const
{
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_2D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);

  // always write out as a double
  if (sizeof(REAL) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    for (int x = 0; x < _totalCells; x++)
      dataDouble[x] = _data[x];

    fwrite((void*)dataDouble, sizeof(double), _totalCells, file);
    delete[] dataDouble;
    fclose(file);
  }
  else
    fwrite((void*)_data, sizeof(REAL), _totalCells, file);
  fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::read(string filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_2D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  // read dimensions
  fread((void*)&_xRes, sizeof(int), 1, file);
  fread((void*)&_yRes, sizeof(int), 1, file);
  _totalCells = _xRes * _yRes;
  if (_data) delete[] _data;
  _data = new REAL[_totalCells];

  // always read in as a double
  if (sizeof(REAL) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    fread((void*)dataDouble, sizeof(double), _totalCells, file);

    for (int x = 0; x < _totalCells; x++)
      _data[x] = dataDouble[x];

    delete[] dataDouble;
  }
  else
    fread((void*)_data, sizeof(REAL), _totalCells, file);
  fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::writePPM(string filename)
{
  FILE *fp;
  unsigned char* pixels = new unsigned char[3 * _totalCells];

  for (int x = 0; x < _totalCells; x++)
  {
    pixels[3 * x] = 255 * _data[x];
    pixels[3 * x + 1] = 255 * _data[x];
    pixels[3 * x + 2] = 255 * _data[x];
  }

  fp = fopen(filename.c_str(), "wb");
  fprintf(fp, "P6\n%d %d\n255\n", _xRes, _yRes);
  fwrite(pixels, 1, _totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::writeMatlab(string filename, string variableName) const
{
  FILE* file;
  file = fopen(filename.c_str(), "w");
  fprintf(file, "%s = [", variableName.c_str());
  for (int y = 0; y < _yRes; y++)
  {
    for (int x = 0; x < _xRes; x++)
      fprintf(file, "%f ", (*this)(x,y));
    fprintf(file, "; ");
  }
  fprintf(file, "];\n");

  fclose(file);
}

// to minimize the number of dependencies in HOBAK (even FFTW!),
// commenting this out for now
#if 0
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::FFT(FIELD_2D& real, FIELD_2D& im)
{
  fftw_complex* forward = NULL;
  fftw_complex* backward= NULL;
  fftw_plan forwardPlan;

  // if it's the first time, create the FFT vars
  forward  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * _totalCells);
  backward = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * _totalCells);

  // resize the fields if needed
  real.resizeAndWipe(_xRes, _yRes);
  im.resizeAndWipe(_xRes, _yRes);

  forwardPlan = fftw_plan_dft_2d(_xRes, _yRes, forward, backward, FFTW_FORWARD, FFTW_ESTIMATE);  

  // populate the forward with the height field
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      int index = x + _xRes * y;
      forward[index][0] = _data[index];
      forward[index][1] = 0.0f;
    }

  // run the FFT
  fftw_execute(forwardPlan);

  // copy to displayable arrays
  REAL* realData = real.data();
  REAL* imData = im.data();
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      int index = x + _xRes * y;
      realData[index] = backward[index][0];
      imData[index] = backward[index][1];
    }

  real.shiftFFT();
  im.shiftFFT();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::shiftFFT()
{
  assert(_xRes % 2 == 1);
  assert(_yRes % 2 == 1);
  REAL* scratch = new REAL[_xRes * _yRes];

  int xHalf = _xRes / 2;
  int yHalf = _yRes / 2;

  int xMod = _xRes % 2;
  int yMod = _yRes % 2;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < xHalf + xMod; x++)
    {
      int index = x + y * _xRes;
      scratch[index] = _data[index + xHalf + xMod];
    }
  for (int y = 0; y < _yRes; y++)
    for (int x = xHalf; x < _xRes; x++)
    {
      int index = x + y * _xRes;
      scratch[index] = _data[index - xHalf];
    }

  for (int y = 0; y < yHalf + yMod; y++)
    for (int x = 0; x < _xRes; x++)
    {
      int original = x + y * _xRes;
      int copy = x + (y + yHalf) * _xRes;
      _data[copy] = scratch[original];
    }

  for (int y = yHalf; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      int original = x + y * _xRes;
      int copy = x + (y - yHalf - yMod) * _xRes;
      _data[copy] = scratch[original];
    }

  delete[] scratch;
}
#endif

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::normalize()
{
  REAL maxFound = 0.0;
  REAL minFound = _data[0];
  for (int x = 0; x < _totalCells; x++)
  {
    maxFound = (_data[x] > maxFound) ? _data[x] : maxFound;
    minFound = (_data[x] < minFound) ? _data[x] : minFound;
  }

  float range = 1.0 / (maxFound - minFound);
  for (int x = 0; x < _totalCells; x++)
    _data[x] = (_data[x] - minFound) * range;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::abs()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = fabs(_data[x]);

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::resizeAndWipe(int xRes, int yRes)
{
  if (_xRes == xRes && _yRes == yRes)
  {
    clear();
    return;
  }

  if (_data)
    delete[] _data;

  _xRes = xRes;
  _yRes = yRes;
  _totalCells = _xRes * _yRes;

  _data = new REAL[_xRes * _yRes];
}

///////////////////////////////////////////////////////////////////////
// set this field to the result of convolving filter and input
///////////////////////////////////////////////////////////////////////
void FIELD_2D::convolve(const FIELD_2D& filter, const FIELD_2D& input)
{
  TIMER functionTimer(__FUNCTION__);

  int filterWidth = filter.xRes();

  assert(filterWidth == filter.yRes());
  assert(filterWidth % 2 == 1);
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);

  int iwidth = input.xRes();
  int iheight = input.yRes();

  // periodic boundaries
  int halfWidth = filterWidth / 2;
  for (int ix = 0; ix < iwidth; ix++)
    for (int iy = 0; iy < iheight; iy++)
    {
      int index = ix + iwidth * iy;
      float vd = 0;
      for (int iix = -halfWidth; iix <= halfWidth; iix++)
        for (int iiy = -halfWidth; iiy <= halfWidth; iiy++)
        {
          int xIndex = ix + iix;
          int yIndex = iy + iiy;

          if (xIndex < 0)
            xIndex += iwidth;
          if (yIndex < 0)
            yIndex += iheight;

          if (xIndex >= iwidth)
            xIndex -= iwidth;
          if (yIndex >= iheight)
            yIndex -= iheight;
          
          vd += filter(iix + halfWidth,iiy + halfWidth) * input(xIndex, yIndex);
        }
      _data[index] = vd;
    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator=(const REAL& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator*=(const REAL& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator/=(const REAL& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] /= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator+=(const REAL& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator-=(const FIELD_2D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator+=(const FIELD_2D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  for (int x = 0; x < _totalCells; x++)
    _data[x] += input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator*=(const FIELD_2D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);

  for (int x = 0; x < _totalCells; x++)
    _data[x] *= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator/=(const FIELD_2D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);

  for (int x = 0; x < _totalCells; x++)
    if (fabs(input[x]) > 1e-6)
      _data[x] /= input[x];
    else
      _data[x] = 0;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator*(const FIELD_2D& A, const REAL alpha)
{
  FIELD_2D result(A);
  result *= alpha;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator/(const FIELD_2D& A, const REAL alpha)
{
  FIELD_2D result(A);
  result /= alpha;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator+(const FIELD_2D& A, const FIELD_2D& B)
{
  FIELD_2D result(A);
  result += B;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator-(const FIELD_2D& A, const FIELD_2D& B)
{
  FIELD_2D result(A);
  result -= B;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator+(const FIELD_2D& A, const REAL alpha)
{
  FIELD_2D result(A);
  result += alpha;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator*(const REAL alpha, const FIELD_2D& A)
{
  return A * alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator+(const REAL alpha, const FIELD_2D& A)
{
  return A + alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator=(const FIELD_2D& A)
{
  resizeAndWipe(A.xRes(), A.yRes());

  for (int x = 0; x < _totalCells; x++)
    _data[x] = A[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
// sum of all entries
///////////////////////////////////////////////////////////////////////
REAL FIELD_2D::sum()
{
  REAL total = 0;
  for (int x = 0; x < _totalCells; x++)
    total += _data[x];

  return total;
}

///////////////////////////////////////////////////////////////////////
// take the log
///////////////////////////////////////////////////////////////////////
void FIELD_2D::log(REAL base)
{
  REAL scale = 1.0 / std::log(base);
  for (int x = 0; x < _totalCells; x++)
    _data[x] = std::log(_data[x]) * scale;
}

///////////////////////////////////////////////////////////////////////
// Compute the elements of the vertical derivative convolution kernel
///////////////////////////////////////////////////////////////////////
void FIELD_2D::verticalDerivativeKernel(double kMax, double dk, double sigma, double L)
{
  assert(_xRes % 2);
  assert(_xRes == _yRes);

  double norm = 0;

  for (double k = 0; k < kMax; k += dk)
    norm += k * k * exp(-sigma * k * k);

  int halfWidth = _xRes / 2;

  for (int i = -halfWidth; i <= halfWidth; i++)
    for (int j = -halfWidth; j <= halfWidth; j++)
    {
      double r = sqrt((float)(i * i + j * j));
      double kern = 0;
      for (double k = 0; k < kMax; k += dk)
        kern += k * k * exp(-sigma * k * k) * j0(r * k);

      (*this)(i + halfWidth, j + halfWidth) = kern / norm;
    }
}

///////////////////////////////////////////////////////////////////////
// Compute a radial Bessel function
///////////////////////////////////////////////////////////////////////
void FIELD_2D::radialBessel()
{
  assert(_xRes % 2);
  assert(_xRes == _yRes);

  int halfWidth = _xRes / 2;
  int kMin = 0;
  int kMax = 5;
  double dk = 0.01;

  for (int i = -halfWidth; i <= halfWidth; i++)
    for (int j = -halfWidth; j <= halfWidth; j++)
    {
      double r = sqrt((float)(i * i + j * j)) / (float)_xRes * 20;
      double kern = 0.0;
      for (double k = kMin; k < kMax; k += dk)
        kern += j0(r * k);

      (*this)(i + halfWidth, j + halfWidth) = kern;
    }
}

///////////////////////////////////////////////////////////////////////
// set to a bessel function
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToBessel(float k)
{
  assert(_xRes % 2);
  assert(_xRes == _yRes);

  int halfWidth = _xRes / 2;

  for (int i = -halfWidth; i <= halfWidth; i++)
    for (int j = -halfWidth; j <= halfWidth; j++)
    {
      double r = sqrt((float)(i * i + j * j)) / (float)_xRes * 20;
      double kern = 0.0;
      kern += j0(r * k);

      (*this)(i + halfWidth, j + halfWidth) = kern;
    }
}

///////////////////////////////////////////////////////////////////////
// upsample the texture by a certain factor
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::nearestNeighborUpsample(int factor)
{
  FIELD_2D result(_xRes * factor, _yRes * factor);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      double value = (*this)(x,y);
      for (int i = 0; i < factor; i++)
        for (int j = 0; j < factor; j++)
          result(factor * x + i, factor * y + j) = value;
    }
  return result;
}

///////////////////////////////////////////////////////////////////////
// get the min of the field
///////////////////////////////////////////////////////////////////////
REAL FIELD_2D::min()
{
  assert(_xRes > 0);
  assert(_yRes > 0);
  REAL result = _data[0];

  for (int i = 0; i < _xRes * _yRes; i++)
    result = (_data[i] < result) ? _data[i] : result;

  return result;
}

///////////////////////////////////////////////////////////////////////
// get the max of the field
///////////////////////////////////////////////////////////////////////
REAL FIELD_2D::max()
{
  assert(_xRes > 0);
  assert(_yRes > 0);
  REAL result = _data[0];

  for (int i = 0; i < _xRes * _yRes; i++)
    result = (_data[i] > result) ? _data[i] : result;

  return result;
}

///////////////////////////////////////////////////////////////////////
// set to a checkboard for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToCheckerboard(int xChecks, int yChecks)
{
  int index = 0;
  for (int x = 0; x < _xRes; x++)
    for (int y = 0; y < _yRes; y++, index++)
    {
      int xMod = (x / (_xRes / xChecks)) % 2;
      int yMod = (y / (_yRes / yChecks)) % 2;

      if ((xMod && yMod) || (!xMod && !yMod))
        _data[index] = 1;
    }
}

///////////////////////////////////////////////////////////////////////
// set to a checkboard for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToRampedCheckerboard(int xChecks, int yChecks)
{
  int index = 0;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++, index++)
    {
      int xMod = (x / (_xRes / xChecks)) % 2;
      int yMod = (y / (_yRes / yChecks)) % 2;

      if ((xMod && yMod) || (!xMod && !yMod))
        _data[index] = (float)y / _yRes;
    }
}

///////////////////////////////////////////////////////////////////////
// pass a field to fieldViewer2D
///////////////////////////////////////////////////////////////////////
void FIELD_2D::fieldViewer(const FIELD_2D& field, string name)
{
  field.write("temp.field");
  string execute("./bin/fieldViewer temp.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// set to a ramp for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToRampX()
{
  int index = 0;
  for (int x = 0; x < _xRes; x++)
    for (int y = 0; y < _yRes; y++, index++)
      _data[index] = (float)x / _xRes;
}

///////////////////////////////////////////////////////////////////////
// set to a ramp for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToRampY()
{
  int index = 0;
  for (int x = 0; x < _xRes; x++)
    for (int y = 0; y < _yRes; y++, index++)
      _data[index] = (float)y / _yRes;
}

///////////////////////////////////////////////////////////////////////
// get the projection of the field in the x direction
///////////////////////////////////////////////////////////////////////
VECTOR FIELD_2D::projectionX()
{
  VECTOR result(_yRes);

  int index = 0;
  for (int x = 0; x < _xRes; x++)
    for (int y = 0; y < _yRes; y++, index++)
      result[y] += _data[index];

  return result; 
}

///////////////////////////////////////////////////////////////////////
// return a field for the Laplacian of this field
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::laplacian()
{
  FIELD_2D result(_xRes, _yRes);

  // assume a unit field
  REAL dx = 1.0 / _xRes;
  REAL dy = 1.0 / _yRes;
  REAL dx2 = dx * dx;
  REAL dy2 = dy * dy;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      REAL yPlus  = safe(x, y+1); 
      REAL yMinus = safe(x, y-1); 
      REAL xPlus  = safe(x+1, y); 
      REAL xMinus = safe(x-1, y);

      result(x,y) = (-2.0 * (*this)(x,y) + xPlus + xMinus) / dx2 + (-2.0 * (*this)(x,y) + yPlus + yMinus) / dy2; 
    }
  return result;
}

///////////////////////////////////////////////////////////////////////
// return a field for the Laplacian of this field
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::laplacian4th()
{
  FIELD_2D result(_xRes, _yRes);

  // assume a unit field
  REAL dx = 1.0 / _xRes;
  REAL dy = 1.0 / _yRes;
  REAL dx2 = dx * dx;
  REAL dy2 = dy * dy;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      REAL yPlus  = safe(x, y+1); 
      REAL yPlusPlus  = safe(x, y+2); 
      REAL yMinus = safe(x, y-1); 
      REAL yMinusMinus = safe(x, y-2); 
      REAL xPlus  = safe(x+1, y); 
      REAL xPlusPlus  = safe(x+2, y); 
      REAL xMinus = safe(x-1, y);
      REAL xMinusMinus = safe(x-2, y);

      result(x,y) = ((-1.0 / 12.0) * xMinusMinus + (4.0 / 3.0) * xMinus + (-5.0 / 2.0) * (*this)(x,y) + (4.0 / 3.0) * xPlus + (-1.0 / 12.0) * xPlusPlus) / dx2;
      result(x,y) += ((-1.0 / 12.0) * yMinusMinus + (4.0 / 3.0) * yMinus + (-5.0 / 2.0) * (*this)(x,y) + (4.0 / 3.0) * yPlus + (-1.0 / 12.0) * yPlusPlus) / dy2;
    }
  return result;
}

///////////////////////////////////////////////////////////////////////
// a safe, toroidal data accessor -- does all bounds checking for you
///////////////////////////////////////////////////////////////////////
REAL FIELD_2D::safe(int x, int y)
{
  while (x < 0)
    x += _xRes;
  while (y < 0)
    y += _yRes;

  while (x > _xRes - 1)
    x -= _xRes;
  while (y > _yRes - 1)
    y -= _yRes;

  return (*this)(x,y);
}

///////////////////////////////////////////////////////////////////////
// return a field for the gradient of this field
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::gradient()
{
  FIELD_2D result(_xRes, _yRes);

  // assume a unit field
  REAL dx = 1.0 / _xRes;
  REAL dy = 1.0 / _yRes;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      REAL yPlus  = safe(x, y+1); 
      REAL yMinus = safe(x, y-1); 
      REAL xPlus  = safe(x+1, y); 
      REAL xMinus = safe(x-1, y);

      result(x,y) = (xPlus - xMinus) / dx + (yPlus - yMinus) / dy;
    }
  return result;
}

///////////////////////////////////////////////////////////////////////
// return the transpose (flip x and y)
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::transpose() const
{
  FIELD_2D result(_yRes, _xRes);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      result(y,x) = (*this)(x,y);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3 FIELD_2D::maxIndex()
{
  REAL maxFound = _data[0];

  VECTOR3 maxFoundIndex;
  int index = 0;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++, index++)
      if (_data[index] > maxFound)
      {
        maxFound = _data[index];

        maxFoundIndex[0] = x;
        maxFoundIndex[1] = y;
      }

  return maxFoundIndex;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3 FIELD_2D::minIndex()
{
  REAL minFound = _data[0];

  VECTOR3 minFoundIndex;
  int index = 0;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++, index++)
      if (_data[index] < minFound)
      {
        minFound = _data[index];

        minFoundIndex[0] = x;
        minFoundIndex[1] = y;
      }

  return minFoundIndex;
}

} // HOBAK
