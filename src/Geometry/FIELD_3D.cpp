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
#include "FIELD_3D.h"
#include "TIMER.h"
#include <omp.h>
#include <iostream>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace std;

namespace HOBAK {

int FIELD_3D::_quinticClamps = 0;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D::FIELD_3D(const int& xRes, const int& yRes, const int& zRes,
    const VECTOR3& center, const VECTOR3& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
 
  try {
    _data = new REAL[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(REAL);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _outside = maxRes() * maxRes();
}

FIELD_3D::FIELD_3D(const double* data, const int& xRes, const int& yRes, const int& zRes,
    const VECTOR3& center, const VECTOR3& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  try {
    _data = new REAL[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(REAL);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _outside = maxRes() * maxRes();

  for (int x = 0; x < _totalCells; x++)
    _data[x] = data[x];
}

FIELD_3D::FIELD_3D(const FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths())
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  try {
    _data = new REAL[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(REAL);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];

  _outside = maxRes() * maxRes();
}

FIELD_3D::FIELD_3D() :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL)
{
}

FIELD_3D::FIELD_3D(const vector<FIELD_2D>& slices)
{
  assert(slices.size() > 0);

  _xRes = slices[0].xRes();
  _yRes = slices[0].yRes();
  _zRes = slices.size();

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  
  try {
    _data = new REAL[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(REAL);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        (*this)(x,y,z) = slices[z](x,y);

  _outside = maxRes() * maxRes();
}

FIELD_3D::FIELD_3D(const char* filename) :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL)
{
  cout << " Reading file " << filename << endl;
  //int size = string(filename).size();
  /*
  if (filename[size - 1] == 'z' && filename[size - 2] == 'g')
    readGz(filename);
  else
    read(filename);
    */
  read(filename);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D::~FIELD_3D()
{
  if (_data)
    delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clear()
{
  TIMER functionTimer(__FUNCTION__);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        _data[x + y * _xRes + z * _slabSize] = 0.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::write(string filename) const
{
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  cout << " Writing file " << filename.c_str() << " ... "; flush(cout);

  // write to the stream
  write(file);

  // close the stream
  fclose(file);

  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::read(string filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  // read from the stream
  read(file);

  // close the file
  fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::resizeAndWipe(int xRes, int yRes, int zRes, const VECTOR3& center, const VECTOR3& lengths)
{
  if (_xRes == xRes && _yRes == yRes && _zRes == zRes)
  {
    _center = center;
    _lengths = lengths;
    clear();

    _dx = _lengths[0] / _xRes;
    _dy = _lengths[1] / _yRes;
    _dz = _lengths[2] / _zRes;
    _invDx = 1.0 / _dx;
    _invDy = 1.0 / _dy;
    _invDz = 1.0 / _dz;
    return;
  }

  if (_data) delete[] _data;

  _xRes = xRes;
  _yRes = yRes;
  _zRes = zRes;
  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  _outside = maxRes() * maxRes();
  _center = center;
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  try {
    _data = new REAL[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(REAL);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator=(const REAL& alpha)
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] = alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator*=(const REAL& alpha)
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator/=(const REAL& alpha)
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] /= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator+=(const REAL& alpha)
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator-=(const REAL& alpha)
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator-=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator+=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] += input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator*=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator^(const FIELD_3D& A, const REAL alpha)
{
  FIELD_3D result(A);

  for (int x = 0; x < result.totalCells(); x++)
    result[x] = pow(result[x], alpha);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const FIELD_3D& A, const REAL alpha)
{
  FIELD_3D result(A);
  result *= alpha;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator/(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D result(A);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < A.totalCells(); x++)
    result[x] *= 1.0 / B[x];

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator/(const FIELD_3D& A, const REAL alpha)
{
  FIELD_3D result(A);
  result /= alpha;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator-(const FIELD_3D& A, const FIELD_3D& B)
{
  assert(A.xRes() == B.xRes());
  assert(A.yRes() == B.yRes());
  assert(A.zRes() == B.zRes());

  FIELD_3D result(A);
  result -= B;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D result(A);
  result += B;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D result(A);
  result *= B;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const FIELD_3D& A, const REAL alpha)
{
  FIELD_3D result(A);
  result += alpha;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const REAL alpha, const FIELD_3D& A)
{
  return A * alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const REAL alpha, const FIELD_3D& A)
{
  return A + alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator=(const FIELD_3D& A)
{
  resizeAndWipe(A.xRes(), A.yRes(), A.zRes(), A.center(), A.lengths());

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] = A[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
// return a slice in the form of a FIELD_2D
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_3D::zSlice(int z) const
{
  FIELD_2D result(_xRes, _yRes);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      result(x,y) = (*this)(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// real-valued cell center coordinates
///////////////////////////////////////////////////////////////////////
VECTOR3 FIELD_3D::cellCenter(int x, int y, int z) const
{
  VECTOR3 halfLengths = (REAL)0.5 * _lengths;

  // set it to the lower corner
  VECTOR3 result = _center - halfLengths;

  // displace to the NNN corner
  result[0] += x * _dx;
  result[1] += y * _dy;
  result[2] += z * _dz;

  // displace it to the cell center
  result[0] += _dx * 0.5;
  result[1] += _dy * 0.5;
  result[2] += _dz * 0.5;

  return result;
}

///////////////////////////////////////////////////////////////////////
// do a union with 'field', assuming both this and field are signed i
// distance fields
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::signedDistanceUnion(const FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);

  FIELD_3D result = field;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if ((*this)(x,y,z) > 0.0 && field(x,y,z) > 0.0)
          result(x,y,z) = (*this)(x,y,z) < field(x,y,z) ? (*this)(x,y,z) : field(x,y,z);
        else if ((*this)(x,y,z) < 0.0 && field(x,y,z) < 0.0)
          result(x,y,z) = (*this)(x,y,z) < field(x,y,z) ? (*this)(x,y,z) : field(x,y,z);
        else if ((*this)(x,y,z) < 0.0 && field(x,y,z) > 0.0)
          result(x,y,z) = (*this)(x,y,z);
        else 
          result(x,y,z) = field(x,y,z);
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// what's the maximum resolution in any direction?
///////////////////////////////////////////////////////////////////////
int FIELD_3D::maxRes()
{
  int result = _xRes;
  if (_yRes > result) result = _yRes;
  if (_zRes > result) result = _zRes;
  return result;
}

///////////////////////////////////////////////////////////////////////
// assuming that SURFACE.initializeSignedDistanceField() has been 
// called on this field, do the fast marching method
///////////////////////////////////////////////////////////////////////
void FIELD_3D::fastMarchingMethod()
{
  MIN_HEAP minHeap;

  // clear the retired nodes
  _retired.clear();
 
  // insert the ones for the forward marching
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        REAL distance = (*this)(x,y,z);
        if (distance < _outside)
        {
          if (distance > 0.0f)
          {
            HEAP_ENTRY entry;
            entry.distance = distance;
            entry.index = index;
            minHeap.insert(entry);
          }
        }
      }

  // march forward
  marchOneway(true, minHeap);

  // insert the ones for the backward marching
  minHeap.clear();
  index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        REAL distance = (*this)(x,y,z);
        if (distance <= 0.0f)
        {
          HEAP_ENTRY entry;
          entry.distance = distance;
          entry.index = index;
          minHeap.insert(entry);
        }
        if (distance >= _outside)
          (*this)(x,y,z) *= -1.0;
      }
  
  // march backward
  marchOneway(false, minHeap);

  // stomp the retired nodes hash when we're done
  _retired.clear();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::fastMarchingNegativeOnly()
{
  MIN_HEAP minHeap;

  // clear the retired nodes
  _retired.clear();
 
  // insert the ones for the backward marching
  minHeap.clear();
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        REAL distance = (*this)(x,y,z);
        if (distance <= 0.0f)
        {
          HEAP_ENTRY entry;
          entry.distance = distance;
          entry.index = index;
          minHeap.insert(entry);
        }
        if (distance >= _outside)
          (*this)(x,y,z) *= -1.0;
      }
  
  // march backward
  marchOneway(false, minHeap);

  // stomp the retired nodes hash when we're done
  _retired.clear();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::fastMarchingPositiveOnly()
{
  MIN_HEAP minHeap;

  // clear the retired nodes
  _retired.clear();
 
  // insert the ones for the backward marching
  minHeap.clear();
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        REAL distance = (*this)(x,y,z);
        if (distance > 0.0f && distance < _outside)
        {
          // make sure one of the neighbors is unmarched
          bool found = false;
          if (x != 0 && (*this)(x-1,y,z) >= _outside)
            found = true;
          if (x != _xRes - 1 && (*this)(x+1,y,z) >= _outside)
            found = true;
          if (y != 0 && (*this)(x,y-1,z) >= _outside)
            found = true;
          if (y != _yRes - 1 && (*this)(x,y+1,z) >= _outside)
            found = true;
          if (z != 0 && (*this)(x,y,z-1) >= _outside)
            found = true;
          if (z != _zRes - 1 && (*this)(x,y,z+1) >= _outside)
            found = true;

          if (!found) continue;

          cout << " Adding " << x << " " << y << " " << z << endl;

          HEAP_ENTRY entry;
          entry.distance = distance;
          entry.index = index;
          minHeap.insert(entry);
        }
      }

  cout << " positive heap size: " << minHeap.size() << endl;
  
  // march backward
  marchOneway(true, minHeap);

  // stomp the retired nodes hash when we're done
  _retired.clear();
}

///////////////////////////////////////////////////////////////////////
// quicksort comparator
///////////////////////////////////////////////////////////////////////
int compare(const void *arg1, const void *arg2)
{
  return *(REAL*)arg1 > *(REAL*)arg2;
}

///////////////////////////////////////////////////////////////////////
// do fast marching in one direction
///////////////////////////////////////////////////////////////////////
void FIELD_3D::marchOneway(bool forward, MIN_HEAP& minHeap)
{
  int pops = 0;

  // do the forward marching
  while (!minHeap.empty())
  {
    // pop off the top and retire it
    HEAP_ENTRY popped = minHeap.popMin();
    _retired[popped.index] = true;

    pops++;

    // popped was just retired, calculate a new distance value for all its neighbors
    int candidate = popped.index;
    int zIndex = candidate / _slabSize;
    int yIndex = (candidate - zIndex * _slabSize) / _xRes;
    int xIndex = (candidate - zIndex * _slabSize - yIndex * _xRes);

    for (int y = 0; y < 6; y++)
    {
      candidate = popped.index;

      if (y == 0)
      {
        if (xIndex < _xRes - 1) candidate++;
        else continue;
      }

      if (y == 1) 
      {
        if (xIndex > 0) candidate--;
        else continue;
      }

      if (y == 2)
      {
        if (yIndex < _yRes - 1) candidate += _xRes;
        else continue;
      }

      if (y == 3) 
      {
        if (yIndex > 0) candidate -= _xRes;
        else continue;
      }

      if (y == 4) 
      {
        if (zIndex < _zRes - 1) candidate += _slabSize;
        else continue;
      }

      if (y == 5)
      {
        if (zIndex > 0) candidate -= _slabSize;
        else continue;
      }

      // account for marching direction
      bool distanceTest = _data[candidate] >= 0.0 ? true : false;
      if (!forward) distanceTest = !distanceTest;

      if (distanceTest && (_retired.find(candidate) == _retired.end()))
      {
        // get the distance values
        REAL distances[3];
        REAL neighbors[6];

        int zIndex = candidate / _slabSize;
        int yIndex = (candidate - zIndex * _slabSize) / _xRes;
        int xIndex = (candidate - zIndex * _slabSize - yIndex * _xRes);

        // x plus and minus
        neighbors[0] = (xIndex < _xRes - 1) ? _data[candidate + 1] : _data[candidate];
        neighbors[1] = (xIndex > 0)         ? _data[candidate - 1] : _data[candidate];

        // y plus and minus
        neighbors[2] = (yIndex < _yRes - 1) ? _data[candidate + _xRes] : _data[candidate];
        neighbors[3] = (yIndex > 0)         ? _data[candidate - _xRes] : _data[candidate];

        // z plus and minus
        neighbors[4] = (zIndex < _zRes - 1) ? _data[candidate + _slabSize] : _data[candidate];
        neighbors[5] = (zIndex > 0)         ? _data[candidate - _slabSize] : _data[candidate];
      
        // find the store the upwind direction along each axis
        //float sign = (forward) ? 1.0f : -1.0f;
        for (int x = 0; x < 3; x++)
        {
          // store the positive direction
          distances[x] = neighbors[2 * x];
          
          // clamp out the wrong side
          if (forward && distances[x] < 0.0f) distances[x] = _outside;
          if (!forward && distances[x] > 0.0f) distances[x] = -_outside;
          
          // see if the negative direction is more upwind
          if (forward)
          {
            if (neighbors[2 * x + 1] < distances[x] && neighbors[2 * x + 1] > 0.0f)
              distances[x] = neighbors[2 * x + 1];
          }
          else if (neighbors[2 * x + 1] > distances[x] && neighbors[2 * x + 1] < 0.0f)
            distances[x] = neighbors[2 * x + 1];
        }

        // clamp out values from the wrong side
        for (int x = 0; x < 3; x++)
          if (forward)
            distances[x] = (distances[x] >= 0.0f) ? distances[x] : _outside;
          else
            distances[x] = (distances[x] <= 0.0f) ? -distances[x] : _outside;
        
        // sort the distances
        qsort((void*)distances, 3, sizeof(REAL), compare);

        // set up the quadratics David's way
        REAL b[3];
        REAL c[3];
        b[0] = distances[0];
        c[0] = distances[0] * distances[0] - 1.0f;
        for (int x = 1; x < 3; x++)
        {
          b[x] = distances[x] + b[x-1];
          c[x] = distances[x] * distances[x] + c[x-1];
        }

        // solve for the right one
        int i = 2;
        REAL discrim = b[i] * b[i] - (REAL)(i + 1) * c[i];
        REAL newDist = (b[i] + sqrtf(discrim)) / (REAL)(i + 1);
        while ((discrim < 0.0f || newDist < distances[i]) && i != -1)
        {
          i--;
          discrim = b[i] * b[i] - (REAL)(i + 1) * c[i];
          
          if (discrim > 0.0f)
            newDist = (b[i] + sqrtf(discrim)) / (REAL)(i + 1);
        }
        if (i == -1)
          cout << __FILE__ << " " << __LINE__ 
               << " Couldn't solve the distance field quadratic!" << endl;
       
        // try and insert the new distance value
        if (newDist < _outside)
        {
          REAL minus = forward ? 1.0f : -1.0f;

          // if it has never been on the heap
          if (_data[candidate] == minus * _outside)
          {
            _data[candidate] = minus * newDist;
            HEAP_ENTRY entry;
            entry.distance = _data[candidate];
            entry.index = candidate;
            minHeap.insert(entry);
          }
          else if (minus * _data[candidate] > newDist)
          {
            _data[candidate] = minus * newDist;
            minHeap.decreaseKey(candidate, _data[candidate]);
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////
// struct for second-order marching
///////////////////////////////////////////////////////////////////////
struct TRIPLET {
  REAL distance;
  REAL secondDistance;
  bool secondFound;
};

///////////////////////////////////////////////////////////////////////
// quicksort comparator for second order fast marching
///////////////////////////////////////////////////////////////////////
int compareTriplets(const void *arg1, const void *arg2) { 
  return ((TRIPLET*)arg1)->distance > ((TRIPLET*)arg2)->distance;
}

//////////////////////////////////////////////////////////////////////
// do fast marching in one direction, second order 
//////////////////////////////////////////////////////////////////////
void FIELD_3D::marchOneway2ndOrder(bool forward, MIN_HEAP& minHeap)
{
  int pops = 0;

  // do the forward marching
  while (!minHeap.empty())
  {
    // pop off the top and retire it
    HEAP_ENTRY popped = minHeap.popMin();
    _retired[popped.index] = true;

    pops++;
  
    // popped was just retired, calculate a new distance value for all its neighbors
    int candidate = popped.index;
    int zIndex = candidate / _slabSize;
    int yIndex = (candidate - zIndex * _slabSize) / _xRes;
    int xIndex = (candidate - zIndex * _slabSize - yIndex * _xRes);

    for (int y = 0; y < 6; y++)
    {
      candidate = popped.index;

      if (y == 0)
      {
        if (xIndex < _xRes - 1) candidate++;
        else continue;
      }

      if (y == 1) 
      {
        if (xIndex > 0) candidate--;
        else continue;
      }

      if (y == 2)
      {
        if (yIndex < _yRes - 1) candidate += _xRes;
        else continue;
      }

      if (y == 3) 
      {
        if (yIndex > 0) candidate -= _xRes;
        else continue;
      }

      if (y == 4) 
      {
        if (zIndex < _zRes - 1) candidate += _slabSize;
        else continue;
      }

      if (y == 5)
      {
        if (zIndex > 0) candidate -= _slabSize;
        else continue;
      }

      // account for marching direction
      bool distanceTest = _data[candidate] >= 0.0 ? true : false;
      if (!forward) distanceTest = !distanceTest;

      if (distanceTest && (_retired.find(candidate) == _retired.end()))
      {
        // triplet of distances and second order distances
        TRIPLET triplets[3];
        
        // update index breakdown for new candidate
        zIndex = candidate / _slabSize;
        yIndex = (candidate - zIndex * _slabSize) / _xRes;
        xIndex = (candidate - zIndex * _slabSize - yIndex * _xRes);

        // for each cardinal direction, add to the quadratic coefficients
        for (int x = 0; x < 3; x++)
        {
          // get the indices for the plus and minus directions
          int plusIndex = 0;
          int minusIndex = 0;
          if (x == 0) 
          { 
            plusIndex  = (xIndex < _xRes - 1) ?  1 : 0;
            minusIndex = (xIndex > 0)         ? -1 : 0; 
          }
          if (x == 1) 
          { 
            plusIndex  = (yIndex < _yRes - 1) ?  _xRes : 0;
            minusIndex = (yIndex > 0)         ? -_xRes : 0;
          }
          if (x == 2) 
          { 
            plusIndex  = (zIndex < _zRes - 1) ?  _slabSize : 0; 
            minusIndex = (zIndex > 0)         ? -_slabSize : 0;
          }

          // get the distance values
          REAL sign = (forward) ? 1.0f : -1.0f;
          REAL plusDistance = (plusIndex != 0) ? _data[candidate + plusIndex] : sign * _outside;
          REAL minusDistance = (minusIndex != 0) ? _data[candidate + minusIndex] : sign * _outside;

          // clamp out points from the wrong side
          if (sign * plusDistance < 0.0f)  plusDistance  = sign * _outside;
          if (sign * minusDistance < 0.0f) minusDistance = sign * _outside;

          // see which one is more upwind
          REAL distance = (fabs(plusDistance) < fabs(minusDistance)) ? plusDistance : minusDistance;
          triplets[x].distance = distance;

          // look for a second order point
          triplets[x].secondFound = false;
          int secondIndex = (fabs(plusDistance) < fabs(minusDistance)) ? 
                            candidate + 2 * plusIndex : 
                            candidate + 2 * minusIndex;
          int zSecondIndex = secondIndex / _slabSize;
          int ySecondIndex = (secondIndex - zSecondIndex * _slabSize) / _xRes;
          int xSecondIndex = (secondIndex - zSecondIndex * _slabSize - ySecondIndex * _xRes);
         
          bool found = true; 
          if (x == 0) 
          {
            if (xSecondIndex < 0)         found = false;
            if (xSecondIndex > _xRes - 1) found = false;
          }
          if (x == 1) 
          {
            if (ySecondIndex < 0)         found = false;
            if (ySecondIndex > _yRes - 1) found = false;
          }
          if (x == 2) 
          {
            if (zSecondIndex < 0)         found = false;
            if (zSecondIndex > _zRes - 1) found = false;
          }
          if (found == false)
          {
            triplets[x].secondDistance = sign * _outside; 
            continue;
          }

          // store the point and check if its valid
          triplets[x].secondDistance = _data[secondIndex];
          bool decreasing = (sign * _data[secondIndex]) <= (sign * distance);

          // make sure it's upwind and finalized
          if ((_retired.find(secondIndex) != _retired.end()) && decreasing)
            triplets[x].secondFound = true;
        }

        // for backwards marching, invert all the values
        if (!forward)
          for (int x = 0; x < 3; x++)
          {
            triplets[x].distance = -triplets[x].distance;
            triplets[x].secondDistance = -triplets[x].secondDistance;
          }

        // sort the distances
        qsort((void*)triplets, 3, sizeof(TRIPLET), compareTriplets);

        // calculate possible discriminants
        REAL a[3];
        REAL b[3];
        REAL c[3];
        REAL discrims[3];
        if (triplets[0].secondFound)
        {
          a[0] = 9.0f / 4.0f;
          b[0] = -6.0f * triplets[0].distance + 1.5f * triplets[0].secondDistance;
          c[0] = 4.0f * triplets[0].distance * triplets[0].distance - 
                 2.0f * triplets[0].distance * triplets[0].secondDistance + 
                 0.25f * triplets[0].secondDistance * triplets[0].secondDistance - 1.0f;
        }
        else
        {
          a[0] = 1.0f;
          b[0] = -2.0f * triplets[0].distance;
          c[0] = triplets[0].distance * triplets[0].distance - 1.0f;
        }
        discrims[0] = b[0] * b[0] - 4.0f * a[0] * c[0];
        for (int x = 1; x < 3; x++)
        {
          if (triplets[x].secondFound)
          {
            a[x] = 9.0f / 4.0f + a[x-1];
            b[x] = -6.0f * triplets[x].distance + 1.5f * triplets[x].secondDistance + b[x-1];
            c[x] = 4.0f * triplets[x].distance * triplets[x].distance - 
                   2.0f * triplets[x].distance * triplets[x].secondDistance + 
                   0.25f * triplets[x].secondDistance * triplets[x].secondDistance + c[x-1];
          }
          else
          {
            a[x] = 1.0f + a[x-1];
            b[x] = -2.0f * triplets[x].distance + b[x-1];
            c[x] = triplets[x].distance * triplets[x].distance + c[x-1];
          }
          discrims[x] = b[x] * b[x] - 4.0f * a[x] * c[x];
        }
     
        // find the first valid discriminant
        int i = 2;
        while (i != -1.0f && discrims[i] <= 0.0f) i--;
        if (i == -1)
        {
          cout << __FILE__ << " " << __LINE__ << " Second Order Fast Marching: no valid discriminant found! " << endl;
          continue;
        }
        
        // solve the quadratic
        REAL newDist = (-b[i] + sqrtf(discrims[i])) / (2.0f * a[i]);

        // try and insert the new distance value
        if (newDist < _outside)
        {
          REAL minus = forward ? 1.0f : -1.0f;

          // if it has never been on the heap
          if (_data[candidate] == minus * _outside)
          {
            _data[candidate] = minus * newDist;
            HEAP_ENTRY entry;
            entry.distance = _data[candidate];
            entry.index = candidate;
            minHeap.insert(entry);
          }
          else if (minus * _data[candidate] > newDist)
          {
            _data[candidate] = minus * newDist;
            minHeap.decreaseKey(candidate, _data[candidate]);
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////
// return the indices of the grid points insides a world-space bounding box
///////////////////////////////////////////////////////////////////////
void FIELD_3D::boundingBoxIndices(const VECTOR3& mins, const VECTOR3& maxs, VECTOR3I& iMins, VECTOR3I& iMaxs)
{
  VECTOR3 ds(_dx, _dy, _dz);

  VECTOR3 fracMins = (mins - (ds * (REAL)0.5) - _center + (_lengths * (REAL)0.5));
  VECTOR3 fracMaxs = (maxs - (ds * (REAL)0.5) - _center + (_lengths * (REAL)0.5));

  // fatten things to account for the integer cast
  iMins[0] = fracMins[0] / _dx - 2;
  iMins[1] = fracMins[1] / _dy - 2;
  iMins[2] = fracMins[2] / _dz - 2;

  iMaxs[0] = fracMaxs[0] / _dx + 2;
  iMaxs[1] = fracMaxs[1] / _dy + 2;
  iMaxs[2] = fracMaxs[2] / _dz + 2;

  // clamp them
  for (int x = 0; x < 3; x++)
  {
    if (iMins[x] < 0) iMins[x] = 0;
    if (iMaxs[x] < 0) iMaxs[x] = 0;
  }
  if (iMins[0] > _xRes - 1) iMins[0] = _xRes - 1;
  if (iMins[1] > _yRes - 1) iMins[1] = _yRes - 1;
  if (iMins[2] > _zRes - 1) iMins[2] = _zRes - 1;

  if (iMaxs[0] > _xRes - 1) iMaxs[0] = _xRes - 1;
  if (iMaxs[1] > _yRes - 1) iMaxs[1] = _yRes - 1;
  if (iMaxs[2] > _zRes - 1) iMaxs[2] = _zRes - 1;
}

///////////////////////////////////////////////////////////////////////
// copy out the boundary
///////////////////////////////////////////////////////////////////////
void FIELD_3D::copyBorderAll()
{
  int index;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      // front slab
      index = x + y * _xRes;
      _data[index] = _data[index + _slabSize];

      // back slab
      index += _totalCells - _slabSize;
      _data[index] = _data[index - _slabSize];
    }

  for (int z = 0; z < _zRes; z++)
    for (int x = 0; x < _xRes; x++)
    {
      // bottom slab
      index = x + z * _slabSize;
      _data[index] = _data[index + _xRes];

      // top slab
      index += _slabSize - _xRes;
      _data[index] = _data[index - _xRes];
    }

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
    {
      // left slab
      index = y * _xRes + z * _slabSize;
      _data[index] = _data[index + 1];

      // right slab
      index += _xRes - 1;
      _data[index] = _data[index - 1];
    }
}

///////////////////////////////////////////////////////////////////////
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
const REAL FIELD_3D::operator()(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  VECTOR3 corner = _center - (REAL)0.5 * _lengths;
  VECTOR3 dxs(_dx, _dy, _dz);
  corner += (REAL)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  // get interpolation weights
  const REAL s1 = positionCopy[0]- x0;
  const REAL s0 = 1.0f - s1;
  const REAL t1 = positionCopy[1]- y0;
  const REAL t0 = 1.0f - t1;
  const REAL u1 = positionCopy[2]- z0;
  const REAL u0 = 1.0f - u1;

  const int i000 = x0 + y0 * _xRes + z0 * _slabSize;
  const int i010 = x0 + y1 * _xRes + z0 * _slabSize;
  const int i100 = x1 + y0 * _xRes + z0 * _slabSize;
  const int i110 = x1 + y1 * _xRes + z0 * _slabSize;
  const int i001 = x0 + y0 * _xRes + z1 * _slabSize;
  const int i011 = x0 + y1 * _xRes + z1 * _slabSize;
  const int i101 = x1 + y0 * _xRes + z1 * _slabSize;
  const int i111 = x1 + y1 * _xRes + z1 * _slabSize;

  // interpolate
  // (indices could be computed once)
  return u0 * (s0 * (t0 * _data[i000] + t1 * _data[i010]) +
               s1 * (t0 * _data[i100] + t1 * _data[i110])) +
         u1 * (s0 * (t0 * _data[i001] + t1 * _data[i011]) +
               s1 * (t0 * _data[i101] + t1 * _data[i111]));
}

///////////////////////////////////////////////////////////////////////
// summed squared entries
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::sumSq()
{
  REAL result = 0;
  for (int i = 0; i < _totalCells; i++)
    result += _data[i] * _data[i];

  return result;
}

///////////////////////////////////////////////////////////////////////
// maximum entry
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::max()
{
  REAL result = _data[0];

  for (int i = 0; i < _totalCells; i++)
    if (_data[i] > result)
      result = _data[i];

  return result;
}

///////////////////////////////////////////////////////////////////////
// maximum entry
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::absMax()
{
  REAL result = fabs(_data[0]);

  for (int i = 0; i < _totalCells; i++)
    if (fabs(_data[i]) > result)
      result = fabs(_data[i]);

  return result;
}

///////////////////////////////////////////////////////////////////////
// build a const field with the given dims
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::constField(const FIELD_3D& dims, REAL value)
{
  FIELD_3D result(dims);
  result = value;
  return result;
}

///////////////////////////////////////////////////////////////////////
// check if any entry is a nan
///////////////////////////////////////////////////////////////////////
bool FIELD_3D::isNan()
{
  int i = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, i++)
        if (isnan(_data[i]))
        {
          cout << " Nan found at: " << x << ", " << y << ", " << z << endl;
          return true;
        }

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute the inverse
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::inverse()
{
  FIELD_3D result(*this);

  for (int x = 0; x < _totalCells; x++)
    result[x] = 1.0 / _data[x];

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::printNeighborhood(int index) const
{
  int z = index / _slabSize;
  int y = (index - z * _slabSize) / _xRes;
  int x = (index - z * _slabSize - y * _xRes);

  printNeighborhood(x,y,z);
}

///////////////////////////////////////////////////////////////////////
// print the neighborhood of a cell for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_3D::printNeighborhood(int x, int y, int z) const
{
  cout << " Neighborhood of (" << x << ", " << y << ", " << z << ")" << endl;
  for (int k = z - 1; k <= z + 1; k++)
  {
    cout.precision(8);
    cout << "[ " << endl;
    for (int j = y - 1; j <= y + 1; j++)
    {
      for (int i = x - 1; i <= x + 1; i++)
      {
        int index = i + j * _xRes + k * _slabSize;
        cout << _data[index] << " ";
      }
      cout << ";" << endl;
    }
    cout << "]" << endl;
  }
}

///////////////////////////////////////////////////////////////////////
// clamp nans to some specified value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clampNans(REAL value)
{
  for (int x = 0; x < _totalCells; x++)
    if (isnan(_data[x]))
      _data[x] = value;
}

///////////////////////////////////////////////////////////////////////
// clamp nans to some specified value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clampInfs(REAL value)
{
  for (int x = 0; x < _totalCells; x++)
    if (isinf(_data[x]))
      _data[x] = value;
}

///////////////////////////////////////////////////////////////////////
// clamp infinities to values in this field
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clampInfs(FIELD_3D& clampField)
{
  for (int x = 0; x < _totalCells; x++)
    if (isinf(_data[x]))
      _data[x] = clampField[x];
}

///////////////////////////////////////////////////////////////////////
// clamp nans to values in this field
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clampNans(FIELD_3D& clampField)
{
  for (int x = 0; x < _totalCells; x++)
    if (isnan(_data[x]))
      _data[x] = clampField[x];
}

///////////////////////////////////////////////////////////////////////
// extend some scalar quantity off of a front, given a signed distance function
///////////////////////////////////////////////////////////////////////
void FIELD_3D::fastExtension(const FIELD_3D& signedDistance)
{
  cout << " Extending scalars ... "; flush(cout);
  // assume that initializeExtensionScalars has been called elsewhere

  // insert the front for marching onto the heap 
  MIN_HEAP minHeap;

  // clear the retired nodes
  _retired.clear();

  FIELD_3D distanceCopy(signedDistance);
  insertFront(true, distanceCopy, minHeap);

  // march forward
  extendOneway(true, distanceCopy, minHeap);

  insertFront(false, distanceCopy, minHeap);

  extendOneway(false, distanceCopy, minHeap);

  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// cell index of a real-valued position
///////////////////////////////////////////////////////////////////////
int FIELD_3D::cellIndex(VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  VECTOR3 corner = _center - (REAL)0.5 * _lengths;
  VECTOR3 dxs(_dx, _dy, _dz);
  corner += (REAL)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  // get interpolation weights
  const REAL s1 = positionCopy[0]- x0;
  const REAL s0 = 1.0f - s1;
  const REAL t1 = positionCopy[1]- y0;
  const REAL t0 = 1.0f - t1;
  const REAL u1 = positionCopy[2]- z0;
  const REAL u0 = 1.0f - u1;

  int xFinal = (s0 > s1) ? x0 : x1;
  int yFinal = (t0 > t1) ? y0 : y1;
  int zFinal = (u0 > u1) ? z0 : z1;

  return xFinal + yFinal * _xRes + zFinal * _slabSize;
}

///////////////////////////////////////////////////////////////////////
// insert the front in preparation for reinitialization or extension
///////////////////////////////////////////////////////////////////////
void FIELD_3D::insertFront(const bool forward, FIELD_3D& distance, MIN_HEAP& minHeap)
{
  // insert the ones for the forward marching
  minHeap.clear();
  for (int i = 0; i < _totalCells; i++)
  {
    bool compare = (forward) ? distance[i] >= 0.0f && distance[i] < _outside
                             : distance[i] <= 0.0f;

    REAL sum = 0.0f;
    if (compare)
    {
      int zIndex = i / _slabSize;
      int yIndex = (i - zIndex * _slabSize) / _xRes;
      int xIndex = (i - zIndex * _slabSize - yIndex * _xRes);
      REAL center = distance[i];

      // x plus and minus
      REAL xPlus  = (xIndex < _xRes - 1) ? distance[i + 1] : _outside;
      REAL xMinus = (xIndex > 0)         ? distance[i - 1] : _outside;

      // y plus and minus
      REAL yPlus  = (yIndex < _yRes - 1) ? distance[i + _xRes] : _outside;
      REAL yMinus = (yIndex > 0)         ? distance[i - _xRes] : _outside;

      // z plus and minus
      REAL zPlus  = (zIndex < _zRes - 1) ? distance[i + _slabSize] : _outside;
      REAL zMinus = (zIndex > 0)         ? distance[i - _slabSize] : _outside;

      REAL interpolate;
      compare = (forward) ? true : yPlus < _outside;
      if (yPlus * center <= 0.0f && compare)
      {
        interpolate = center / (center - yPlus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : yMinus < _outside;
      if (yMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - yMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : xMinus < _outside;
      if (xMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - xMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : xPlus < _outside;
      if (xPlus * center <= 0.0f && compare)
      {
        interpolate = center / (center - xPlus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : zMinus < _outside;
      if (zMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - zMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : zPlus < _outside;
      if (zPlus * center <= 0.0f && compare)
      {
        interpolate = center / (center - zPlus);
        sum += 1.0f / (interpolate * interpolate);
      }

      REAL sign = (forward) ? 1.0f : -1.0f;

      if (sum > 0.0f)
      {
        REAL finalDistance = sign / sqrtf(sum);
       
        HEAP_ENTRY entry;
        entry.distance = finalDistance;
        entry.index = i;
        distance[i] = finalDistance;
        minHeap.insert(entry);
      }
      else
        distance[i] = sign * _outside;
    }
    if (!forward && distance[i] >= _outside)
      distance[i] = -_outside;
  }
}

//////////////////////////////////////////////////////////////////////
// do extension in one direction
//////////////////////////////////////////////////////////////////////
void FIELD_3D::extendOneway(bool forward, FIELD_3D& distance, MIN_HEAP& minHeap)
{
  int pops = 0;
  
  // do the forward marching
  while (!minHeap.empty())
  {
    // pop off the top and retire it
    HEAP_ENTRY popped = minHeap.popMin();
    _retired[popped.index] = true;

    pops++;

    // popped was just retired, calculate a new distance value for all its neighbors
    int candidate = popped.index;
    int zIndex = candidate / _slabSize;
    int yIndex = (candidate - zIndex * _slabSize) / _xRes;
    int xIndex = (candidate - zIndex * _slabSize - yIndex * _xRes);
  
    for (int y = 0; y < 6; y++)
    {
      candidate = popped.index;

      if (y == 0)
      {
        if (xIndex < _xRes - 1) candidate++;
        else continue;
      }

      if (y == 1) 
      {
        if (xIndex > 0) candidate--;
        else continue;
      }

      if (y == 2)
      {
        if (yIndex < _yRes - 1) candidate += _xRes;
        else continue;
      }

      if (y == 3) 
      {
        if (yIndex > 0) candidate -= _xRes;
        else continue;
      }

      if (y == 4) 
      {
        if (zIndex < _zRes - 1) candidate += _slabSize;
        else continue;
      }

      if (y == 5)
      {
        if (zIndex > 0) candidate -= _slabSize;
        else continue;
      }

      // account for marching direction
      bool distanceTest = distance[candidate] >= 0.0 ? true : false;
      if (!forward) distanceTest = !distanceTest;

      if (distanceTest && (_retired.find(candidate) == _retired.end()))
      {
        // get the distance values
        REAL distances[3];
        REAL extensions[3];
        REAL neighbors[6];
        REAL extensionNeighbors[6];

        int zIndex = candidate / _slabSize;
        int yIndex = (candidate - zIndex * _slabSize) / _xRes;
        int xIndex = (candidate - zIndex * _slabSize - yIndex * _xRes);

        // x plus and minus
        extensionNeighbors[0] = (xIndex < _xRes - 1) ? _data[candidate + 1] : _data[candidate];
        extensionNeighbors[1] = (xIndex > 0)         ? _data[candidate - 1] : _data[candidate];
        neighbors[0] = (xIndex < _xRes - 1) ? distance[candidate + 1] : distance[candidate];
        neighbors[1] = (xIndex > 0)         ? distance[candidate - 1] : distance[candidate];

        // y plus and minus
        extensionNeighbors[2] = (yIndex < _yRes - 1) ? _data[candidate + _xRes] : _data[candidate];
        extensionNeighbors[3] = (yIndex > 0)         ? _data[candidate - _xRes] : _data[candidate];
        neighbors[2] = (yIndex < _yRes - 1) ? distance[candidate + _xRes] : distance[candidate];
        neighbors[3] = (yIndex > 0)         ? distance[candidate - _xRes] : distance[candidate];

        // z plus and minus
        extensionNeighbors[4] = (zIndex < _zRes - 1) ? _data[candidate + _slabSize] : _data[candidate];
        extensionNeighbors[5] = (zIndex > 0)         ? _data[candidate - _slabSize] : _data[candidate];
        neighbors[4] = (zIndex < _zRes - 1) ? distance[candidate + _slabSize] : distance[candidate];
        neighbors[5] = (zIndex > 0)         ? distance[candidate - _slabSize] : distance[candidate];
      
        for (int x = 0; x < 3; x++)
        {
          // store the positive direction
          distances[x] = neighbors[2 * x];
          extensions[x] = extensionNeighbors[2 * x];
          
          // clamp out the wrong side
          if (forward && distances[x] < 0.0f)  distances[x] = _outside;
          if (!forward && distances[x] > 0.0f) distances[x] = -_outside;
          
          // see if the negative direction is more upwind
          if (forward)
          {
            if (neighbors[2 * x + 1] < distances[x] && neighbors[2 * x + 1] > 0.0f)
            {
              distances[x] = neighbors[2 * x + 1];
              extensions[x] = extensionNeighbors[2 * x + 1];
            }
          }
          else if (neighbors[2 * x + 1] > distances[x] && neighbors[2 * x + 1] < 0.0f)
          {
            distances[x] = neighbors[2 * x + 1];
            extensions[x] = extensionNeighbors[2 * x + 1];
          }
        }

        // clamp out values from the wrong side
        for (int x = 0; x < 3; x++)
          if (forward)
            distances[x] = (distances[x] >= 0.0f) ? distances[x] : _outside;
          else
            distances[x] = (distances[x] <= 0.0f) ? -distances[x] : _outside;
      
        // do a 3-element bubble sort 
        REAL temp;
        REAL tempPoint;
        if (distances[0] > distances[1])
        {
          temp = distances[0];
          distances[0] = distances[1];
          distances[1] = temp;
          tempPoint = extensions[0];
          extensions[0] = extensions[1];
          extensions[1] = tempPoint;
        }
        if (distances[1] > distances[2])
        {
          temp = distances[1];
          distances[1] = distances[2];
          distances[2] = temp;
          tempPoint = extensions[1];
          extensions[1] = extensions[2];
          extensions[2] = tempPoint;
        }
        if (distances[0] > distances[1])
        {
          temp = distances[0];
          distances[0] = distances[1];
          distances[1] = temp;
          tempPoint = extensions[0];
          extensions[0] = extensions[1];
          extensions[1] = tempPoint;
        }
        
        // set up the quadratics David's way
        REAL b[3];
        REAL c[3];
        b[0] = distances[0];
        c[0] = distances[0] * distances[0] - 1.0f;
        for (int x = 1; x < 3; x++)
        {
          b[x] = distances[x] + b[x-1];
          c[x] = distances[x] * distances[x] + c[x-1];
        }

        // solve for the right one
        int i = 2;
        REAL discrim = b[i] * b[i] - (REAL)(i + 1) * c[i];
        REAL newDist = (b[i] + sqrtf(discrim)) / (REAL)(i + 1);
        while ((discrim < 0.0f || newDist < distances[i]) && i != -1)
        {
          i--;
          discrim = b[i] * b[i] - (REAL)(i + 1) * c[i];
          
          if (discrim > 0.0f)
            newDist = (b[i] + sqrtf(discrim)) / (REAL)(i + 1);
        }
        if (i == -1)
          cout << __FILE__ << " " << __LINE__ 
               << " Couldn't solve the distance field quadratic!" << endl;
       
        // get the extension
        if (newDist < fabs(distance[candidate]))
        {
          if (i == 2)
          {
            REAL interpolate = 1.0f / (3.0f * newDist - distances[2] - distances[1] - distances[0]);
            REAL diffs[] = {(newDist - distances[0]), (newDist - distances[1]), (newDist - distances[2])};
            _data[candidate] = (diffs[0] * extensions[0] +
                                diffs[1] * extensions[1] +
                                diffs[2] * extensions[2]) * interpolate;
          }
          else if (i == 1)
          {
            REAL interpolate = 1.0f / (2.0f * newDist - distances[1] - distances[0]);
            REAL diffs[] = {(newDist - distances[0]), (newDist - distances[1])};
            _data[candidate] = (diffs[0] * extensions[0] +
                                diffs[1] * extensions[1]) * interpolate;
          }
          else
            _data[candidate] = extensions[0];
        }

        // try and insert the new distance value
        if (newDist < _outside)
        {
          REAL minus = forward ? 1.0f : -1.0f;

          // if it has never been on the heap
          if (fabs(distance[candidate]) >= _outside)
          {
            distance[candidate] = minus * newDist;
            HEAP_ENTRY entry;
            entry.distance = distance[candidate];
            entry.index = candidate;
            minHeap.insert(entry);
          }
          else if (minus * distance[candidate] > newDist)
          {
            distance[candidate] = minus * newDist;
            minHeap.decreaseKey(candidate, distance[candidate]);
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////
// pass a field to fieldViewer3D
///////////////////////////////////////////////////////////////////////
void FIELD_3D::fieldViewer(const FIELD_3D& field, string name)
{
  field.write("temp3d.field");
  string execute("./fieldViewer3D temp3d.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// pass a field to fieldViewer3D
///////////////////////////////////////////////////////////////////////
void FIELD_3D::fieldViewerYZ(const FIELD_3D& field, string name)
{
  field.write("temp3d.field");
  string execute("./fieldViewer3D temp3d.field \"");
  execute = execute + name + string("\" yz &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// pass a field with a distance field (inside/outside) overlay
// to a overlayFieldViewer3D
///////////////////////////////////////////////////////////////////////
void FIELD_3D::overlayFieldViewerYZ(const FIELD_3D& field, const FIELD_3D& distance, string name)
{
  field.write("temp3d.field");
  distance.write("temp3d.distance.field");
  string execute("./overlayFieldViewer3D temp3d.field temp3d.distance.field \"");
  execute = execute + name + string("\" yz &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// get the integer indices of a spatial position
///////////////////////////////////////////////////////////////////////
void FIELD_3D::indices(const VECTOR3& position, int* x)
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  VECTOR3 corner = _center - (REAL)0.5 * _lengths;
  VECTOR3 dxs(_dx, _dy, _dz);
  corner += (REAL)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  x[0] = (int)positionCopy[0];
  x[1] = (int)positionCopy[1];
  x[2] = (int)positionCopy[2];
}

///////////////////////////////////////////////////////////////////////
// load a PhysBAM level set
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readPhysBAM(const char* filename)
{
  if (_data) delete[] _data;

  // level set contains
  //
  // counts (TV_INT) vector<int, 3>
  // domain (RANGE) probably 6 floats
  // mac_offset (T) single floating point, 0 or 0.5
  //
  // then is reads in a scalar array, which contains
  //
  // length2 - an int, which seems to always equal 1
  // domain (RANGE<TV>) not TV_INT, float ranges, again?
  // the entries - scalars, 
  // everything appears to be single precision

  FILE* file = fopen(filename, "rb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open " << filename << "!" << endl;
    exit(0);
  }

  int xRes, yRes, zRes;

  // counts - grid resolutions without padding
  fread((void*)&xRes, sizeof(int), 1, file);
  fread((void*)&yRes, sizeof(int), 1, file);
  fread((void*)&zRes, sizeof(int), 1, file);

  // domain
  float xMin, yMin, zMin;
  float xMax, yMax, zMax;
  fread((void*)&xMin, sizeof(float), 1, file);
  fread((void*)&yMin, sizeof(float), 1, file);
  fread((void*)&zMin, sizeof(float), 1, file);
  fread((void*)&xMax, sizeof(float), 1, file);
  fread((void*)&yMax, sizeof(float), 1, file);
  fread((void*)&zMax, sizeof(float), 1, file);

  // MAC offset
  float macOffset;
  fread((void*)&macOffset, sizeof(float), 1, file);

  // length2
  int length2;
  fread((void*)&length2, sizeof(int), 1, file);

  // domain (grid resolutions with padding)
  int xPaddedMin, xPaddedMax;
  int yPaddedMin, yPaddedMax;
  int zPaddedMin, zPaddedMax;
  fread((void*)&xPaddedMin, sizeof(int), 1, file);
  fread((void*)&xPaddedMax, sizeof(int), 1, file);
  fread((void*)&yPaddedMin, sizeof(int), 1, file);
  fread((void*)&yPaddedMax, sizeof(int), 1, file);
  fread((void*)&zPaddedMin, sizeof(int), 1, file);
  fread((void*)&zPaddedMax, sizeof(int), 1, file);

  _xRes = (xPaddedMax - xPaddedMin) + 1;
  _yRes = (yPaddedMax - yPaddedMin) + 1;
  _zRes = (zPaddedMax - zPaddedMin) + 1;

  int temp = _xRes;
  _xRes = _zRes;
  _zRes = temp;

  // set longest dimension to 1
  int biggest = (_xRes > _yRes) ? _xRes : _yRes;
  biggest = (_zRes > biggest) ? _zRes : biggest;
  _lengths[0] = (REAL)_xRes / biggest;
  _lengths[1] = (REAL)_yRes / biggest;
  _lengths[2] = (REAL)_zRes / biggest;

  cout << " lengths: " << _lengths << endl;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new REAL[_totalCells];

  if (sizeof(REAL) == sizeof(float))
    fread((void*)_data, sizeof(float), _totalCells, file);
  else
    for (int x = 0; x < _totalCells; x++)
    {
      float single;
      fread((void*)&single, sizeof(float), 1, file);
      _data[x] = single;
    }
  fclose(file);

  // find the biggest length
  REAL maxLength = (_lengths[1] > _lengths[0]) ? _lengths[1] : _lengths[0];
  maxLength = (_lengths[2] > maxLength) ? _lengths[2] : maxLength;

  (*this) *= 1.0 / maxLength;
}

///////////////////////////////////////////////////////////////////////
// set to a checkerboard solid texture
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setToSolidCheckboard(int xChecks, int yChecks, int zChecks)
{
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        int xMod = (x / (_xRes / xChecks)) % 2;
        int yMod = (y / (_yRes / yChecks)) % 2;
        int zMod = (z / (_zRes / zChecks)) % 2;

        if (((xMod && yMod) || (!xMod && !yMod)) && zMod)
          _data[index] = 1;
        
        if (!((xMod && yMod) || (!xMod && !yMod)) && !zMod)
          _data[index] = 1;
      }
}

///////////////////////////////////////////////////////////////////////
// set to a checkerboard solid texture
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setToGrayCheckerboard(int xChecks, int yChecks, int zChecks)
{
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        int xMod = (x / (_xRes / xChecks)) % 2;
        int yMod = (y / (_yRes / yChecks)) % 2;
        int zMod = (z / (_zRes / zChecks)) % 2;

        if (((xMod && yMod) || (!xMod && !yMod)) && zMod)
          _data[index] = 0.25;
        else if (!((xMod && yMod) || (!xMod && !yMod)) && !zMod)
          _data[index] = 0.25;
        else
          _data[index] = -0.25;
      }
}

///////////////////////////////////////////////////////////////////////
// triqunitic interpolation lookup
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::quinticLookup(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  VECTOR3 corner = _center - (REAL)0.5 * _lengths;
  VECTOR3 dxs(_dx, _dy, _dz);
  corner += (REAL)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const REAL xInterp = positionCopy[0]- x1;
  const REAL yInterp = positionCopy[1]- y1;
  const REAL zInterp = positionCopy[2]- z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const REAL p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const REAL p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const REAL p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const REAL p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const REAL p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const REAL p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const REAL p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const REAL p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const REAL p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const REAL p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const REAL p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const REAL p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const REAL p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const REAL p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const REAL p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const REAL p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const REAL z0Points[] = {quinticInterp(xInterp, p0), quinticInterp(xInterp, p1), quinticInterp(xInterp, p2), quinticInterp(xInterp, p3)};
  const REAL z1Points[] = {quinticInterp(xInterp, p4), quinticInterp(xInterp, p5), quinticInterp(xInterp, p6), quinticInterp(xInterp, p7)};
  const REAL z2Points[] = {quinticInterp(xInterp, p8), quinticInterp(xInterp, p9), quinticInterp(xInterp, p10), quinticInterp(xInterp, p11)};
  const REAL z3Points[] = {quinticInterp(xInterp, p12), quinticInterp(xInterp, p13), quinticInterp(xInterp, p14), quinticInterp(xInterp, p15)};

  const REAL finalPoints[] = {quinticInterp(yInterp, z0Points), quinticInterp(yInterp, z1Points), quinticInterp(yInterp, z2Points), quinticInterp(yInterp, z3Points)};

  return quinticInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// tricubic interpolation lookup
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::cubicLookup(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  const VECTOR3 corner = _center - (REAL)0.5 * _lengths + (REAL)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const REAL xInterp = positionCopy[0] - x1;
  const REAL yInterp = positionCopy[1] - y1;
  const REAL zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const REAL p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const REAL p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const REAL p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const REAL p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const REAL p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const REAL p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const REAL p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const REAL p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const REAL p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const REAL p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const REAL p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const REAL p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const REAL p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const REAL p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const REAL p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const REAL p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const REAL z0Points[] = {cubicInterp(xInterp, p0), cubicInterp(xInterp, p1), cubicInterp(xInterp, p2), cubicInterp(xInterp, p3)};
  const REAL z1Points[] = {cubicInterp(xInterp, p4), cubicInterp(xInterp, p5), cubicInterp(xInterp, p6), cubicInterp(xInterp, p7)};
  const REAL z2Points[] = {cubicInterp(xInterp, p8), cubicInterp(xInterp, p9), cubicInterp(xInterp, p10), cubicInterp(xInterp, p11)};
  const REAL z3Points[] = {cubicInterp(xInterp, p12), cubicInterp(xInterp, p13), cubicInterp(xInterp, p14), cubicInterp(xInterp, p15)};

  const REAL finalPoints[] = {cubicInterp(yInterp, z0Points), cubicInterp(yInterp, z1Points), cubicInterp(yInterp, z2Points), cubicInterp(yInterp, z3Points)};

  return cubicInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// a debugging version of linear interpolation -- allows
// a watch point to be inserted for a grid index
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::lerpDebug(const VECTOR3& position, int x, int y, int z) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  VECTOR3 corner = _center - (REAL)0.5 * _lengths;
  VECTOR3 dxs(_dx, _dy, _dz);
  corner += (REAL)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  // get interpolation weights
  const REAL s1 = positionCopy[0]- x0;
  const REAL s0 = 1.0f - s1;
  const REAL t1 = positionCopy[1]- y0;
  const REAL t0 = 1.0f - t1;
  const REAL u1 = positionCopy[2]- z0;
  const REAL u0 = 1.0f - u1;

  if ((x == 78 && y == 76 && z == 94) || (x == 73 && y == 76 && z == 94))
  {
    cout << " xyz: " << x << " " << y << " " << z << endl;
    cout << " s: " << s0 << " " << s1 << endl;
    cout << " t: " << t0 << " " << t1 << endl;
    cout << " u: " << u0 << " " << u1 << endl;
  }

  const int i000 = x0 + y0 * _xRes + z0 * _slabSize;
  const int i010 = x0 + y1 * _xRes + z0 * _slabSize;
  const int i100 = x1 + y0 * _xRes + z0 * _slabSize;
  const int i110 = x1 + y1 * _xRes + z0 * _slabSize;
  const int i001 = x0 + y0 * _xRes + z1 * _slabSize;
  const int i011 = x0 + y1 * _xRes + z1 * _slabSize;
  const int i101 = x1 + y0 * _xRes + z1 * _slabSize;
  const int i111 = x1 + y1 * _xRes + z1 * _slabSize;

  // interpolate
  // (indices could be computed once)
  return u0 * (s0 * (t0 * _data[i000] + t1 * _data[i010]) +
               s1 * (t0 * _data[i100] + t1 * _data[i110])) +
         u1 * (s0 * (t0 * _data[i001] + t1 * _data[i011]) +
               s1 * (t0 * _data[i101] + t1 * _data[i111]));
}

///////////////////////////////////////////////////////////////////////
// triquartic interpolation lookup
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::quarticLookup(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  const VECTOR3 corner = _center - (REAL)0.5 * _lengths + (REAL)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const REAL xInterp = positionCopy[0] - x1;
  const REAL yInterp = positionCopy[1] - y1;
  const REAL zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const REAL p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const REAL p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const REAL p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const REAL p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const REAL p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const REAL p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const REAL p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const REAL p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const REAL p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const REAL p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const REAL p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const REAL p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const REAL p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const REAL p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const REAL p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const REAL p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const REAL z0Points[] = {quarticInterp(xInterp, p0), quarticInterp(xInterp, p1), quarticInterp(xInterp, p2), quarticInterp(xInterp, p3)};
  const REAL z1Points[] = {quarticInterp(xInterp, p4), quarticInterp(xInterp, p5), quarticInterp(xInterp, p6), quarticInterp(xInterp, p7)};
  const REAL z2Points[] = {quarticInterp(xInterp, p8), quarticInterp(xInterp, p9), quarticInterp(xInterp, p10), quarticInterp(xInterp, p11)};
  const REAL z3Points[] = {quarticInterp(xInterp, p12), quarticInterp(xInterp, p13), quarticInterp(xInterp, p14), quarticInterp(xInterp, p15)};

  const REAL finalPoints[] = {quarticInterp(yInterp, z0Points), quarticInterp(yInterp, z1Points), quarticInterp(yInterp, z2Points), quarticInterp(yInterp, z3Points)};

  return quarticInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// clamp to the nearest neighbor
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::nearestNeighborLookup(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  const VECTOR3 corner = _center - (REAL)0.5 * _lengths + (REAL)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1 = x0 + 1;

  int y0 = (int)positionCopy[1];
  int y1 = y0 + 1;
  
  int z0 = (int)positionCopy[2];
  int z1 = z0 + 1;

  const REAL xInterp = positionCopy[0] - x0;
  const REAL yInterp = positionCopy[1] - y0;
  const REAL zInterp = positionCopy[2] - z0;

  x0 = (x0 < 0) ? 0 : x0;
  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  x1 = (x1 < 0) ? 0 : x1;
  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;

  y0 = (y0 < 0) ? 0 : y0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  y1 = (y1 < 0) ? 0 : y1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;

  z0 = (z0 < 0) ? 0 : z0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;
  z1 = (z1 < 0) ? 0 : z1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  int x = (xInterp < 0.5) ? x0 : x1;
  int y = (yInterp < 0.5) ? y0 : y1;
  int z = (zInterp < 0.5) ? z0 : z1;

  return (*this)(x,y,z);
}

#define QUARTIC(NAME, INDEX, COEFF, POINTS) \
  const REAL fim1##NAME##INDEX = POINTS[0];\
  const REAL fi##NAME##INDEX   = POINTS[1];\
  const REAL fip1##NAME##INDEX = POINTS[2];\
  const REAL fip2##NAME##INDEX = POINTS[3];\
  const REAL x##NAME##INDEX = COEFF; \
  const REAL p1##NAME##INDEX = fi##NAME##INDEX + ((fip1##NAME##INDEX - fim1##NAME##INDEX) + (fip1##NAME##INDEX - 2.0 * fi##NAME##INDEX + fim1##NAME##INDEX) * x##NAME##INDEX) * 0.5 * x##NAME##INDEX;\
  const REAL p2##NAME##INDEX = fi##NAME##INDEX + ((-fip2##NAME##INDEX + 4.0 * fip1##NAME##INDEX - 3 * fi##NAME##INDEX) + (fip2##NAME##INDEX - 2.0 * fip1##NAME##INDEX + fi##NAME##INDEX) * x##NAME##INDEX) * 0.5 * x##NAME##INDEX;\
  const REAL C1##NAME##INDEX = (2 - x##NAME##INDEX) / 3.0;\
  const REAL C2##NAME##INDEX = (x##NAME##INDEX + 1) / 3.0;\
  const REAL middle##NAME##INDEX = -76 * fip1##NAME##INDEX * fi##NAME##INDEX;\
  const REAL fip1Sq##NAME##INDEX = fip1##NAME##INDEX * fip1##NAME##INDEX;\
  const REAL fiSq##NAME##INDEX = fi##NAME##INDEX * fi##NAME##INDEX;\
  const REAL IS1##NAME##INDEX = (26 * fip1##NAME##INDEX * fim1##NAME##INDEX - 52 * fi##NAME##INDEX * fim1##NAME##INDEX + middle##NAME##INDEX + 25 * fip1Sq##NAME##INDEX + 64 * fiSq##NAME##INDEX + 13 * fim1##NAME##INDEX * fim1##NAME##INDEX) / 12.0 + 1e-6;\
  const REAL IS2##NAME##INDEX = (26 * fip2##NAME##INDEX * fi##NAME##INDEX - 52 * fip2##NAME##INDEX * fip1##NAME##INDEX + middle##NAME##INDEX + 25 * fiSq##NAME##INDEX + 64 * fip1Sq##NAME##INDEX + 13 * fip2##NAME##INDEX * fip2##NAME##INDEX) / 12.0 + 1e-6;\
  const REAL alpha1##NAME##INDEX = C1##NAME##INDEX / (IS1##NAME##INDEX * IS1##NAME##INDEX);\
  const REAL alpha2##NAME##INDEX = C2##NAME##INDEX / (IS2##NAME##INDEX * IS2##NAME##INDEX);\
  const REAL sum##NAME##INDEX = 1.0 / (alpha1##NAME##INDEX + alpha2##NAME##INDEX);\
  const REAL w1##NAME##INDEX = alpha1##NAME##INDEX * sum##NAME##INDEX;\
  const REAL w2##NAME##INDEX = alpha2##NAME##INDEX * sum##NAME##INDEX;\
  NAME[INDEX] = w1##NAME##INDEX * p1##NAME##INDEX + w2##NAME##INDEX * p2##NAME##INDEX;

///////////////////////////////////////////////////////////////////////
// triquartic interpolation lookup with aggressive inlining
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::quarticLookupInlined(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  const VECTOR3 corner = _center - (REAL)0.5 * _lengths + (REAL)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const REAL xInterp = positionCopy[0] - x1;
  const REAL yInterp = positionCopy[1] - y1;
  const REAL zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const REAL p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const REAL p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const REAL p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const REAL p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const REAL p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const REAL p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const REAL p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const REAL p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const REAL p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const REAL p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const REAL p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const REAL p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const REAL p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const REAL p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const REAL p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const REAL p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  //REAL z0Points[] = {quarticInterp(xInterp, p0), quarticInterp(xInterp, p1), quarticInterp(xInterp, p2), quarticInterp(xInterp, p3)};
  REAL z0Points[4];
  QUARTIC(z0Points, 0, xInterp, p0);
  QUARTIC(z0Points, 1, xInterp, p1);
  QUARTIC(z0Points, 2, xInterp, p2);
  QUARTIC(z0Points, 3, xInterp, p3);

  //const REAL z1Points[] = {quarticInterp(xInterp, p4), quarticInterp(xInterp, p5), quarticInterp(xInterp, p6), quarticInterp(xInterp, p7)};
  REAL z1Points[4];
  QUARTIC(z1Points, 0, xInterp, p4);
  QUARTIC(z1Points, 1, xInterp, p5);
  QUARTIC(z1Points, 2, xInterp, p6);
  QUARTIC(z1Points, 3, xInterp, p7);

  //const REAL z2Points[] = {quarticInterp(xInterp, p8), quarticInterp(xInterp, p9), quarticInterp(xInterp, p10), quarticInterp(xInterp, p11)};
  REAL z2Points[4];
  QUARTIC(z2Points, 0, xInterp, p8);
  QUARTIC(z2Points, 1, xInterp, p9);
  QUARTIC(z2Points, 2, xInterp, p10);
  QUARTIC(z2Points, 3, xInterp, p11);
 
  //const REAL z3Points[] = {quarticInterp(xInterp, p12), quarticInterp(xInterp, p13), quarticInterp(xInterp, p14), quarticInterp(xInterp, p15)};
  REAL z3Points[4];
  QUARTIC(z3Points, 0, xInterp, p12);
  QUARTIC(z3Points, 1, xInterp, p13);
  QUARTIC(z3Points, 2, xInterp, p14);
  QUARTIC(z3Points, 3, xInterp, p15);

  //const REAL finalPoints[] = {quarticInterp(yInterp, z0Points), quarticInterp(yInterp, z1Points), quarticInterp(yInterp, z2Points), quarticInterp(yInterp, z3Points)};
  REAL finalPoints[4];
  QUARTIC(finalPoints, 0, yInterp, z0Points);
  QUARTIC(finalPoints, 1, yInterp, z1Points);
  QUARTIC(finalPoints, 2, yInterp, z2Points);
  QUARTIC(finalPoints, 3, yInterp, z3Points);

  const REAL fim1 = finalPoints[0];
  const REAL fi   = finalPoints[1];
  const REAL fip1 = finalPoints[2];
  const REAL fip2 = finalPoints[3];
  const REAL x = zInterp;

  const REAL p1Final = fi + ((fip1 - fim1) + (fip1 - 2.0 * fi + fim1) * x) * 0.5 * x;
  const REAL p2Final = fi + ((-fip2 + 4.0 * fip1 - 3 * fi) + (fip2 - 2.0 * fip1 + fi) * x) * 0.5 * x;

  const REAL C1 = (2 - x) * (1.0 / 3.0);
  const REAL C2 = (x + 1) * (1.0 / 3.0);

  const REAL middle = -76 * fip1 * fi;
  const REAL fip1Sq = fip1 * fip1;
  const REAL fiSq = fi * fi;

  const REAL IS1 = (26 * fip1 * fim1 - 52 * fi * fim1 + middle + 25 * fip1Sq + 64 * fiSq + 13 * fim1 * fim1) / 12.0 + 1e-6;
  const REAL IS2 = (26 * fip2 * fi - 52 * fip2 * fip1 + middle + 25 * fiSq + 64 * fip1Sq + 13 * fip2 * fip2) / 12.0 + 1e-6;

  const REAL alpha1 = C1 / (IS1 * IS1);
  const REAL alpha2 = C2 / (IS2 * IS2);

  const REAL sum = 1.0 / (alpha1 + alpha2);
  const REAL w1 = alpha1 * sum;
  const REAL w2 = alpha2 * sum;

  const REAL result = w1 * p1Final + w2 * p2Final;

  return result;
}

///////////////////////////////////////////////////////////////////////
// triquartic interpolation lookup, with clamping
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::quarticLookupClamped(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  const VECTOR3 corner = _center - (REAL)0.5 * _lengths + (REAL)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const REAL xInterp = positionCopy[0] - x1;
  const REAL yInterp = positionCopy[1] - y1;
  const REAL zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const REAL p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const REAL p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const REAL p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const REAL p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const REAL p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const REAL p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const REAL p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const REAL p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const REAL p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const REAL p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const REAL p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const REAL p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const REAL p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const REAL p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const REAL p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const REAL p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const REAL z0Points[] = {quarticInterpClamped(xInterp, p0), quarticInterpClamped(xInterp, p1), quarticInterpClamped(xInterp, p2), quarticInterpClamped(xInterp, p3)};
  const REAL z1Points[] = {quarticInterpClamped(xInterp, p4), quarticInterpClamped(xInterp, p5), quarticInterpClamped(xInterp, p6), quarticInterpClamped(xInterp, p7)};
  const REAL z2Points[] = {quarticInterpClamped(xInterp, p8), quarticInterpClamped(xInterp, p9), quarticInterpClamped(xInterp, p10), quarticInterpClamped(xInterp, p11)};
  const REAL z3Points[] = {quarticInterpClamped(xInterp, p12), quarticInterpClamped(xInterp, p13), quarticInterpClamped(xInterp, p14), quarticInterpClamped(xInterp, p15)};

  const REAL finalPoints[] = {quarticInterpClamped(yInterp, z0Points), quarticInterpClamped(yInterp, z1Points), quarticInterpClamped(yInterp, z2Points), quarticInterpClamped(yInterp, z3Points)};

  return quarticInterpClamped(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// trisextic interpolation lookup
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::sexticLookup(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  const VECTOR3 corner = _center - (REAL)0.5 * _lengths + (REAL)0.5 * dxs();

  // recenter position
  positionCopy -= corner;
  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x2 = (int)positionCopy[0];
  const int x5    = x2 + 3;
  const int x4    = x2 + 2;
  const int x3    = x2 + 1;
  const int x1    = x2 - 1;
  const int x0    = x2 - 2;

  const int y2 = (int)positionCopy[1];
  const int y5    = y2 + 3;
  const int y4    = y2 + 2;
  const int y3    = y2 + 1;
  const int y1    = y2 - 1;
  const int y0    = y2 - 2;

  const int z2 = (int)positionCopy[2];
  const int z5    = z2 + 3;
  const int z4    = z2 + 2;
  const int z3    = z2 + 1;
  const int z1    = z2 - 1;
  const int z0    = z2 - 2;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x5 >= _xRes || y5 >= _yRes || z5 >= _zRes)
    return (*this)(position);

  const REAL xInterp = positionCopy[0] - x1;
  const REAL yInterp = positionCopy[1] - y1;
  const REAL zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  const int z4Slab = z4 * _slabSize;
  const int z5Slab = z5 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;
  const int y4x = y4 * _xRes;
  const int y5x = y5 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;
  const int y4z0 = y4x + z0Slab;
  const int y5z0 = y5x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;
  const int y4z1 = y4x + z1Slab;
  const int y5z1 = y5x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;
  const int y4z2 = y4x + z2Slab;
  const int y5z2 = y5x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;
  const int y4z3 = y4x + z3Slab;
  const int y5z3 = y5x + z3Slab;

  const int y0z4 = y0x + z4Slab;
  const int y1z4 = y1x + z4Slab;
  const int y2z4 = y2x + z4Slab;
  const int y3z4 = y3x + z4Slab;
  const int y4z4 = y4x + z4Slab;
  const int y5z4 = y5x + z4Slab;

  const int y0z5 = y0x + z5Slab;
  const int y1z5 = y1x + z5Slab;
  const int y2z5 = y2x + z5Slab;
  const int y3z5 = y3x + z5Slab;
  const int y4z5 = y4x + z5Slab;
  const int y5z5 = y5x + z5Slab;

  // do the z0 slice
  const REAL p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0], _data[x4 + y0z0], _data[x5 + y0z0]};
  const REAL p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0], _data[x4 + y1z0], _data[x5 + y1z0]};
  const REAL p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0], _data[x4 + y2z0], _data[x5 + y2z0]};
  const REAL p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0], _data[x4 + y3z0], _data[x5 + y3z0]};
  const REAL p4[] = {_data[x0 + y4z0], _data[x1 + y4z0], _data[x2 + y4z0], _data[x3 + y4z0], _data[x4 + y4z0], _data[x5 + y4z0]};
  const REAL p5[] = {_data[x0 + y5z0], _data[x1 + y5z0], _data[x2 + y5z0], _data[x3 + y5z0], _data[x4 + y5z0], _data[x5 + y5z0]};

  // do the z1 slice
  const REAL p6[]  = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1], _data[x4 + y0z1], _data[x5 + y0z1]};
  const REAL p7[]  = {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1], _data[x4 + y1z1], _data[x5 + y1z1]};
  const REAL p8[]  = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1], _data[x4 + y2z1], _data[x5 + y2z1]};
  const REAL p9[]  = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1], _data[x4 + y3z1], _data[x5 + y3z1]};
  const REAL p10[] = {_data[x0 + y4z1], _data[x1 + y4z1], _data[x2 + y4z1], _data[x3 + y4z1], _data[x4 + y4z1], _data[x5 + y4z1]};
  const REAL p11[] = {_data[x0 + y5z1], _data[x1 + y5z1], _data[x2 + y5z1], _data[x3 + y5z1], _data[x4 + y5z1], _data[x5 + y5z1]};

  // do the z2 slice
  const REAL p12[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2], _data[x4 + y0z2], _data[x5 + y0z2]};
  const REAL p13[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2], _data[x4 + y1z2], _data[x5 + y1z2]};
  const REAL p14[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2], _data[x4 + y2z2], _data[x5 + y2z2]};
  const REAL p15[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2], _data[x4 + y3z2], _data[x5 + y3z2]};
  const REAL p16[] = {_data[x0 + y4z2], _data[x1 + y4z2], _data[x2 + y4z2], _data[x3 + y4z2], _data[x4 + y4z2], _data[x5 + y4z2]};
  const REAL p17[] = {_data[x0 + y5z2], _data[x1 + y5z2], _data[x2 + y5z2], _data[x3 + y5z2], _data[x4 + y5z2], _data[x5 + y5z2]};

  // do the z3 slice
  const REAL p18[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3], _data[x4 + y0z3], _data[x5 + y0z3]};
  const REAL p19[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3], _data[x4 + y1z3], _data[x5 + y1z3]};
  const REAL p20[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3], _data[x4 + y2z3], _data[x5 + y2z3]};
  const REAL p21[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3], _data[x4 + y3z3], _data[x5 + y3z3]};
  const REAL p22[] = {_data[x0 + y4z3], _data[x1 + y4z3], _data[x2 + y4z3], _data[x3 + y4z3], _data[x4 + y4z3], _data[x5 + y4z3]};
  const REAL p23[] = {_data[x0 + y5z3], _data[x1 + y5z3], _data[x2 + y5z3], _data[x3 + y5z3], _data[x4 + y5z3], _data[x5 + y5z3]};

  // do the z4 slice
  const REAL p24[] = {_data[x0 + y0z4], _data[x1 + y0z4], _data[x2 + y0z4], _data[x3 + y0z4], _data[x4 + y0z4], _data[x5 + y0z4]};
  const REAL p25[] = {_data[x0 + y1z4], _data[x1 + y1z4], _data[x2 + y1z4], _data[x3 + y1z4], _data[x4 + y1z4], _data[x5 + y1z4]};
  const REAL p26[] = {_data[x0 + y2z4], _data[x1 + y2z4], _data[x2 + y2z4], _data[x3 + y2z4], _data[x4 + y2z4], _data[x5 + y2z4]};
  const REAL p27[] = {_data[x0 + y3z4], _data[x1 + y3z4], _data[x2 + y3z4], _data[x3 + y3z4], _data[x4 + y3z4], _data[x5 + y3z4]};
  const REAL p28[] = {_data[x0 + y4z4], _data[x1 + y4z4], _data[x2 + y4z4], _data[x3 + y4z4], _data[x4 + y4z4], _data[x5 + y4z4]};
  const REAL p29[] = {_data[x0 + y5z4], _data[x1 + y5z4], _data[x2 + y5z4], _data[x3 + y5z4], _data[x4 + y5z4], _data[x5 + y5z4]};

  // do the z5 slice
  const REAL p30[] = {_data[x0 + y0z5], _data[x1 + y0z5], _data[x2 + y0z5], _data[x3 + y0z5], _data[x4 + y0z5], _data[x5 + y0z5]};
  const REAL p31[] = {_data[x0 + y1z5], _data[x1 + y1z5], _data[x2 + y1z5], _data[x3 + y1z5], _data[x4 + y1z5], _data[x5 + y1z5]};
  const REAL p32[] = {_data[x0 + y2z5], _data[x1 + y2z5], _data[x2 + y2z5], _data[x3 + y2z5], _data[x4 + y2z5], _data[x5 + y2z5]};
  const REAL p33[] = {_data[x0 + y3z5], _data[x1 + y3z5], _data[x2 + y3z5], _data[x3 + y3z5], _data[x4 + y3z5], _data[x5 + y3z5]};
  const REAL p34[] = {_data[x0 + y4z5], _data[x1 + y4z5], _data[x2 + y4z5], _data[x3 + y4z5], _data[x4 + y4z5], _data[x5 + y4z5]};
  const REAL p35[] = {_data[x0 + y5z5], _data[x1 + y5z5], _data[x2 + y5z5], _data[x3 + y5z5], _data[x4 + y5z5], _data[x5 + y5z5]};

  const REAL z0Points[] = {sexticInterp(xInterp, p0), 
                           sexticInterp(xInterp, p1), 
                           sexticInterp(xInterp, p2), 
                           sexticInterp(xInterp, p3),
                           sexticInterp(xInterp, p4),
                           sexticInterp(xInterp, p5)};

  const REAL z1Points[] = {sexticInterp(xInterp, p6), 
                           sexticInterp(xInterp, p7), 
                           sexticInterp(xInterp, p8), 
                           sexticInterp(xInterp, p9),
                           sexticInterp(xInterp, p10),
                           sexticInterp(xInterp, p11)};

  const REAL z2Points[] = {sexticInterp(xInterp, p12), 
                           sexticInterp(xInterp, p13), 
                           sexticInterp(xInterp, p14), 
                           sexticInterp(xInterp, p15),
                           sexticInterp(xInterp, p16),
                           sexticInterp(xInterp, p17)};

  const REAL z3Points[] = {sexticInterp(xInterp, p18), 
                           sexticInterp(xInterp, p19), 
                           sexticInterp(xInterp, p20), 
                           sexticInterp(xInterp, p21),
                           sexticInterp(xInterp, p22),
                           sexticInterp(xInterp, p23)};

  const REAL z4Points[] = {sexticInterp(xInterp, p24), 
                           sexticInterp(xInterp, p25), 
                           sexticInterp(xInterp, p26), 
                           sexticInterp(xInterp, p27),
                           sexticInterp(xInterp, p28),
                           sexticInterp(xInterp, p29)};

  const REAL z5Points[] = {sexticInterp(xInterp, p30), 
                           sexticInterp(xInterp, p31), 
                           sexticInterp(xInterp, p32), 
                           sexticInterp(xInterp, p33),
                           sexticInterp(xInterp, p34),
                           sexticInterp(xInterp, p35)};


  const REAL finalPoints[] = {sexticInterp(yInterp, z0Points), 
                              sexticInterp(yInterp, z1Points), 
                              sexticInterp(yInterp, z2Points), 
                              sexticInterp(yInterp, z3Points),
                              sexticInterp(yInterp, z4Points),
                              sexticInterp(yInterp, z5Points)};

  return sexticInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// trisextic interpolation lookup, with clamping
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::sexticLookupClamped(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  const VECTOR3 corner = _center - (REAL)0.5 * _lengths + (REAL)0.5 * dxs();

  // recenter position
  positionCopy -= corner;
  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x2 = (int)positionCopy[0];
  const int x5    = x2 + 3;
  const int x4    = x2 + 2;
  const int x3    = x2 + 1;
  const int x1    = x2 - 1;
  const int x0    = x2 - 2;

  const int y2 = (int)positionCopy[1];
  const int y5    = y2 + 3;
  const int y4    = y2 + 2;
  const int y3    = y2 + 1;
  const int y1    = y2 - 1;
  const int y0    = y2 - 2;

  const int z2 = (int)positionCopy[2];
  const int z5    = z2 + 3;
  const int z4    = z2 + 2;
  const int z3    = z2 + 1;
  const int z1    = z2 - 1;
  const int z0    = z2 - 2;


  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x5 >= _xRes || y5 >= _yRes || z5 >= _zRes)
    return (*this)(position);

  const REAL xInterp = positionCopy[0] - x1;
  const REAL yInterp = positionCopy[1] - y1;
  const REAL zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  const int z4Slab = z4 * _slabSize;
  const int z5Slab = z5 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;
  const int y4x = y4 * _xRes;
  const int y5x = y5 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;
  const int y4z0 = y4x + z0Slab;
  const int y5z0 = y5x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;
  const int y4z1 = y4x + z1Slab;
  const int y5z1 = y5x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;
  const int y4z2 = y4x + z2Slab;
  const int y5z2 = y5x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;
  const int y4z3 = y4x + z3Slab;
  const int y5z3 = y5x + z3Slab;

  const int y0z4 = y0x + z4Slab;
  const int y1z4 = y1x + z4Slab;
  const int y2z4 = y2x + z4Slab;
  const int y3z4 = y3x + z4Slab;
  const int y4z4 = y4x + z4Slab;
  const int y5z4 = y5x + z4Slab;

  const int y0z5 = y0x + z5Slab;
  const int y1z5 = y1x + z5Slab;
  const int y2z5 = y2x + z5Slab;
  const int y3z5 = y3x + z5Slab;
  const int y4z5 = y4x + z5Slab;
  const int y5z5 = y5x + z5Slab;

  // do the z0 slice
  const REAL p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0], _data[x4 + y0z0], _data[x5 + y0z0]};
  const REAL p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0], _data[x4 + y1z0], _data[x5 + y1z0]};
  const REAL p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0], _data[x4 + y2z0], _data[x5 + y2z0]};
  const REAL p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0], _data[x4 + y3z0], _data[x5 + y3z0]};
  const REAL p4[] = {_data[x0 + y4z0], _data[x1 + y4z0], _data[x2 + y4z0], _data[x3 + y4z0], _data[x4 + y4z0], _data[x5 + y4z0]};
  const REAL p5[] = {_data[x0 + y5z0], _data[x1 + y5z0], _data[x2 + y5z0], _data[x3 + y5z0], _data[x4 + y5z0], _data[x5 + y5z0]};

  // do the z1 slice
  const REAL p6[]  = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1], _data[x4 + y0z1], _data[x5 + y0z1]};
  const REAL p7[]  = {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1], _data[x4 + y1z1], _data[x5 + y1z1]};
  const REAL p8[]  = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1], _data[x4 + y2z1], _data[x5 + y2z1]};
  const REAL p9[]  = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1], _data[x4 + y3z1], _data[x5 + y3z1]};
  const REAL p10[] = {_data[x0 + y4z1], _data[x1 + y4z1], _data[x2 + y4z1], _data[x3 + y4z1], _data[x4 + y4z1], _data[x5 + y4z1]};
  const REAL p11[] = {_data[x0 + y5z1], _data[x1 + y5z1], _data[x2 + y5z1], _data[x3 + y5z1], _data[x4 + y5z1], _data[x5 + y5z1]};

  // do the z2 slice
  const REAL p12[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2], _data[x4 + y0z2], _data[x5 + y0z2]};
  const REAL p13[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2], _data[x4 + y1z2], _data[x5 + y1z2]};
  const REAL p14[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2], _data[x4 + y2z2], _data[x5 + y2z2]};
  const REAL p15[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2], _data[x4 + y3z2], _data[x5 + y3z2]};
  const REAL p16[] = {_data[x0 + y4z2], _data[x1 + y4z2], _data[x2 + y4z2], _data[x3 + y4z2], _data[x4 + y4z2], _data[x5 + y4z2]};
  const REAL p17[] = {_data[x0 + y5z2], _data[x1 + y5z2], _data[x2 + y5z2], _data[x3 + y5z2], _data[x4 + y5z2], _data[x5 + y5z2]};

  // do the z3 slice
  const REAL p18[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3], _data[x4 + y0z3], _data[x5 + y0z3]};
  const REAL p19[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3], _data[x4 + y1z3], _data[x5 + y1z3]};
  const REAL p20[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3], _data[x4 + y2z3], _data[x5 + y2z3]};
  const REAL p21[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3], _data[x4 + y3z3], _data[x5 + y3z3]};
  const REAL p22[] = {_data[x0 + y4z3], _data[x1 + y4z3], _data[x2 + y4z3], _data[x3 + y4z3], _data[x4 + y4z3], _data[x5 + y4z3]};
  const REAL p23[] = {_data[x0 + y5z3], _data[x1 + y5z3], _data[x2 + y5z3], _data[x3 + y5z3], _data[x4 + y5z3], _data[x5 + y5z3]};

  // do the z4 slice
  const REAL p24[] = {_data[x0 + y0z4], _data[x1 + y0z4], _data[x2 + y0z4], _data[x3 + y0z4], _data[x4 + y0z4], _data[x5 + y0z4]};
  const REAL p25[] = {_data[x0 + y1z4], _data[x1 + y1z4], _data[x2 + y1z4], _data[x3 + y1z4], _data[x4 + y1z4], _data[x5 + y1z4]};
  const REAL p26[] = {_data[x0 + y2z4], _data[x1 + y2z4], _data[x2 + y2z4], _data[x3 + y2z4], _data[x4 + y2z4], _data[x5 + y2z4]};
  const REAL p27[] = {_data[x0 + y3z4], _data[x1 + y3z4], _data[x2 + y3z4], _data[x3 + y3z4], _data[x4 + y3z4], _data[x5 + y3z4]};
  const REAL p28[] = {_data[x0 + y4z4], _data[x1 + y4z4], _data[x2 + y4z4], _data[x3 + y4z4], _data[x4 + y4z4], _data[x5 + y4z4]};
  const REAL p29[] = {_data[x0 + y5z4], _data[x1 + y5z4], _data[x2 + y5z4], _data[x3 + y5z4], _data[x4 + y5z4], _data[x5 + y5z4]};

  // do the z5 slice
  const REAL p30[] = {_data[x0 + y0z5], _data[x1 + y0z5], _data[x2 + y0z5], _data[x3 + y0z5], _data[x4 + y0z5], _data[x5 + y0z5]};
  const REAL p31[] = {_data[x0 + y1z5], _data[x1 + y1z5], _data[x2 + y1z5], _data[x3 + y1z5], _data[x4 + y1z5], _data[x5 + y1z5]};
  const REAL p32[] = {_data[x0 + y2z5], _data[x1 + y2z5], _data[x2 + y2z5], _data[x3 + y2z5], _data[x4 + y2z5], _data[x5 + y2z5]};
  const REAL p33[] = {_data[x0 + y3z5], _data[x1 + y3z5], _data[x2 + y3z5], _data[x3 + y3z5], _data[x4 + y3z5], _data[x5 + y3z5]};
  const REAL p34[] = {_data[x0 + y4z5], _data[x1 + y4z5], _data[x2 + y4z5], _data[x3 + y4z5], _data[x4 + y4z5], _data[x5 + y4z5]};
  const REAL p35[] = {_data[x0 + y5z5], _data[x1 + y5z5], _data[x2 + y5z5], _data[x3 + y5z5], _data[x4 + y5z5], _data[x5 + y5z5]};

  const REAL z0Points[] = {sexticInterpClamped(xInterp, p0), 
                           sexticInterpClamped(xInterp, p1), 
                           sexticInterpClamped(xInterp, p2), 
                           sexticInterpClamped(xInterp, p3),
                           sexticInterpClamped(xInterp, p4),
                           sexticInterpClamped(xInterp, p5)};

  const REAL z1Points[] = {sexticInterpClamped(xInterp, p6), 
                           sexticInterpClamped(xInterp, p7), 
                           sexticInterpClamped(xInterp, p8), 
                           sexticInterpClamped(xInterp, p9),
                           sexticInterpClamped(xInterp, p10),
                           sexticInterpClamped(xInterp, p11)};

  const REAL z2Points[] = {sexticInterpClamped(xInterp, p12), 
                           sexticInterpClamped(xInterp, p13), 
                           sexticInterpClamped(xInterp, p14), 
                           sexticInterpClamped(xInterp, p15),
                           sexticInterpClamped(xInterp, p16),
                           sexticInterpClamped(xInterp, p17)};

  const REAL z3Points[] = {sexticInterpClamped(xInterp, p18), 
                           sexticInterpClamped(xInterp, p19), 
                           sexticInterpClamped(xInterp, p20), 
                           sexticInterpClamped(xInterp, p21),
                           sexticInterpClamped(xInterp, p22),
                           sexticInterpClamped(xInterp, p23)};

  const REAL z4Points[] = {sexticInterpClamped(xInterp, p24), 
                           sexticInterpClamped(xInterp, p25), 
                           sexticInterpClamped(xInterp, p26), 
                           sexticInterpClamped(xInterp, p27),
                           sexticInterpClamped(xInterp, p28),
                           sexticInterpClamped(xInterp, p29)};

  const REAL z5Points[] = {sexticInterpClamped(xInterp, p30), 
                           sexticInterpClamped(xInterp, p31), 
                           sexticInterpClamped(xInterp, p32), 
                           sexticInterpClamped(xInterp, p33),
                           sexticInterpClamped(xInterp, p34),
                           sexticInterpClamped(xInterp, p35)};


  const REAL finalPoints[] = {sexticInterpClamped(yInterp, z0Points), 
                              sexticInterpClamped(yInterp, z1Points), 
                              sexticInterpClamped(yInterp, z2Points), 
                              sexticInterpClamped(yInterp, z3Points),
                              sexticInterpClamped(yInterp, z4Points),
                              sexticInterpClamped(yInterp, z5Points)};

  return sexticInterpClamped(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// tricubic interpolation lookup
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::cubicLookupUnclamped(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  const VECTOR3 corner = _center - (REAL)0.5 * _lengths + (REAL)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const REAL xInterp = positionCopy[0] - x1;
  const REAL yInterp = positionCopy[1] - y1;
  const REAL zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const REAL p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const REAL p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const REAL p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const REAL p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const REAL p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const REAL p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const REAL p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const REAL p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const REAL p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const REAL p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const REAL p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const REAL p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const REAL p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const REAL p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const REAL p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const REAL p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const REAL z0Points[] = {cubicInterpUnclamped(xInterp, p0), cubicInterp(xInterp, p1), cubicInterp(xInterp, p2), cubicInterp(xInterp, p3)};
  const REAL z1Points[] = {cubicInterpUnclamped(xInterp, p4), cubicInterp(xInterp, p5), cubicInterp(xInterp, p6), cubicInterp(xInterp, p7)};
  const REAL z2Points[] = {cubicInterpUnclamped(xInterp, p8), cubicInterp(xInterp, p9), cubicInterp(xInterp, p10), cubicInterp(xInterp, p11)};
  const REAL z3Points[] = {cubicInterpUnclamped(xInterp, p12), cubicInterp(xInterp, p13), cubicInterp(xInterp, p14), cubicInterp(xInterp, p15)};

  const REAL finalPoints[] = {cubicInterpUnclamped(yInterp, z0Points), cubicInterp(yInterp, z1Points), cubicInterp(yInterp, z2Points), cubicInterp(yInterp, z3Points)};

  return cubicInterpUnclamped(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// tricubic interpolation lookup
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::cubicNewtonLookup(const VECTOR3& position) const
{
  VECTOR3 positionCopy = position;

  // get the lower corner position
  VECTOR3 corner = _center - (REAL)0.5 * _lengths;
  VECTOR3 dxs(_dx, _dy, _dz);
  corner += (REAL)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x1 = (int)positionCopy[0];
  int x2    = x1 + 1;
  int x3    = x1 + 2;
  int x0    = x1 - 1;

  int y1 = (int)positionCopy[1];
  int y2    = y1 + 1;
  int y3    = y1 + 2;
  int y0    = y1 - 1;
  
  int z1 = (int)positionCopy[2];
  int z2    = z1 + 1;
  int z3    = z1 + 2;
  int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  REAL points[4];
  const REAL xInterp = positionCopy[0]- x1;
  const REAL yInterp = positionCopy[1]- y1;
  const REAL zInterp = positionCopy[2]- z1;

  REAL finalPoints[4];

  // do the z0 slice
  REAL z0Points[4];
  points[0] = _data[x0 + y0 * _xRes + z0 * _slabSize];
  points[1] = _data[x1 + y0 * _xRes + z0 * _slabSize];
  points[2] = _data[x2 + y0 * _xRes + z0 * _slabSize];
  points[3] = _data[x3 + y0 * _xRes + z0 * _slabSize];
  z0Points[0] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y1 * _xRes + z0 * _slabSize];
  points[1] = _data[x1 + y1 * _xRes + z0 * _slabSize];
  points[2] = _data[x2 + y1 * _xRes + z0 * _slabSize];
  points[3] = _data[x3 + y1 * _xRes + z0 * _slabSize];
  z0Points[1] = cubicNewtonInterp(xInterp, points);
  
  points[0] = _data[x0 + y2 * _xRes + z0 * _slabSize];
  points[1] = _data[x1 + y2 * _xRes + z0 * _slabSize];
  points[2] = _data[x2 + y2 * _xRes + z0 * _slabSize];
  points[3] = _data[x3 + y2 * _xRes + z0 * _slabSize];
  z0Points[2] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y3 * _xRes + z0 * _slabSize];
  points[1] = _data[x1 + y3 * _xRes + z0 * _slabSize];
  points[2] = _data[x2 + y3 * _xRes + z0 * _slabSize];
  points[3] = _data[x3 + y3 * _xRes + z0 * _slabSize];
  z0Points[3] = cubicNewtonInterp(xInterp, points);
  finalPoints[0] = cubicNewtonInterp(yInterp, z0Points);

  // do the z1 slice
  REAL z1Points[4];
  points[0] = _data[x0 + y0 * _xRes + z1 * _slabSize];
  points[1] = _data[x1 + y0 * _xRes + z1 * _slabSize];
  points[2] = _data[x2 + y0 * _xRes + z1 * _slabSize];
  points[3] = _data[x3 + y0 * _xRes + z1 * _slabSize];
  z1Points[0] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y1 * _xRes + z1 * _slabSize];
  points[1] = _data[x1 + y1 * _xRes + z1 * _slabSize];
  points[2] = _data[x2 + y1 * _xRes + z1 * _slabSize];
  points[3] = _data[x3 + y1 * _xRes + z1 * _slabSize];
  z1Points[1] = cubicNewtonInterp(xInterp, points);
  
  points[0] = _data[x0 + y2 * _xRes + z1 * _slabSize];
  points[1] = _data[x1 + y2 * _xRes + z1 * _slabSize];
  points[2] = _data[x2 + y2 * _xRes + z1 * _slabSize];
  points[3] = _data[x3 + y2 * _xRes + z1 * _slabSize];
  z1Points[2] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y3 * _xRes + z1 * _slabSize];
  points[1] = _data[x1 + y3 * _xRes + z1 * _slabSize];
  points[2] = _data[x2 + y3 * _xRes + z1 * _slabSize];
  points[3] = _data[x3 + y3 * _xRes + z1 * _slabSize];
  z1Points[3] = cubicNewtonInterp(xInterp, points);
  finalPoints[1] = cubicNewtonInterp(yInterp, z1Points);

  // do the z2 slice
  REAL z2Points[4];
  points[0] = _data[x0 + y0 * _xRes + z2 * _slabSize];
  points[1] = _data[x1 + y0 * _xRes + z2 * _slabSize];
  points[2] = _data[x2 + y0 * _xRes + z2 * _slabSize];
  points[3] = _data[x3 + y0 * _xRes + z2 * _slabSize];
  z2Points[0] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y1 * _xRes + z2 * _slabSize];
  points[1] = _data[x1 + y1 * _xRes + z2 * _slabSize];
  points[2] = _data[x2 + y1 * _xRes + z2 * _slabSize];
  points[3] = _data[x3 + y1 * _xRes + z2 * _slabSize];
  z2Points[1] = cubicNewtonInterp(xInterp, points);
  
  points[0] = _data[x0 + y2 * _xRes + z2 * _slabSize];
  points[1] = _data[x1 + y2 * _xRes + z2 * _slabSize];
  points[2] = _data[x2 + y2 * _xRes + z2 * _slabSize];
  points[3] = _data[x3 + y2 * _xRes + z2 * _slabSize];
  z2Points[2] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y3 * _xRes + z2 * _slabSize];
  points[1] = _data[x1 + y3 * _xRes + z2 * _slabSize];
  points[2] = _data[x2 + y3 * _xRes + z2 * _slabSize];
  points[3] = _data[x3 + y3 * _xRes + z2 * _slabSize];
  z2Points[3] = cubicNewtonInterp(xInterp, points);
  finalPoints[2] = cubicNewtonInterp(yInterp, z2Points);

  // do the z3 slice
  REAL z3Points[4];
  points[0] = _data[x0 + y0 * _xRes + z3 * _slabSize];
  points[1] = _data[x1 + y0 * _xRes + z3 * _slabSize];
  points[2] = _data[x2 + y0 * _xRes + z3 * _slabSize];
  points[3] = _data[x3 + y0 * _xRes + z3 * _slabSize];
  z3Points[0] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y1 * _xRes + z3 * _slabSize];
  points[1] = _data[x1 + y1 * _xRes + z3 * _slabSize];
  points[2] = _data[x2 + y1 * _xRes + z3 * _slabSize];
  points[3] = _data[x3 + y1 * _xRes + z3 * _slabSize];
  z3Points[1] = cubicNewtonInterp(xInterp, points);
  
  points[0] = _data[x0 + y2 * _xRes + z3 * _slabSize];
  points[1] = _data[x1 + y2 * _xRes + z3 * _slabSize];
  points[2] = _data[x2 + y2 * _xRes + z3 * _slabSize];
  points[3] = _data[x3 + y2 * _xRes + z3 * _slabSize];
  z3Points[2] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y3 * _xRes + z3 * _slabSize];
  points[1] = _data[x1 + y3 * _xRes + z3 * _slabSize];
  points[2] = _data[x2 + y3 * _xRes + z3 * _slabSize];
  points[3] = _data[x3 + y3 * _xRes + z3 * _slabSize];
  z3Points[3] = cubicNewtonInterp(xInterp, points);
  finalPoints[3] = cubicNewtonInterp(yInterp, z3Points);

  return cubicNewtonInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite that clamps to the immediate neighborhood
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::cubicInterpClamped(const REAL interp, const REAL* points)
{
  REAL d0 = (points[2] - points[0]) * 0.5;
  REAL d1 = (points[3] - points[1]) * 0.5;

  const REAL deltak = (points[2] - points[1]);

  // do monotonic interpolation
  if (deltak * d0 < 0.0)
    d0 = 0;
  if (deltak * d1 < 0.0)
    d1 = 0;

  const REAL a0 = points[1];
  const REAL a1 = d0;
  const REAL a2 = 3.0 * deltak - 2.0 * d0 - d1;
  const REAL a3 = -2.0 * deltak + d0 + d1;

  const REAL squared = interp * interp;
  const REAL cubed = squared * interp;
  REAL trial = a3 * cubed + a2 * squared + a1 * interp + a0;

  const REAL intervalMax = (points[1] > points[2]) ? points[1] : points[2];
  const REAL intervalMin = (points[1] > points[2]) ? points[2] : points[1];

  trial = (trial > intervalMax) ? intervalMax : trial;
  trial = (trial < intervalMin) ? intervalMin : trial;

  return trial;
}

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite interpolation
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::cubicInterp(const REAL interp, const REAL* points)
{
  REAL d0 = (points[2] - points[0]) * 0.5;
  REAL d1 = (points[3] - points[1]) * 0.5;

  REAL deltak = (points[2] - points[1]);

  // do monotonic interpolation
  if (deltak * d0 < 0.0)
    d0 = 0;
  if (deltak * d1 < 0.0)
    d1 = 0;

  REAL a0 = points[1];
  REAL a1 = d0;
  REAL a2 = 3.0 * deltak - 2.0 * d0 - d1;
  REAL a3 = -2.0 * deltak + d0 + d1;

  REAL squared = interp * interp;
  REAL cubed = squared * interp;
  return a3 * cubed + a2 * squared + a1 * interp + a0;
}

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite interpolation
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::cubicInterpUnclamped(const REAL interp, const REAL* points)
{
  REAL d0 = (points[2] - points[0]) * 0.5;
  REAL d1 = (points[3] - points[1]) * 0.5;

  REAL deltak = (points[2] - points[1]);

  REAL a0 = points[1];
  REAL a1 = d0;
  REAL a2 = 3.0 * deltak - 2.0 * d0 - d1;
  REAL a3 = -2.0 * deltak + d0 + d1;

  REAL squared = interp * interp;
  REAL cubed = squared * interp;
  return a3 * cubed + a2 * squared + a1 * interp + a0;
}

///////////////////////////////////////////////////////////////////////
// do a quintic Hermite interpolation
//
// The constraint matrix is:
// A = [ 0 0 0 0 0 1;0 0 0 0 1 0;0 0 0 2 0 0 ; 1 1 1 1 1 1 ;5 4 3 2 1 0; 20 12 6 2 0 0 ]
//
// where b = [y0 dy0 ddy0 d1 dy1 ddy1]';
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::quinticInterp(REAL interp, const REAL* points)
{
  // do monotone increasing
  if (points[1] > points[2])
  {
    const REAL newPts[] = {points[3], points[2], points[1], points[0]};
    return quinticInterp(1 - interp, newPts); 
  }

  const REAL dy0 = (points[2] - points[0]) * 0.5;
  const REAL dy1= (points[3] - points[1]) * 0.5;
  const REAL ddy0 = (points[2] - 2.0 * points[1] + points[0]);
  const REAL ddy1 = (points[3] - 2.0 * points[2] + points[1]);
  const REAL y0 = points[1];
  const REAL y1 = points[2];

  // make it monotonic -- the pow calls cause a 10x slowdown!
  const REAL A = 2 * dy0;
  const REAL B = 2 * dy1;
  const REAL C = ddy0;
  const REAL D = ddy1;

  const REAL alpha = 4.0 * (B - D) / (pow(A, (REAL)0.25) * pow(B, (REAL)0.75));
  const REAL beta  = 6.0 * (D - C - 2.0 * B - 2.0 * A + 5.0) / (sqrt(A) * sqrt(B));
  const REAL gamma = 4.0 * (A + C) / (pow(A, (REAL)0.75) * pow(B, (REAL)0.25));

  const REAL test = (beta > 6.0) ? -2.0 * sqrt(beta - 2.0) : -(beta + 2.0) * 0.5f;
  if (!(alpha > test && gamma > test))
  {
    _quinticClamps++;
    return cubicInterp(interp, points);
  }

  // the quintic interpolation matrix
  const REAL a =  -6. * y0 - 3. * dy0 - 0.5 * ddy0 + 6.  * y1 - 3. * dy1 + 0.5 * ddy1;
  const REAL b =  15. * y0 + 8. * dy0 + 1.5 * ddy0 - 15. * y1 + 7. * dy1 - ddy1;
  const REAL c = -10. * y0 - 6. * dy0 - 1.5 * ddy0 + 10. * y1 - 4. * dy1 + 0.5 * ddy1;
  const REAL d =  0.5 *ddy0;
  const REAL e = dy0; 
  const REAL f = y0;

  const REAL x2 = interp * interp;
  const REAL x3 = x2 * interp;
  const REAL x4 = x2 * x2;
  const REAL x5 = x3 * x2;

  return a * x5 + b * x4 + c * x3 + d * x2 + e * interp + f;
}

///////////////////////////////////////////////////////////////////////
// do a quartic WENO interpolation
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::quarticInterp(const REAL interp, const REAL* points)
{
  const REAL& fim1 = points[0];
  const REAL& fi   = points[1];
  const REAL& fip1 = points[2];
  const REAL& fip2 = points[3];
  const REAL& x = interp;
  const REAL xHalf = interp * 0.5;

  const REAL p1 = fi + ((fip1 - fim1) + (fip1 - 2.0 * fi + fim1) * x) * xHalf;
  const REAL p2 = fi + ((-fip2 + 4.0 * fip1 - 3.0 * fi) + (fip2 - 2.0 * fip1 + fi) * x) * xHalf;

  const REAL third = 1.0 / 3.0;
  const REAL C1 = (2.0 - x) * third;
  const REAL C2 = (x + 1.0) * third;

  const REAL middle = -76.0 * fip1 * fi;
  const REAL fip1Sq = fip1 * fip1;
  const REAL fiSq = fi * fi;

  const REAL twelfth = 1.0 / 12.0;
  const REAL IS1 = (26.0 * fip1 * fim1 - 52.0 * fi * fim1 + middle + 25.0 * fip1Sq + 64.0 * fiSq + 13.0 * fim1 * fim1) * twelfth + 1e-6f;
  const REAL IS2 = (26.0 * fip2 * fi - 52.0 * fip2 * fip1 + middle + 25.0 * fiSq + 64.0 * fip1Sq + 13.0 * fip2 * fip2) * twelfth + 1e-6f;

  const REAL alpha1 = C1 / (IS1 * IS1);
  const REAL alpha2 = C2 / (IS2 * IS2);

  return (alpha1 * p1 + alpha2 * p2) / (alpha1 + alpha2);
}

///////////////////////////////////////////////////////////////////////
// do a quartic WENO interpolation
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::quarticInterpClamped(const REAL interp, const REAL* points)
{
  const REAL fim1 = points[0];
  const REAL fi   = points[1];
  const REAL fip1 = points[2];
  const REAL fip2 = points[3];
  const REAL x = interp;

  const REAL p1 = fi + ((fip1 - fim1) + (fip1 - 2.0 * fi + fim1) * x) * 0.5 * x;
  const REAL p2 = fi + ((-fip2 + 4.0 * fip1 - 3 * fi) + (fip2 - 2.0 * fip1 + fi) * x) * 0.5 * x;

  const REAL C1 = (2 - x) / 3.0;
  const REAL C2 = (x + 1) / 3.0;

  const REAL middle = -76 * fip1 * fi;
  const REAL fip1Sq = fip1 * fip1;
  const REAL fiSq = fi * fi;

  const REAL eps = 1e-6;
  const REAL IS1 = (26.0 * fip1 * fim1 - 52 * fi * fim1 + middle + 25 * fip1Sq + 64 * fiSq + 13 * fim1 * fim1) / 12.0 + eps;
  const REAL IS2 = (26.0 * fip2 * fi - 52 * fip2 * fip1 + middle + 25 * fiSq + 64 * fip1Sq + 13 * fip2 * fip2) / 12.0 + eps;

  const REAL alpha1 = C1 / (IS1 * IS1);
  const REAL alpha2 = C2 / (IS2 * IS2);

  const REAL sum = alpha1 + alpha2;
  const REAL w1 = alpha1 / sum;
  const REAL w2 = alpha2 / sum;

  const REAL result = w1 * p1 + w2 * p2;

  const REAL intervalMax = (points[1] > points[2]) ? points[1] : points[2];
  const REAL intervalMin = (points[1] > points[2]) ? points[2] : points[1];

  return (result < intervalMax) ? ((result > intervalMin) ? result : intervalMin) : intervalMax;

  // falling back to linear is probably a bad idea -- it introduced both a C0 and C1
  // discontinuity.
  //const REAL linear = x * points[2] + (1.0 - x) * points[1];
  //return (result < intervalMax) ? ((result > intervalMin) ? result : linear) : linear;
}

///////////////////////////////////////////////////////////////////////
// do a cubic Newton interpolation
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::cubicNewtonInterp(REAL interp, REAL* points)
{
  REAL x = interp;
  const REAL x1 = -1;
  const REAL x2 = 0;
  const REAL x3 = 1;
  const REAL x4 = 2;

  REAL y1 = points[0];
  REAL y2 = points[1];
  REAL y3 = points[2];
  REAL y4 = points[3];

  REAL c1 = y1;
  REAL c2 = (y2 - c1) / (x2 - x1);
  REAL c3 = (y3 - (c1 + c2 * (x3 - x1))) / ((x3 - x1) * (x3 - x2));
  REAL c4 = (y4 - ((c1 + c2 * (x4 - x1) + c3 * (x4 - x1) * (x4 - x2)))) / ((x4 - x1) * (x4 - x2) * (x4 - x3));

  return c1 + c2 * (x - x1) + c3 * (x - x1) * (x - x2) + c4 * (x - x1) * (x - x2) * (x - x3);
}

///////////////////////////////////////////////////////////////////////
// Compute the elements of the vertical derivative convolution kernel
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setToVerticalDerivativeKernel(double kMax, double dk, double sigma, double L)
{
  assert(_xRes % 2);
  assert(_xRes == _yRes);
  assert(_xRes == _zRes);

  double norm = 0;

  for (double k = 0; k < kMax; k += dk)
    norm += k * k * exp(-sigma * k * k);

  int halfWidth = _xRes / 2;

  for (int h = -halfWidth; h <= halfWidth; h++)
    for (int i = -halfWidth; i <= halfWidth; i++)
      for (int j = -halfWidth; j <= halfWidth; j++)
      {
        double r = sqrt((REAL)(i * i + j * j + h * h));
        double kern = 0;
        for (double k = 0; k < kMax; k += dk)
          kern += k * k * (sqrt(1.0 + k * k * L * L)) * exp(-sigma * k * k) * j0(r * k);
          //kern += k * k * (sqrt(k * k * L * L)) * exp(-sigma * k * k) * j0(r * k);

        (*this)(i + halfWidth, j + halfWidth, h + halfWidth) = kern / norm;
      }
}

///////////////////////////////////////////////////////////////////////
// set whole field to a Gaussian
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setToGaussian(REAL amplitude, VECTOR3 sigmas)
{
  REAL dx = 1.0 / (_xRes - 1);
  REAL dy = 1.0 / (_yRes - 1);
  REAL dz = 1.0 / (_zRes - 1);

  REAL x0 = 0.5;
  REAL y0 = 0.5;
  REAL z0 = 0.5;
  REAL xREAL = 0;
  REAL yREAL = 0;
  REAL zREAL = 0;

  for (int z = 0; z < _zRes; z++)
  {
    for (int y = 0; y < _yRes; y++)
    {
      for (int x = 0; x < _xRes; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        REAL xValue = xREAL;
        REAL yValue = yREAL;
        REAL zValue = zREAL;

        xValue = xValue - x0;
        xValue *= xValue;
        xValue *= 1.0 / (2.0 * sigmas[0] * sigmas[0]);

        yValue = yValue - y0;
        yValue *= yValue;
        yValue *= 1.0 / (2.0 * sigmas[1] * sigmas[1]);
        
        zValue = zValue - z0;
        zValue *= zValue;
        zValue *= 1.0 / (2.0 * sigmas[2] * sigmas[2]);

        _data[index] = exp(-(xValue + yValue + zValue));

        xREAL += dx;
      }
      xREAL = 0;
      yREAL += dy;
    }
    yREAL = 0;
    zREAL += dz;
  }

  // normalize
  (*this) *= amplitude / sum();
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolve(const FIELD_3D& filter)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(*this);

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = zHalf; z < _zRes - zHalf; z++)
    for (int y = yHalf; y < _yRes - yHalf; y++)
      for (int x = xHalf; x < _xRes - xHalf; x++)
      {
        result(x,y,z) = 0;
        for (int fz = 0; fz < filter.zRes(); fz++)
          for (int fy = 0; fy < filter.yRes(); fy++)
            for (int fx = 0; fx < filter.xRes(); fx++)
              result(x,y,z) += filter(fx,fy,fz) * (*this)(x + (fx - xHalf), 
                                                         y + (fy - yHalf), 
                                                         z + (fz - zHalf)); 
      }

  return result;
}


#define FILTERENTRY(X) const int filterIndex##X = filterPartial + X; \
                       const int dataIndex##X   = dataPartial + X; \
                       resultEntry += filter[filterIndex##X] * _data[dataIndex##X];

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveFast15(const FIELD_3D& filter)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(*this);

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = zHalf; z < _zRes - zHalf; z++)
    for (int y = yHalf; y < _yRes - yHalf; y++)
      for (int x = xHalf; x < _xRes - xHalf; x++)
      {
        // compute all the partial index stuff possible
        const int dataPartialTop = (z - zHalf) * _slabSize - xHalf + (y - yHalf) * _xRes + x;
        REAL resultEntry = 0;
        for (int fz = 0; fz < filter.zRes(); fz++)
        {
          // compute the z-varying index components
          const int filterPartialZ = fz * 225;
          const int dataPartialZ = fz * _slabSize  + dataPartialTop;
          for (int fy = 0; fy < filter.yRes(); fy++)
          {
            // compute the x-varying index components
            const int filterPartial = fy * 15 + filterPartialZ;
            const int dataPartial = fy * _xRes + dataPartialZ;
            FILTERENTRY(0);
            FILTERENTRY(1);
            FILTERENTRY(2);
            FILTERENTRY(3);
            FILTERENTRY(4);
            FILTERENTRY(5);
            FILTERENTRY(6);
            FILTERENTRY(7);
            FILTERENTRY(8);
            FILTERENTRY(9);
            FILTERENTRY(10);
            FILTERENTRY(11);
            FILTERENTRY(12);
            FILTERENTRY(13);
            FILTERENTRY(14);
          }
        }

        const int resultIndex = x + y * _xRes + z * _slabSize;
        result[resultIndex] = resultEntry;
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveToroidal(const FIELD_3D& filter)
{
  FIELD_3D result(*this);

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

#pragma omp parallel
  {
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        result(x,y,z) = 0;
        for (int fz = 0; fz < filter.zRes(); fz++)
          for (int fy = 0; fy < filter.yRes(); fy++)
            for (int fx = 0; fx < filter.xRes(); fx++)
            {
              int ffx = x + fx - xHalf;
              int ffy = y + fy - yHalf;
              int ffz = z + fz - zHalf;

              if (ffx < 0) ffx += _xRes;
              if (ffy < 0) ffy += _yRes;
              if (ffz < 0) ffz += _zRes;

              ffx = ffx % _xRes;
              ffy = ffy % _yRes;
              ffz = ffz % _zRes;

              result(x,y,z) += filter(fx,fy,fz) * (*this)(ffx, ffy, ffz);
            }
      }
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveNarrowBand(const FIELD_3D& filter, const FIELD_3D& distance, int maxCells)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(*this);
  result = 0;

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  // for the time being, assume the fields are the same red
  assert(_xRes == distance.xRes());
  assert(_yRes == distance.yRes());
  assert(_zRes == distance.zRes());

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

  REAL invDx = 1.0 / distance.dx();

  // find which cells are inside the band
  vector<int> xs;
  vector<int> ys;
  vector<int> zs;
  for (int z = zHalf; z < _zRes - zHalf; z++)
    for (int y = yHalf; y < _yRes - yHalf; y++)
      for (int x = xHalf; x < _xRes - xHalf; x++)
      {
        REAL currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance < maxCells)
        {
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
        }
      }

  int size = xs.size();
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int i = 0; i < size; i++)
  {
    int x = xs[i];
    int y = ys[i];
    int z = zs[i];
    for (int fz = 0; fz < filter.zRes(); fz++)
      for (int fy = 0; fy < filter.yRes(); fy++)
        for (int fx = 0; fx < filter.xRes(); fx++)
          result(x,y,z) += filter(fx,fy,fz) * (*this)(x + (fx - xHalf), 
                                                     y + (fy - yHalf), 
                                                     z + (fz - zHalf)); 
  }
  cout << endl;

  return result;
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveNarrowBand(const FIELD_3D& filter, const vector<int>& narrowBand)
{
  FIELD_3D result(*this);
  result = 0;

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

  // find which cells are inside the band
  vector<int> filteredNarrowBand;
  for (unsigned int index = 0; index < narrowBand.size(); index += 3)
  {
    int x = narrowBand[index];
    int y = narrowBand[index + 1];
    int z = narrowBand[index + 2];

    if (x < xHalf) continue;
    if (x > _xRes - xHalf - 1) continue; 
    if (y < yHalf) continue;
    if (y > _yRes - yHalf - 1) continue; 
    if (z < zHalf) continue;
    if (z > _zRes - zHalf - 1) continue;

    filteredNarrowBand.push_back(x);
    filteredNarrowBand.push_back(y);
    filteredNarrowBand.push_back(z);
  }

  int size = filteredNarrowBand.size();
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int i = 0; i < size; i += 3)
  {
    int x = filteredNarrowBand[i];
    int y = filteredNarrowBand[i + 1];
    int z = filteredNarrowBand[i + 2];

    for (int fz = 0; fz < filter.zRes(); fz++)
      for (int fy = 0; fy < filter.yRes(); fy++)
        for (int fx = 0; fx < filter.xRes(); fx++)
          result(x,y,z) += filter(fx,fy,fz) * (*this)(x + (fx - xHalf), 
                                                      y + (fy - yHalf), 
                                                      z + (fz - zHalf)); 
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveNarrowBandFast15(const FIELD_3D& filter, const vector<int>& narrowBand)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(*this);
  result = 0;

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

  // find which cells are inside the band
  vector<int> filteredNarrowBand;
  for (unsigned int index = 0; index < narrowBand.size(); index += 3)
  {
    int x = narrowBand[index];
    int y = narrowBand[index + 1];
    int z = narrowBand[index + 2];

    if (x < xHalf) continue;
    if (y < yHalf) continue;
    if (z < zHalf) continue;
    if (x > _xRes - xHalf - 1) continue; 
    if (y > _yRes - yHalf - 1) continue; 
    if (z > _zRes - zHalf - 1) continue;

    filteredNarrowBand.push_back(x);
    filteredNarrowBand.push_back(y);
    filteredNarrowBand.push_back(z);
  }

  int size = filteredNarrowBand.size();
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int i = 0; i < size; i += 3)
  {
    const int x = filteredNarrowBand[i];
    const int y = filteredNarrowBand[i + 1];
    const int z = filteredNarrowBand[i + 2];

    REAL resultEntry = 0;
    const int dataPartialTop = (z - zHalf) * _slabSize - xHalf + (y - yHalf) * _xRes + x;

    for (int fz = 0; fz < filter.zRes(); fz++)
    {
      const int filterPartialZ = fz * 225;
      const int dataPartialZ = fz * _slabSize  + dataPartialTop;
      for (int fy = 0; fy < filter.yRes(); fy++)
      {
        const int filterPartial = fy * 15 + filterPartialZ;
        const int dataPartial = fy * _xRes + dataPartialZ;

        FILTERENTRY(0);
        FILTERENTRY(1);
        FILTERENTRY(2);
        FILTERENTRY(3);
        FILTERENTRY(4);
        FILTERENTRY(5);
        FILTERENTRY(6);
        FILTERENTRY(7);
        FILTERENTRY(8);
        FILTERENTRY(9);
        FILTERENTRY(10);
        FILTERENTRY(11);
        FILTERENTRY(12);
        FILTERENTRY(13);
        FILTERENTRY(14);
      }
    }

    const int resultIndex = x + y * _xRes + z * _slabSize;
    result[resultIndex] = resultEntry;
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveNarrowBandFast15(const FIELD_3D& filter, const FIELD_3D& distance, const REAL maxRadius)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(*this);
  result = 0;
  const REAL invDx = 1.0 / distance.dx();

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  const int xHalf = filter.xRes() / 2;
  const int yHalf = filter.yRes() / 2;
  const int zHalf = filter.zRes() / 2;

  const int xRes = _xRes;
  const int yRes = _yRes;
  const int zRes = _zRes;

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        REAL currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance >= maxRadius) continue;
        if (x < xHalf || y < yHalf || z < zHalf || 
            x > _xRes - xHalf - 1 || y > _yRes - yHalf - 1 || z > _zRes - zHalf - 1) continue;

        REAL resultEntry = 0;
        const int dataPartialTop = (z - zHalf) * _slabSize - xHalf + (y - yHalf) * _xRes + x;

        for (int fz = 0; fz < filter.zRes(); fz++)
        {
          const int filterPartialZ = fz * 225;
          const int dataPartialZ = fz * _slabSize  + dataPartialTop;
          for (int fy = 0; fy < filter.yRes(); fy++)
          {
            const int filterPartial = fy * 15 + filterPartialZ;
            const int dataPartial = fy * _xRes + dataPartialZ;

            FILTERENTRY(0);
            FILTERENTRY(1);
            FILTERENTRY(2);
            FILTERENTRY(3);
            FILTERENTRY(4);
            FILTERENTRY(5);
            FILTERENTRY(6);
            FILTERENTRY(7);
            FILTERENTRY(8);
            FILTERENTRY(9);
            FILTERENTRY(10);
            FILTERENTRY(11);
            FILTERENTRY(12);
            FILTERENTRY(13);
            FILTERENTRY(14);
          }
        }

        const int resultIndex = x + y * _xRes + z * _slabSize;
        result[resultIndex] = resultEntry;
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// stomp all the value outside a narrow band to zero in order to 
// boost compression
///////////////////////////////////////////////////////////////////////
void FIELD_3D::stompOutsideNarrowBand(const FIELD_3D& distance, const int maxCells)
{
  TIMER functionTimer(__FUNCTION__);
  const REAL invDx = 1.0 / distance.dx();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        REAL currentDistance = fabs(distance(x,y,z) * invDx);

        if (currentDistance > maxCells)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// get the sum of the field
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::sum()
{
  REAL result = 0;
  for (int x = 0; x < _totalCells; x++)
    result += _data[x];

  return result;
}

///////////////////////////////////////////////////////////////////////
// insert a Gaussian at a specific point
///////////////////////////////////////////////////////////////////////
void FIELD_3D::insertGaussian(const VECTOR3& center, const REAL amplitude, const VECTOR3 sigmas)
{
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        VECTOR3 cell = cellCenter(x,y,z);

        REAL xValue = cell[0] - center[0];
        REAL yValue = cell[1] - center[1];
        REAL zValue = cell[2] - center[2];

        xValue *= xValue;
        xValue *= 1.0 / (2.0 * sigmas[0] * sigmas[0]);

        yValue *= yValue;
        yValue *= 1.0 / (2.0 * sigmas[1] * sigmas[1]);
        
        zValue *= zValue;
        zValue *= 1.0 / (2.0 * sigmas[2] * sigmas[2]);

        int index = x + y * _xRes + z * _slabSize;
        _data[index] += exp(-(xValue + yValue + zValue));
      }

  (*this) *= amplitude;
}

///////////////////////////////////////////////////////////////////////
// reset dimensions
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setLengths(const VECTOR3& lengths)
{
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;
}

///////////////////////////////////////////////////////////////////////
// determine how many non-zero entries are in the filter
///////////////////////////////////////////////////////////////////////
int FIELD_3D::nonZeroEntries()
{
  int nonZero = 0;

  for (int x = 0; x < _totalCells; x++)
    if (fabs(_data[x]) > 0)
      nonZero++;

  return nonZero;
}

///////////////////////////////////////////////////////////////////////
// return the projection of the filter in Z direction
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_3D::zProjection()
{
  FIELD_2D result(_xRes, _yRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y) += (*this)(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// return the projection of the filter in Y direction
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_3D::yProjection()
{
  FIELD_2D result(_xRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,z) += (*this)(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// return the projection of the filter in X direction
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_3D::xProjection()
{
  FIELD_2D result(_yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(y,z) += (*this)(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// set a given z slice to the given 2D field
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setSliceZ(const int& z, const FIELD_2D& slice)
{
  assert(slice.xRes() == _xRes);
  assert(slice.yRes() == _yRes);
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      (*this)(x,y,z) = slice(x,y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::fieldMax()
{
  assert(_totalCells > 0);

  REAL result = _data[0];

  for (int x = 0; x < _totalCells; x++)
    if (_data[x] > result)
      result = _data[x];

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3 FIELD_3D::maxIndex()
{
  REAL maxFound = _data[0];

  VECTOR3 maxFoundIndex;
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        if (_data[index] > maxFound)
        {
          maxFound = _data[index];

          maxFoundIndex[0] = x;
          maxFoundIndex[1] = y;
          maxFoundIndex[2] = z;
        }

  return maxFoundIndex;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::fieldMin()
{
  assert(_totalCells > 0);

  REAL result = _data[0];

  for (int x = 0; x < _totalCells; x++)
    if (_data[x] < result)
      result = _data[x];

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3 FIELD_3D::minIndex()
{
  REAL minFound = _data[0];

  VECTOR3 minFoundIndex;
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        if (_data[index] < minFound)
        {
          minFound = _data[index];

          minFoundIndex[0] = x;
          minFoundIndex[1] = y;
          minFoundIndex[2] = z;
        }

  return minFoundIndex;
}

///////////////////////////////////////////////////////////////////////
// normalize the data
///////////////////////////////////////////////////////////////////////
void FIELD_3D::normalize()
{
  REAL totalMin = fieldMin();
  REAL totalMax = fieldMax();

  for (int x = 0; x < _totalCells; x++)
    _data[x] = (_data[x] - totalMin) / (totalMax - totalMin);
}

///////////////////////////////////////////////////////////////////////
// flip the z and y coordinates
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::flipZY() const
{
  FIELD_3D result(_xRes, _zRes, _yRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,z,y) = (*this)(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// flip the x and y coordinates
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::flipXY() const
{
  FIELD_3D result(_yRes, _xRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(y,x,z) = (*this)(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// flip the x and z coordinates
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::flipXZ() const
{
  FIELD_3D result(_zRes, _yRes, _xRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(z,y,x) = (*this)(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// create a mirror image along Z
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::mirrorZ() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = (*this)(x,y,(_zRes - 1) - z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// create a mirror image along Y
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::mirrorY() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = (*this)(x,(_yRes - 1) - y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// create a mirror image along Z
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::mirrorX() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = (*this)((_xRes - 1) - x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudini(string filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  bool volumeFound = false;
  char buffer[256];
  string search("Volume");
  while (!volumeFound)
  {
    fscanf(file, "%s", buffer);
    if (search.compare(buffer) == 0)
      volumeFound = true;
  }

  // read off the collision, vel.x, vel.y, and vel.z volumes
  readHoudiniField(file, true);
  readHoudiniField(file, false);
  readHoudiniField(file, false);
  readHoudiniField(file, false);
  
  // get the distance field
  readHoudiniField(file, false);

  fclose(file);
}
	  
///////////////////////////////////////////////////////////////////////
// read in a triplet of fields for velocity
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudiniVel(const string filename, FIELD_3D& xVelocity, FIELD_3D& yVelocity, FIELD_3D& zVelocity)
{
  FILE* file;
  file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
	  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
	  cout << " FIELD_3D read failed! " << endl;
	  cout << " Could not open file " << filename.c_str() << endl;
	  exit(0);
  }
  
  bool volumeFound = false;
  char buffer[256];
  string search("Volume");
  while (!volumeFound)
  {
	  fscanf(file, "%s", buffer);
	  if (search.compare(buffer) == 0)
		  volumeFound = true;
  }

  xVelocity = FIELD_3D::readHoudiniField(file);
  yVelocity = FIELD_3D::readHoudiniField(file);
  zVelocity = FIELD_3D::readHoudiniField(file);

  fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudiniSurf(string filename)
{
  cout << " Reading Houdini distance file ..."; flush(cout);
  FILE* file;
  file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
	  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
	  cout << " FIELD_3D read failed! " << endl;
	  cout << " Could not open file " << filename.c_str() << endl;
	  exit(0);
  }
  
  bool volumeFound = false;
  char buffer[256];
  string search("Volume");
  while (!volumeFound)
  {
	  fscanf(file, "%s", buffer);
	  if (search.compare(buffer) == 0)
		  volumeFound = true;
  }
  
  // get the distance field
  readHoudiniField(file, true);
  
  fclose(file);
  cout << " done." << endl;
}

///////////////////////////////////////////////////////////////////////
// This is hard-coded for the paddle example
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudini12Surf(string filename)
{
  cout << " Reading Houdini 12 distance file " << filename.c_str() << " ..."; flush(cout);
  int size = filename.size();
  if (filename[size - 1] == 'z' && filename[size - 2] == 'g')
  {
    cout << " File is Gzipped! " << endl;
    exit(0);
  }
  
  // MAGIC NUMBER - just for the paddle example
	_lengths[0] = 2.5;
	_lengths[1] = 2.5;
	_lengths[2] = 2.5;

  if (_data) delete[] _data;

  // MAGIC NUMBER - just for the paddle example
  _xRes = 100; 
  _yRes = 100; 
  _zRes = 100;
  
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new REAL[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _outside = maxRes() * maxRes();

  FILE* file;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open filename " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  char buffer[512];
  bool tileFound = false;

  while (!tileFound)
  {
    fscanf(file, "%s", buffer);
    string search(buffer); 

    // look for the one that says "tiles"
    size_t findVolume = search.find("\"tiles\"");
    if (findVolume != string::npos)
      tileFound = true;
  }

  int xTotalTiles = (_xRes / 16);
  int yTotalTiles = (_yRes / 16);
  int zTotalTiles = (_zRes / 16);

  if (xTotalTiles * 16 != _xRes) xTotalTiles++;
  if (yTotalTiles * 16 != _yRes) yTotalTiles++;
  if (zTotalTiles * 16 != _zRes) zTotalTiles++;

  for (int zTile = 0; zTile < zTotalTiles; zTile++)
    for (int yTile = 0; yTile < yTotalTiles; yTile++)
      for (int xTile = 0; xTile < xTotalTiles; xTile++)
      {
        // peel off '['
        fscanf(file, "%s", buffer);

        if (buffer[0] != '[')
          break;

        // peel off 'compression'
        fscanf(file, "%s", buffer);

        // peel off whitespace and "
        char single;
        fscanf(file, " %c", &single);
        double data;
        fscanf(file, "data\",[");

        bool closingBracketFound = false;

        // read in a block
        vector<double> blockData;
        while (!closingBracketFound)
        {
          fscanf(file, "%lf", &data);
          blockData.push_back(data);

          // peek at the next token
          fscanf(file, "%c", &single);

          if (single == ']')
            closingBracketFound = true;
        }

        int xTileSize = 16;
        int yTileSize = 16;
        int zTileSize = 16;

        if (xTile == xTotalTiles - 1)
          xTileSize = _xRes % 16;
        if (yTile == yTotalTiles - 1)
          yTileSize = _yRes % 16;
        if (zTile == zTotalTiles - 1)
          zTileSize = _zRes % 16;

        int xTileOffset = xTile * 16;
        int yTileOffset = yTile * 16;
        int zTileOffset = zTile * 16;

        int index = 0;
        for (int z = 0; z < zTileSize; z++)
          for (int y = 0; y < yTileSize; y++)
            for (int x = 0; x < xTileSize; x++, index++)
              (*this)(xTileOffset + x, yTileOffset + y,zTileOffset + z) = blockData[index];

        fscanf(file, "%s", buffer);
      }
  fclose(file);
}

///////////////////////////////////////////////////////////////////////
// read a Houdini field off from a file stream -- it is assumed that the
// file is already advanced to the beginning of the field
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudiniField(FILE* file, bool storeValues)
{
  // peg it to a double
  typedef Eigen::Matrix<double, 3, 1> VECTOR3D;
  VECTOR3D v0, v1, v2;

  char buffer[256];
  int index;
  int xRes, yRes, zRes;
  REAL negativeTwo;

  fscanf(file, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i %i %i", &index, 
      &v0[0], &v0[1], &v0[2],
      &v1[0], &v1[1], &v1[2],
      &v2[0], &v2[1], &v2[2],
      &negativeTwo,
      &xRes, &yRes, &zRes);
  
	_lengths[0] = v0[0];
	_lengths[1] = v1[1];
	_lengths[2] = v2[2];	

  // initialize things if we aren't going to throw the field away
  if (storeValues)
  {
    if (_data) delete[] _data;
    _xRes = xRes; 
    _yRes = yRes; 
    _zRes = zRes; 

    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;
    _data = new REAL[_totalCells];

    _dx = _lengths[0] / _xRes;
    _dy = _lengths[1] / _yRes;
    _dz = _lengths[2] / _zRes;
    _invDx = 1.0 / _dx;
    _invDy = 1.0 / _dy;
    _invDz = 1.0 / _dz;

    _outside = maxRes() * maxRes();
  }

  // read in "streak"
  fscanf(file, "%s", buffer);

  // read in two zeros
  int zero1, zero2;
  fscanf(file, "%i %i", &zero1, &zero2);

  // read in "iso" or "invisible"
  fscanf(file, "%s", buffer);

  // read in zero and ten
  fscanf(file, "%i %i", &zero1, &zero2);

  if (storeValues)
  {
	  for (int y = 0; y < xRes * yRes * zRes; y++)
    {
      double entry;
      fscanf(file, "%lf", &entry);
      _data[y] = (REAL)entry;
	  }
  }
  else
  {
    double data;
	  for (int y = 0; y < xRes * yRes * zRes; y++){
      fscanf(file, "%lf", &data);
	  }
  }

  // chomp [0 4]
  fscanf(file, "%s", buffer);
  fscanf(file, "%s", buffer);
}

///////////////////////////////////////////////////////////////////////
// read a Houdini field off from a file stream -- it is assumed that the
// file is already advanced to the beginning of the field
//
// a static version for the velocity field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::readHoudiniField(FILE* file)
{
  // peg it to a double
  typedef Eigen::Matrix<double, 3, 1> VECTOR3D;
  VECTOR3D v0, v1, v2;

  char buffer[256];
  int index;
  int xRes, yRes, zRes;
  REAL negativeTwo;
  fscanf(file, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i %i %i", &index, 
      &v0[0], &v0[1], &v0[2],
      &v1[0], &v1[1], &v1[2],
      &v2[0], &v2[1], &v2[2],
      &negativeTwo,
      &xRes, &yRes, &zRes);

  VECTOR3 lengths(v0[0], v1[1], v2[2]);
  VECTOR3 center;

  // initialize things if we aren't going to throw the field away
  FIELD_3D result(xRes, yRes, zRes, center, lengths);

  // read in "streak"
  fscanf(file, "%s", buffer);

  // read in two zeros
  int zero1, zero2;
  fscanf(file, "%i %i", &zero1, &zero2);

  // read in "iso" or "invisible"
  fscanf(file, "%s", buffer);

  // read in zero and ten
  fscanf(file, "%i %i", &zero1, &zero2);

  REAL* data = result.data();
	for (int y = 0; y < xRes * yRes * zRes; y++)
  {
    double entry;
    fscanf(file, "%lf", &entry);
    data[y] = (REAL)entry;
  }

  // chomp [0 4]
  fscanf(file, "%s", buffer);
  fscanf(file, "%s", buffer);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::Dx(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const REAL right = (x < _xRes - 1) ? _data[index + 1] : _data[index];
  const REAL left  = (x > 0)         ? _data[index - 1] : _data[index];
  const REAL denom = (x > 0 && x < _xRes -1) ? 1.0 / (2.0 * _dx) : 1.0 / _dx;
  return (right - left) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::DDx(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const REAL right = (x < _xRes - 1) ? _data[index + 1] : _data[index];
  const REAL left  = (x > 0)         ? _data[index - 1] : _data[index];
  const REAL denom = 1.0 / (_dx * _dx);
  return (right - 2.0 * _data[index] + left) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::Dy(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const REAL up   = (y < _yRes - 1) ? _data[index + _xRes] : _data[index];
  const REAL down = (y > 0)         ? _data[index - _xRes] : _data[index];
  const REAL denom = (y > 0 && y < _yRes -1) ? 1.0 / (2.0 * _dy) : 1.0 / _dy;
  return (up - down) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::DDy(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const REAL up   = (y < _yRes - 1) ? _data[index + _xRes] : _data[index];
  const REAL down = (y > 0)         ? _data[index - _xRes] : _data[index];
  const REAL denom = 1.0 / (_dy * _dy);
  return (up - 2.0 * _data[index] + down) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::Dz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const REAL in  = (z < _zRes - 1) ? _data[index + _slabSize] : _data[index];
  const REAL out = (z > 0)         ? _data[index - _slabSize] : _data[index];
  const REAL denom = (z > 0 && z < _zRes -1) ? 1.0 / (2.0 * _dz) : 1.0 / _dz;
  return (in - out) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::DDz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const REAL in  = (z < _zRes - 1) ? _data[index + _slabSize] : _data[index];
  const REAL out = (z > 0)         ? _data[index - _slabSize] : _data[index];
  const REAL denom = 1.0 / (_dz * _dz);
  return (in - out) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::DDxy(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  // if it's on the border, just use the center cell
  const int xPlus  = (x < _xRes - 1) ? 1 : 0;
  const int xMinus = (x > 0) ? -1 : 0;
  const int yPlus  = (y < _yRes - 1) ? _xRes : 0;
  const int yMinus = (y > 0) ? -_xRes : 0;

  const REAL plusPlus   = _data[index + xPlus + yPlus];
  const REAL plusMinus  = _data[index + xPlus + yMinus];
  const REAL minusPlus  = _data[index + xMinus + yPlus];
  const REAL minusMinus = _data[index + xMinus + yMinus];

  return (plusPlus - plusMinus - minusPlus + minusMinus) / (_dx * _dy);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::DDxz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  // if it's on the border, just use the center cell
  const int xPlus  = (x < _xRes - 1) ? 1 : 0;
  const int xMinus = (x > 0) ? -1 : 0;
  const int zPlus  = (z < _zRes - 1) ? _slabSize : 0;
  const int zMinus = (z > 0) ? -_slabSize : 0;

  const REAL plusPlus   = _data[index + xPlus + zPlus];
  const REAL plusMinus  = _data[index + xPlus + zMinus];
  const REAL minusPlus  = _data[index + xMinus + zPlus];
  const REAL minusMinus = _data[index + xMinus + zMinus];

  return (plusPlus - plusMinus - minusPlus + minusMinus) / (_dx * _dz);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::DDyz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  // if it's on the border, just use the center cell
  const int yPlus  = (y < _yRes - 1) ? _xRes : 0;
  const int yMinus = (y > 0) ? -_xRes : 0;
  const int zPlus  = (z < _zRes - 1) ? _slabSize : 0;
  const int zMinus = (z > 0) ? -_slabSize : 0;

  const REAL plusPlus   = _data[index + yPlus + zPlus];
  const REAL plusMinus  = _data[index + yPlus + zMinus];
  const REAL minusPlus  = _data[index + yMinus + zPlus];
  const REAL minusMinus = _data[index + yMinus + zMinus];

  return (plusPlus - plusMinus - minusPlus + minusMinus) / (_dy * _dz);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::Dx() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = Dx(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::Dy() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = Dy(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::Dz() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = Dz(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDx() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = DDx(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDy() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = DDy(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDz() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = DDz(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDxy() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = DDxy(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDxz() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = DDxz(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDyz() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x,y,z) = DDyz(x,y,z);

  return result;
}

///////////////////////////////////////////////////////////////////////
// get the mean curvature
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::meanCurvature() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        REAL px = Dx(x,y,z);
        REAL py = Dy(x,y,z);
        REAL pz = Dz(x,y,z);
        
        REAL pxx = DDx(x,y,z);
        REAL pyy = DDy(x,y,z);
        REAL pzz = DDz(x,y,z);
        
        REAL pxy = DDxy(x,y,z);
        REAL pxz = DDxz(x,y,z);
        REAL pyz = DDyz(x,y,z);

        result(x,y,z) = (pyy + pzz) * px * px +
                       (pxx + pzz) * py * py +
                       (pxx + pyy) * pz * pz - 
                       2.0 * (px * py * pxy + px * pz * pxz + py * pz * pyz);

        REAL denom = pow(px + py + pz, (REAL)1.5);
        denom = (fabs(denom) > 1e-6) ? 1.0 / denom : 0;

        result(x,y,z) *= denom;
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// get the Gaussian curvature
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::gaussianCurvature() const
{
  FIELD_3D result(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        REAL px = Dx(x,y,z);
        REAL py = Dy(x,y,z);
        REAL pz = Dz(x,y,z);
        
        REAL pxx = DDx(x,y,z);
        REAL pyy = DDy(x,y,z);
        REAL pzz = DDz(x,y,z);
        
        REAL pxy = DDxy(x,y,z);
        REAL pxz = DDxz(x,y,z);
        REAL pyz = DDyz(x,y,z);

        result(x,y,z) = (pyy + pzz - pyz * pyz) * px * px +
                        (pxx + pzz - pxz * pxz) * py * py +
                        (pxx + pyy - pxy * pxy) * pz * pz + 
                        2.0 * (px * py * (pxz * pyz - pxy * pzz) + 
                               py * pz * (pxy * pxz - pyz * pxx) +
                               px * pz * (pxy * pyz - pxz * pyy));

        REAL denom = pow(px + py + pz, (REAL)2.0);
        denom = (fabs(denom) > 1e-6) ? 1.0 / denom : 0;

        result(x,y,z) *= denom;
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// mask out any values past a certain distance
///////////////////////////////////////////////////////////////////////
void FIELD_3D::maskByDistance(const FIELD_3D& distanceField, const REAL distance)
{
  assert(_xRes == distanceField.xRes());
  assert(_yRes == distanceField.yRes());
  assert(_zRes == distanceField.zRes());

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (fabs(distanceField(x,y,z)) > distance)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// clamp the field to a min and max
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clamp(const REAL minValue, const REAL maxValue)
{
  TIMER functionTimer(__FUNCTION__);
  assert(minValue <= maxValue);

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        bool greater = (_data[index] > maxValue);
        bool lesser = (_data[index] < minValue);

        if (greater)
        {
          _data[index] = maxValue;
          continue;
        }

        if (lesser)
          _data[index] = minValue;
      }
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleSextic(int xRes, int yRes, int zRes) const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(xRes, yRes, zRes, _center, _lengths);

  int index = 0;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++)
      {
        VECTOR3 center = result.cellCenter(x,y,z);
        result[index] = this->sexticLookupClamped(center);
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleCubic(int xRes, int yRes, int zRes) const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(xRes, yRes, zRes, _center, _lengths);

  int index = 0;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++)
      {
        VECTOR3 center = result.cellCenter(x,y,z);
        result[index] = this->cubicLookup(center);
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleCubicUnclamped(int xRes, int yRes, int zRes) const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(xRes, yRes, zRes, _center, _lengths);
  const int slabSize = xRes * yRes;

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        VECTOR3 center = result.cellCenter(x,y,z);
        result[index] = this->cubicLookupUnclamped(center);
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleCubicUnclampedNarrowBand(int upResFactor, const FIELD_3D& distanceField, const int maxRadius) const
{
  TIMER functionTimer(__FUNCTION__);

  int xRes = upResFactor * _xRes;
  int yRes = upResFactor * _yRes;
  int zRes = upResFactor * _zRes;

  int maxRes = (xRes > yRes) ? xRes : yRes;
  maxRes = (zRes > maxRes) ? zRes : maxRes;

  FIELD_3D result(xRes, yRes, zRes, _center, _lengths);
  //const int slabSize = xRes * yRes;
  const REAL invDx = 1.0 / distanceField.dx();
  const float outside = sqrt(2.0 * maxRes * maxRes); 

  result = outside;

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        const REAL currentDistance = fabs(distanceField(x,y,z) * invDx);
        if (currentDistance >= maxRadius) 
          continue;

        const int xStart = x * upResFactor;
        const int yStart = y * upResFactor;
        const int zStart = z * upResFactor;
        for (int k = 0; k < upResFactor; k++)
        {
          int zUpres = zStart + k;
          for (int j = 0; j < upResFactor; j++)
          {
            int yUpres = yStart + j;
            for (int i = 0; i < upResFactor; i++)
            {
              int xUpres = xStart + i;

              VECTOR3 center = result.cellCenter(xUpres,yUpres,zUpres);
              result(xUpres, yUpres, zUpres) = this->cubicLookupUnclamped(center);
            }
          }
        }
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleQuintic(int xRes, int yRes, int zRes) const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(xRes, yRes, zRes, _center, _lengths);

  int index = 0;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++)
      {
        VECTOR3 center = result.cellCenter(x,y,z);
        result[index] = this->quinticLookup(center);
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// take the square root of the field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::squareRoot()
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(*this);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        result[index] = std::sqrt(_data[index]);

  return result;
}

///////////////////////////////////////////////////////////////////////
// from "Level Set Surface Editing Operators", Museth et al. 2002
// "Geometric Surface Processing via Normal Maps", Tasdizen 2003
///////////////////////////////////////////////////////////////////////
void FIELD_3D::principalCurvatures(FIELD_3D& minCurvature, FIELD_3D& maxCurvature) const
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        MATRIX3 N = Dnormal(x,y,z);
        VECTOR3 n = normal(x,y,z);
        MATRIX3 outer = n * n.transpose();

        MATRIX3 B = N * (MATRIX3::Identity() - outer);

        REAL D = sqrt(B.squaredNorm());
        REAL H = B.trace() * 0.5;

        REAL discrim = D * D * 0.5 - H * H;

        if (discrim < 0.0)
        {
          minCurvature[index] = 0;
          maxCurvature[index] = 0;
          continue;
        }

        REAL root = sqrt(discrim);
        REAL k1 = H + root;
        REAL k2 = H - root;

        maxCurvature[index] = (fabs(k1) > fabs(k2)) ? k1 : k2;
        minCurvature[index] = (fabs(k1) > fabs(k2)) ? k2 : k1;
      }
}

///////////////////////////////////////////////////////////////////////
// from "Level Set Surface Editing Operators", Museth et al. 2002
// "Geometric Surface Processing via Normal Maps", Tasdizen 2003
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::principalCurvature() const
{
  FIELD_3D result(*this);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        MATRIX3 N = Dnormal(x,y,z);
        VECTOR3 n = normal(x,y,z);
        MATRIX3 outer = n * n.transpose();

        MATRIX3 B = N * (MATRIX3::Identity() - outer);

        REAL D = sqrt(B.squaredNorm());
        REAL H = B.trace() * 0.5;

        REAL discrim = D * D * 0.5 - H * H;

        if (discrim < 0.0)
        {
          result[index] = 0;
          continue;
        }

        REAL root = sqrt(discrim);
        REAL k1 = H + root;
        REAL k2 = H - root;

        result[index] = (fabs(k1) > fabs(k2)) ? k1 : k2;
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// get the normal at a point
///////////////////////////////////////////////////////////////////////
VECTOR3 FIELD_3D::normal(int x, int y, int z) const
{
  assert (x >= 0 && x < _xRes);
  assert (y >= 0 && y < _yRes);
  assert (z >= 0 && z < _zRes);

  VECTOR3 result;
  result[0] = Dx(x,y,z);
  result[1] = Dy(x,y,z);
  result[2] = Dz(x,y,z);

  result.normalize();

  return result;
}

///////////////////////////////////////////////////////////////////////
// get the normal at a point
///////////////////////////////////////////////////////////////////////
MATRIX3 FIELD_3D::Dnormal(int x, int y, int z) const
{
  const VECTOR3 left  = (x == 0)         ? normal(x,y,z) : normal(x-1,y,z);
  const VECTOR3 right = (x == _xRes - 1) ? normal(x,y,z) : normal(x+1,y,z);
  const REAL dx     = (x == 0 || x == _xRes - 1) ? 1.0 / _dx : 0.5 / _dx;

  const VECTOR3 down = (y == 0)         ? normal(x,y,z) : normal(x,y-1,z);
  const VECTOR3 up   = (y == _yRes - 1) ? normal(x,y,z) : normal(x,y+1,z);
  const REAL dy    = (y == 0 || y == _yRes - 1) ? 1.0 / _dy : 0.5 / _dy;

  const VECTOR3 out = (z == 0)         ? normal(x,y,z) : normal(x,y,z-1);
  const VECTOR3 in  = (z == _zRes - 1) ? normal(x,y,z) : normal(x,y,z+1);
  const REAL dz   = (z == 0 || z == _zRes - 1) ? 1.0 / _dz : 0.5 / _dz;

  MATRIX3 result;
  result.col(0) = (right - left) * dx;
  result.col(1) = (up - down) * dy;
  result.col(2) = (in - out) * dz;

  return result;
}

///////////////////////////////////////////////////////////////////////
// compute laplace everywhere (not the most efficient way for now)
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::computeLaplace(REAL factor) const
{
  FIELD_3D result( *this );

  { 
    // get the fine-res indices to loop over
    const int zStart = 1;
    const int yStart = 1;
    const int xStart = 1;

    int zEnd   = this->zRes()-1;
    int yEnd   = this->yRes()-1;
    int xEnd   = this->xRes()-1;

    for (int z = zStart; z < zEnd; z++)
      for (int y = yStart; y < yEnd; y++)
        for (int x = xStart; x < xEnd; x++) {
          // directly downscale, to prevent tiny accelerations for the wave eq.
          result(x,y,z) = this->laplace(x,y,z) * factor;
		}
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::laplace(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  const int xp = (x < _xRes-1) ? x+1 : x;
  const int xm = (x > 0      ) ? x-1 : x;
  const int yp = (y < _yRes-1) ? y+1 : y;
  const int ym = (y > 0      ) ? y-1 : y;
  const int zp = (z < _zRes-1) ? z+1 : z;
  const int zm = (z > 0      ) ? z-1 : z;

  REAL value = 0.;
  value += (*this)(xm,y,z) + (*this)(xp,y,z);
  value += (*this)(x,ym,z) + (*this)(x,yp,z);
  value += (*this)(x,y,zm) + (*this)(x,y,zp);
  value +=  -6.* (*this)(x,y,z+0);

  // assume dx==dy==dz
  const REAL denom = (_dx*_dx);

  return value / denom;
}

///////////////////////////////////////////////////////////////////////
// BLAS-like interface, output += alpha * input
///////////////////////////////////////////////////////////////////////
void FIELD_3D::axpy(const REAL& alpha, const FIELD_3D& input, FIELD_3D& output)
{
  assert(input.xRes() == output.xRes());
  assert(input.yRes() == output.yRes());
  assert(input.zRes() == output.zRes());

  int totalCells = input.totalCells();
  for (int x = 0; x < totalCells; x++)
    output[x] += alpha * input[x];
}

///////////////////////////////////////////////////////////////////////
// clamp the field to a min and max
///////////////////////////////////////////////////////////////////////
void FIELD_3D::bandPass(const REAL minValue, const REAL maxValue)
{
  assert(minValue < maxValue);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if ((*this)(x,y,z) < minValue)
          (*this)(x,y,z) = 0;
        if ((*this)(x,y,z) > maxValue)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// isolate values near the current value, within a certain width
///////////////////////////////////////////////////////////////////////
void FIELD_3D::isolateBand(const REAL target, const REAL width)
{
  int unstomped = 0;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        REAL value = fabs((*this)(x,y,z));

        if (value > target + width || value < target - width)
          (*this)(x,y,z) = 0;
        else
          unstomped++;
      }
}

///////////////////////////////////////////////////////////////////////
// single explicit diffusion step
///////////////////////////////////////////////////////////////////////
void FIELD_3D::blur(const REAL dt)
{
  FIELD_3D temp(*this);

  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _slabSize;
        REAL mean = -6.0 * _data[index] + 
                    _data[index + 1] + _data[index - 1] +
                    _data[index + _xRes] + _data[index - _xRes] +
                    _data[index + _slabSize] + _data[index - _slabSize];

        temp(x,y,z) = _data[index] + dt * mean;
      }

  *this = temp;
}

///////////////////////////////////////////////////////////////////////
// get band-limited curvature, with some diffusion
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::bandLimitedCurvature(const REAL target, const REAL width)
{
  FIELD_3D result(*this);

  result = this->principalCurvature();

  result.maskByDistance(*this, _dx * 5);

  result.isolateBand(target, width);
 
  return result;
}

///////////////////////////////////////////////////////////////////////
// set to the absolute value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::absoluteValue()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = fabs(_data[x]);
}

///////////////////////////////////////////////////////////////////////
// do a soft bandpass where there's a gradual cubic falloff
///////////////////////////////////////////////////////////////////////
void FIELD_3D::softBandPass(const REAL band, const REAL falloff)
{
  for (int x = 0; x < _totalCells; x++)
    if (fabs(fabs(_data[x]) - band) > falloff)
      _data[x] = 0;
    else
    {
      REAL interp = fabs(fabs(_data[x]) - band) / falloff;
      REAL squared = interp * interp;
      REAL scaling = 2 * squared * interp - 3 * squared + 1;
      _data[x] *= scaling;
    }
}

///////////////////////////////////////////////////////////////////////
// write a VECTOR3 to a file stream
///////////////////////////////////////////////////////////////////////
static void writeVector3(const VECTOR3& v, FILE* file)
{
  for (int x = 0; x < 3; x++)
    fwrite((void*)&v[x], sizeof(REAL), 1, file);
}

///////////////////////////////////////////////////////////////////////
// write out a field to a file stream
///////////////////////////////////////////////////////////////////////
void FIELD_3D::write(FILE* file) const
{
  assert(_xRes > 0);
  assert(_yRes > 0);
  assert(_zRes > 0);

  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);
  fwrite((void*)&_zRes, sizeof(int), 1, file);
  //_center.write(file);
  //_lengths.write(file);
  writeVector3(_center, file);
  writeVector3(_lengths, file);

  // always write out as a double
  if (sizeof(REAL) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    for (int x = 0; x < _totalCells; x++)
      dataDouble[x] = _data[x];

    // if the field sizes get too large, the file system starts complaining, so only write out 
    // one slab at a time
    for (int z = 0; z < _zRes; z++)
      fwrite((void*)(&dataDouble[z * _slabSize]), sizeof(double), _slabSize, file);
    
    delete[] dataDouble;
  }
  else
  {
    for (int z = 0; z < _zRes; z++)
      fwrite((void*)(&_data[z * _slabSize]), sizeof(double), _slabSize, file);
  }
}


///////////////////////////////////////////////////////////////////////
// read a VECTOR3 from a file stream
///////////////////////////////////////////////////////////////////////
static void readVector3(VECTOR3& v, FILE* file)
{
  for (int x = 0; x < 3; x++)
    fread((void*)&v[x], sizeof(REAL), 1, file);
}

///////////////////////////////////////////////////////////////////////
// read in a field from a file stream
///////////////////////////////////////////////////////////////////////
void FIELD_3D::read(FILE* file)
{
  // read dimensions
  fread((void*)&_xRes, sizeof(int), 1, file);
  fread((void*)&_yRes, sizeof(int), 1, file);
  fread((void*)&_zRes, sizeof(int), 1, file);
  readVector3(_center, file);
  readVector3(_lengths, file);
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;
  _outside = maxRes() * maxRes();

  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;

  assert(_xRes > 0);
  assert(_yRes > 0);
  assert(_zRes > 0);

  if (_data) delete[] _data;
  try {
    _data = new REAL[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(REAL);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  // always read in as a double
  if (sizeof(REAL) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];

    // if the field sizes get too large, the file system starts complaining, so only write out 
    // one slab at a time
    for (int z = 0; z < _zRes; z++)
      fread((void*)(&dataDouble[z * _slabSize]), sizeof(double), _slabSize, file);

    for (int x = 0; x < _totalCells; x++)
      _data[x] = dataDouble[x];

    delete[] dataDouble;
  }
  else
  {
    for (int z = 0; z < _zRes; z++)
      fread((void*)(&_data[z * _slabSize]), sizeof(double), _slabSize, file);
  }
}

///////////////////////////////////////////////////////////////////////
// see if the field has any data in it yet
///////////////////////////////////////////////////////////////////////
const bool FIELD_3D::initialized() const
{
  if (_xRes < 0 || _yRes < 0 || _zRes < 0 || _totalCells < 0 || _data == NULL)
    return false;

  return true;
}

///////////////////////////////////////////////////////////////////////
// copy values out into the border, assuming that "borderSize" is the 
// width of the grid padding
///////////////////////////////////////////////////////////////////////
void FIELD_3D::copyIntoBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x == borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x - i,y,z) = value;
        }
        if (x == _xRes - 1 - borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x + i,y,z) = value;
        }              
        if (y == borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y - i,z) = value;
        }
        if (y == _yRes - 1 - borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y+i,z) = value;
        }              
        if (z == borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z - i) = value;
        }
        if (z == _zRes - 1 - borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z+i) = value;
        }

        if (x == borderSize && y == borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y - j,z) = value;
        }
        if (x == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y + j,z) = value;
        }

        // handle the corners
        if (x == borderSize && z == borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y,z-j) = value;
        }
        if (x == _xRes - 1 - borderSize && z == _zRes - 1 - borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y,z+j) = value;
        }

        if (z == borderSize && y == borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y - j,z -i) = value;
        }
        if (z == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          REAL value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y + j,z+i) = value;
        }

      }
}

///////////////////////////////////////////////////////////////////////
// pass back a field with a new padding of size "paddingSize"
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::withAddedPadding(int paddingSize) const
{
  // new length, with padding
  VECTOR3 newLengths = _lengths;
  newLengths[0] += paddingSize * 2 * _dx;
  newLengths[1] += paddingSize * 2 * _dx;
  newLengths[2] += paddingSize * 2 * _dx;

  FIELD_3D result(_xRes + 2 * paddingSize, 
                 _yRes + 2 * paddingSize, 
                 _zRes + 2 * paddingSize, _center, newLengths);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        result(x + paddingSize,
              y + paddingSize,
              z + paddingSize) = (*this)(x,y,z);

  result.copyIntoBorder(paddingSize);

  return result;
}

///////////////////////////////////////////////////////////////////////
// stomp the border to zero
///////////////////////////////////////////////////////////////////////
void FIELD_3D::stompBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x < borderSize || y < borderSize || z < borderSize ||
            x >= _xRes - borderSize || y >= _yRes - borderSize || z >= _zRes - borderSize)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// set the border to value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setBorder(int borderSize, REAL value)
{
  TIMER functionTimer(__FUNCTION__);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x < borderSize || y < borderSize || z < borderSize ||
            x >= _xRes - borderSize || y >= _yRes - borderSize || z >= _zRes - borderSize)
          (*this)(x,y,z) = value;
      }
}

///////////////////////////////////////////////////////////////////////
// compute the narrow band indices for this object, which is assumed 
// to be a SDF
///////////////////////////////////////////////////////////////////////
vector<int> FIELD_3D::computeNarrowBand(REAL maxCellDistance) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<int> triplets;
  REAL invDx = 1.0 / _dx;

  // find which cells are inside the band
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        REAL currentDistance = fabs( (*this)(x,y,z) * invDx);
        if (currentDistance < maxCellDistance)
        {
          triplets.push_back(x);
          triplets.push_back(y);
          triplets.push_back(z);
        }
      }

  return triplets;
}

///////////////////////////////////////////////////////////////////////
// do a diff between two fields, but only along a narrow band
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::narrowBandDiff(const FIELD_3D& rhs, const vector<int>& narrowBand) const
{
  assert(rhs.xRes() == _xRes);
  assert(rhs.yRes() == _yRes);
  assert(rhs.zRes() == _zRes);

  FIELD_3D result(rhs);
  result.clear();

  for (unsigned int index = 0; index < narrowBand.size(); index += 3)
  {
    int x = narrowBand[index];
    int y = narrowBand[index + 1];
    int z = narrowBand[index + 2];

    result(x,y,z) = (*this)(x,y,z) - rhs(x,y,z);
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
// do a diff between two fields, but only along a narrow band, and 
// not along borders where the filter skipped
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::filteredNarrowBandDiff(const FIELD_3D& rhs, const vector<int>& narrowBand, const FIELD_3D& filter)
{
  assert(rhs.xRes() == _xRes);
  assert(rhs.yRes() == _yRes);
  assert(rhs.zRes() == _zRes);

  FIELD_3D result(rhs);
  result.clear();

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

  // find which cells are inside the band
  vector<int> filteredNarrowBand;
  for (unsigned int index = 0; index < narrowBand.size(); index += 3)
  {
    int x = narrowBand[index];
    int y = narrowBand[index + 1];
    int z = narrowBand[index + 2];

    if (x < xHalf) continue;
    if (x > _xRes - xHalf - 1) continue; 
    if (y < yHalf) continue;
    if (y > _yRes - yHalf - 1) continue; 
    if (z < zHalf) continue;
    if (z > _zRes - zHalf - 1) continue;

    filteredNarrowBand.push_back(x);
    filteredNarrowBand.push_back(y);
    filteredNarrowBand.push_back(z);
  }

  for (unsigned int index = 0; index < filteredNarrowBand.size(); index += 3)
  {
    int x = filteredNarrowBand[index];
    int y = filteredNarrowBand[index + 1];
    int z = filteredNarrowBand[index + 2];

    result(x,y,z) = (*this)(x,y,z) - rhs(x,y,z);
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
// return a field showing the narrow band
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::narrowBandField(int bandwidth) const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D result(*this);
  result.clear();
  REAL invDx = 1.0 / _dx;

  // find which cells are inside the band
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        REAL currentDistance = fabs( (*this)(x,y,z) * invDx);
        if (currentDistance < bandwidth)
          result(x,y,z) = 1;
      }

  return result;
}

///////////////////////////////////////////////////////////////////////
// stomp everything outside of the narrow band
///////////////////////////////////////////////////////////////////////
void FIELD_3D::stompOutsideNarrowBand(const vector<int>& narrowBand)
{
  // build a lookup table
  map<int, bool> bandHash;
  for (unsigned int x = 0; x < narrowBand.size(); x++)
    bandHash[narrowBand[x]] = true;

  // do a lookup before stomping each entry
  for (int x = 0; x < _totalCells; x++)
    if (bandHash.find(x) == bandHash.end())
      _data[x] = 0;
}

///////////////////////////////////////////////////////////////////////
// do a WENO6 interpolation
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::sexticInterp(const REAL interp, const REAL* points)
{
  const REAL fm2 = points[0];
  const REAL fm1 = points[1];
  const REAL f   = points[2];
  const REAL fp1 = points[3];
  const REAL fp2 = points[4];
  const REAL fp3 = points[5];

  const REAL frac = 1.0 / 180.0;

  const REAL fp1Sq = fp1 * fp1;
  const REAL fp2Sq = fp2 * fp2;
  const REAL fSq = f * f;
  const REAL fm1Sq = fm1 * fm1;
  const REAL ffp1 = f * fp1;
  const REAL ffm1 = f * fm1;
  const REAL fp1fm1 = fp1 * fm1;
  const REAL fp1fp2 = fp1 * fp2;
  const REAL eps = 1e-6;

  const REAL IS1 = (814  * fp1Sq 
             + 4326 * fSq 
             + 2976 * fm1Sq 
             + 244  * fm2 * fm2 
             - 3579 * ffp1 
             - 6927 * ffm1 
             + 1854 * f * fm2 
             + 2634 * fp1fm1
             - 683  * fp1 * fm2
             - 1659 * fm1 * fm2) * frac + eps;

  const REAL IS2 = (1986 * fp1Sq 
             + 1986 * fSq
             + 244  * fm1Sq
             + 244  * fp2Sq
             + 1074 * f * fp2
             - 3777 * ffp1
             - 1269 * ffm1
             + 1074 * fp1fm1
             - 1269 * fp1fp2
             - 293  * fp2 * fm1) * frac + eps;

  const REAL IS3 = (814 * fSq
             + 4326 * fp1Sq
             + 2976 * fp2Sq
             + 244  * fp3 * fp3
             - 683  * f * fp3
             + 2634 * f * fp2
             - 3579 * ffp1
             - 6927 * fp1fp2
             + 1854 * fp1 * fp3
             - 1659 * fp2 * fp3) * frac + eps;

  const REAL interpPlusOne = interp + 1;
  const REAL interpPlusTwo = interp + 2;
  const REAL C1 = (2.0 - interp) * (3.0 - interp) * 0.05;
  // note the double negatives have been handled below
  const REAL C2 = (3.0 - interp) * (interp + 2.0) * 0.1;
  const REAL C3 = (interpPlusTwo) * (interpPlusOne)  * 0.05;

  const REAL alpha1 = C1 / (IS1 * IS1);
  const REAL alpha2 = C2 / (IS2 * IS2);
  const REAL alpha3 = C3 / (IS3 * IS3);

  const REAL sum = 1.0 / (alpha1 + alpha2 + alpha3);
  const REAL w1 = alpha1 * sum;
  const REAL w2 = alpha2 * sum;
  const REAL w3 = alpha3 * sum;

  const REAL interpMinusOne = interp - 1;
  const REAL interpMinusTwo = interp - 2;
  const REAL sixth = 1.0 / 6.0;
  const REAL p1 = fm2
            + (fm1 - fm2) * interpPlusTwo
            + (f - 2.0 * fm1 + fm2) * interpPlusTwo * interpPlusOne * 0.5
            + (fp1 - 3.0 * f + 3.0 * fm1 - fm2) * interpPlusTwo * interpPlusOne * interp * sixth;

  const REAL p2 = fm1
            + (f - fm1) * (interpPlusOne)
            + (fp1 - 2.0 * f + fm1) * (interpPlusOne) * interp * 0.5
            + (fp2 - 3.0 * fp1 + 3.0 * f - fm1) * (interpPlusOne) * interp * (interpMinusOne) * sixth;

  const REAL p3 = f
            + (fp1 - f) * interp
            + (fp2 - 2.0 * fp1 + f) * interp * (interpMinusOne) * 0.5
            + (fp3 - 3.0 * fp2 + 3.0 * fp1 - f) * interp * (interpMinusOne) * (interpMinusTwo) * sixth;

  const REAL result = w1 * p1 + w2 * p2 + w3 * p3;

  return result;
}

///////////////////////////////////////////////////////////////////////
// do a WENO6 interpolation
///////////////////////////////////////////////////////////////////////
REAL FIELD_3D::sexticInterpClamped(const REAL interp, const REAL* points)
{
  const REAL fm2 = points[0];
  const REAL fm1 = points[1];
  const REAL f   = points[2];
  const REAL fp1 = points[3];
  const REAL fp2 = points[4];
  const REAL fp3 = points[5];

  const REAL frac = 1.0 / 180.0;

  const REAL fp1Sq = fp1 * fp1;
  const REAL fp2Sq = fp2 * fp2;
  const REAL fSq = f * f;
  const REAL fm1Sq = fm1 * fm1;
  const REAL ffp1 = f * fp1;
  const REAL ffm1 = f * fm1;
  const REAL fp1fm1 = fp1 * fm1;
  const REAL fp1fp2 = fp1 * fp2;
  const REAL eps = 1e-6;

  const REAL IS1 = (814  * fp1Sq 
             + 4326 * fSq 
             + 2976 * fm1Sq 
             + 244  * fm2 * fm2 
             - 3579 * ffp1 
             - 6927 * ffm1 
             + 1854 * f * fm2 
             + 2634 * fp1fm1
             - 683  * fp1 * fm2
             - 1659 * fm1 * fm2) * frac + eps;

  const REAL IS2 = (1986 * fp1Sq 
             + 1986 * fSq
             + 244  * fm1Sq
             + 244  * fp2Sq
             + 1074 * f * fp2
             - 3777 * ffp1
             - 1269 * ffm1
             + 1074 * fp1fm1
             - 1269 * fp1fp2
             - 293  * fp2 * fm1) * frac + eps;

  const REAL IS3 = (814 * fSq
             + 4326 * fp1Sq
             + 2976 * fp2Sq
             + 244  * fp3 * fp3
             - 683  * f * fp3
             + 2634 * f * fp2
             - 3579 * ffp1
             - 6927 * fp1fp2
             + 1854 * fp1 * fp3
             - 1659 * fp2 * fp3) * frac + eps;

  const REAL interpPlusOne = interp + 1;
  const REAL interpPlusTwo = interp + 2;
  const REAL C1 = (2.0 - interp) * (3.0 - interp) * 0.05;
  // note the double negatives have been handled below
  const REAL C2 = (3.0 - interp) * (interp + 2.0) * 0.1;
  const REAL C3 = (interpPlusTwo) * (interpPlusOne)  * 0.05;

  const REAL alpha1 = C1 / (IS1 * IS1);
  const REAL alpha2 = C2 / (IS2 * IS2);
  const REAL alpha3 = C3 / (IS3 * IS3);

  const REAL sum = 1.0 / (alpha1 + alpha2 + alpha3);
  const REAL w1 = alpha1 * sum;
  const REAL w2 = alpha2 * sum;
  const REAL w3 = alpha3 * sum;

  const REAL interpMinusOne = interp - 1;
  const REAL interpMinusTwo = interp - 2;
  const REAL sixth = 1.0 / 6.0;
  const REAL p1 = fm2
            + (fm1 - fm2) * interpPlusTwo
            + (f - 2.0 * fm1 + fm2) * interpPlusTwo * interpPlusOne * 0.5
            + (fp1 - 3.0 * f + 3.0 * fm1 - fm2) * interpPlusTwo * interpPlusOne * interp * sixth;

  const REAL p2 = fm1
            + (f - fm1) * (interpPlusOne)
            + (fp1 - 2.0 * f + fm1) * (interpPlusOne) * interp * 0.5
            + (fp2 - 3.0 * fp1 + 3.0 * f - fm1) * (interpPlusOne) * interp * (interpMinusOne) * sixth;

  const REAL p3 = f
            + (fp1 - f) * interp
            + (fp2 - 2.0 * fp1 + f) * interp * (interpMinusOne) * 0.5
            + (fp3 - 3.0 * fp2 + 3.0 * fp1 - f) * interp * (interpMinusOne) * (interpMinusTwo) * sixth;

  const REAL result = w1 * p1 + w2 * p2 + w3 * p3;

  const REAL intervalMax = (points[3] > points[2]) ? points[3] : points[2];
  const REAL intervalMin = (points[3] > points[2]) ? points[2] : points[3];

  return (result < intervalMax) ? ((result > intervalMin) ? result : intervalMin) : intervalMax;
}

///////////////////////////////////////////////////////////////////////
// stomp everything within the grid "corners" assuming a padding 
// thickness of "bordersize"
///////////////////////////////////////////////////////////////////////
void FIELD_3D::stompCorners(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
      
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int i = 1; i <= borderSize; i++)
        (*this)(borderSize - i,y,z) = 0;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int i = 1; i <= borderSize; i++)
        (*this)(_xRes - 1 - borderSize + i,y,z) = 0;

  for (int z = 0; z < _zRes; z++)
    for (int x = 0; x < _xRes; x++)
      for (int i = 1; i <= borderSize; i++)
        (*this)(x, borderSize - i,z) = 0;

  for (int z = 0; z < _zRes; z++)
    for (int x = 0; x < _xRes; x++)
      for (int i = 1; i <= borderSize; i++)
        (*this)(x,_yRes - 1 - borderSize + i,z) = 0;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      for (int i = 1; i <= borderSize; i++)
        (*this)(x,y,borderSize - i) = 0;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      for (int i = 1; i <= borderSize; i++)
        (*this)(x,y,_zRes - 1 - borderSize +i) = 0;
}

///////////////////////////////////////////////////////////////////////
// set everything in the specified interval to a given value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setInterval(const int xMin, const int xMax, const int yMin, const int yMax, const int zMin, const int zMax, const REAL value)
{
  assert(xMin >= 0);
  assert(xMin <= _xRes);
  assert(xMax >= 0);
  assert(xMax <= _xRes);

  assert(yMin >= 0);
  assert(yMin <= _yRes);
  assert(yMax >= 0);
  assert(yMax <= _yRes);

  assert(zMin >= 0);
  assert(zMin <= _zRes);
  assert(zMax >= 0);
  assert(zMax <= _zRes);

  for (int x = xMin; x < xMax; x++)
    for (int y = yMin; y < yMax; y++)
      for (int z = zMin; z < zMax; z++)
        (*this)(x,y,z) = value;
}

} // HOBAK

