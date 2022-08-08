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
#include "VECTOR3_FIELD_3D.h"
#include <omp.h>

#include <set>

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const int& xRes, const int& yRes, const int& zRes,
    const VECTOR3& center, const VECTOR3& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VECTOR3[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(double* data, const int& xRes, const int& yRes, const int& zRes,
    const VECTOR3& center, const VECTOR3& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VECTOR3[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] = data[3 * x];
    _data[x][1] = data[3 * x + 1];
    _data[x][2] = data[3 * x + 2];
  }
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const VECTOR3_FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths()), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VECTOR3[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths()), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VECTOR3[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D() :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL), _initialized(false)
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D::~VECTOR3_FIELD_3D()
{
  delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clear()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x].setZero();
}

///////////////////////////////////////////////////////////////////////
// reset the lengths to something else, and recompute all the 
// dimesions as well
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setLengths(const VECTOR3& lengths)
{
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

///////////////////////////////////////////////////////////////////////
// create a field of the grid positions of the passed in grid
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::cellCenters(const FIELD_3D& input)
{
  int xRes = input.xRes();
  int yRes = input.yRes();
  int zRes = input.zRes();
  const VECTOR3& center = input.center();
  const VECTOR3& lengths = input.lengths();
  VECTOR3_FIELD_3D final(xRes, yRes, zRes, center, lengths);

  int index = 0;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++)
        final[index] = input.cellCenter(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the gradient of a scalar field
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::gradient(const FIELD_3D& input)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = input.xRes();
  const int yRes = input.yRes();
  const int zRes = input.zRes();
  const int slabSize = xRes * yRes;
  const int totalCells = xRes * yRes * zRes;
  const VECTOR3& center = input.center();
  const VECTOR3& lengths = input.lengths();
  VECTOR3_FIELD_3D final(xRes, yRes, zRes, center, lengths);
 
  const REAL dxHalfInv = 0.5 / input.dx();
  const REAL dyHalfInv = 0.5 / input.dy();
  const REAL dzHalfInv = 0.5 / input.dz();

  // do the x middles
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 1; x < xRes - 1; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][0] = (input[index + 1] - input[index - 1]) * dxHalfInv;
      }

  // do the y middles
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 1; y < yRes - 1; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][1] = (input[index + xRes] - input[index - xRes]) * dyHalfInv;
      }

  // do the z middles
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 1; z < zRes - 1; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][2] = (input[index + slabSize] - input[index - slabSize]) * dzHalfInv;
      }

  // reset dx's to a single cell
  REAL dxInv = 1.0 / input.dx();
  REAL dyInv = 1.0 / input.dy();
  REAL dzInv = 1.0 / input.dz();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
    {
      // front slab
      int index = x + y * xRes;
      final[index][2] = (input[index + slabSize] - input[index]) * dzInv;

      // back slab
      index += totalCells - slabSize;
      final[index][2] = (input[index] - input[index - slabSize]) * dzInv;
    }

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int x = 0; x < xRes; x++)
    {
      // bottom slab
      int index = x + z * slabSize;
      final[index][1] = (input[index + xRes] - input[index]) * dyInv;

      // top slab
      index += slabSize - xRes;
      final[index][1] = (input[index] - input[index - xRes]) * dyInv;
    }

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
    {
      // left slab
      int index = y * xRes + z * slabSize;
      final[index][0] = (input[index + 1] - input[index]) * dxInv;

      // right slab
      index += xRes - 1;
      final[index][0] = (input[index] - input[index - 1]) * dxInv;
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::scalarField(int component) const
{
  assert(component >= 0 && component < 3);

  FIELD_3D final(_xRes, _yRes, _zRes, _center, _lengths);

  for (int x = 0; x < _totalCells; x++)
    final[x] = _data[x][component];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::magnitudeField() const
{
  FIELD_3D final(_xRes, _yRes, _zRes, _center, _lengths);

  for (int x = 0; x < _totalCells; x++)
    final[x] = _data[x].norm();

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x].transpose() * v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator*(const FIELD_3D& u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(v);

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] * v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// sum two vector fields
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator+(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());
  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] + v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// diff two vector fields
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator-(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());
  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] - v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// return a grid of values at the given spatial positions
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::compose(const FIELD_3D& values, const VECTOR3_FIELD_3D& positions)
{
  assert(positions.xRes() == values.xRes());
  assert(positions.yRes() == values.yRes());
  assert(positions.zRes() == values.zRes());

  // intialize to the same dims as the input
  FIELD_3D final(values);

  int i = 0;
  for (int z = 0; z < values.zRes(); z++)
    for (int y = 0; y < values.yRes(); y++)
      for (int x = 0; x < values.xRes(); x++, i++)
        final[i] = values(positions[i]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// return a grid of values at the given spatial positions
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::compose(const VECTOR3_FIELD_3D& values, const VECTOR3_FIELD_3D& positions)
{
  assert(positions.xRes() == values.xRes());
  assert(positions.yRes() == values.yRes());
  assert(positions.zRes() == values.zRes());

  // intialize to the same dims as the input
  VECTOR3_FIELD_3D final(values);

  int i = 0;
  for (int z = 0; z < values.zRes(); z++)
    for (int y = 0; y < values.yRes(); y++)
      for (int x = 0; x < values.xRes(); x++, i++)
        final[i] = values(positions[i]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
const VECTOR3 VECTOR3_FIELD_3D::operator()(const VECTOR3& position) const
{
  /*
  VECTOR3 positionCopy = position;

  // get the lower corner position
  VECTOR3 corner = _center - (REAL)0.5 * _lengths;
  VECTOR3 dxs(_dx, _dy, _dz);
  corner += (REAL)0.5 * dxs;

  // recenter position
  positionCopy -= corner;
  */
  VECTOR3 positionCopy = position - _center + (REAL)0.5 * _lengths - (REAL)0.5 * dxs();

  positionCopy[0] *= 1.0 / _dx;
  positionCopy[1] *= 1.0 / _dy;
  positionCopy[2] *= 1.0 / _dz;

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
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
VECTOR3 VECTOR3_FIELD_3D::debugPositionOperator(const VECTOR3& position) const
{
  /*
  VECTOR3 positionCopy = position;

  // get the lower corner position
  VECTOR3 corner = _center - (REAL)0.5 * _lengths;
  VECTOR3 dxs(_dx, _dy, _dz);
  corner += (REAL)0.5 * dxs;

  // recenter position
  positionCopy -= corner;
  */
  VECTOR3 positionCopy = position - _center + (REAL)0.5 * _lengths - (REAL)0.5 * dxs();

  positionCopy[0] *= 1.0 / _dx;
  positionCopy[1] *= 1.0 / _dy;
  positionCopy[2] *= 1.0 / _dz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  cout << " recomputed fine: " << x0 << " " << y0 << " " << z0 << endl;
  cout << " recomputed dx: " << _lengths[0] / _xRes << endl;
  cout << " length: " << _lengths << endl;
  cout << " center: " << _center << endl;

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
// norms
///////////////////////////////////////////////////////////////////////
REAL VECTOR3_FIELD_3D::sumMagnitudes()
{
  REAL final = 0;
  for (int i = 0; i < _totalCells; i++)
    final += _data[i].norm();

  return final;
}

///////////////////////////////////////////////////////////////////////
// norms
///////////////////////////////////////////////////////////////////////
REAL VECTOR3_FIELD_3D::maxMagnitudes()
{
  REAL final = _data[0].norm();
  for (int i = 0; i < _totalCells; i++)
    if (_data[i].norm() > final)
      final = _data[i].norm();

  return final;
}

///////////////////////////////////////////////////////////////////////
// check if any entry is a nan
///////////////////////////////////////////////////////////////////////
bool VECTOR3_FIELD_3D::isNan()
{
  for (int x = 0; x < _totalCells; x++)
    for (int y = 0; y < 3; y++)
      if (isnan(_data[x][y]))
        return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator-=(const VECTOR3_FIELD_3D& input)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator+=(const VECTOR3_FIELD_3D& input)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator*=(const REAL& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator=(const REAL& value)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = VECTOR3::Constant(value);

  return *this;
}

///////////////////////////////////////////////////////////////////////
// extend some vector quantity off of a front, given a signed distance function
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::fastExtension(const FIELD_3D& signedDistance)
{
  cout << " Extending vectors ... "; flush(cout);
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
// insert the front in preparation for reinitialization or extension
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::insertFront(const bool forward, FIELD_3D& distance, MIN_HEAP& minHeap)
{
  REAL outside = distance.outside();

  // insert the ones for the forward marching
  minHeap.clear();
  for (int i = 0; i < _totalCells; i++)
  {
    bool compare = (forward) ? distance[i] >= 0.0f && distance[i] < outside
                             : distance[i] <= 0.0f;

    REAL sum = 0.0f;
    if (compare)
    {
      int zIndex = i / _slabSize;
      int yIndex = (i - zIndex * _slabSize) / _xRes;
      int xIndex = (i - zIndex * _slabSize - yIndex * _xRes);
      REAL center = distance[i];

      // x plus and minus
      REAL xPlus  = (xIndex < _xRes - 1) ? distance[i + 1] : outside;
      REAL xMinus = (xIndex > 0)         ? distance[i - 1] : outside;

      // y plus and minus
      REAL yPlus  = (yIndex < _yRes - 1) ? distance[i + _xRes] : outside;
      REAL yMinus = (yIndex > 0)         ? distance[i - _xRes] : outside;

      // z plus and minus
      REAL zPlus  = (zIndex < _zRes - 1) ? distance[i + _slabSize] : outside;
      REAL zMinus = (zIndex > 0)         ? distance[i - _slabSize] : outside;

      REAL interpolate;
      compare = (forward) ? true : yPlus < outside;
      if (yPlus * center <= 0.0f && compare)
      {
        interpolate = center / (center - yPlus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : yMinus < outside;
      if (yMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - yMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : xMinus < outside;
      if (xMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - xMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : xPlus < outside;
      if (xPlus * center <= 0.0f && compare)
      {
        interpolate = center / (center - xPlus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : zMinus < outside;
      if (zMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - zMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : zPlus < outside;
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
        distance[i] = sign * outside;
    }
    if (!forward && distance[i] >= outside)
      distance[i] = -outside;
  }
}

//////////////////////////////////////////////////////////////////////
// do extension in one direction
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::extendOneway(bool forward, FIELD_3D& distance, MIN_HEAP& minHeap)
{
  int pops = 0;
  REAL outside = distance.outside();
  
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
        VECTOR3 extensions[3];
        REAL neighbors[6];
        VECTOR3 extensionNeighbors[6];

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
          if (forward && distances[x] < 0.0f)  distances[x] = outside;
          if (!forward && distances[x] > 0.0f) distances[x] = -outside;
          
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
            distances[x] = (distances[x] >= 0.0f) ? distances[x] : outside;
          else
            distances[x] = (distances[x] <= 0.0f) ? -distances[x] : outside;
      
        // do a 3-element bubble sort 
        REAL temp;
        VECTOR3 tempPoint;
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
        float b[3];
        float c[3];
        b[0] = distances[0];
        c[0] = distances[0] * distances[0] - 1.0f;
        for (int x = 1; x < 3; x++)
        {
          b[x] = distances[x] + b[x-1];
          c[x] = distances[x] * distances[x] + c[x-1];
        }

        // solve for the right one
        int i = 2;
        float discrim = b[i] * b[i] - (float)(i + 1) * c[i];
        float newDist = (b[i] + sqrtf(discrim)) / (float)(i + 1);
        while ((discrim < 0.0f || newDist < distances[i]) && i != -1)
        {
          i--;
          discrim = b[i] * b[i] - (float)(i + 1) * c[i];
          
          if (discrim > 0.0f)
            newDist = (b[i] + sqrtf(discrim)) / (float)(i + 1);
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
        if (newDist < outside)
        {
          float minus = forward ? 1.0f : -1.0f;

          // if it has never been on the heap
          if (fabs(distance[candidate]) >= outside)
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

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator=(const VECTOR3_FIELD_3D& input)
{
  if (input.xRes() != _xRes ||
      input.yRes() != _yRes ||
      input.zRes() != _zRes)
  {
    delete[] _data;

    _xRes = input.xRes();
    _yRes = input.yRes();
    _zRes = input.zRes();

    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;
    _data = new VECTOR3[_totalCells];
  }

  _center = input.center();
  _lengths = input.lengths();

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] = input[x];

  _initialized = input._initialized;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// set the values in the field to the values at the closest points
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setToClosestPointValues(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(input);

  //int before = input.quinticClamps();
  cout << " Using quartic ... "; flush(cout);
  //cout << " Using quartic clamped ... "; flush(cout);
  //cout << " Using sextic ... "; flush(cout);
  //cout << " Using sextic clamped ... "; flush(cout);
  //cout << " Using quinitic ... "; flush(cout);
  //cout << " Using cubic ... "; flush(cout);
  //cout << " Using linear ... "; flush(cout);

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < closestPoints.zRes(); z++)
    for (int y = 0; y < closestPoints.yRes(); y++)
      for (int x = 0; x < closestPoints.xRes(); x++)
        //final(x,y,z) = input.lerpDebug(closestPoints(x,y,z), x,y,z);
        //final(x,y,z) = input(closestPoints(x,y,z));
        //final(x,y,z) = input.cubicLookup(closestPoints(x,y,z));
        //final(x,y,z) = input.quinticLookup(closestPoints(x,y,z));
        final(x,y,z) = input.quarticLookup(closestPoints(x,y,z));
        //final(x,y,z) = input.quarticLookupClamped(closestPoints(x,y,z));
        //final(x,y,z) = input.sexticLookup(closestPoints(x,y,z));
        //final(x,y,z) = input.sexticLookupClamped(closestPoints(x,y,z));
        //final(x,y,z) = input.cubicNewtonLookup(closestPoints(x,y,z));

  /*
  int totalCalls = 21 * closestPoints.totalCells();
  int after = input.quinticClamps();
  int clamps = after - before;

  cout << " Quintic clamps: " << clamps << " of " << totalCalls << " " << 100.0 * clamps / (REAL)totalCalls << "%" << endl;
  */

  return final;
}

//////////////////////////////////////////////////////////////////////
// set the values in the field to the values at the closest points
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBand(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int maxCells)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Setting to closest points, narrow banded ..."; flush(cout);
  FIELD_3D final(input);
  REAL invDx = 1.0 / distance.dx();

  const int xRes = closestPoints.xRes();
  const int yRes = closestPoints.yRes();
  const int zRes = closestPoints.zRes();

  // find which cells are inside the band
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        REAL currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance >= maxCells) continue;

        final(x,y,z) = input.quarticLookup(closestPoints(x,y,z));
      }
  cout << " done." << endl;

  /*
  // find which cells are inside the band
  vector<int> xs;
  vector<int> ys;
  vector<int> zs;
  for (int z = 0; z < closestPoints.zRes(); z++)
    for (int y = 0; y < closestPoints.yRes(); y++)
      for (int x = 0; x < closestPoints.xRes(); x++)
      {
        REAL currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance < maxCells)
        {
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
        }
      }

  const int size = xs.size();
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int i = 0; i < size; i++)
  {
    int x = xs[i];
    int y = ys[i];
    int z = zs[i];
    //final(x,y,z) = input(closestPoints(x,y,z));
    //final(x,y,z) = input.cubicLookup(closestPoints(x,y,z));
    final(x,y,z) = input.quarticLookup(closestPoints(x,y,z));
    //final(x,y,z) = input.quarticLookupInlined(closestPoints(x,y,z));
    //final(x,y,z) = input.quinticLookup(closestPoints(x,y,z));
    //final(x,y,z) = input.cubicNewtonLookup(closestPoints(x,y,z));
  }
  cout << " done." << endl;
  */

  return final;
}

//////////////////////////////////////////////////////////////////////
// set the values in the field to the values at the closest points,
// but exclude a central core of cells, so it's a hose, not a band
//
// Might be better to called it a "cored band"
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int maxRadius, int coreRadius)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Setting to closest points, narrow band, frozen core..."; flush(cout);
  FIELD_3D final(input);
  //const REAL invDx = 1.0 / distance.dx();
  const REAL invDx = 1.0 / input.dx();

  const int xRes = closestPoints.xRes();
  const int yRes = closestPoints.yRes();
  const int zRes = closestPoints.zRes();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        REAL currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance >= maxRadius || currentDistance <= coreRadius) continue;
        
        final(x,y,z) = input.quarticLookup(closestPoints(x,y,z));
        //final(x,y,z) = input.sexticLookup(closestPoints(x,y,z));
      }
  cout << " done." << endl;

  return final;
  /*
  TIMER functionTimer(__FUNCTION__);
  cout << " Setting to closest points, narrow band, frozen core..."; flush(cout);
  FIELD_3D final(input);
  REAL invDx = 1.0 / distance.dx();

  // find which cells are inside the band
  vector<int> xs;
  vector<int> ys;
  vector<int> zs;
  for (int z = 0; z < closestPoints.zRes(); z++)
    for (int y = 0; y < closestPoints.yRes(); y++)
      for (int x = 0; x < closestPoints.xRes(); x++)
      {
        REAL currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance < maxRadius && currentDistance > coreRadius)
        {
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
        }
      }

  int size = xs.size();
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int i = 0; i < size; i++)
  {
    int x = xs[i];
    int y = ys[i];
    int z = zs[i];
    //final(x,y,z) = input(closestPoints(x,y,z));
    //final(x,y,z) = input.cubicLookup(closestPoints(x,y,z));
    final(x,y,z) = input.quarticLookup(closestPoints(x,y,z));
    //final(x,y,z) = input.quarticLookupInlined(closestPoints(x,y,z));
    //final(x,y,z) = input.quinticLookup(closestPoints(x,y,z));
    //final(x,y,z) = input.cubicNewtonLookup(closestPoints(x,y,z));
  }
  cout << " done." << endl;

  return final;
  */
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setToClosestPointValuesFrozenCore(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int coreRadius)
{
  TIMER frozenCore("Building frozen core");

  cout << " Setting to closest points, narrow band, frozen core radius " << coreRadius << " ..."; flush(cout);
  FIELD_3D final(input);
  REAL invDx = 1.0 / distance.dx();

  // find which cells are inside the band
  vector<int> xs;
  vector<int> ys;
  vector<int> zs;
  for (int z = 0; z < closestPoints.zRes(); z++)
    for (int y = 0; y < closestPoints.yRes(); y++)
      for (int x = 0; x < closestPoints.xRes(); x++)
      {
        REAL currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance > coreRadius)
        {
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
        }
      }

  TIMER functionTimer(__FUNCTION__);
  int size = xs.size();
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int i = 0; i < size; i++)
  {
    int x = xs[i];
    int y = ys[i];
    int z = zs[i];
    //final(x,y,z) = input(closestPoints(x,y,z));
    //final(x,y,z) = input.cubicLookup(closestPoints(x,y,z));
    final(x,y,z) = input.quarticLookup(closestPoints(x,y,z));
    //final(x,y,z) = input.quarticLookupInlined(closestPoints(x,y,z));
    //final(x,y,z) = input.quinticLookup(closestPoints(x,y,z));
    //final(x,y,z) = input.cubicNewtonLookup(closestPoints(x,y,z));
  }
  cout << " done." << endl;

  return final;
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::write(const string filename) const
{
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " VECTOR3_FIELD_3D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  write(file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::read(FILE* file, VECTOR3& v)
{
  double data[3];
  fread((void*)data, sizeof(double), 3, file);
  v[0] = data[0];
  v[1] = data[1];
  v[2] = data[2];
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::write(FILE* file, const VECTOR3& v) const
{
  double data[3] = {v[0], v[1], v[2]};
  fwrite((void*)data, sizeof(double), 3, file);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::write(FILE* file) const
{
  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);
  fwrite((void*)&_zRes, sizeof(int), 1, file);
  //_center.write(file);
  //_lengths.write(file);
  write(file, _center);
  write(file, _lengths);

  double* dataDouble = new double[3 * _totalCells];
  for (int x = 0; x < _totalCells; x++)
  {
    dataDouble[3 * x] = _data[x][0];
    dataDouble[3 * x + 1] = _data[x][1];
    dataDouble[3 * x + 2] = _data[x][2];
  }

  fwrite((void*)dataDouble, sizeof(double), 3 * _totalCells, file);
  delete[] dataDouble;
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::read(const string filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " VECTOR3_FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  read(file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::read(FILE* file)
{
  // read dimensions
  fread((void*)&_xRes, sizeof(int), 1, file);
  fread((void*)&_yRes, sizeof(int), 1, file);
  fread((void*)&_zRes, sizeof(int), 1, file);
  //_center.read(file);
  //_lengths.read(file);
  read(file, _center);
  read(file, _lengths);
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  if (_data) delete[] _data;
  _data = new VECTOR3[_totalCells];

  double* dataDouble = new double[3 * _totalCells];
  fread((void*)dataDouble, sizeof(double), 3 * _totalCells, file);
  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0]= dataDouble[3 * x];
    _data[x][1]= dataDouble[3 * x + 1];
    _data[x][2]= dataDouble[3 * x + 2];
  }
  delete[] dataDouble;
  _initialized = true;
}

/*
//////////////////////////////////////////////////////////////////////
// draw to GL
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::draw(const VECTOR3_FIELD_3D& origins) const
{
  int stride = 6;

  glPointSize(1);

  glBegin(GL_POINTS);
  glColor4f(1,0,0,1);
  for (int z = 0; z < _zRes; z+=stride)
    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VECTOR3& origin = origins(x,y,z);
        glVertex3f(origin[0], origin[1], origin[2]);
      }
  glEnd();

  glBegin(GL_POINTS);
  glColor4f(0,0,1,1);
  for (int z = 0; z < _zRes; z+=stride)
    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VECTOR3& endpoint = (*this)(x,y,z);
        glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
      }
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// draw to GL
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::drawZSlice(const VECTOR3_FIELD_3D& origins, const int zSlice, const REAL scale, const int stride) const
{
  float radius = 0.001;
  
  glBegin(GL_LINES);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VECTOR3& origin = origins(x,y,z);
        VECTOR3 endpoint = origin + scale * (*this)(x,y,z);

        glColor4f(1,1,1,1);
        glVertex3f(origin[0], origin[1], origin[2]);
        glColor4f(0,0,0,0);
        glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
      }
  }
  glEnd();

  glColor4f(1,0,0,1);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VECTOR3& origin = origins(x,y,z);
        glPushMatrix();
        glTranslatef(origin[0], origin[1], origin[2]);
        glutSolidSphere(radius, 10,10);
        glPopMatrix();
      }
  }
}

//////////////////////////////////////////////////////////////////////
// draw to GL, but assume the field represents closest points,
// so don't add the origin to the endpoint
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::drawClosestPointZSlice(const VECTOR3_FIELD_3D& origins, const int zSlice, const REAL scale, const int stride) const
{
  float radius = 0.00125;
  
  glBegin(GL_LINES);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VECTOR3& origin = origins(x,y,z);
        VECTOR3 endpoint = (*this)(x,y,z);

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

  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VECTOR3& origin = origins(x,y,z);
        const VECTOR3 endpoint = (*this)(x,y,z);

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
*/

//////////////////////////////////////////////////////////////////////
// dump to a viewer
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::fieldViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name)
{
  field.write("temp3d.vector.field");
  distanceField.write("temp3d.field");
  string execute("./vectorFieldViewer3D temp3d.vector.field temp3d.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

//////////////////////////////////////////////////////////////////////
// dump to a viewer
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::overlayFieldViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name)
{
  field.write("temp3d.vector.field");
  distanceField.write("temp3d.field");
  string execute("./overlayVectorFieldViewer3D temp3d.vector.field temp3d.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

//////////////////////////////////////////////////////////////////////
// dump to a viewer
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::closestPointViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name)
{
  field.write("temp3d.vector.field");
  distanceField.write("temp3d.field");
  //string execute("./bin/closestPointViewer3D temp3d.vector.field temp3d.field \"");
  string execute("./closestPointViewer3D temp3d.vector.field temp3d.field \"");
  execute = execute + name + string("\" &");
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// flip the x and y coordinates
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::flipXY() const
{
  VECTOR3_FIELD_3D final(_yRes, _xRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(y,x,z) = (*this)(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// flip the x and y coordinates
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::flipXZ() const
{
  VECTOR3_FIELD_3D final(_zRes, _yRes, _xRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        final(z,y,x) = (*this)(x,y,z);

        REAL temp = final(z,y,x)[0];
        final(z,y,x)[0] = final(z,y,x)[2];
        final(z,y,x)[2] = temp;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// flip the z and y coordinates
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::flipZY() const
{
  VECTOR3_FIELD_3D final(_xRes, _zRes, _yRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,z,y) = (*this)(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advect(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();
 
 /* 
  const VECTOR3 corner = oldField.center() - (REAL)0.5 * oldField.lengths() + (REAL)0.5 * oldField.dxs();
  const VECTOR3 dxs = oldField.dxs();

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        VECTOR3 position(x - dt * velocity[0], y - dt * velocity[1], z - dt * velocity[2]);
        position[0] *= dxs[0];
        position[1] *= dxs[1];
        position[2] *= dxs[2];
        position += corner;
        newField[index] = oldField.quarticLookup(position);
      }
      */

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        //const VECTOR3 velocity = velocityGrid[index];
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        // get interpolation weights
        const REAL s1 = xTrace - x0;
        const REAL s0 = 1.0f - s1;
        const REAL t1 = yTrace - y0;
        const REAL t0 = 1.0f - t1;
        const REAL u1 = zTrace - z0;
        const REAL u0 = 1.0f - u1;

        const int i000 = x0 + y0 * xRes + z0 * slabSize;
        const int i010 = x0 + y1 * xRes + z0 * slabSize;
        const int i100 = x1 + y0 * xRes + z0 * slabSize;
        const int i110 = x1 + y1 * xRes + z0 * slabSize;
        const int i001 = x0 + y0 * xRes + z1 * slabSize;
        const int i011 = x0 + y1 * xRes + z1 * slabSize;
        const int i101 = x1 + y0 * xRes + z1 * slabSize;
        const int i111 = x1 + y1 * xRes + z1 * slabSize;

        // interpolate
        // (indices could be computed once)
        newField[index] = u0 * (s0 * (t0 * oldField[i000] +
                                      t1 * oldField[i010]) +
                                s1 * (t0 * oldField[i100] +
                                      t1 * oldField[i110])) +
                          u1 * (s0 * (t0 * oldField[i001] +
                                      t1 * oldField[i011]) +
                                s1 * (t0 * oldField[i101] +
                                      t1 * oldField[i111]));
      }
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectNarrowBand(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField, const FIELD_3D& distance, const int maxCells)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();
  const REAL invDx = 1.0 / oldField.dx();
 
  const VECTOR3 corner = oldField.center() - (REAL)0.5 * oldField.lengths() + (REAL)0.5 * oldField.dxs();
  const VECTOR3 dxs = oldField.dxs();

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distance(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        VECTOR3 position(x - dt * velocity[0], y - dt * velocity[1], z - dt * velocity[2]);
        position[0] *= dxs[0];
        position[1] *= dxs[1];
        position[2] *= dxs[2];
        position += corner;
        newField[index] = oldField.quarticLookup(position);
      }

  /*
  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distance(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        // get interpolation weights
        const REAL s1 = xTrace - x0;
        const REAL s0 = 1.0f - s1;
        const REAL t1 = yTrace - y0;
        const REAL t0 = 1.0f - t1;
        const REAL u1 = zTrace - z0;
        const REAL u0 = 1.0f - u1;

        const int i000 = x0 + y0 * xRes + z0 * slabSize;
        const int i010 = x0 + y1 * xRes + z0 * slabSize;
        const int i100 = x1 + y0 * xRes + z0 * slabSize;
        const int i110 = x1 + y1 * xRes + z0 * slabSize;
        const int i001 = x0 + y0 * xRes + z1 * slabSize;
        const int i011 = x0 + y1 * xRes + z1 * slabSize;
        const int i101 = x1 + y0 * xRes + z1 * slabSize;
        const int i111 = x1 + y1 * xRes + z1 * slabSize;

        // interpolate
        // (indices could be computed once)
        newField[index] = u0 * (s0 * (t0 * oldField[i000] +
                                      t1 * oldField[i010]) +
                                s1 * (t0 * oldField[i100] +
                                      t1 * oldField[i110])) +
                          u1 * (s0 * (t0 * oldField[i001] +
                                      t1 * oldField[i011]) +
                                s1 * (t0 * oldField[i101] +
                                      t1 * oldField[i111]));
      }
      */
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectNarrowBandLinear(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField, const FIELD_3D& distance, const int maxCells)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();
  const REAL invDx = 1.0 / oldField.dx();
 
  const VECTOR3 corner = oldField.center() - (REAL)0.5 * oldField.lengths() + (REAL)0.5 * oldField.dxs();
  const VECTOR3 dxs = oldField.dxs();

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distance(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        VECTOR3 position(x - dt * velocity[0], y - dt * velocity[1], z - dt * velocity[2]);
        position[0] *= dxs[0];
        position[1] *= dxs[1];
        position[2] *= dxs[2];
        position += corner;
        //newField[index] = oldField.quarticLookup(position);
        newField[index] = oldField(position);
      }

  /*
  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distance(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        // get interpolation weights
        const REAL s1 = xTrace - x0;
        const REAL s0 = 1.0f - s1;
        const REAL t1 = yTrace - y0;
        const REAL t0 = 1.0f - t1;
        const REAL u1 = zTrace - z0;
        const REAL u0 = 1.0f - u1;

        const int i000 = x0 + y0 * xRes + z0 * slabSize;
        const int i010 = x0 + y1 * xRes + z0 * slabSize;
        const int i100 = x1 + y0 * xRes + z0 * slabSize;
        const int i110 = x1 + y1 * xRes + z0 * slabSize;
        const int i001 = x0 + y0 * xRes + z1 * slabSize;
        const int i011 = x0 + y1 * xRes + z1 * slabSize;
        const int i101 = x1 + y0 * xRes + z1 * slabSize;
        const int i111 = x1 + y1 * xRes + z1 * slabSize;

        // interpolate
        // (indices could be computed once)
        newField[index] = u0 * (s0 * (t0 * oldField[i000] +
                                      t1 * oldField[i010]) +
                                s1 * (t0 * oldField[i100] +
                                      t1 * oldField[i110])) +
                          u1 * (s0 * (t0 * oldField[i001] +
                                      t1 * oldField[i011]) +
                                s1 * (t0 * oldField[i101] +
                                      t1 * oldField[i111]));
      }
      */
}

///////////////////////////////////////////////////////////////////////
// do an advection of a subset of the grid
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::advectSubset(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, const vector<int>& subset)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D newField(oldField);
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();

#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < subset.size(); i++)
  {
    const int index = subset[i];
    const int z = index / slabSize;
    const int y = (index - z * slabSize) / xRes;
    const int x = index - z * slabSize - y * xRes;
    
    // backtrace
    const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
    REAL xTrace = x - dt * velocity[0];
    REAL yTrace = y - dt * velocity[1];
    REAL zTrace = z - dt * velocity[2];

    // clamp backtrace to grid boundaries
    xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
    xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
    yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
    yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
    zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
    zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

    // locate neighbors to interpolate
    const int x0 = (int)xTrace;
    const int x1 = x0 + 1;
    const int y0 = (int)yTrace;
    const int y1 = y0 + 1;
    const int z0 = (int)zTrace;
    const int z1 = z0 + 1;

    // get interpolation weights
    const REAL s1 = xTrace - x0;
    const REAL s0 = 1.0f - s1;
    const REAL t1 = yTrace - y0;
    const REAL t0 = 1.0f - t1;
    const REAL u1 = zTrace - z0;
    const REAL u0 = 1.0f - u1;

    const int z0slabSize = z0 * slabSize;
    const int z1slabSize = z1 * slabSize;
    const int y0xRes = y0 * xRes;
    const int y1xRes = y1 * xRes;

    const int i000 = x0 + y0xRes + z0slabSize;
    const int i010 = x0 + y1xRes + z0slabSize;
    const int i100 = x1 + y0xRes + z0slabSize;
    const int i110 = x1 + y1xRes + z0slabSize;
    const int i001 = x0 + y0xRes + z1slabSize;
    const int i011 = x0 + y1xRes + z1slabSize;
    const int i101 = x1 + y0xRes + z1slabSize;
    const int i111 = x1 + y1xRes + z1slabSize;

    // interpolate
    // (indices could be computed once)
    newField[index] = u0 * (s0 * (t0 * oldField[i000] +
                                  t1 * oldField[i010]) +
                            s1 * (t0 * oldField[i100] +
                                  t1 * oldField[i110])) +
                      u1 * (s0 * (t0 * oldField[i001] +
                                  t1 * oldField[i011]) +
                            s1 * (t0 * oldField[i101] +
                                  t1 * oldField[i111]));
  }
  return newField;
}

///////////////////////////////////////////////////////////////////////
// get a single closest point
///////////////////////////////////////////////////////////////////////
REAL VECTOR3_FIELD_3D::getClosestPointValue(const VECTOR3_FIELD_3D& targetGradient, const VECTOR3& startPosition, const FIELD_3D& distanceField, const FIELD_3D& toExtend)
{
  TIMER functionTimer(__FUNCTION__);

  int maxSteps = 100;
  REAL fictionalDt = 1.0;
  
  REAL diff = 1;
  int steps = 0;
  VECTOR3 position = startPosition;
  VECTOR3 targetDelta = targetGradient(position);

  while (diff > 1e-6 && steps < maxSteps)
  {
    REAL targetDistance = distanceField(position);
    diff = fabs(targetDistance);

    // if the gradient is zero, stop trying, because it probably
    // wandered off the grid
    //if (norm2(targetDelta) < 1e-6) break;
    if (targetDelta.norm() < 1e-6) break;

    // go ahead and always do first -- second gives the occasional wacky value
    // that throws off the stability of the simulation
    targetDelta.normalize();
    VECTOR3 move = (targetDistance) * targetDelta * fictionalDt;

    position = position - move;

    // update the direction
    targetDelta = targetGradient(position);
    steps++;
  }

  // do the extension here
  return toExtend.quarticLookup(position);
}

///////////////////////////////////////////////////////////////////////
// compute the extension field for some subset of points
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::computeExtensionFieldSubset(const FIELD_3D& distanceField, const FIELD_3D& toExtend, const vector<int>& subset)
{
  TIMER functionTimer(__FUNCTION__);

  // get the cell centers 
  VECTOR3_FIELD_3D cellCenters = VECTOR3_FIELD_3D::cellCenters(toExtend);

  // build the gradient field
  VECTOR3_FIELD_3D targetGradient = VECTOR3_FIELD_3D::gradient(distanceField);

  // perform Nacelle on all the tagged extension cells
  FIELD_3D extendedOld(toExtend);
  int maxSteps = 100;
  REAL fictionalDt = 1.0;
  cout << " Computing extension field ..."; flush(cout);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < subset.size(); i++)
  {
    REAL diff = 1;
    int steps = 0;
    VECTOR3& position = cellCenters[subset[i]];
    VECTOR3 targetDelta = targetGradient(position);

    while (diff > 1e-6 && steps < maxSteps)
    {
      REAL targetDistance = distanceField(position);
      diff = fabs(targetDistance);

      // if the gradient is zero, stop trying, because it probably
      // wandered off the grid
      //if (norm2(targetDelta) < 1e-6) break;
      if (targetDelta.norm() < 1e-6) break;

      // go ahead and always do first -- second gives the occasional wacky value
      // that throws off the stability of the simulation
      targetDelta.normalize();
      VECTOR3 move = (targetDistance) * targetDelta * fictionalDt;

      position = position - move;

      // update the direction
      targetDelta = targetGradient(position);
      steps++;
    }

    // do the extension here
    extendedOld[subset[i]] = toExtend.quarticLookup(position);
  }

  return extendedOld;
}

///////////////////////////////////////////////////////////////////////
// compute the extension field for a masked set of points
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::computeExtensionFieldMasked(const FIELD_3D& distanceField, const FIELD_3D& toExtend, const FIELD_3D& mask)
{
  TIMER functionTimer(__FUNCTION__);

  // get the cell centers 
  VECTOR3_FIELD_3D cellCenters = VECTOR3_FIELD_3D::cellCenters(toExtend);

  // build the gradient field
  VECTOR3_FIELD_3D targetGradient = VECTOR3_FIELD_3D::gradient(distanceField);

  const int xRes = mask.xRes();
  const int yRes = mask.yRes();
  const int zRes = mask.zRes();

  //const REAL invDx = 1.0 / distanceField.dx();
  const REAL invDx = 1.0 / distanceField.dx();

  // perform Nacelle on all the tagged extension cells
  FIELD_3D extendedOld(toExtend);
  int maxSteps = 100;
  REAL fictionalDt = 1.0;
  cout << " Computing masked extension field ..."; flush(cout);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        VECTOR3& position = cellCenters(x,y,z);
        REAL currentDistance = fabs(distanceField(position) * invDx);

        //if (mask(x,y,z) < 0.5) continue;
        if (mask(x,y,z) < 0.5 || currentDistance <= 1.0) continue;

        REAL diff = 1;
        int steps = 0;
        VECTOR3 targetDelta = targetGradient(position);

        while (diff > 1e-6 && steps < maxSteps)
        {
          const REAL targetDistance = distanceField(position);
          diff = fabs(targetDistance);

          // if the gradient is zero, stop trying, because it probably
          // wandered off the grid
          //if (norm2(targetDelta) < 1e-6) break;
          if (targetDelta.norm() < 1e-6) break;

          // go ahead and always do first -- second gives the occasional wacky value
          // that throws off the stability of the simulation
          targetDelta.normalize();
          const VECTOR3 move = (targetDistance) * targetDelta * fictionalDt;

          position = position - move;

          // update the direction
          targetDelta = targetGradient(position);
          steps++;
        }

        // do the extension here
        //extendedOld(x,y,z) = toExtend.quarticLookup(position);
        extendedOld(x,y,z) = toExtend.nearestNeighborLookup(position);
      }

  return extendedOld;
}

///////////////////////////////////////////////////////////////////////
// compute the extension field for a masked set of points
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::computeExtensionFieldMaskedBadAdvect(const FIELD_3D& distanceField, const FIELD_3D& toExtend, const FIELD_3D& mask)
{
  TIMER functionTimer(__FUNCTION__);

  // get the cell centers 
  VECTOR3_FIELD_3D cellCenters = VECTOR3_FIELD_3D::cellCenters(toExtend);

  // build the gradient field
  VECTOR3_FIELD_3D targetGradient = VECTOR3_FIELD_3D::gradient(distanceField);

  const int xRes = mask.xRes();
  const int yRes = mask.yRes();
  const int zRes = mask.zRes();

  //const REAL invDx = 1.0 / distanceField.dx();

  // perform Nacelle on all the tagged extension cells
  FIELD_3D extendedOld(toExtend);
  int maxSteps = 100;
  REAL fictionalDt = 1.0;
  cout << " Computing bad masked extension field ..."; flush(cout);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        VECTOR3& position = cellCenters(x,y,z);
        //REAL currentDistance = fabs(distanceField(position) * invDx);

        //if (mask(x,y,z) < 0.5) continue;
        if (mask(x,y,z) < 0.5) continue;

        REAL diff = 1;
        int steps = 0;
        VECTOR3 targetDelta = targetGradient(position);

        while (diff > 1e-6 && steps < maxSteps)
        {
          const REAL targetDistance = distanceField(position);
          diff = fabs(targetDistance);

          // if the gradient is zero, stop trying, because it probably
          // wandered off the grid
          //if (norm2(targetDelta) < 1e-6) break;
          if (targetDelta.norm() < 1e-6) break;

          // go ahead and always do first -- second gives the occasional wacky value
          // that throws off the stability of the simulation
          targetDelta.normalize();
          const VECTOR3 move = (targetDistance) * targetDelta * fictionalDt;

          position = position - move;

          // update the direction
          targetDelta = targetGradient(position);
          steps++;
        }

        // do the extension here
        extendedOld(x,y,z) = toExtend.quarticLookup(position);
        //extendedOld(x,y,z) = toExtend.nearestNeighborLookup(position);
      }

  return extendedOld;
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectLazy(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField, FIELD_3D& newField)
{
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();

  // build the gradient field
  REAL invDx = 1.0 / oldField.dx();
  
  // which cells need extension values?
  map<int, bool> extensionMap;

  // which cells will appear in the final advected field?
  vector<int> finalIndices;

  // build a list of cells to populate the extension for
  cout << " Computing extension band ..."; flush(cout);
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // add it to the list of final
        int index = x + y * xRes + z * slabSize;
        finalIndices.push_back(index);

        // backtrace, store those cells for extension later.
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        const int i000 = x0 + y0 * xRes + z0 * slabSize;
        const int i010 = x0 + y1 * xRes + z0 * slabSize;
        const int i100 = x1 + y0 * xRes + z0 * slabSize;
        const int i110 = x1 + y1 * xRes + z0 * slabSize;
        const int i001 = x0 + y0 * xRes + z1 * slabSize;
        const int i011 = x0 + y1 * xRes + z1 * slabSize;
        const int i101 = x1 + y0 * xRes + z1 * slabSize;
        const int i111 = x1 + y1 * xRes + z1 * slabSize;

        extensionMap[i000] = true;
        extensionMap[i001] = true;
        extensionMap[i010] = true;
        extensionMap[i011] = true;
        extensionMap[i100] = true;
        extensionMap[i101] = true;
        extensionMap[i110] = true;
        extensionMap[i111] = true;
      }

  // convert extension map to a vector
  vector<int> extensionVector;
  map<int,bool>::const_iterator iter;
  for (iter = extensionMap.begin(); iter != extensionMap.end(); iter++)
    extensionVector.push_back(iter->first);

  // perform Nacelle on all the tagged extension cells
  FIELD_3D extendedOld = computeExtensionFieldSubset(distanceOld, oldField, extensionVector); 

  cout << " Computing final advection ..."; flush(cout);
  newField = advectSubset(dt, velocityGrid, extendedOld, finalIndices);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectMacCormack(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2)
{
  TIMER functionTimer(__FUNCTION__);
	FIELD_3D& phiHatN  = temp1;
	FIELD_3D& phiHatN1 = temp2;

	const int sx = oldField.xRes();
	const int sy = oldField.yRes();
	const int sz = oldField.zRes();

	for (int x = 0; x < sx * sy * sz; x++)
		phiHatN[x] = phiHatN1[x] = oldField[x];

	FIELD_3D& phiN    = oldField;
	FIELD_3D& phiN1   = newField;

	// phiHatN1 = A(phiN)
  cout << " forward ... ";flush(cout);
	advect(dt, velocityGrid, phiN, phiHatN1);

	// phiHatN = A^R(phiHatN1)
  cout << " backward ... ";flush(cout);
	advect(-1.0 * dt, velocityGrid, phiHatN1, phiHatN);

  phiN1.clear();
  phiN1 += phiN;
  phiN1 -= phiHatN;
  phiN1 *= 0.5;
  phiN1 += phiHatN1;

	// clamp any newly created extrema
  cout << " clamping extrema ... ";flush(cout);
	clampExtrema(dt, velocityGrid, oldField, newField);

	// if the error estimate was bad, revert to first order
  cout << " clamping rays ... ";flush(cout);
	clampOutsideRays(dt, velocityGrid, oldField, phiHatN1, newField);
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// set the closes point values needed for MacCormack
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setMacCormackClosestPoints(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField)
{
  TIMER functionTimer(__FUNCTION__);
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();
  REAL invDx = 1.0 / oldField.dx();

  FIELD_3D closestPointValues(oldField);
  VECTOR3_FIELD_3D targetGradient = VECTOR3_FIELD_3D::gradient(distanceOld);

  // build a list of cells to populate the extension for
  TIMER forwardTimer("Forward band");
  cout << " Computing forward extension band ..."; flush(cout);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // add it to the list of final
        //int index = x + y * xRes + z * slabSize;

        // backtrace, store those cells for extension later.
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

        closestPointValues[i000] = getClosestPointValue(targetGradient, oldField.cellCenter(x0, y0, z0), distanceOld, oldField); 
        closestPointValues[i001] = getClosestPointValue(targetGradient, oldField.cellCenter(x0, y0, z1), distanceOld, oldField); 
        closestPointValues[i010] = getClosestPointValue(targetGradient, oldField.cellCenter(x0, y1, z0), distanceOld, oldField); 
        closestPointValues[i011] = getClosestPointValue(targetGradient, oldField.cellCenter(x0, y1, z1), distanceOld, oldField); 
        closestPointValues[i100] = getClosestPointValue(targetGradient, oldField.cellCenter(x1, y0, z0), distanceOld, oldField); 
        closestPointValues[i101] = getClosestPointValue(targetGradient, oldField.cellCenter(x1, y0, z1), distanceOld, oldField); 
        closestPointValues[i110] = getClosestPointValue(targetGradient, oldField.cellCenter(x1, y1, z0), distanceOld, oldField); 
        closestPointValues[i111] = getClosestPointValue(targetGradient, oldField.cellCenter(x1, y1, z1), distanceOld, oldField); 
      }

  // see what cells show up in the backward advection
  cout << " Computing backward extension band ..."; flush(cout);
  TIMER backwardTimer("Backward band");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // add it to the list of final
        //int index = x + y * xRes + z * slabSize;

        // backtrace, store those cells for extension later.
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x + dt * velocity[0];
        REAL yTrace = y + dt * velocity[1];
        REAL zTrace = z + dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

        closestPointValues[i000] = getClosestPointValue(targetGradient, oldField.cellCenter(x0, y0, z0), distanceOld, oldField); 
        closestPointValues[i001] = getClosestPointValue(targetGradient, oldField.cellCenter(x0, y0, z1), distanceOld, oldField); 
        closestPointValues[i010] = getClosestPointValue(targetGradient, oldField.cellCenter(x0, y1, z0), distanceOld, oldField); 
        closestPointValues[i011] = getClosestPointValue(targetGradient, oldField.cellCenter(x0, y1, z1), distanceOld, oldField); 
        closestPointValues[i100] = getClosestPointValue(targetGradient, oldField.cellCenter(x1, y0, z0), distanceOld, oldField); 
        closestPointValues[i101] = getClosestPointValue(targetGradient, oldField.cellCenter(x1, y0, z1), distanceOld, oldField); 
        closestPointValues[i110] = getClosestPointValue(targetGradient, oldField.cellCenter(x1, y1, z0), distanceOld, oldField); 
        closestPointValues[i111] = getClosestPointValue(targetGradient, oldField.cellCenter(x1, y1, z1), distanceOld, oldField); 
      }

  /*
  // what cells need to be extended so that correct values appear 
  // upon backwards advection?
  TIMER backwardExtensionTimer("Backward extension band");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < backwardExtensionVector.size(); i++)
  {
    int index = backwardExtensionVector[i];
    const int z = index / slabSize;
    const int y = (index - z * slabSize) / xRes;
    const int x = index - z * slabSize - y * xRes;

    // backtrace, store those cells for extension later.
    const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
    REAL xTrace = x - dt * velocity[0];
    REAL yTrace = y - dt * velocity[1];
    REAL zTrace = z - dt * velocity[2];

    // clamp backtrace to grid boundaries
    xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
    xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
    yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
    yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
    zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
    zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

    // locate neighbors to interpolate
    const int x0 = (int)xTrace;
    const int x1 = x0 + 1;
    const int y0 = (int)yTrace;
    const int y1 = y0 + 1;
    const int z0 = (int)zTrace;
    const int z1 = z0 + 1;

    const int z0slabSize = z0 * slabSize;
    const int z1slabSize = z1 * slabSize;
    const int y0xRes = y0 * xRes;
    const int y1xRes = y1 * xRes;

    const int i000 = x0 + y0xRes + z0slabSize;
    const int i010 = x0 + y1xRes + z0slabSize;
    const int i100 = x1 + y0xRes + z0slabSize;
    const int i110 = x1 + y1xRes + z0slabSize;
    const int i001 = x0 + y0xRes + z1slabSize;
    const int i011 = x0 + y1xRes + z1slabSize;
    const int i101 = x1 + y0xRes + z1slabSize;
    const int i111 = x1 + y1xRes + z1slabSize;

    int threadID = omp_get_thread_num();
    forwardExtensionMapPerThread[threadID].insert(i000);
    forwardExtensionMapPerThread[threadID].insert(i001);
    forwardExtensionMapPerThread[threadID].insert(i010);
    forwardExtensionMapPerThread[threadID].insert(i011);
    forwardExtensionMapPerThread[threadID].insert(i100);
    forwardExtensionMapPerThread[threadID].insert(i101);
    forwardExtensionMapPerThread[threadID].insert(i110);
    forwardExtensionMapPerThread[threadID].insert(i111);

    finalMapPerThread[threadID].insert(i000);
    finalMapPerThread[threadID].insert(i001);
    finalMapPerThread[threadID].insert(i010);
    finalMapPerThread[threadID].insert(i011);
    finalMapPerThread[threadID].insert(i100);
    finalMapPerThread[threadID].insert(i101);
    finalMapPerThread[threadID].insert(i110);
    finalMapPerThread[threadID].insert(i111);
  }
  */

  return closestPointValues;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectMacCormackNarrowBand(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2)
{
  TIMER functionTimer(__FUNCTION__);
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();
  REAL invDx = 1.0 / oldField.dx();

  // get the number of OMP threads
  //static int totalThreads = omp_get_max_threads();

  // track what cell to do the closest point transform on
  FIELD_3D& forwardExtend = temp1;
  forwardExtend.clear();
  REAL* forwardExtendData = forwardExtend.data();

  FIELD_3D& backwardExtend = temp2;
  backwardExtend.clear();
  REAL* backwardExtendData = backwardExtend.data();

  // build a list of cells to populate the extension for
  TIMER forwardTimer("Forward band");
  cout << " Computing forward extension band ..."; flush(cout);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // add it to the list of final
        //int index = x + y * xRes + z * slabSize;

        // backtrace, store those cells for extension later.
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

        //int threadID = omp_get_thread_num();

#pragma omp atomic
        forwardExtendData[i000] += 1;
#pragma omp atomic
        forwardExtendData[i001] += 1;
#pragma omp atomic
        forwardExtendData[i010] += 1;
#pragma omp atomic
        forwardExtendData[i011] += 1;
#pragma omp atomic
        forwardExtendData[i100] += 1;
#pragma omp atomic
        forwardExtendData[i101] += 1;
#pragma omp atomic
        forwardExtendData[i110] += 1;
#pragma omp atomic
        forwardExtendData[i111] += 1;
      }

  // see what cells show up in the backward advection
  TIMER backwardTimer("Backward band");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // add it to the list of final
        //int index = x + y * xRes + z * slabSize;

        // backtrace, store those cells for extension later.
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x + dt * velocity[0];
        REAL yTrace = y + dt * velocity[1];
        REAL zTrace = z + dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

        //int threadID = omp_get_thread_num();
        
#pragma omp atomic
        backwardExtendData[i000] += 1;
#pragma omp atomic
        backwardExtendData[i001] += 1;
#pragma omp atomic
        backwardExtendData[i010] += 1;
#pragma omp atomic
        backwardExtendData[i011] += 1;
#pragma omp atomic
        backwardExtendData[i100] += 1;
#pragma omp atomic
        backwardExtendData[i101] += 1;
#pragma omp atomic
        backwardExtendData[i110] += 1;
#pragma omp atomic
        backwardExtendData[i111] += 1;
      }

  // what cells need to be extended so that correct values appear 
  // upon backwards advection?
  TIMER backwardExtensionTimer("Backward extension band");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // add it to the list of final
        int index = x + y * xRes + z * slabSize;

        if (backwardExtendData[index] < 0.5) 
          continue;

        // backtrace, store those cells for extension later.
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

#pragma omp atomic
        forwardExtendData[i000] += 1;
#pragma omp atomic
        forwardExtendData[i001] += 1;
#pragma omp atomic
        forwardExtendData[i010] += 1;
#pragma omp atomic
        forwardExtendData[i011] += 1;
#pragma omp atomic
        forwardExtendData[i100] += 1;
#pragma omp atomic
        forwardExtendData[i101] += 1;
#pragma omp atomic
        forwardExtendData[i110] += 1;
#pragma omp atomic
        forwardExtendData[i111] += 1;
      }

  // add the backward extends to the forward
  forwardExtend += backwardExtend;

  // perform Nacelle on all the tagged extension cells
  FIELD_3D extendedOld = computeExtensionFieldMasked(distanceOld, oldField, forwardExtend); 

	FIELD_3D& phiN    = extendedOld;
	FIELD_3D& phiN1   = newField;
	FIELD_3D& phiHatN  = temp1;
	FIELD_3D& phiHatN1 = temp2;

  cout << " advecting forward ... ";flush(cout);
  advect(dt, velocityGrid, phiN, phiHatN1);

  cout << " advecting backward ... ";flush(cout);
	advect(-1.0 * dt, velocityGrid, phiHatN1, phiHatN);

  TIMER arithmeticTimer("Arithmetic field ops");
  phiN1.clear();
  phiN1 += phiN;
  phiN1 -= phiHatN;
  phiN1 *= 0.5;
  phiN1 += phiHatN1;

  cout << " clamping ... ";flush(cout);

	// clamp any newly created extrema
	clampExtrema(dt, velocityGrid, extendedOld, newField);

	// if the error estimate was bad, revert to first order
	clampOutsideRays(dt, velocityGrid, extendedOld, phiHatN1, newField);

  cout << " done." << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectMacCormackLazy(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2)
{
  TIMER functionTimer(__FUNCTION__);
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();
  REAL invDx = 1.0 / oldField.dx();

  // get the number of OMP threads
  static int totalThreads = omp_get_max_threads();

  // This block is very slow on Snow Leopard for some reason. Is OpenMP doing something under the hood
  // that sets things up to be parallel in a clever way?
  TIMER setTimer("Set allocation");
  //vector<set<int> > forwardExtensionMapPerThread(totalThreads);
  //vector<set<int> > finalMapPerThread(totalThreads);
  set<int>* forwardExtensionMapPerThread = new set<int>[totalThreads];
  set<int>* finalMapPerThread = new set<int>[totalThreads];

  // build a list of cells to populate the extension for
  TIMER forwardTimer("Forward band");
  cout << " Computing forward extension band ..."; flush(cout);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // add it to the list of final
        int index = x + y * xRes + z * slabSize;

        // backtrace, store those cells for extension later.
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

        int threadID = omp_get_thread_num();

        finalMapPerThread[threadID].insert(index);
        forwardExtensionMapPerThread[threadID].insert(i000);
        forwardExtensionMapPerThread[threadID].insert(i001);
        forwardExtensionMapPerThread[threadID].insert(i010);
        forwardExtensionMapPerThread[threadID].insert(i011);
        forwardExtensionMapPerThread[threadID].insert(i100);
        forwardExtensionMapPerThread[threadID].insert(i101);
        forwardExtensionMapPerThread[threadID].insert(i110);
        forwardExtensionMapPerThread[threadID].insert(i111);
      }

  // which cells are touched on the backward extension?
  //vector<set<int> > backwardExtensionMapPerThread(totalThreads);
  set<int>* backwardExtensionMapPerThread = new set<int>[totalThreads];

  // see what cells show up in the backward advection
  TIMER backwardTimer("Backward band");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        REAL currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // add it to the list of final
        int index = x + y * xRes + z * slabSize;

        // backtrace, store those cells for extension later.
        const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
        REAL xTrace = x + dt * velocity[0];
        REAL yTrace = y + dt * velocity[1];
        REAL zTrace = z + dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

        int threadID = omp_get_thread_num();
        
        finalMapPerThread[threadID].insert(index);
        backwardExtensionMapPerThread[threadID].insert(i000);
        backwardExtensionMapPerThread[threadID].insert(i001);
        backwardExtensionMapPerThread[threadID].insert(i010);
        backwardExtensionMapPerThread[threadID].insert(i011);
        backwardExtensionMapPerThread[threadID].insert(i100);
        backwardExtensionMapPerThread[threadID].insert(i101);
        backwardExtensionMapPerThread[threadID].insert(i110);
        backwardExtensionMapPerThread[threadID].insert(i111);
      }

  // convert to a vector so OpenMP will support it
  vector<int> backwardExtensionVector;
  //map<int, bool>::const_iterator iter;
  set<int>::const_iterator iter;
  for (int x = 0; x < totalThreads; x++)
    for (iter = backwardExtensionMapPerThread[x].begin(); iter != backwardExtensionMapPerThread[x].end(); iter++)
      backwardExtensionVector.push_back(*iter);

  // what cells need to be extended so that correct values appear 
  // upon backwards advection?
  TIMER backwardExtensionTimer("Backward extension band");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < backwardExtensionVector.size(); i++)
  {
    int index = backwardExtensionVector[i];
    const int z = index / slabSize;
    const int y = (index - z * slabSize) / xRes;
    const int x = index - z * slabSize - y * xRes;

    // backtrace, store those cells for extension later.
    const VECTOR3 velocity = velocityGrid(oldField.cellCenter(x,y,z));
    REAL xTrace = x - dt * velocity[0];
    REAL yTrace = y - dt * velocity[1];
    REAL zTrace = z - dt * velocity[2];

    // clamp backtrace to grid boundaries
    xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
    xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
    yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
    yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
    zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
    zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

    // locate neighbors to interpolate
    const int x0 = (int)xTrace;
    const int x1 = x0 + 1;
    const int y0 = (int)yTrace;
    const int y1 = y0 + 1;
    const int z0 = (int)zTrace;
    const int z1 = z0 + 1;

    const int z0slabSize = z0 * slabSize;
    const int z1slabSize = z1 * slabSize;
    const int y0xRes = y0 * xRes;
    const int y1xRes = y1 * xRes;

    const int i000 = x0 + y0xRes + z0slabSize;
    const int i010 = x0 + y1xRes + z0slabSize;
    const int i100 = x1 + y0xRes + z0slabSize;
    const int i110 = x1 + y1xRes + z0slabSize;
    const int i001 = x0 + y0xRes + z1slabSize;
    const int i011 = x0 + y1xRes + z1slabSize;
    const int i101 = x1 + y0xRes + z1slabSize;
    const int i111 = x1 + y1xRes + z1slabSize;

    int threadID = omp_get_thread_num();
    forwardExtensionMapPerThread[threadID].insert(i000);
    forwardExtensionMapPerThread[threadID].insert(i001);
    forwardExtensionMapPerThread[threadID].insert(i010);
    forwardExtensionMapPerThread[threadID].insert(i011);
    forwardExtensionMapPerThread[threadID].insert(i100);
    forwardExtensionMapPerThread[threadID].insert(i101);
    forwardExtensionMapPerThread[threadID].insert(i110);
    forwardExtensionMapPerThread[threadID].insert(i111);

    finalMapPerThread[threadID].insert(i000);
    finalMapPerThread[threadID].insert(i001);
    finalMapPerThread[threadID].insert(i010);
    finalMapPerThread[threadID].insert(i011);
    finalMapPerThread[threadID].insert(i100);
    finalMapPerThread[threadID].insert(i101);
    finalMapPerThread[threadID].insert(i110);
    finalMapPerThread[threadID].insert(i111);
  }

  // start registering the rest under the current function name
  TIMER mapMergeTimer("Map merge");

  // join together the extension maps
  //map<int, bool> extensionMap;
  set<int> extensionMap;
  for (int x = 0; x < totalThreads; x++)
    for (iter = backwardExtensionMapPerThread[x].begin(); iter != backwardExtensionMapPerThread[x].end(); iter++)
      extensionMap.insert(*iter);

  for (int x = 0; x < totalThreads; x++)
    for (iter = forwardExtensionMapPerThread[x].begin(); iter != forwardExtensionMapPerThread[x].end(); iter++)
      extensionMap.insert(*iter);

  // convert extension maps to a vectors
  vector<int> extensionVector;
  for (iter = extensionMap.begin(); iter != extensionMap.end(); iter++)
    extensionVector.push_back(*iter);

  vector<int> finalIndices;
  for (int x = 0; x < totalThreads; x++)
    for (iter = finalMapPerThread[x].begin(); iter != finalMapPerThread[x].end(); iter++)
      finalIndices.push_back(*iter);

  TIMER functionTimer2(__FUNCTION__);

  // perform Nacelle on all the tagged extension cells
  FIELD_3D extendedOld = computeExtensionFieldSubset(distanceOld, oldField, extensionVector); 

	FIELD_3D& phiN    = extendedOld;
	FIELD_3D& phiN1   = newField;
	FIELD_3D& phiHatN  = temp1;
	FIELD_3D& phiHatN1 = temp2;

  // TODO: don't do a copy here?

  cout << " advecting forward ... ";flush(cout);
  phiHatN1 = advectSubset(dt, velocityGrid, phiN, finalIndices);

  cout << " advecting backward ... ";flush(cout);
	phiHatN  = advectSubset(-1.0 * dt, velocityGrid, phiHatN1, finalIndices);

  TIMER arithmeticTimer("Arithmetic field ops");
  phiN1.clear();
  phiN1 += phiN;
  phiN1 -= phiHatN;
  phiN1 *= 0.5;
  phiN1 += phiHatN1;

  cout << " clamping ... ";flush(cout);

	// clamp any newly created extrema
	clampExtremaSubset(dt, velocityGrid, extendedOld, newField, finalIndices);

	// if the error estimate was bad, revert to first order
	clampOutsideRaysSubset(dt, velocityGrid, extendedOld, phiHatN1, newField, finalIndices);

  cout << " done." << endl;

  // cleanup
  TIMER cleanupTimer("Cleanup timer");

  for (int x = 0; x < totalThreads; x++)
  {
    forwardExtensionMapPerThread[x].clear();
    finalMapPerThread[x].clear();
    backwardExtensionMapPerThread[x].clear();
  }
  delete[] forwardExtensionMapPerThread;
  delete[] finalMapPerThread;
  delete[] backwardExtensionMapPerThread;
}

///////////////////////////////////////////////////////////////////////
// Clamp the extrema generated by the BFECC error correction
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampExtrema(const REAL dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, FIELD_3D& newField)
{
  TIMER functionTimer(__FUNCTION__);
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  //const int slabSize = oldField.slabSize();

#pragma omp parallel
#pragma omp for schedule(dynamic)
	for (int z = 1; z < zRes-1; z++)
		for (int y = 1; y < yRes-1; y++)
			for (int x = 1; x < xRes-1; x++)
			{
				//const int index = x + y * xRes + z * slabSize;
				// backtrace
        const VECTOR3 velocity = velocityField(oldField.cellCenter(x,y,z));
        //if (norm2(velocity) < 1e-6) continue;
        if (velocity.norm() < 1e-6) continue;

        REAL xTrace = x - dt * velocity[0];
        REAL yTrace = y - dt * velocity[1];
        REAL zTrace = z - dt * velocity[2];

				// clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

        REAL minField = oldField(x0,y0,z0);
        REAL maxField = minField;

        const REAL x0y1z0 = oldField(x0,y1,z0);
        const REAL x1y0z0 = oldField(x1,y0,z0);
        const REAL x1y1z0 = oldField(x1,y1,z0);
        const REAL x0y0z1 = oldField(x0,y0,z1);
        const REAL x0y1z1 = oldField(x0,y1,z1);
        const REAL x1y0z1 = oldField(x1,y0,z1);
        const REAL x1y1z1 = oldField(x1,y1,z1);

        minField = (x0y1z0 < minField) ? x0y1z0 : minField;
        maxField = (x0y1z0 > maxField) ? x0y1z0 : maxField;

        minField = (x1y0z0 < minField) ? x1y0z0 : minField;
        maxField = (x1y0z0 > maxField) ? x1y0z0 : maxField;

        minField = (x1y1z0 < minField) ? x1y1z0 : minField;
        maxField = (x1y1z0 > maxField) ? x1y1z0 : maxField;

        minField = (x0y0z1 < minField) ? x0y0z1 : minField;
        maxField = (x0y0z1 > maxField) ? x0y0z1 : maxField;

        minField = (x0y1z1 < minField) ? x0y1z1 : minField;
        maxField = (x0y1z1 > maxField) ? x0y1z1 : maxField;

        minField = (x1y0z1 < minField) ? x1y0z1 : minField;
        maxField = (x1y0z1 > maxField) ? x1y0z1 : maxField;

        minField = (x1y1z1 < minField) ? x1y1z1 : minField;
        maxField = (x1y1z1 > maxField) ? x1y1z1 : maxField;

        const REAL newValue = newField(x,y,z);

        newField(x,y,z) = (newValue > maxField) ? maxField : newValue;
        newField(x,y,z) = (newValue < minField) ? minField : newValue;
      }
}

///////////////////////////////////////////////////////////////////////
// Clamp the extrema generated by the BFECC error correction
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampExtremaSubset(const REAL dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, FIELD_3D& newField, const vector<int>& subset)
{
  TIMER functionTimer(__FUNCTION__);
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();

#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (unsigned int i = 0; i < subset.size(); i++)
  {
    const int index = subset[i];
    const int z = index / slabSize;
    const int y = (index - z * slabSize) / xRes;
    const int x = index - z * slabSize - y * xRes;

    // backtrace
    const VECTOR3 velocity = velocityField(oldField.cellCenter(x,y,z));
    //if (norm2(velocity) < 1e-6) continue;
    if (velocity.norm() < 1e-6) continue;

    REAL xTrace = x - dt * velocity[0];
    REAL yTrace = y - dt * velocity[1];
    REAL zTrace = z - dt * velocity[2];

    // clamp backtrace to grid boundaries
    xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
    xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
    yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
    yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
    zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
    zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

    // locate neighbors to interpolate
    const int x0 = (int)xTrace;
    const int x1 = x0 + 1;
    const int y0 = (int)yTrace;
    const int y1 = y0 + 1;
    const int z0 = (int)zTrace;
    const int z1 = z0 + 1;

    REAL minField = oldField(x0,y0,z0);
    REAL maxField = minField;

    const REAL x0y1z0 = oldField(x0,y1,z0);
    const REAL x1y0z0 = oldField(x1,y0,z0);
    const REAL x1y1z0 = oldField(x1,y1,z0);
    const REAL x0y0z1 = oldField(x0,y0,z1);
    const REAL x0y1z1 = oldField(x0,y1,z1);
    const REAL x1y0z1 = oldField(x1,y0,z1);
    const REAL x1y1z1 = oldField(x1,y1,z1);

    minField = (x0y1z0 < minField) ? x0y1z0 : minField;
    maxField = (x0y1z0 > maxField) ? x0y1z0 : maxField;

    minField = (x1y0z0 < minField) ? x1y0z0 : minField;
    maxField = (x1y0z0 > maxField) ? x1y0z0 : maxField;

    minField = (x1y1z0 < minField) ? x1y1z0 : minField;
    maxField = (x1y1z0 > maxField) ? x1y1z0 : maxField;

    minField = (x0y0z1 < minField) ? x0y0z1 : minField;
    maxField = (x0y0z1 > maxField) ? x0y0z1 : maxField;

    minField = (x0y1z1 < minField) ? x0y1z1 : minField;
    maxField = (x0y1z1 > maxField) ? x0y1z1 : maxField;

    minField = (x1y0z1 < minField) ? x1y0z1 : minField;
    maxField = (x1y0z1 > maxField) ? x1y0z1 : maxField;

    minField = (x1y1z1 < minField) ? x1y1z1 : minField;
    maxField = (x1y1z1 > maxField) ? x1y1z1 : maxField;

    const REAL newValue = newField(x,y,z);

    newField(x,y,z) = (newValue > maxField) ? maxField : newValue;
    newField(x,y,z) = (newValue < minField) ? minField : newValue;
  }
}

//////////////////////////////////////////////////////////////////////
// Reverts any backtraces that go into boundaries back to first 
// order -- in this case the error correction term was totally
// incorrect
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampOutsideRays(const REAL dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, const FIELD_3D& oldAdvection, FIELD_3D& newField)
{
  TIMER functionTimer(__FUNCTION__);
  // this actually is correct -- should be using the res of the
  // low res grid
  const int sx = velocityField.xRes();
  const int sy = velocityField.yRes();
  const int sz = velocityField.zRes();
  //const int slabSize = oldField.slabSize();
	//for (int z = 1; z < sz - 1; z++)
	//	for (int y = 1; y < sy - 1; y++)
	//		for (int x = 1; x < sx - 1; x++)
#pragma omp parallel
#pragma omp for schedule(dynamic)
	for (int z = 0; z < oldField.zRes(); z++)
		for (int y = 0; y < oldField.yRes(); y++)
			for (int x = 0; x < oldField.xRes(); x++)
			{
				//const int index = x + y * sx + z * slabSize;
				// backtrace
        VECTOR3 velocity = velocityField(oldField.cellCenter(x,y,z));

        // this line seems to *significantly* reduce banding
        //if (norm2(velocity) < 1e-6) continue;
        if (velocity.norm() < 1e-6) continue;

        velocity *= dt;

				float xBackward = x + velocity[0];
				float yBackward = y + velocity[1];
				float zBackward = z + velocity[2];

				float xTrace    = x - velocity[0];
				float yTrace    = x - velocity[1];
				float zTrace    = x - velocity[2];

				// see if it goes outside the boundaries
				bool hasObstacle = 
					(zTrace < 1.0f)    || (zTrace > sz - 2.0f) ||
					(yTrace < 1.0f)    || (yTrace > sy - 2.0f) ||
					(xTrace < 1.0f)    || (xTrace > sx - 2.0f) ||
					(zBackward < 1.0f) || (zBackward > sz - 2.0f) ||
					(yBackward < 1.0f) || (yBackward > sy - 2.0f) ||
					(xBackward < 1.0f) || (xBackward > sx - 2.0f);

				// reuse old advection instead of doing another one...
				//if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }
				if(hasObstacle) { newField(x,y,z) = oldAdvection(x,y,z); continue; }

        /*
				// clamp to prevent an out of bounds access when looking into
				// the _obstacles array
				zTrace = (zTrace < 0.5f) ? 0.5f : zTrace;
				zTrace = (zTrace > sz - 1.5f) ? sz - 1.5f : zTrace;
				yTrace = (yTrace < 0.5f) ? 0.5f : yTrace;
				yTrace = (yTrace > sy - 1.5f) ? sy - 1.5f : yTrace;
				xTrace = (xTrace < 0.5f) ? 0.5f : xTrace;
				xTrace = (xTrace > sx - 1.5f) ? sx - 1.5f : xTrace;

				// locate neighbors to interpolate,
				// do backward first since we will use the forward indices if a
				// reversion is actually necessary
				zBackward = (zBackward < 0.5f) ? 0.5f : zBackward;
				zBackward = (zBackward > sz - 1.5f) ? sz - 1.5f : zBackward;
				yBackward = (yBackward < 0.5f) ? 0.5f : yBackward;
				yBackward = (yBackward > sy - 1.5f) ? sy - 1.5f : yBackward;
				xBackward = (xBackward < 0.5f) ? 0.5f : xBackward;
				xBackward = (xBackward > sx - 1.5f) ? sx - 1.5f : xBackward;

				int x0 = (int)xBackward;
				int x1 = x0 + 1;
				int y0 = (int)yBackward;
				int y1 = y0 + 1;
				int z0 = (int)zBackward;
				int z1 = z0 + 1;
				if(obstacles && !hasObstacle) {
					hasObstacle = hasObstacle || 
						obstacles[x0 + y0 * sx + z0*slabSize] ||
						obstacles[x0 + y1 * sx + z0*slabSize] ||
						obstacles[x1 + y0 * sx + z0*slabSize] ||
						obstacles[x1 + y1 * sx + z0*slabSize] ||
						obstacles[x0 + y0 * sx + z1*slabSize] ||
						obstacles[x0 + y1 * sx + z1*slabSize] ||
						obstacles[x1 + y0 * sx + z1*slabSize] ||
						obstacles[x1 + y1 * sx + z1*slabSize] ;
				}
				// reuse old advection instead of doing another one...
				if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }

				x0 = (int)xTrace;
				x1 = x0 + 1;
				y0 = (int)yTrace;
				y1 = y0 + 1;
				z0 = (int)zTrace;
				z1 = z0 + 1;
				if(obstacles && !hasObstacle) {
					hasObstacle = hasObstacle || 
						obstacles[x0 + y0 * sx + z0*slabSize] ||
						obstacles[x0 + y1 * sx + z0*slabSize] ||
						obstacles[x1 + y0 * sx + z0*slabSize] ||
						obstacles[x1 + y1 * sx + z0*slabSize] ||
						obstacles[x0 + y0 * sx + z1*slabSize] ||
						obstacles[x0 + y1 * sx + z1*slabSize] ||
						obstacles[x1 + y0 * sx + z1*slabSize] ||
						obstacles[x1 + y1 * sx + z1*slabSize] ;
				} // obstacle array
				// reuse old advection instead of doing another one...
				if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }

				// see if either the forward or backward ray went into
				// a boundary
				if (hasObstacle) {
					// get interpolation weights
					float s1 = xTrace - x0;
					float s0 = 1.0f - s1;
					float t1 = yTrace - y0;
					float t0 = 1.0f - t1;
					float u1 = zTrace - z0;
					float u0 = 1.0f - u1;

					const int i000 = x0 + y0 * sx + z0 * slabSize;
					const int i010 = x0 + y1 * sx + z0 * slabSize;
					const int i100 = x1 + y0 * sx + z0 * slabSize;
					const int i110 = x1 + y1 * sx + z0 * slabSize;
					const int i001 = x0 + y0 * sx + z1 * slabSize;
					const int i011 = x0 + y1 * sx + z1 * slabSize;
					const int i101 = x1 + y0 * sx + z1 * slabSize;
					const int i111 = x1 + y1 * sx + z1 * slabSize;

					// interpolate, (indices could be computed once)
					newField[index] = u0 * (s0 * (
								t0 * oldField[i000] +
								t1 * oldField[i010]) +
							s1 * (t0 * oldField[i100] +
								t1 * oldField[i110])) +
						u1 * (s0 * (t0 * oldField[i001] +
									t1 * oldField[i011]) +
								s1 * (t0 * oldField[i101] +
									t1 * oldField[i111])); 
				}
        */
			} // xyz
}

//////////////////////////////////////////////////////////////////////
// Reverts any backtraces that go into boundaries back to first 
// order -- in this case the error correction term was totally
// incorrect
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampOutsideRaysSubset(const REAL dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, const FIELD_3D& oldAdvection, FIELD_3D& newField, const vector<int>& subset)
{
  TIMER functionTimer(__FUNCTION__);
  // this actually is correct -- should be using the res of the
  // low res grid
  const int sx = velocityField.xRes();
  const int sy = velocityField.yRes();
  const int sz = velocityField.zRes();
  const int xRes = oldField.xRes();
  const int slabSize = oldField.slabSize();

  for (unsigned int i = 0; i < subset.size(); i++)
  {
    const int index = subset[i];
    const int z = index / slabSize;
    const int y = (index - z * slabSize) / xRes;
    const int x = index - z * slabSize - y * xRes;
    
    // backtrace
    VECTOR3 velocity = velocityField(oldField.cellCenter(x,y,z));

    // this line seems to *significantly* reduce banding
    //if (norm2(velocity) < 1e-6) continue;
    if (velocity.norm() < 1e-6) continue;

    velocity *= dt;

    float xBackward = x + velocity[0];
    float yBackward = y + velocity[1];
    float zBackward = z + velocity[2];

    float xTrace    = x - velocity[0];
    float yTrace    = x - velocity[1];
    float zTrace    = x - velocity[2];

    // see if it goes outside the boundaries
    bool hasObstacle = 
      (zTrace < 1.0f)    || (zTrace > sz - 2.0f) ||
      (yTrace < 1.0f)    || (yTrace > sy - 2.0f) ||
      (xTrace < 1.0f)    || (xTrace > sx - 2.0f) ||
      (zBackward < 1.0f) || (zBackward > sz - 2.0f) ||
      (yBackward < 1.0f) || (yBackward > sy - 2.0f) ||
      (xBackward < 1.0f) || (xBackward > sx - 2.0f);

    // reuse old advection instead of doing another one...
    //if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }
    if(hasObstacle) { newField(x,y,z) = oldAdvection(x,y,z); continue; }

    /*
    // clamp to prevent an out of bounds access when looking into
    // the _obstacles array
    zTrace = (zTrace < 0.5f) ? 0.5f : zTrace;
    zTrace = (zTrace > sz - 1.5f) ? sz - 1.5f : zTrace;
    yTrace = (yTrace < 0.5f) ? 0.5f : yTrace;
    yTrace = (yTrace > sy - 1.5f) ? sy - 1.5f : yTrace;
    xTrace = (xTrace < 0.5f) ? 0.5f : xTrace;
    xTrace = (xTrace > sx - 1.5f) ? sx - 1.5f : xTrace;

    // locate neighbors to interpolate,
    // do backward first since we will use the forward indices if a
    // reversion is actually necessary
    zBackward = (zBackward < 0.5f) ? 0.5f : zBackward;
    zBackward = (zBackward > sz - 1.5f) ? sz - 1.5f : zBackward;
    yBackward = (yBackward < 0.5f) ? 0.5f : yBackward;
    yBackward = (yBackward > sy - 1.5f) ? sy - 1.5f : yBackward;
    xBackward = (xBackward < 0.5f) ? 0.5f : xBackward;
    xBackward = (xBackward > sx - 1.5f) ? sx - 1.5f : xBackward;

    int x0 = (int)xBackward;
    int x1 = x0 + 1;
    int y0 = (int)yBackward;
    int y1 = y0 + 1;
    int z0 = (int)zBackward;
    int z1 = z0 + 1;
    if(obstacles && !hasObstacle) {
      hasObstacle = hasObstacle || 
        obstacles[x0 + y0 * sx + z0*slabSize] ||
        obstacles[x0 + y1 * sx + z0*slabSize] ||
        obstacles[x1 + y0 * sx + z0*slabSize] ||
        obstacles[x1 + y1 * sx + z0*slabSize] ||
        obstacles[x0 + y0 * sx + z1*slabSize] ||
        obstacles[x0 + y1 * sx + z1*slabSize] ||
        obstacles[x1 + y0 * sx + z1*slabSize] ||
        obstacles[x1 + y1 * sx + z1*slabSize] ;
    }
    // reuse old advection instead of doing another one...
    if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }

    x0 = (int)xTrace;
    x1 = x0 + 1;
    y0 = (int)yTrace;
    y1 = y0 + 1;
    z0 = (int)zTrace;
    z1 = z0 + 1;
    if(obstacles && !hasObstacle) {
      hasObstacle = hasObstacle || 
        obstacles[x0 + y0 * sx + z0*slabSize] ||
        obstacles[x0 + y1 * sx + z0*slabSize] ||
        obstacles[x1 + y0 * sx + z0*slabSize] ||
        obstacles[x1 + y1 * sx + z0*slabSize] ||
        obstacles[x0 + y0 * sx + z1*slabSize] ||
        obstacles[x0 + y1 * sx + z1*slabSize] ||
        obstacles[x1 + y0 * sx + z1*slabSize] ||
        obstacles[x1 + y1 * sx + z1*slabSize] ;
    } // obstacle array
    // reuse old advection instead of doing another one...
    if(hasObstacle) { newField[index] = oldAdvection[index]; continue; }

    // see if either the forward or backward ray went into
    // a boundary
    if (hasObstacle) {
      // get interpolation weights
      float s1 = xTrace - x0;
      float s0 = 1.0f - s1;
      float t1 = yTrace - y0;
      float t0 = 1.0f - t1;
      float u1 = zTrace - z0;
      float u0 = 1.0f - u1;

      const int i000 = x0 + y0 * sx + z0 * slabSize;
      const int i010 = x0 + y1 * sx + z0 * slabSize;
      const int i100 = x1 + y0 * sx + z0 * slabSize;
      const int i110 = x1 + y1 * sx + z0 * slabSize;
      const int i001 = x0 + y0 * sx + z1 * slabSize;
      const int i011 = x0 + y1 * sx + z1 * slabSize;
      const int i101 = x1 + y0 * sx + z1 * slabSize;
      const int i111 = x1 + y1 * sx + z1 * slabSize;

      // interpolate, (indices could be computed once)
      newField[index] = u0 * (s0 * (
            t0 * oldField[i000] +
            t1 * oldField[i010]) +
          s1 * (t0 * oldField[i100] +
            t1 * oldField[i110])) +
        u1 * (s0 * (t0 * oldField[i001] +
              t1 * oldField[i011]) +
            s1 * (t0 * oldField[i101] +
              t1 * oldField[i111])); 
    }
    */
  } // xyz
}

///////////////////////////////////////////////////////////////////////
// real-valued cell center coordinates
///////////////////////////////////////////////////////////////////////
VECTOR3 VECTOR3_FIELD_3D::cellCenter(int x, int y, int z) const
{
  VECTOR3 halfLengths = (REAL)0.5 * _lengths;

  // set it to the lower corner
  VECTOR3 final = _center - halfLengths;

  // displace to the NNN corner
  final[0] += x * _dx;
  final[1] += y * _dy;
  final[2] += z * _dz;

  // displace it to the cell center
  final[0] += _dx * 0.5;
  final[1] += _dy * 0.5;
  final[2] += _dz * 0.5;

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::resizeAndWipe(int xRes, int yRes, int zRes, const VECTOR3& center, const VECTOR3& lengths)
{
  if (_xRes == xRes && _yRes == yRes && _zRes == zRes)
  {
    _center = center;
    _lengths = lengths;
    clear();

    _dx = _lengths[0] / _xRes;
    _dy = _lengths[1] / _yRes;
    _dz = _lengths[2] / _zRes;
    return;
  }

  if (_data) delete[] _data;

  _xRes = xRes;
  _yRes = yRes;
  _zRes = zRes;
  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  _lengths = lengths;
  _center = center;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  try {
    _data = new VECTOR3[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " VECTOR3_FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(VECTOR3);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _initialized = true;
}

///////////////////////////////////////////////////////////////////////
// normalize all the vectors in the field
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::normalize()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x].normalize();
}

///////////////////////////////////////////////////////////////////////
// copy values out into the border, assuming that "borderSize" is the 
// width of the grid padding
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyIntoBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x == borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x - i,y,z) = value;
        }
        if (x == _xRes - 1 - borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x + i,y,z) = value;
        }              
        if (y == borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y - i,z) = value;
        }
        if (y == _yRes - 1 - borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y+i,z) = value;
        }              
        if (z == borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z - i) = value;
        }
        if (z == _zRes - 1 - borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z+i) = value;
        }

        // handle the corners
        if (x == borderSize && z == borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y,z-j) = value;
        }
        if (x == _xRes - 1 - borderSize && z == _zRes - 1 - borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y,z+j) = value;
        }

        if (z == borderSize && y == borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y - j,z -i) = value;
        }
        if (z == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          VECTOR3 value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y + j,z+i) = value;
        }

      }
}

///////////////////////////////////////////////////////////////////////
// pass back a field with a new padding of size "paddingSize"
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::withAddedPadding(int paddingSize) const
{
  // new length, with padding
  VECTOR3 newLengths = _lengths;
  newLengths[0] += paddingSize * 2 * _dx;
  newLengths[1] += paddingSize * 2 * _dx;
  newLengths[2] += paddingSize * 2 * _dx;

  VECTOR3_FIELD_3D final(_xRes + 2 * paddingSize, 
                         _yRes + 2 * paddingSize, 
                         _zRes + 2 * paddingSize, _center, newLengths);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x + paddingSize,
              y + paddingSize,
              z + paddingSize) = (*this)(x,y,z);

  final.copyIntoBorder(paddingSize);

  return final;
}

//////////////////////////////////////////////////////////////////////
// compute closest point field
//
// For each cell center in input, compute the position of the closest point
// in the surface field
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::computeClosestPoints(const FIELD_3D& input, const FIELD_3D& surfaceField, bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  VECTOR3_FIELD_3D final = VECTOR3_FIELD_3D::cellCenters(input);
  VECTOR3_FIELD_3D targetGradient = VECTOR3_FIELD_3D::gradient(surfaceField);
  
  // DEBUG
  //int maxSteps = 1000;
  int maxSteps = 100;
  //int maxSteps = 10;
  //int maxSeen = 0;

  int xRes = input.xRes();
  int yRes = input.yRes();
  int zRes = input.zRes();

  REAL fictionalDt = 1.0;
  //REAL fictionalDt = 0.95;
  //REAL fictionalDt = 0.5;

  int maxOuts = 0;
  int totalIterations = 0;

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        REAL diff = 1;
        int steps = 0;
        VECTOR3& position = final(x,y,z);
        VECTOR3 targetDelta = targetGradient(position);

        while (diff > 1e-6 && steps < maxSteps)
        {
          REAL targetDistance = surfaceField(position);
          diff = fabs(targetDistance);

          // if the gradient is zero, stop trying, because it probably
          // wandered off the grid
          //if (norm2(targetDelta) < 1e-6) break;
          if (targetDelta.norm() < 1e-6) break;

          // go ahead and always do first -- second gives the occasional wacky value
          // that throws off the stability of the simulation
          targetDelta.normalize();
          //VECTOR3 move = (targetDistance) * targetDelta;
          VECTOR3 move = (targetDistance) * targetDelta * fictionalDt;

          position = position - move;

          // update the direction
          targetDelta = targetGradient(position);
          steps++;
        }

        /*
        if (steps == maxSteps)
        {
#pragma omp atomic 
          maxOuts++;
        }
        */
#pragma omp atomic
        totalIterations += steps; 
      }

  cout << " Closest point iteration (full grid) maxed out " << 100.0f * (REAL)maxOuts / (xRes * yRes * zRes) << "% of the time" << endl;
  cout << " Mean number of iterations: " << (REAL)totalIterations / (xRes * yRes * zRes) << endl;

  return final;
}

//////////////////////////////////////////////////////////////////////
// compute closest point field
//
// For each cell center in input, compute the position of the closest point
// in the surface field
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::computeClosestPointsNarrowBand(const FIELD_3D& input, const FIELD_3D& surfaceField, const int maxCells, VECTOR3_FIELD_3D& final)
{
  TIMER functionTimer(__FUNCTION__);
  const VECTOR3_FIELD_3D targetGradient = VECTOR3_FIELD_3D::gradient(surfaceField);
  
  const REAL invDx = 1.0 / surfaceField.dx();

  /*
  // find which cells are inside the band
  vector<int> xs;
  vector<int> ys;
  vector<int> zs;
  for (int z = 0; z < surfaceField.zRes(); z++)
    for (int y = 0; y < surfaceField.yRes(); y++)
      for (int x = 0; x < surfaceField.xRes(); x++)
      {
        REAL currentDistance = fabs(surfaceField(x,y,z) * invDx);
        if (currentDistance < maxCells)
        {
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
        }
      }
      */

  // DEBUG
  const int maxSteps = 100;
  const int xRes = input.xRes();
  const int yRes = input.yRes();
  const int zRes = input.zRes();

  //REAL fictionalDt = 1;
  //REAL fictionalDt = 0.5;

  const VECTOR3 lowerCorner = input.center() - (REAL)0.5 * input.lengths(); 
  const VECTOR3 upperCorner = input.center() + (REAL)0.5 * input.lengths(); 

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        REAL currentDistance = fabs(surfaceField(x,y,z) * invDx);
        if (currentDistance >= maxCells) continue;

        REAL diff = 1;
        int steps = 0;
        VECTOR3& position = final(x,y,z);
        position = input.cellCenter(x,y,z);
        VECTOR3 targetDelta = targetGradient(position);

        while (diff > 1e-6 && steps < maxSteps)
        {
          REAL targetDistance = surfaceField(position);
          diff = fabs(targetDistance);

          // if the gradient is zero, stop trying, because it probably
          // wandered off the grid
          //if (norm2(targetDelta) < 1e-6) break;

          // go ahead and always do first -- second gives the occasional wacky value
          // that throws off the stability of the simulation
          targetDelta.normalize();
          VECTOR3 move = (targetDistance) * targetDelta;
          //VECTOR3 move = (targetDistance) * targetDelta * fictionalDt;

          position = position - move;

          // update the direction
          targetDelta = targetGradient(position);
          steps++;
        }

        // clamp position to something reasonable
        position[0] = (position[0] < lowerCorner[0]) ? lowerCorner[0] : position[0]; 
        position[1] = (position[1] < lowerCorner[1]) ? lowerCorner[1] : position[1]; 
        position[2] = (position[2] < lowerCorner[2]) ? lowerCorner[2] : position[2]; 

        position[0] = (position[0] > upperCorner[0]) ? upperCorner[0] : position[0]; 
        position[1] = (position[1] > upperCorner[1]) ? upperCorner[1] : position[1]; 
        position[2] = (position[2] > upperCorner[2]) ? upperCorner[2] : position[2]; 
      }
}

//////////////////////////////////////////////////////////////////////
// getting nans and infs from reading in the Houdini file -- stomp them here.
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::fixNanInfs()
{
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++)
          if (isinf((*this)(x,y,z)[i]) || isnan((*this)(x,y,z)[i]))
            (*this)(x,y,z)[i] = 0;
}

} // HOBAK
