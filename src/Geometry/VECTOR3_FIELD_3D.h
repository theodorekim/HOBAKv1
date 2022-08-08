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
#ifndef VECTOR3_FIELD_3D_H
#define VECTOR3_FIELD_3D_H

#include <cmath>
#include <string>
#include <map>
#include <iostream>

#include "SETTINGS.h"
#include "FIELD_3D.h"
#include "util/TIMER.h"

#ifndef VECTORFIELDVIEW3D
#define VECTORFIELDVIEW3D(x,y) VECTOR3_FIELD_3D::fieldViewer(x, y, VARNAME(x)); sleep(1);
#endif

#ifndef OVERLAYVECTORFIELDVIEW3D
#define OVERLAYVECTORFIELDVIEW3D(x,y) VECTOR3_FIELD_3D::overlayFieldViewer(x, y, VARNAME(x)); sleep(1);
#endif

#ifndef CLOSESTPOINTVIEW3D
#define CLOSESTPOINTVIEW3D(x,y) VECTOR3_FIELD_3D::closestPointViewer(x, y, VARNAME(x)); sleep(1);
#endif

namespace HOBAK {

using namespace std;

class VECTOR3_FIELD_3D {
public:
  VECTOR3_FIELD_3D();
  VECTOR3_FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VECTOR3& center = VECTOR3(0,0,0), const VECTOR3& lengths = VECTOR3(1,1,1));
  VECTOR3_FIELD_3D(double* data, const int& xRes, const int& yRes, const int& zRes, const VECTOR3& center = VECTOR3(0,0,0), const VECTOR3& lengths = VECTOR3(1,1,1));
  VECTOR3_FIELD_3D(const VECTOR3_FIELD_3D& m);
  VECTOR3_FIELD_3D(const FIELD_3D& m);
  ~VECTOR3_FIELD_3D();

  // reallocate everything
  void resizeAndWipe(int xRes, int yRes, int zRes, const VECTOR3& center = VECTOR3(0,0,0), const VECTOR3& lengths = VECTOR3(1,1,1));

  // accessors
  inline VECTOR3& operator()(int x, int y, int z) { 
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  const VECTOR3 operator()(int x, int y, int z) const 
  { 
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  const VECTOR3 operator()(const VECTOR3& position) const;
  inline VECTOR3& operator[](int x) { 
    assert(x >= 0);
    assert(x < _totalCells);
    return _data[x]; 
  };
  const VECTOR3 operator[](int x) const { 
    assert(x >= 0);
    assert(x < _totalCells);
    return _data[x]; 
  };
  VECTOR3* data() { return _data; };
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int zRes() const { return _zRes; };
  const VECTOR3 dims() const { return VECTOR3(_xRes, _yRes, _zRes); };
  const REAL dx() const { return _dx; };
  const REAL dy() const { return _dy; };
  const REAL dz() const { return _dz; };
  VECTOR3 dxs() const { return VECTOR3(_dx, _dy, _dz); };
  const int slabSize() const { return _slabSize; };
  const VECTOR3 center() const { return _center; };
  const VECTOR3 lengths() const { return _lengths; };
  const int totalCells() const { return _totalCells; };
  const bool initialized() const { return _initialized; };
  const VECTOR3 constEntry(int index) const { return _data[index]; };
  void clear();
 
  // reset the lengths to something else, and recompute all the dimesions as well
  void setLengths(const VECTOR3& lengths);
  
  void setCenter(const VECTOR3& center) { _center = center; };

  // what's the maximum resolution in any direction?
  int maxRes();

  // create a field of the grid positions of the passed in grid
  static VECTOR3_FIELD_3D cellCenters(const FIELD_3D& input);
  
  // take the gradient of a scalar field
  static VECTOR3_FIELD_3D gradient(const FIELD_3D& input);

  // return a grid of values at the given spatial positions
  static FIELD_3D compose(const FIELD_3D& values, const VECTOR3_FIELD_3D& positions);
  static VECTOR3_FIELD_3D compose(const VECTOR3_FIELD_3D& values, const VECTOR3_FIELD_3D& positions);

  // retrieve the components
  FIELD_3D scalarField(int component) const;
  FIELD_3D magnitudeField() const;

  // norms
  REAL sumMagnitudes();
  REAL maxMagnitudes();

  // check if any entry is a nan
  bool isNan();

  // overloaded operators
  VECTOR3_FIELD_3D& operator-=(const VECTOR3_FIELD_3D& input);
  VECTOR3_FIELD_3D& operator+=(const VECTOR3_FIELD_3D& input);
  VECTOR3_FIELD_3D& operator*=(const REAL& value);
  VECTOR3_FIELD_3D& operator=(const REAL& value);
  VECTOR3_FIELD_3D& operator=(const VECTOR3_FIELD_3D& input);

  // extend some vector quantity off of a front, given a signed distance function
  void fastExtension(const FIELD_3D& signedDistance);

  // set the values in the field to the values at the closest points
  static FIELD_3D setToClosestPointValues(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints);
  static FIELD_3D setToClosestPointValuesNarrowBand(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int maxCells);
  static FIELD_3D setToClosestPointValuesNarrowBandFrozenCore(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int maxRadius, int coreRadius = 1);
  static FIELD_3D setToClosestPointValuesFrozenCore(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int coreRadius = 1);

  // real-valued cell center coordinates
  VECTOR3 cellCenter(int x, int y, int z) const;

  // file IO
  void write(const string filename) const;
  void read(const string filename);
  void write(FILE* file) const;
  void read(FILE* file);
  void write(FILE* file, const VECTOR3& v) const;
  void read(FILE* file, VECTOR3& v);

  /*
  // draw to GL, with the cell centers of 'field' as the vector origins
  void draw(const VECTOR3_FIELD_3D& origins) const;
  void drawZSlice(const VECTOR3_FIELD_3D& origins, const int zSlice, const REAL scale = 1, const int stride = 4) const;
  void drawClosestPointZSlice(const VECTOR3_FIELD_3D& origins, const int zSlice, const REAL scale = 1, const int stride = 4) const;
  */

  // dump to a viewer
  static void fieldViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name);
  static void overlayFieldViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name);
  static void closestPointViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name);

  // flip the z and y coordinates
  VECTOR3_FIELD_3D flipZY() const;
  
  // flip the x and y coordinates
  VECTOR3_FIELD_3D flipXY() const;
  
  // flip the x and z coordinates
  VECTOR3_FIELD_3D flipXZ() const;

  // advect using first order semi-Lagrangian
  static void advect(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField);
  static void advectNarrowBand(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField, const FIELD_3D& distance, const int maxCells);
  static void advectNarrowBandLinear(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField, const FIELD_3D& distance, const int maxCells);
  static void advectLazy(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField, FIELD_3D& newField);
  
  // advect with second order MacCormack
  static void advectMacCormack(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2);
  static void advectMacCormackLazy(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2);
  static void advectMacCormackNarrowBand(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2);

  // version of the () operator for debugging purposes
  VECTOR3 debugPositionOperator(const VECTOR3& position) const;

  // normalize all the vectors in the field
  void normalize();

  // copy values out into the border, assuming that "borderSize" is the width of the grid padding
  void copyIntoBorder(int borderSize);

  // pass back a field with a new padding of size "paddingSize"
  VECTOR3_FIELD_3D withAddedPadding(int paddingSize) const;

  // compute closest point field
  //
  // For each cell center in input, compute the position of the closest point
  // in the surface field
  static VECTOR3_FIELD_3D computeClosestPoints(const FIELD_3D& input, const FIELD_3D& surfaceField, bool verbose = true);
  static void computeClosestPointsNarrowBand(const FIELD_3D& input, const FIELD_3D& surfaceField, const int maxCells, VECTOR3_FIELD_3D& final);

  // compute the extension field for some subset of points
  static FIELD_3D computeExtensionFieldSubset(const FIELD_3D& distanceField, const FIELD_3D& toExtend, const vector<int>& subset);
  
  // compute the extension field for some subset of points
  static FIELD_3D computeExtensionFieldMasked(const FIELD_3D& distanceField, const FIELD_3D& toExtend, const FIELD_3D& mask);
  static FIELD_3D computeExtensionFieldMaskedBadAdvect(const FIELD_3D& distanceField, const FIELD_3D& toExtend, const FIELD_3D& mask);

  // do an advection of a subset of the grid
  static FIELD_3D advectSubset(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, const vector<int>& subset);

  // set the closes point values needed for MacCormack
  static FIELD_3D setMacCormackClosestPoints(const REAL dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField);

  // Clamp the extrema generated by the BFECC error correction
  static void clampExtrema(const REAL dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, FIELD_3D& newField);
  static void clampOutsideRays(const REAL dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, const FIELD_3D& oldAdvection, FIELD_3D& newField);
  static void clampExtremaSubset(const REAL dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, FIELD_3D& newField, const vector<int>& subset);
  static void clampOutsideRaysSubset(const REAL dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, const FIELD_3D& oldAdvection, FIELD_3D& newField, const vector<int>& subset);

private:
  int _xRes;
  int _yRes;
  int _zRes;

  int _slabSize;
  int _totalCells;

  VECTOR3* _data;

  // center position of the grid
  VECTOR3 _center;

  // lengths of the x,y,z dimensions of the grid
  VECTOR3 _lengths;

  // physical lengths
  REAL _dx;
  REAL _dy;
  REAL _dz;

  // has this field been allocated?
  bool _initialized;
  
  // retirement hash for fast extension
  // note that entry existence is the important things, not whether
  // the value is true or false. If the entry even exists, it was retired
  map<int, bool> _retired;
  
  // insert the front in preparation for reinitialization or extension
  void insertFront(const bool forward, FIELD_3D& distance, MIN_HEAP& minHeap);
  
  // do fast extension in one direction
  void extendOneway(bool forward, FIELD_3D& distance, MIN_HEAP& minHeap);

  // get a single closest point
  static REAL getClosestPointValue(const VECTOR3_FIELD_3D& targetGradient, const VECTOR3& startPosition, const FIELD_3D& distanceField, const FIELD_3D& toExtend);

  // getting nans and infs from reading in the Houdini file -- stomp them here.
  void fixNanInfs();
};

// take the field dot product
FIELD_3D operator*(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v);
VECTOR3_FIELD_3D operator*(const FIELD_3D& u, const VECTOR3_FIELD_3D& v);
VECTOR3_FIELD_3D operator+(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v);

// diff two vector fields
VECTOR3_FIELD_3D operator-(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v);

} // HOBAK

#endif
