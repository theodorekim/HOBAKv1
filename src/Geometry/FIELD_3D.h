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
#ifndef FIELD_3D_H
#define FIELD_3D_H

#include <cmath>
#include <string>
#include <map>
#include <unistd.h>

#include "FIELD_2D.h"
#include "MIN_HEAP.h"

#ifndef VARNAME
#define VARNAME(x) #x
#endif

// macro to debug 3D fields, only the first one has been ported
// into HOBAK so far
#ifndef FIELDVIEW3D
#define FIELDVIEW3D(x) FIELD_3D::fieldViewer(x, VARNAME(x)); sleep(1);
#endif
#ifndef OVERLAYFIELDVIEW3D
#define OVERLAYFIELDVIEW3D(x,y) FIELD_3D::overlayFieldViewer(x, y, VARNAME(x)); sleep(1);
#endif
#ifndef FIELDVIEW3DYZ
#define FIELDVIEW3DYZ(x) FIELD_3D::fieldViewerYZ(x, VARNAME(x)); sleep(1);
#endif
#ifndef OVERLAYFIELDVIEW3DYZ
#define OVERLAYFIELDVIEW3DYZ(x,y) FIELD_3D::overlayFieldViewerYZ(x, y, VARNAME(x)); sleep(1);
#endif

namespace HOBAK {

class FIELD_3D 
{
public:
  FIELD_3D();
  FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VECTOR3& center = VECTOR3(0,0,0), const VECTOR3& lengths = VECTOR3(1,1,1));
  FIELD_3D(const double* data, const int& xRes, const int& yRes, const int& zRes, const VECTOR3& center = VECTOR3(0,0,0), const VECTOR3& lengths = VECTOR3(1,1,1));
  FIELD_3D(const FIELD_3D& m);
  FIELD_3D(const std::vector<FIELD_2D>& slices);
  FIELD_3D(const char* filename);
  ~FIELD_3D();

  // accessors
  inline REAL& operator()(int x, int y, int z) { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  const REAL operator()(int x, int y, int z) const { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  const REAL operator()(const VECTOR3& position) const; 

  inline REAL& operator[](int x) { return _data[x]; };
  const REAL operator[](int x) const { return _data[x]; };
  REAL* data() { return _data; };
  const REAL* dataConst() const { return _data; };
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int zRes() const { return _zRes; };
  VECTOR3 dxs() const { return VECTOR3(_dx, _dy, _dz); };
  VECTOR3 invDxs() const { return VECTOR3(_invDx, _invDy, _invDz); };
  const REAL dx() const { return _dx; };
  const REAL dy() const { return _dy; };
  const REAL dz() const { return _dz; };
  const REAL invDx() const { return _invDx; };
  const REAL invDy() const { return _invDy; };
  const REAL invDz() const { return _invDz; };
  const int slabSize() const { return _slabSize; };
  const VECTOR3 center() const { return _center; };
  const VECTOR3 lengths() const { return _lengths; };
  const int totalCells() const { return _totalCells; };
  const int outside() const { return _outside; };
  const int& quinticClamps() const { return _quinticClamps; };
  const bool initialized() const;
  const int maxRes() const { return (_xRes > _yRes) ? ((_zRes > _xRes) ? _zRes : _xRes) : ((_zRes > _yRes) ? _zRes : _yRes); };

  // hack to make Houdini work
  void setDx(const REAL newDx) { _dx = newDx; };

  // reset dimensions
  void setCenter(const VECTOR3& center) { _center = center; };
  void setLengths(const VECTOR3& lengths);
  
  // various order lookup functions
  REAL sexticLookup(const VECTOR3& position) const;
  REAL sexticLookupClamped(const VECTOR3& position) const;
  REAL quinticLookup(const VECTOR3& position) const;
  REAL quarticLookup(const VECTOR3& position) const;
  REAL quarticLookupInlined(const VECTOR3& position) const;
  REAL quarticLookupClamped(const VECTOR3& position) const;
  REAL cubicLookup(const VECTOR3& position) const;
  REAL cubicLookupUnclamped(const VECTOR3& position) const;
  REAL cubicNewtonLookup(const VECTOR3& position) const;
  REAL lerpDebug(const VECTOR3& position, int x, int y, int z) const;
  REAL nearestNeighborLookup(const VECTOR3& position) const;

  void clear();
  
  // real-valued cell center coordinates
  VECTOR3 cellCenter(int x, int y, int z) const;
  
  // cell index of a real-valued position
  int cellIndex(VECTOR3& position) const;

  // overloaded operators
  FIELD_3D& operator=(const REAL& alpha);
  FIELD_3D& operator=(const FIELD_3D& A);
  FIELD_3D& operator*=(const REAL& alpha);
  FIELD_3D& operator/=(const REAL& alpha);
  FIELD_3D& operator+=(const REAL& alpha);
  FIELD_3D& operator-=(const REAL& alpha);
  FIELD_3D& operator-=(const FIELD_3D& input);
  FIELD_3D& operator+=(const FIELD_3D& input);
  FIELD_3D& operator*=(const FIELD_3D& input);

  // BLAS-like interface, output += alpha * input
  static void axpy(const REAL& alpha, const FIELD_3D& input, FIELD_3D& output);

  // IO functions
  void write(std::string filename) const;
  void read(std::string filename);
  void readHoudini(std::string filename);
  void readHoudiniSurf(std::string filename);
  static void readHoudiniVel(const std::string filename, FIELD_3D& xVelocity, FIELD_3D& yVelocity, FIELD_3D& zVelocity);

  // This is hard-coded for the paddle example
  void readHoudini12Surf(std::string filename);

  // streaming IO functions
  void write(FILE* file) const;
  void read(FILE* file);
	
  // load a PhysBAM level set
  void readPhysBAM(const char* filename);

  void resizeAndWipe(int xRes, int yRes, int zRes, const VECTOR3& center = VECTOR3(0,0,0), const VECTOR3& lengths = VECTOR3(1,1,1));

  // return a slice in the form of a FIELD_2D
  FIELD_2D zSlice(int z) const;
  
  // return the projection of the field in different directions
  FIELD_2D zProjection();
  FIELD_2D yProjection();
  FIELD_2D xProjection();
  
  // do a union with 'field', assuming both this and field are signed distance fields
  FIELD_3D signedDistanceUnion(const FIELD_3D& field);

  // what's the maximum resolution in any direction?
  int maxRes();

  // assuming that SURFACE.initializeSignedDistanceField() has been called on this field,
  // do the fast marching method
  void fastMarchingMethod();

  // only fast march the inside (negative SDF) of the fluid to handle immersed obstacles
  void fastMarchingNegativeOnly();
  
  // only fast march the outside (positive SDF) of the fluid to handle immersed obstacles
  void fastMarchingPositiveOnly();

  // extend some scalar quantity off of a front, given a signed distance function
  void fastExtension(const FIELD_3D& signedDistance);
  
  // return the indices of the grid points insides a world-space bounding box
  void boundingBoxIndices(const VECTOR3& mins, const VECTOR3& maxs, VECTOR3I& iMins, VECTOR3I& iMaxs);

  // norms
  REAL sumSq();
  REAL max();
  REAL absMax();

  // build a const field with the given dims
  static FIELD_3D constField(const FIELD_3D& dims, REAL value);

  // check if any entry is a nan
  bool isNan();

  // compute the inverse
  FIELD_3D inverse();

  // print the neighborhood of a cell for debugging
  void printNeighborhood(int x, int y, int z) const;
  void printNeighborhood(int index) const;

  // clamp nans to some specified value
  void clampNans(REAL value = 0);

  // clamp nansto values in this field
  void clampNans(FIELD_3D& clampField);

  // clamp infinities to some specified value
  void clampInfs(REAL value = 0);
  
  // clamp infinities to values in this field
  void clampInfs(FIELD_3D& clampField);

  // dump to a viewer
  static void fieldViewer(const FIELD_3D& field, std::string name);
  static void fieldViewerYZ(const FIELD_3D& field, std::string name);
  
  static void overlayFieldViewer(const FIELD_3D& field, const FIELD_3D& distance, std::string name);
  static void overlayFieldViewerYZ(const FIELD_3D& field, const FIELD_3D& distance, std::string name);

  // get the integer indices of a spatial position
  void indices(const VECTOR3& position, int* x);

  // set to a checkerboard solid texture
  void setToSolidCheckboard(int xChecks = 10, int yChecks = 10, int zChecks = 10);
  void setToGrayCheckerboard(int xChecks = 10, int yChecks = 10, int zChecks = 10);

  // set to vertical derivative kernel
  void setToVerticalDerivativeKernel(double kMax = 10, double dk = 0.01, double sigma = 1.0, double L = 1);
  
  // set whole field to a Gaussian
  void setToGaussian(REAL amplitude = 1.0, VECTOR3 sigmas = VECTOR3(0.1, 0.1, 0.1));
  
  // set to a Gaussian
  void insertGaussian(const VECTOR3& center, const REAL amplitude = 1.0, const VECTOR3 sigmas = VECTOR3(0.1, 0.1, 0.1));

  // convolve this field with a smaller field
  FIELD_3D convolve(const FIELD_3D& filter);
  FIELD_3D convolveFast15(const FIELD_3D& filter);
  FIELD_3D convolveToroidal(const FIELD_3D& filter);
  FIELD_3D convolveNarrowBand(const FIELD_3D& filter, const FIELD_3D& distance, int maxCells);
  FIELD_3D convolveNarrowBand(const FIELD_3D& filter, const std::vector<int>& narrowBand);
  FIELD_3D convolveNarrowBandFast15(const FIELD_3D& filter, const std::vector<int>& narrowBand);
  FIELD_3D convolveNarrowBandFast15(const FIELD_3D& filter, const FIELD_3D& distance, const REAL maxRadius);

  // stomp all the value outside a narrow band to zero in order to boost compression
  void stompOutsideNarrowBand(const std::vector<int>& narrowBand);

  // compute the narrow band indices for this object, which is assumed to be a SDF
  // the returned vector is (x,y,z) triplets
  std::vector<int> computeNarrowBand(REAL maxCellDistance) const;

  // sum of the entire field
  REAL sum();

  // vector of the field's dimesions
  VECTOR3 dims() const { return VECTOR3(_xRes, _yRes, _zRes); };

  // determine how many non-zero entries are in the filter
  int nonZeroEntries();

  // set a given z slice to the given 2D field
  void setSliceZ(const int& z, const FIELD_2D& slice);

  // normalize the data
  void normalize();

  // field minimum
  REAL fieldMin();

  // field maximum
  REAL fieldMax();
  
  // field maximum cell index
  VECTOR3 maxIndex();

  // field minimum cell index
  VECTOR3 minIndex();

  // flip the z and y coordinates
  FIELD_3D flipZY() const;
  
  // flip the x and y coordinates
  FIELD_3D flipXY() const;

  // flip the x and z coordinates
  FIELD_3D flipXZ() const;
  
  // create a mirror image along different axes
  FIELD_3D mirrorZ() const;
  FIELD_3D mirrorY() const;
  FIELD_3D mirrorX() const;

  // copy out the boundary
  void copyBorderAll();

  // copy values out into the border, assuming that "borderSize" is the width of the grid padding
  void copyIntoBorder(int borderSize);

  // first order spatial derivatives
  // on the border, difference falls back to first order (not centered) difference
  inline REAL Dx(int x, int y, int z) const; 
  inline REAL Dy(int x, int y, int z) const; 
  inline REAL Dz(int x, int y, int z) const; 
  
  // second order spatial derivatives
  // on the border, center cell is copied to outside, and centered difference
  // is still taken
  inline REAL DDx(int x, int y, int z) const; 
  inline REAL DDy(int x, int y, int z) const; 
  inline REAL DDz(int x, int y, int z) const; 

  // mixed derivatives
  // on the border, center cell is copied to outside, and centered difference
  // is still taken
  inline REAL DDxy(int x, int y, int z) const;
  inline REAL DDxz(int x, int y, int z) const; 
  inline REAL DDyz(int x, int y, int z) const;

  // evaluate laplace operator at cell
  REAL laplace(int x, int y, int z) const; 

  // get a field of the entire derivative
  FIELD_3D Dx() const;
  FIELD_3D Dy() const;
  FIELD_3D Dz() const;
  FIELD_3D DDx() const;
  FIELD_3D DDy() const;
  FIELD_3D DDz() const;
  FIELD_3D DDxy() const;
  FIELD_3D DDxz() const;
  FIELD_3D DDyz() const;

  // get the curvature
  FIELD_3D meanCurvature() const;
  FIELD_3D gaussianCurvature() const;
  FIELD_3D principalCurvature() const;
  void principalCurvatures(FIELD_3D& minCurvature, FIELD_3D& maxCurvature) const;

  // compute laplace everywhere (default factor scales to similar wave prop. speed like iwave for wave equation)
  FIELD_3D computeLaplace(REAL factor = (0.001 * -0.25) ) const;

  // mask out any values past a certain distance
  void maskByDistance(const FIELD_3D& distanceField, const REAL distance);

  // clamp the field to a min and max
  void clamp(const REAL minValue, const REAL maxValue);

  // get a resampled version
  FIELD_3D resampleCubic(int xRes, int yRes, int zRes) const;
  FIELD_3D resampleCubicUnclamped(int xRes, int yRes, int zRes) const;
  FIELD_3D resampleCubicUnclampedNarrowBand(int upResFactor, const FIELD_3D& distanceField, const int maxRadius) const;
  FIELD_3D resampleQuintic(int xRes, int yRes, int zRes) const;
  FIELD_3D resampleSextic(int xRes, int yRes, int zRes) const;

  // take the square root of the field
  FIELD_3D squareRoot();

  // get the normal at a point
  VECTOR3 normal(int x, int y, int z) const;

  // get the derivative of the normal at a point
  MATRIX3 Dnormal(int x, int y, int z) const;

  // do a cubic Hermite interpolation
  static REAL cubicInterp(const REAL interp, const REAL* points);
  
  // do a cubic Hermite that clamps to the immediate neighborhood
  static REAL cubicInterpClamped(const REAL interp, const REAL* points);

  // do a cubic Hermite interpolation, but turn off monotonic clamping
  static REAL cubicInterpUnclamped(const REAL interp, const REAL* points);
 
  // do a quartic WENO interpolation
  static REAL quarticInterp(const REAL interp, const REAL* points);
  static REAL quarticInterpClamped(const REAL interp, const REAL* points);

  // do a quintic Hermite interpolation
  static REAL quinticInterp(const REAL interp, const REAL* points);
 
  // do a WENO6 interpolation
  static REAL sexticInterp(const REAL interp, const REAL* points);
  static REAL sexticInterpClamped(const REAL interp, const REAL* points);

  // stomp to zero anything in the field that is not between the min and max
  void bandPass(const REAL minValue, const REAL maxValue);
  
  // isolate values near the current value, within a certain width
  void isolateBand(const REAL target, const REAL width);

  // single explicit diffusion step
  void blur(REAL dt = 1.0);

  // get band-limited curvature, with some diffusion
  FIELD_3D bandLimitedCurvature(const REAL target, const REAL width);

  // set to the absolute value
  void absoluteValue();

  // do a soft bandpass where there's a gradual falloff
  void softBandPass(const REAL band, const REAL falloff);

  // pass back a field with a new padding of size "paddingSize"
  FIELD_3D withAddedPadding(int paddingSize) const;

  // stomp the border to zero
  void stompBorder(int borderSize);

  // set the border to value
  void setBorder(int borderSize, REAL value);

  // stomp corners to zero
  void stompCorners(int borderSize);

  // do a diff between two fields, but only along a narrow band
  FIELD_3D narrowBandDiff(const FIELD_3D& rhs, const std::vector<int>& narrowBand) const;

  // do a diff between two fields, but only along a narrow band, and not along
  // borders where the filter skipped
  FIELD_3D filteredNarrowBandDiff(const FIELD_3D& rhs, const std::vector<int>& narrowBand, const FIELD_3D& filter);

  // return a field showing the narrow band
  FIELD_3D narrowBandField(int bandwidth) const;

  // set everything in the specified interval to a given value
  void setInterval(const int xMin, const int xMax, const int yMin, const int yMax, const int zMin, const int zMax, const REAL value = 0);

  // stomp everything outside a narrow band to zero
  void stompOutsideNarrowBand(const FIELD_3D& distance, const int maxRadius);

private:
  int _xRes;
  int _yRes;
  int _zRes;

  int _slabSize;
  int _totalCells;

  REAL* _data;

  // what fast marching considers "outside"
  int _outside;

  // center position of the grid
  VECTOR3 _center;

  // lengths of the x,y,z dimensions of the grid
  VECTOR3 _lengths;

  // physical lengths
  REAL _dx;
  REAL _dy;
  REAL _dz;

  REAL _invDx;
  REAL _invDy;
  REAL _invDz;

  // retirement hash for fast marching
  // note that entry existence is the important things, not whether
  // the value is true or false. If the entry even exists, it was retired
  std::map<int, bool> _retired;

  // track how many times the quintic is clamped
  static int _quinticClamps;

  // do fast marching in one direction
  void marchOneway(bool forward, MIN_HEAP& minHeap);
  
  void marchOneway2ndOrder(bool forward, MIN_HEAP& minHeap);

  // insert the front in preparation for reinitialization or extension
  void insertFront(const bool forward, FIELD_3D& distance, MIN_HEAP& minHeap);
  
  // do fast extension in one direction
  void extendOneway(bool forward, FIELD_3D& distance, MIN_HEAP& minHeap);

  // do a cubic Newton interpolation
  static REAL cubicNewtonInterp(REAL interp, REAL* points);

  // read a Houdini field off from a file stream -- it is assumed that the
  // file is already advanced to the beginning of the field
  void readHoudiniField(FILE* file, bool storeValues);

  // a static version for reading in velocity
  static FIELD_3D readHoudiniField(FILE* file);
};

FIELD_3D operator^(const FIELD_3D& A, const REAL alpha);
FIELD_3D operator*(const FIELD_3D& A, const REAL alpha);
FIELD_3D operator/(const FIELD_3D& A, const REAL alpha);
FIELD_3D operator/(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator+(const FIELD_3D& A, const REAL alpha);
FIELD_3D operator*(const REAL alpha, const FIELD_3D& A);
FIELD_3D operator+(const REAL alpha, const FIELD_3D& A);
FIELD_3D operator-(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator+(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator*(const FIELD_3D& A, const FIELD_3D& B);

} // HOBAK

#endif
