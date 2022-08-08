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
#include "COLLISION_UTIL.h"
#include <iostream>

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
// gradient of the cross product used to compute the triangle normal,
// vertex-face case
///////////////////////////////////////////////////////////////////////
MATRIX3x12 crossGradientVF(const std::vector<VECTOR3>& e)
{
  MATRIX3x12 crossMatrix;

  const REAL e0x = e[0][0];
  const REAL e0y = e[0][1];
  const REAL e0z = e[0][2];

  const REAL e2x = e[2][0];
  const REAL e2y = e[2][1];
  const REAL e2z = e[2][2];

  crossMatrix.col(0) = VECTOR3(0, 0, 0); 
  crossMatrix.col(1) = VECTOR3(0, 0, 0); 
  crossMatrix.col(2) = VECTOR3(0, 0, 0); 
  crossMatrix.col(3) = VECTOR3(0, -e0z, e0y); 
  crossMatrix.col(4) = VECTOR3(e0z, 0, -e0x);
  crossMatrix.col(5) = VECTOR3(-e0y, e0x, 0);
  crossMatrix.col(6) = VECTOR3(0, (e0z - e2z), (-e0y + e2y));
  crossMatrix.col(7) = VECTOR3((-e0z + e2z), 0, (e0x - e2x));
  crossMatrix.col(8) = VECTOR3((e0y - e2y), (-e0x + e2x), 0);
  crossMatrix.col(9) = VECTOR3(0, e2z, -e2y);
  crossMatrix.col(10) = VECTOR3(-e2z, 0, e2x);
  crossMatrix.col(11) = VECTOR3(e2y, -e2x, 0);

  return crossMatrix;
}

///////////////////////////////////////////////////////////////////////
// gradient of the cross product used to compute the normal,
// edge-edge case
///////////////////////////////////////////////////////////////////////
MATRIX3x12 crossGradientEE(const std::vector<VECTOR3>& e)
{
  MATRIX3x12 crossMatrix;

  const REAL e0x = e[0][0];
  const REAL e0y = e[0][1];
  const REAL e0z = e[0][2];

  const REAL e1x = e[1][0];
  const REAL e1y = e[1][1];
  const REAL e1z = e[1][2];

  crossMatrix.col(0) = VECTOR3(0, -e1z, e1y);
  crossMatrix.col(1) = VECTOR3(e1z, 0, -e1x);
  crossMatrix.col(2) = VECTOR3(-e1y, e1x, 0);

  crossMatrix.col(3) = VECTOR3(0, e1z, -e1y);
  crossMatrix.col(4) = VECTOR3(-e1z, 0, e1x);
  crossMatrix.col(5) = VECTOR3(e1y, -e1x, 0);

  crossMatrix.col(6) = VECTOR3(0, e0z, -e0y);
  crossMatrix.col(7) = VECTOR3(-e0z, 0, e0x);
  crossMatrix.col(8) = VECTOR3(e0y, -e0x, 0);

  crossMatrix.col(9)  = VECTOR3(0, -e0z, e0y);
  crossMatrix.col(10) = VECTOR3(e0z, 0, -e0x);
  crossMatrix.col(11) = VECTOR3(-e0y, e0x, 0);

  return crossMatrix;
}

///////////////////////////////////////////////////////////////////////
// gradient of the triangle normal, vertex-face case
///////////////////////////////////////////////////////////////////////
MATRIX3x12 normalGradientVF(const std::vector<VECTOR3>& e)
{
  //crossed = cross(e2, e0);
  VECTOR3 crossed = e[2].cross(e[0]);
  REAL crossNorm = crossed.norm();
  const REAL crossNormCubedInv = 1.0 / pow(crossed.dot(crossed), 1.5);
  MATRIX3x12 crossMatrix = crossGradientVF(e);

  //final = zeros(3,12);
  //for i = 1:12
  //  crossColumn = crossMatrix(:,i);
  //  final(:,i) = (1 / crossNorm) * crossColumn - ((crossed' * crossColumn) / crossNormCubed) * crossed;
  //end
  MATRIX3x12 result;
  for (int i = 0; i < 12; i++)
  {
    const VECTOR3 crossColumn = crossMatrix.col(i);
    result.col(i) = (1.0 / crossNorm) * crossColumn - 
                    ((crossed.dot(crossColumn)) * crossNormCubedInv) * crossed;
  }
  return result;
}

///////////////////////////////////////////////////////////////////////
// gradient of the normal, edge-edge case
///////////////////////////////////////////////////////////////////////
MATRIX3x12 normalGradientEE(const std::vector<VECTOR3>& e)
{
  VECTOR3 crossed = e[1].cross(e[0]);
  const REAL crossNorm = crossed.norm();
  const REAL crossNormInv = (crossNorm > 1e-8) ? 1.0 / crossed.norm() : 0.0;
  const REAL crossNormCubedInv = (crossNorm > 1e-8) ? 1.0 / pow(crossed.dot(crossed), 1.5) : 0.0;
  MATRIX3x12 crossMatrix = crossGradientEE(e);

  MATRIX3x12 result;
  for (int i = 0; i < 12; i++)
  {
    const VECTOR3 crossColumn = crossMatrix.col(i);
    result.col(i) = crossNormInv * crossColumn - 
                    ((crossed.dot(crossColumn)) * crossNormCubedInv) * crossed;
  }
  return result;
}

///////////////////////////////////////////////////////////////////////
// one entry of the rank-3 hessian of the cross product used to compute 
// the triangle normal, vertex-face case
///////////////////////////////////////////////////////////////////////
VECTOR3 crossHessianVF(const int iIn, const int jIn)
{
  int i = iIn;
  int j = jIn;

  if (i > j)
  {
    int temp = j;
    j = i;
    i = temp;
  }

  if ((i == 5 && j == 7)  || (i == 8 && j == 10) || (i == 4 && j == 11))
    return VECTOR3(1, 0, 0);

  if ((i == 6 && j == 11) || (i == 3 && j == 8) || (i == 5 && j == 9))
    return VECTOR3(0, 1, 0);

  if ((i == 4 && j == 6)  || (i == 7 && j == 9) || (i == 3 && j == 10))
    return VECTOR3(0, 0, 1);

  if ((i == 7 && j == 11) || (i == 4 && j == 8) || (i == 5 && j == 10))
    return VECTOR3(-1, 0, 0);

  if ((i == 5 && j == 6)  || (i == 8 && j == 9) || (i == 3 && j == 11))
    return VECTOR3(0, -1, 0);

  if ((i == 6 && j == 10) || (i == 3 && j == 7) || (i == 4 && j == 9))
    return VECTOR3(0, 0, -1);

  return VECTOR3(0, 0, 0);
}

///////////////////////////////////////////////////////////////////////
// one entry of the rank-3 hessian of the cross product used to compute 
// the triangle normal, edge-edge case
///////////////////////////////////////////////////////////////////////
VECTOR3 crossHessianEE(const int iIn, const int jIn)
{
  int i = iIn;
  int j = jIn;

  if (i > j)
  {
    int temp = j;
    j = i;
    i = temp;
  }

  if ((i == 1 && j == 11)  || (i == 2 && j == 7) || (i == 4 && j == 8) || (i == 5 && j == 10))
    return VECTOR3(1, 0, 0);

  if ((i == 0 && j == 8) || (i == 2 && j == 9) || (i == 3 && j == 11) || (i == 5 && j == 6))
    return VECTOR3(0, 1, 0);

  if ((i == 0 && j == 10)  || (i == 1 && j == 6) || (i == 3 && j == 7) || (i == 4 && j == 9))
    return VECTOR3(0, 0, 1);

  if ((i == 1 && j == 8) || (i == 2 && j == 10) || (i == 4 && j == 11) || (i == 5 && j == 7))
    return VECTOR3(-1, 0, 0);

  if ((i == 0 && j == 11) || (i == 2 && j == 6) || (i == 3 && j == 8) || (i == 5 && j == 9))
    return VECTOR3(0, -1, 0);

  if ((i == 0 && j == 7) || (i == 1 && j == 9) || (i == 3 && j == 10) || (i == 4 && j == 6))
    return VECTOR3(0, 0, -1);

  return VECTOR3(0, 0, 0);
}

///////////////////////////////////////////////////////////////////////
// hessian of the triangle normal, vertex-face case
///////////////////////////////////////////////////////////////////////
std::vector<MATRIX12> normalHessianVF(const std::vector<VECTOR3>& e)
{
  using namespace std;

  vector<MATRIX12> H(3);
  for (int i = 0; i < 3; i++)
    H[i].setZero();

  //crossed = cross(e2, e0);
  //crossNorm = norm(crossed);
  //crossGradient = cross_gradient(x);
  VECTOR3 crossed = e[2].cross(e[0]);
  MATRIX3x12 crossGrad = crossGradientVF(e);
  const VECTOR3& z = crossed;
  
  //denom15 = (z' * z) ^ (1.5);
  REAL denom15 = pow(crossed.dot(crossed), 1.5);
  REAL denom25 = pow(crossed.dot(crossed), 2.5);

  for (int j = 0; j < 12; j++)
    for (int i = 0; i < 12; i++)
    {
      VECTOR3 zGradi = crossGrad.col(i);
      VECTOR3 zGradj = crossGrad.col(j);
      VECTOR3 zHessianij = crossHessianVF(i,j);

      // z = cross(e2, e0);
      // zGrad = crossGradientVF(:,i);
      // alpha= (z' * zGrad) / (z' * z) ^ (1.5);
      REAL a = z.dot(crossGrad.col(i)) / denom15;

      // final = (zGradj' * zGradi) / denom15 + (z' * cross_hessian(i,j)) / denom15;
      // final = final - 3 * ((z' * zGradi) / denom25) * (zGradj' * z);
      REAL aGrad = (zGradj.dot(zGradi)) / denom15 + z.dot(crossHessianVF(i,j)) / denom15;
      aGrad -= 3.0 * (z.dot(zGradi) / denom25) * zGradj.dot(z);
      
      //entry = -((zGradj' * z) / denom15) * zGradi + 1 / norm(z) * zHessianij - alpha * zGradj - alphaGradj * z;
      VECTOR3 entry = -((zGradj.dot(z)) / denom15) * zGradi + 1.0 / z.norm() * zHessianij - a * zGradj - aGrad * z;

      H[0](i,j) = entry[0];
      H[1](i,j) = entry[1];
      H[2](i,j) = entry[2];
    }
  return H;
}

///////////////////////////////////////////////////////////////////////
// hessian of the triangle normal, edge-edge case
///////////////////////////////////////////////////////////////////////
std::vector<MATRIX12> normalHessianEE(const std::vector<VECTOR3>& e)
{
  using namespace std;

  vector<MATRIX12> H(3);
  for (int i = 0; i < 3; i++)
    H[i].setZero();

  VECTOR3 crossed = e[1].cross(e[0]);
  MATRIX3x12 crossGrad = crossGradientEE(e);
  const VECTOR3& z = crossed;
  
  //denom15 = (z' * z) ^ (1.5);
  REAL denom15 = pow(crossed.dot(crossed), 1.5);
  REAL denom25 = pow(crossed.dot(crossed), 2.5);

  for (int j = 0; j < 12; j++)
    for (int i = 0; i < 12; i++)
    {
      VECTOR3 zGradi = crossGrad.col(i);
      VECTOR3 zGradj = crossGrad.col(j);
      VECTOR3 zHessianij = crossHessianEE(i,j);

      // z = cross(e2, e0);
      // zGrad = crossGradientVF(:,i);
      // alpha= (z' * zGrad) / (z' * z) ^ (1.5);
      REAL a = z.dot(crossGrad.col(i)) / denom15;

      // final = (zGradj' * zGradi) / denom15 + (z' * cross_hessian(i,j)) / denom15;
      // final = final - 3 * ((z' * zGradi) / denom25) * (zGradj' * z);
      REAL aGrad = (zGradj.dot(zGradi)) / denom15 + 
                   z.dot(crossHessianEE(i,j)) / denom15;
      aGrad -= 3.0 * (z.dot(zGradi) / denom25) * zGradj.dot(z);
      
      //entry = -((zGradj' * z) / denom15) * zGradi + 
      //          1 / norm(z) * zHessianij - 
      //          alpha * zGradj - alphaGradj * z;
      VECTOR3 entry = -((zGradj.dot(z)) / denom15) * zGradi + 
                        1.0 / z.norm() * zHessianij - 
                        a * zGradj - aGrad * z;

      H[0](i,j) = entry[0];
      H[1](i,j) = entry[1];
      H[2](i,j) = entry[2];
    }
  return H;
}

///////////////////////////////////////////////////////////////////////
// get the barycentric coordinate of the projection of v[0] onto the triangle
// formed by v[1], v[2], v[3]
///////////////////////////////////////////////////////////////////////
VECTOR3 getBarycentricCoordinates(const std::vector<VECTOR3>& vertices)
{
  const VECTOR3 v0 = vertices[1];
  const VECTOR3 v1 = vertices[2];
  const VECTOR3 v2 = vertices[3];
    
  const VECTOR3 e1 = v1 - v0;
  const VECTOR3 e2 = v2 - v0;
  const VECTOR3 n = e1.cross(e2);
  const VECTOR3 nHat = n / n.norm();
  const VECTOR3 v = vertices[0] - (nHat.dot(vertices[0] - v0)) * nHat;

  // get the barycentric coordinates
  const VECTOR3 na = (v2 - v1).cross(v - v1);
  const VECTOR3 nb = (v0 - v2).cross(v - v2);
  const VECTOR3 nc = (v1 - v0).cross(v - v0);
  const VECTOR3 barycentric(n.dot(na) / n.squaredNorm(),
                            n.dot(nb) / n.squaredNorm(),
                            n.dot(nc) / n.squaredNorm());

  return barycentric;
}

///////////////////////////////////////////////////////////////////////
// get the barycentric coordinate of the projection of v[0] onto the triangle
// formed by v[1], v[2], v[3]
///////////////////////////////////////////////////////////////////////
VECTOR3 getBarycentricCoordinates(const VECTOR12& vertices)
{
  std::vector<VECTOR3> vs(4);
  for (int x = 0; x < 4; x++)
  {
    vs[x][0] = vertices[3 * x];
    vs[x][1] = vertices[3 * x + 1];
    vs[x][2] = vertices[3 * x + 2];
  }
  return getBarycentricCoordinates(vs);
}

///////////////////////////////////////////////////////////////////////
// find the distance from a line segment (v1, v2) to a point (v0)
///////////////////////////////////////////////////////////////////////
REAL pointLineDistance(const VECTOR3 v0, const VECTOR3& v1, const VECTOR3& v2)
{
  const VECTOR3 e0 = v0 - v1;
  const VECTOR3 e1 = v2 - v1;
  const VECTOR3 e1hat = e1 / e1.norm();
  const REAL projection = e0.dot(e1hat);

  // if it projects onto the line segment, use that length
  if (projection > 0.0 && projection < e1.norm())
  {
    const VECTOR3 normal = e0 - projection * e1hat;
    return normal.norm();
  }

  // if it doesn't, find the point-point distances
  const REAL diff01 = (v0 - v1).norm();
  const REAL diff02 = (v0 - v2).norm();

  return (diff01 < diff02) ? diff01 : diff02;
}

///////////////////////////////////////////////////////////////////////
// find the distance from a line segment (v1, v2) to a point (v0)
///////////////////////////////////////////////////////////////////////
REAL pointLineDistanceDebug(const VECTOR3 v0, const VECTOR3& v1, const VECTOR3& v2)
{
  using namespace std;

  const VECTOR3 e0 = v0 - v1;
  const VECTOR3 e1 = v2 - v1;
  const VECTOR3 e1hat = e1 / e1.norm();
  const REAL projection = e0.dot(e1hat);
  //cout << " projection: " << projection << endl;

  // if it projects onto the line segment, use that length
  if (projection > 0.0 && projection < e1.norm())
  {
    const VECTOR3 normal = e0 - projection * e1hat;
    return normal.norm();
  }

  // if it doesn't, find the point-point distances
  const REAL diff01 = (v0 - v1).norm();
  const REAL diff02 = (v0 - v2).norm();

  return (diff01 < diff02) ? diff01 : diff02;
}

///////////////////////////////////////////////////////////////////////
// get the linear interpolation coordinates from v0 to the line segment
// between v1 and v2
///////////////////////////////////////////////////////////////////////
VECTOR2 getLerp(const VECTOR3 v0, const VECTOR3& v1, const VECTOR3& v2)
{
  const VECTOR3 e0 = v0 - v1;
  const VECTOR3 e1 = v2 - v1;
  const VECTOR3 e1hat = e1 / e1.norm();
  const REAL projection = e0.dot(e1hat);

  if (projection < 0.0)
    return VECTOR2(1.0, 0.0);

  if (projection >= e1.norm())
    return VECTOR2(0.0, 1.0);

  const REAL ratio = projection / e1.norm();
  return VECTOR2(1.0 - ratio, ratio);
}

///////////////////////////////////////////////////////////////////////
// get the linear interpolation coordinates from v0 to the line segment
// between v1 and v2
///////////////////////////////////////////////////////////////////////
VECTOR2 getLerpDebug(const VECTOR3 v0, const VECTOR3& v1, const VECTOR3& v2)
{
  const VECTOR3 e0 = v0 - v1;
  const VECTOR3 e1 = v2 - v1;
  const VECTOR3 e1hat = e1 / e1.norm();
  const REAL projection = e0.dot(e1hat);

  if (projection < 0.0)
  {
    return VECTOR2(1.0, 0.0);
  }

  if (projection >= e1.norm())
  {
    return VECTOR2(0.0, 1.0);
  }

  const REAL ratio = projection / e1.norm();
  return VECTOR2(1.0 - ratio, ratio);
}

///////////////////////////////////////////////////////////////////////
// get the barycentric coordinate of the projection of v[0] onto the triangle
// formed by v[1], v[2], v[3]
//
// but, if the projection is actually outside, project to all of the
// edges and find the closest point that's still inside the triangle
///////////////////////////////////////////////////////////////////////
VECTOR3 getInsideBarycentricCoordinates(const std::vector<VECTOR3>& vertices)
{
  VECTOR3 barycentric = getBarycentricCoordinates(vertices);

  // if it's already inside, we're all done
  if (barycentric[0] >= 0.0 &&
      barycentric[1] >= 0.0 &&
      barycentric[2] >= 0.0)
    return barycentric;

  // find distance to all the line segments
  //
  // there's lots of redundant computation between here and getLerp,
  // but let's get it working and see if it fixes the actual
  // artifact before optimizing
  REAL distance12 = pointLineDistance(vertices[0], vertices[1], vertices[2]);
  REAL distance23 = pointLineDistance(vertices[0], vertices[2], vertices[3]);
  REAL distance31 = pointLineDistance(vertices[0], vertices[3], vertices[1]);

  // less than or equal is important here, otherwise fallthrough breaks
  if (distance12 <= distance23 && distance12 <= distance31)
  {
    VECTOR2 lerp = getLerp(vertices[0], vertices[1], vertices[2]);
    barycentric[0] = lerp[0];
    barycentric[1] = lerp[1];
    barycentric[2] = 0.0;
    return barycentric;
  }
  
  // less than or equal is important here, otherwise fallthrough breaks
  if (distance23 <= distance12 && distance23 <= distance31)
  {
    VECTOR2 lerp = getLerp(vertices[0], vertices[2], vertices[3]);
    barycentric[0] = 0.0;
    barycentric[1] = lerp[0];
    barycentric[2] = lerp[1];
    return barycentric;
  }

  // else it must be the 31 case
  VECTOR2 lerp = getLerp(vertices[0], vertices[3], vertices[1]);
  barycentric[0] = lerp[1];
  barycentric[1] = 0.0;
  barycentric[2] = lerp[0];
  return barycentric;
}

///////////////////////////////////////////////////////////////////////
// get the barycentric coordinate of the projection of v[0] onto the triangle
// formed by v[1], v[2], v[3]
//
// but, if the projection is actually outside, project to all of the
// edges and find the closest point that's still inside the triangle
///////////////////////////////////////////////////////////////////////
VECTOR3 getInsideBarycentricCoordinatesDebug(const std::vector<VECTOR3>& vertices)
{
  using namespace std;

  VECTOR3 barycentric = getBarycentricCoordinates(vertices);
  //cout << " bary: " << barycentric.transpose() << endl;

  // if it's already inside, we're all done
  if (barycentric[0] >= 0.0 &&
      barycentric[1] >= 0.0 &&
      barycentric[2] >= 0.0)
    return barycentric;

  // find distance to all the line segments
  //
  // there's lots of redundant computation between here and getLerp,
  // but let's get it working and see if it fixes the actual
  // artifact before optimizing
  REAL distance12 = pointLineDistanceDebug(vertices[0], vertices[1], vertices[2]);
  REAL distance23 = pointLineDistanceDebug(vertices[0], vertices[2], vertices[3]);
  REAL distance31 = pointLineDistanceDebug(vertices[0], vertices[3], vertices[1]);

  //cout << " distances: " << distance12 << " " << distance23 << " " << distance31 << endl;
  //cout << " 02 distance: " << (vertices[0] - vertices[2]).norm() << endl;

  if (distance12 <= distance23 && distance12 <= distance31)
  {
    VECTOR2 lerp = getLerp(vertices[0], vertices[1], vertices[2]);
    barycentric[0] = lerp[0];
    barycentric[1] = lerp[1];
    barycentric[2] = 0.0;
    return barycentric;
  }
  
  if (distance23 <= distance12 && distance23 <= distance31)
  {
    VECTOR2 lerp = getLerp(vertices[0], vertices[2], vertices[3]);
    barycentric[0] = 0.0;
    barycentric[1] = lerp[0];
    barycentric[2] = lerp[1];
    return barycentric;
  }

  // else it must be the 31 case
  VECTOR2 lerp = getLerpDebug(vertices[0], vertices[3], vertices[1]);
  barycentric[0] = lerp[1];
  barycentric[1] = 0.0;
  barycentric[2] = lerp[0];
  return barycentric;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 flattenVertices(const std::vector<VECTOR3>& v)
{
  VECTOR12 x;
  int entry = 0;
  for (int j = 0; j < 4; j++)
    for (int i = 0; i < 3; i++, entry++)
      x[entry] = v[j][i];
  return x;
}

///////////////////////////////////////////////////////////////////////
// does this face and edge intersect?
///////////////////////////////////////////////////////////////////////
bool faceEdgeIntersection(const std::vector<VECTOR3>& triangleVertices, 
                          const std::vector<VECTOR3>& edgeVertices)
{
  assert(triangleVertices.size() == 3);
  assert(edgeVertices.size() == 2);

  const VECTOR3& a = triangleVertices[0];
  const VECTOR3& b = triangleVertices[1];
  const VECTOR3& c = triangleVertices[2];

  const VECTOR3& origin = edgeVertices[0];
  const VECTOR3& edgeDiff = (edgeVertices[1] - edgeVertices[0]);
  const VECTOR3& direction = edgeDiff.normalized();

  const VECTOR3 geometricNormal = ((b - a).cross(c - a)).normalized();

  const VECTOR3 diff = a - origin;
  REAL denom = direction.dot(geometricNormal);
  if (fabs(denom) <= 0.0) return false;

  REAL t = diff.dot(geometricNormal) / denom;
  if (t < 0) return false;

  VECTOR3 h = origin + direction * t;

  VECTOR3 test = (b - a).cross(h - a);
  if (geometricNormal.dot(test) < 0) return false; 
  test = (c - b).cross(h - b);
  if (geometricNormal.dot(test) < 0) return false; 
  test = (a - c).cross(h - c);
  if (geometricNormal.dot(test) < 0) return false; 

  if (t < edgeDiff.norm())
    return true;

  return false;
}


} // HOBAK
