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
//////////////////////////////////////////////////////////////////////////////
// This is a modified version of the Wild Magic DistVector3Triangle3 class.
// This is what is used to generate the distance grid.
//
// The license info is as below:
//
// Wild Magic Source Code
// David Eberly
// http://www.geometrictools.com
// Copyright (c) 1998-2008
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.  The license is available for reading at
// either of the locations:
//     http://www.gnu.org/copyleft/lgpl.html
//     http://www.geometrictools.com/License/WildMagicLicense.pdf
//
// Version: 4.0.1 (2007/05/06)
//----------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
static REAL pointFaceDistanceSq(const VECTOR3& v0, const VECTOR3& v1, const VECTOR3& v2,
                                const VECTOR3& point)
{
  VECTOR3 kDiff = v0 - point;
  VECTOR3 kEdge0 = v1 - v0;
  VECTOR3 kEdge1 = v2 - v0;
  REAL fA00 = kEdge0.squaredNorm();
  REAL fA01 = kEdge0.dot(kEdge1);
  REAL fA11 = kEdge1.squaredNorm();
  REAL fB0 = kDiff.dot(kEdge0);
  REAL fB1 = kDiff.dot(kEdge1);
  REAL fC = kDiff.squaredNorm();
  REAL fDet = fabs(fA00*fA11-fA01*fA01);
  REAL fS = fA01*fB1-fA11*fB0;
  REAL fT = fA01*fB0-fA00*fB1;
  REAL fSqrDistance;

  if (fS + fT <= fDet) {
      if (fS < (REAL)0.0) {
          if (fT < (REAL)0.0)  // region 4
          {
              if (fB0 < (REAL)0.0) {
                  fT = (REAL)0.0;
                  if (-fB0 >= fA00) {
                      fS = (REAL)1.0;
                      fSqrDistance = fA00+((REAL)2.0)*fB0+fC;
                  }
                  else {
                      fS = -fB0/fA00;
                      fSqrDistance = fB0*fS+fC;
                  }
              }
              else {
                  fS = (REAL)0.0;
                  if (fB1 >= (REAL)0.0) {
                      fT = (REAL)0.0;
                      fSqrDistance = fC;
                  }
                  else if (-fB1 >= fA11) {
                      fT = (REAL)1.0;
                      fSqrDistance = fA11+((REAL)2.0)*fB1+fC;
                  }
                  else {
                      fT = -fB1/fA11;
                      fSqrDistance = fB1*fT+fC;
                  }
              }
          }
          else  // region 3
          {
              fS = (REAL)0.0;
              if (fB1 >= (REAL)0.0) {
                  fT = (REAL)0.0;
                  fSqrDistance = fC;
              }
              else if (-fB1 >= fA11) {
                  fT = (REAL)1.0;
                  fSqrDistance = fA11+((REAL)2.0)*fB1+fC;
              }
              else {
                  fT = -fB1/fA11;
                  fSqrDistance = fB1*fT+fC;
              }
          }
      }
      else if (fT < (REAL)0.0)  // region 5
      {
          fT = (REAL)0.0;
          if (fB0 >= (REAL)0.0) {
              fS = (REAL)0.0;
              fSqrDistance = fC;
          }
          else if (-fB0 >= fA00) {
              fS = (REAL)1.0;
              fSqrDistance = fA00+((REAL)2.0)*fB0+fC;
          }
          else {
              fS = -fB0/fA00;
              fSqrDistance = fB0*fS+fC;
          }
      }
      else  // region 0
      {
          // minimum at interior point
          REAL fInvDet = ((REAL)1.0)/fDet;
          fS *= fInvDet;
          fT *= fInvDet;
          fSqrDistance = fS*(fA00*fS+fA01*fT+((REAL)2.0)*fB0) +
              fT*(fA01*fS+fA11*fT+((REAL)2.0)*fB1)+fC;
      }
  }
  else {
      REAL fTmp0, fTmp1, fNumer, fDenom;

      if (fS < (REAL)0.0)  // region 2
      {
          fTmp0 = fA01 + fB0;
          fTmp1 = fA11 + fB1;
          if (fTmp1 > fTmp0) {
              fNumer = fTmp1 - fTmp0;
              fDenom = fA00-2.0f*fA01+fA11;
              if (fNumer >= fDenom) {
                  fS = (REAL)1.0;
                  fT = (REAL)0.0;
                  fSqrDistance = fA00+((REAL)2.0)*fB0+fC;
              }
              else {
                  fS = fNumer/fDenom;
                  fT = (REAL)1.0 - fS;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                      fT*(fA01*fS+fA11*fT+((REAL)2.0)*fB1)+fC;
              }
          }
          else {
              fS = (REAL)0.0;
              if (fTmp1 <= (REAL)0.0) {
                  fT = (REAL)1.0;
                  fSqrDistance = fA11+((REAL)2.0)*fB1+fC;
              }
              else if (fB1 >= (REAL)0.0) {
                  fT = (REAL)0.0;
                  fSqrDistance = fC;
              }
              else {
                  fT = -fB1/fA11;
                  fSqrDistance = fB1*fT+fC;
              }
          }
      }
      else if (fT < (REAL)0.0)  // region 6
      {
          fTmp0 = fA01 + fB1;
          fTmp1 = fA00 + fB0;
          if (fTmp1 > fTmp0) {
              fNumer = fTmp1 - fTmp0;
              fDenom = fA00-((REAL)2.0)*fA01+fA11;
              if (fNumer >= fDenom) {
                  fT = (REAL)1.0;
                  fS = (REAL)0.0;
                  fSqrDistance = fA11+((REAL)2.0)*fB1+fC;
              }
              else {
                  fT = fNumer/fDenom;
                  fS = (REAL)1.0 - fT;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+((REAL)2.0)*fB0) +
                      fT*(fA01*fS+fA11*fT+((REAL)2.0)*fB1)+fC;
              }
          }
          else {
              fT = (REAL)0.0;
              if (fTmp1 <= (REAL)0.0) {
                  fS = (REAL)1.0;
                  fSqrDistance = fA00+((REAL)2.0)*fB0+fC;
              }
              else if (fB0 >= (REAL)0.0) {
                  fS = (REAL)0.0;
                  fSqrDistance = fC;
              }
              else {
                  fS = -fB0/fA00;
                  fSqrDistance = fB0*fS+fC;
              }
          }
      }
      else  // region 1
      {
          fNumer = fA11 + fB1 - fA01 - fB0;
          if (fNumer <= (REAL)0.0) {
              fS = (REAL)0.0;
              fT = (REAL)1.0;
              fSqrDistance = fA11+((REAL)2.0)*fB1+fC;
          }
          else {
              fDenom = fA00-2.0f*fA01+fA11;
              if (fNumer >= fDenom) {
                  fS = (REAL)1.0;
                  fT = (REAL)0.0;
                  fSqrDistance = fA00+((REAL)2.0)*fB0+fC;
              }
              else {
                  fS = fNumer/fDenom;
                  fT = (REAL)1.0 - fS;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+((REAL)2.0)*fB0) +
                      fT*(fA01*fS+fA11*fT+((REAL)2.0)*fB1)+fC;
              }
          }
      }
  }

  // account for numerical round-off error
  if (fSqrDistance < (REAL)0.0) {
      fSqrDistance = (REAL)0.0;
  }
  return fSqrDistance;
}

//////////////////////////////////////////////////////////////////////////////
// see if vertex-triangle distance test is correct
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Collision detection tests", "[Collision detection]" )
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING TET_MESH collision detection is correct" << endl;
  cout << "=============================================================== " << endl;

  VECTOR3 v0(0,0,0.1);
  VECTOR3 v1(1,0,0.1);
  VECTOR3 v2(0,1,0.1);

  vector<VECTOR3> vs;
  vs.push_back(VECTOR3(0,0,0.123));
  vs.push_back(VECTOR3(0,0,-0.123));  // v0 test

  vs.push_back(VECTOR3(1,0,0.123));
  vs.push_back(VECTOR3(2,0,0.123));
  vs.push_back(VECTOR3(1,0,-0.123));
  vs.push_back(VECTOR3(2,0,-0.123)); // v1 test

  vs.push_back(VECTOR3(0,1,0.123));
  vs.push_back(VECTOR3(0,2,0.123));
  vs.push_back(VECTOR3(0,1,-0.123));
  vs.push_back(VECTOR3(0,2,-0.123)); // v2 test
  
  vs.push_back(VECTOR3(1,1,0.123));
  vs.push_back(VECTOR3(1,1,-0.123));  // edge 1 test
  
  vs.push_back(VECTOR3(0.5,0,0.123));
  vs.push_back(VECTOR3(0.5,0,-0.123));
  vs.push_back(VECTOR3(0.5,-0.1,0.123));
  vs.push_back(VECTOR3(0.5,-0.1,-0.123)); // edge 0 test
  
  vs.push_back(VECTOR3(0,0.5,0.123));
  vs.push_back(VECTOR3(0,0.5,-0.123));
  vs.push_back(VECTOR3(-0.1,0.5,0.123));
  vs.push_back(VECTOR3(-0.1,0.5,-0.123)); // edge 2 test

  // try some structured tests
  bool allPassed = true; 
  for (unsigned int i = 0; i < vs.size(); i++)
  { 
    const REAL distanceSq = pointFaceDistanceSq(v0, v1, v2, vs[i]);
    const REAL trial = TET_MESH::pointTriangleDistance(v0,v1, v2,vs[i]);
    REAL diff = fabs(trial * trial - distanceSq);

    if (diff > 1e-8)
    {
      cout << " diff; " << diff << endl;
      cout << " ground: " << sqrt(distanceSq) << endl;
      cout << " guess:  " << trial << endl;
      allPassed = false;
      break;
    }
  }
  cout << " STRUCTURED TESTS PASSED " << endl;
  REQUIRE(allPassed);

  // try some randomized tests
  for (unsigned int i = 0; i < 10000; i++)
  {
    v0 = randomVector3(10);
    v1 = randomVector3(10);
    v2 = randomVector3(10);
    VECTOR3 v = randomVector3(10);
    
    const REAL distanceSq = pointFaceDistanceSq(v0, v1, v2, v);
    const REAL trial = TET_MESH::pointTriangleDistance(v0,v1, v2,v);
    REAL diff = fabs(trial * trial - distanceSq);

    if (diff > 1e-8)
    {
      cout << " diff; " << diff << endl;
      cout << " ground: " << sqrt(distanceSq) << endl;
      cout << " guess:  " << trial << endl;
      allPassed = false;
      break;
    }
  }
  cout << " RANDOM TESTS PASSED " << endl;
  REQUIRE(allPassed);

  cout << " Point-triangle tests PASSED" << endl;
}

//////////////////////////////////////////////////////////////////////////////
// edge-edge unit tests
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Edge-edge tests", "[Edge-edge testing]" )
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING EDGE-EDGE TESTS" << endl;
  cout << "=============================================================== " << endl;
  //const REAL h = 0.25;
  const VECTOR3 b0(-1,0,0);
  const VECTOR3 b1( 1,0,0);
  REAL h = 1.0;
  VECTOR3 a0(0,h,-1);
  VECTOR3 a1(0,h, 1);

  // test out the intersection points
  VECTOR3 aPoint, bPoint;
  IntersectLineSegments(a0, a1, b0, b1, aPoint, bPoint);

  // we got this back, right?
  const VECTOR3 aPointGround(0,h,0);
  const VECTOR3 bPointGround(0,0,0);
  cout << " a intersect: " << endl << aPoint << endl;
  cout << " b intersect: " << endl << bPoint << endl;
  REQUIRE((aPoint - aPointGround).norm() < 1e-8);
  REQUIRE((bPoint - bPointGround).norm() < 1e-8);

  // test out the interpolation coordinates
  const VECTOR3 e0 = a1 - a0;
  const VECTOR3 e1 = b1 - b0;

  VECTOR2 a,b;
  a[1] = (aPoint - a0).norm() / e0.norm();
  a[0] = 1.0 - a[1];
  b[1] = (bPoint - b0).norm() / e1.norm();
  b[0] = 1.0 - b[1];

  cout << " a coordinates: " << endl << a << endl;
  cout << " b coordinates: " << endl << b << endl;

  const VECTOR2 half(0.5, 0.5);
  REQUIRE((a - half).norm() < 1e-8);
  REQUIRE((b - half).norm() < 1e-8);

  cout << "=============================================================== " << endl;
  cout << " ZERO LENGTH COLLISON GRADIENT TEST" << endl;
  cout << "=============================================================== " << endl;
  // test out the gradient
  vector<VECTOR3> vertices;

  a0 = VECTOR3(0,-h,-1);
  a1 = VECTOR3(0,-h,1);
  vertices.push_back(a0);
  vertices.push_back(a1);
  vertices.push_back(b0);
  vertices.push_back(b1);

  using namespace HOBAK::VOLUME;
  {
    EDGE_COLLISION material(1.0, 1.0);
    VECTOR12 gradient = material.gradient(vertices,a,b);
    cout << " edge gradient: " << endl << gradient << endl;

    // right at the original spring length, so gradient should be nothing
    REQUIRE(gradient.norm() <= 1e-8);

    // make sure that if it's less than epsilon, gradients point the right way 
    cout << "=============================================================== " << endl;
    cout << " LESS THAN EPSILON COLLISON GRADIENT TEST" << endl;
    cout << "=============================================================== " << endl;
    h = 0.5;
    a0 = VECTOR3(0,-h,-1);
    a1 = VECTOR3(0,-h, 1);
    vertices[0] = a0;
    vertices[1] = a1;
    gradient = material.gradient(vertices, a,b);
    cout << " edge gradient: " << endl << gradient << endl;

    VECTOR3 gradient0(gradient[0], gradient[1], gradient[2]);
    VECTOR3 gradient1(gradient[3], gradient[4], gradient[5]);
    VECTOR3 gradient2(gradient[6], gradient[7], gradient[8]);
    VECTOR3 gradient3(gradient[9], gradient[10], gradient[11]);

    VECTOR3 halfGradient(0, -0.5, 0);
    // this flips as the gradient convention flips
    REQUIRE((gradient0 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient1 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient2 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient3 - halfGradient).norm() < 1e-8);
    /*
    REQUIRE((gradient0 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient1 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient2 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient3 + halfGradient).norm() < 1e-8);
    */

    // make sure that if it's greater than epsilon, gradients point the right way 
    cout << "=============================================================== " << endl;
    cout << " GREATER THAN EPSILON COLLISON GRADIENT TEST" << endl;
    cout << "=============================================================== " << endl;
    h = 1.5;
    a0 = VECTOR3(0,-h,-1);
    a1 = VECTOR3(0,-h, 1);
    vertices[0] = a0;
    vertices[1] = a1;
    gradient = material.gradient(vertices, a,b);
    cout << " edge gradient: " << endl << gradient << endl;

    gradient0 = VECTOR3(gradient[0], gradient[1], gradient[2]);
    gradient1 = VECTOR3(gradient[3], gradient[4], gradient[5]);
    gradient2 = VECTOR3(gradient[6], gradient[7], gradient[8]);
    gradient3 = VECTOR3(gradient[9], gradient[10], gradient[11]);

    // this flips as the gradient convention flips
    REQUIRE((gradient0 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient1 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient2 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient3 + halfGradient).norm() < 1e-8);
    /*
    REQUIRE((gradient0 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient1 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient2 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient3 - halfGradient).norm() < 1e-8);
    */
  }
  {
    h = 1.0;
    a0 = VECTOR3(0,-h,-1);
    a1 = VECTOR3(0,-h,1);
    vertices[0] = a0;
    vertices[1] = a1;
    EDGE_SQRT_COLLISION material(1.0, 1.0);
    VECTOR12 gradient = material.EDGE_COLLISION::gradient(vertices,a,b);
    cout << " edge gradient: " << endl << gradient << endl;

    // right at the original spring length, so gradient should be nothing
    REQUIRE(gradient.norm() <= 1e-8);

    // make sure that if it's less than epsilon, gradients point the right way 
    cout << "=============================================================== " << endl;
    cout << " LESS THAN EPSILON COLLISON GRADIENT TEST" << endl;
    cout << "=============================================================== " << endl;
    h = 0.5;
    a0 = VECTOR3(0,-h,-1);
    a1 = VECTOR3(0,-h, 1);
    vertices[0] = a0;
    vertices[1] = a1;
    gradient = material.EDGE_COLLISION::gradient(vertices, a,b);
    cout << " edge gradient: " << endl << gradient << endl;

    VECTOR3 gradient0(gradient[0], gradient[1], gradient[2]);
    VECTOR3 gradient1(gradient[3], gradient[4], gradient[5]);
    VECTOR3 gradient2(gradient[6], gradient[7], gradient[8]);
    VECTOR3 gradient3(gradient[9], gradient[10], gradient[11]);

    VECTOR3 halfGradient(0, -0.5, 0);
    // this flips as the gradient convention flips
    REQUIRE((gradient0 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient1 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient2 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient3 - halfGradient).norm() < 1e-8);
    /*
    REQUIRE((gradient0 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient1 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient2 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient3 + halfGradient).norm() < 1e-8);
    */

    // make sure that if it's greater than epsilon, gradients point the right way 
    cout << "=============================================================== " << endl;
    cout << " GREATER THAN EPSILON COLLISON GRADIENT TEST" << endl;
    cout << "=============================================================== " << endl;
    h = 1.5;
    a0 = VECTOR3(0,-h,-1);
    a1 = VECTOR3(0,-h, 1);
    vertices[0] = a0;
    vertices[1] = a1;
    gradient = material.EDGE_COLLISION::gradient(vertices, a,b);
    cout << " edge gradient: " << endl << gradient << endl;

    gradient0 = VECTOR3(gradient[0], gradient[1], gradient[2]);
    gradient1 = VECTOR3(gradient[3], gradient[4], gradient[5]);
    gradient2 = VECTOR3(gradient[6], gradient[7], gradient[8]);
    gradient3 = VECTOR3(gradient[9], gradient[10], gradient[11]);

    // this flips as the gradient convention flips
    REQUIRE((gradient0 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient1 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient2 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient3 + halfGradient).norm() < 1e-8);
    /*
    REQUIRE((gradient0 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient1 + halfGradient).norm() < 1e-8);
    REQUIRE((gradient2 - halfGradient).norm() < 1e-8);
    REQUIRE((gradient3 - halfGradient).norm() < 1e-8);
    */
  }

  /*
  cout << "=============================================================== " << endl;
  cout << " RUNTIME UNIT TEST " << endl;
  cout << "=============================================================== " << endl;
  vertices[0] = VECTOR3(0, 1, -1);
  vertices[1] = VECTOR3(0, 1, 1);
  vertices[2] = VECTOR3(-1, 1.125, 0);
  vertices[3] = VECTOR3( 1, 1.125, 0);
  material = EDGE_COLLISION(1.0, 0.05);
  gradient = material.gradient(vertices, a,b);
  cout << " a: " << endl << a << endl;
  cout << " b: " << endl << b << endl;
  cout << " edge gradient: " << endl << gradient << endl;
  */

  cout << "=============================================================== " << endl;
  cout << " FACE-EDGE INTERSECTION TEST" << endl;
  cout << "=============================================================== " << endl;

  // pierces right through the center
  vector<VECTOR3> triangleVertices(3);
  triangleVertices[0] = VECTOR3(0,0,0);
  triangleVertices[1] = VECTOR3(1,0,0);
  triangleVertices[2] = VECTOR3(0,1,0);
  vector<VECTOR3> edgeVertices(2);
  edgeVertices[0] = VECTOR3(0.25,0.25,-1);
  edgeVertices[1] = VECTOR3(0.25,0.25,1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == true);
  edgeVertices[0] = VECTOR3(0.25,0.25,1);
  edgeVertices[1] = VECTOR3(0.25,0.25,-1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == true);
  cout << " CORRECT HIT" << endl;

  // starts above the triangle, should miss
  edgeVertices[0] = VECTOR3(0.25,0.25,2);
  edgeVertices[1] = VECTOR3(0.25,0.25,1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == false);
  edgeVertices[0] = VECTOR3(0.25,0.25,1);
  edgeVertices[1] = VECTOR3(0.25,0.25,2);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == false);
  cout << " CORRECT MISS" << endl;
 
  // hits right along the edge 
  edgeVertices[0] = VECTOR3(0.5,0.5,-1);
  edgeVertices[1] = VECTOR3(0.5,0.5,1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == true);
  edgeVertices[0] = VECTOR3(0.5,0.5,1);
  edgeVertices[1] = VECTOR3(0.5,0.5,-1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == true);
  cout << " CORRECT HIT" << endl;
  
  // starts and ends below the triangle, should miss
  edgeVertices[0] = VECTOR3(0.25,0.25,-2);
  edgeVertices[1] = VECTOR3(0.25,0.25,-1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == false);
  edgeVertices[0] = VECTOR3(0.25,0.25,-1);
  edgeVertices[1] = VECTOR3(0.25,0.25,-2);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == false);
  cout << " CORRECT MISS" << endl;

  // hits right on the origin
  edgeVertices[0] = VECTOR3(0,0,-1);
  edgeVertices[1] = VECTOR3(0,0,1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == true);
  edgeVertices[0] = VECTOR3(0,0,1);
  edgeVertices[1] = VECTOR3(0,0,-1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == true);
  cout << " CORRECT HIT" << endl;

  // misses right before the origin
  edgeVertices[0] = VECTOR3(-0.0001,-0.0001,-1);
  edgeVertices[1] = VECTOR3(-0.0001,-0.0001,1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == false);
  edgeVertices[0] = VECTOR3(-0.0001,-0.0001,1);
  edgeVertices[1] = VECTOR3(-0.0001,-0.0001,-1);
  assert(faceEdgeIntersection(triangleVertices, edgeVertices) == false);
  cout << " CORRECT MISS" << endl;
}
