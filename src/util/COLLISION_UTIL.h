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
#ifndef COLLISION_UTIL_H
#define COLLISION_UTIL_H

#include "SETTINGS.h"

namespace HOBAK {
  // gradient of the triangle normal, vertex-face case
  MATRIX3x12 normalGradientVF(const std::vector<VECTOR3>& e);
  
  // gradient of a normal, edge-edge case
  MATRIX3x12 normalGradientEE(const std::vector<VECTOR3>& e);
 
  // hessian of the triangle normal, vertex-face case
  std::vector<MATRIX12> normalHessianVF(const std::vector<VECTOR3>& e);
  
  // hessian of the normal, edge-edge case
  std::vector<MATRIX12> normalHessianEE(const std::vector<VECTOR3>& e);

  // gradient of the cross product used to compute the
  // triangle normal, vertex-face case
  MATRIX3x12 crossGradientVF(const std::vector<VECTOR3>& e);
 
  // gradient of the cross product used to compute the normal, edge-edge case
  MATRIX3x12 crossGradientEE(const std::vector<VECTOR3>& e);

  // one entry of the rank-3 hessian of the cross product used to compute the
  // triangle normal, vertex-face case
  VECTOR3 crossHessianVF(const int i, const int j);
  
  // one entry of the rank-3 hessian of the cross product used to compute the
  // normal, edge-edge case
  VECTOR3 crossHessianEE(const int i, const int j);

  // get the barycentric coordinate of the projection of v[0] onto the triangle
  // formed by v[1], v[2], v[3]
  VECTOR3 getBarycentricCoordinates(const std::vector<VECTOR3>& v);
  VECTOR3 getBarycentricCoordinates(const VECTOR12& vertices);

  // get the barycentric coordinate of the projection of v[0] onto the triangle
  // formed by v[1], v[2], v[3]
  //
  // but, if the projection is actually outside, project to all of the
  // edges and find the closest point that's still inside the triangle
  VECTOR3 getInsideBarycentricCoordinates(const std::vector<VECTOR3>& vertices);
  VECTOR3 getInsideBarycentricCoordinatesDebug(const std::vector<VECTOR3>& vertices);

  // flatten vertices into a vector
  VECTOR12 flattenVertices(const std::vector<VECTOR3>& v);

  // does this face and edge intersect?
  bool faceEdgeIntersection(const std::vector<VECTOR3>& triangleVertices, 
                            const std::vector<VECTOR3>& edgeVertices);
}

#endif
