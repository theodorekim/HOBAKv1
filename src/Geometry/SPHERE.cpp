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
#include "SPHERE.h"
#include <iostream>

using namespace std;

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
// positions are defined using R * S * x + t
///////////////////////////////////////////////////////////////////////
SPHERE::SPHERE(const VECTOR3& center, const REAL& scale)
{
  _scale = MATRIX3::Identity() * scale; 
  _rotation = MATRIX3::Identity();
  _translation = center;
  _scaleInverse = _scale.inverse();

  _name = string("SPHERE");
}

SPHERE::~SPHERE()
{
}

///////////////////////////////////////////////////////////////////////
// is a point inside the sphere?
///////////////////////////////////////////////////////////////////////
bool SPHERE::inside(const VECTOR3& point) const 
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);
  REAL radius = transformed.norm();

  if (radius < 1.0)
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// distance to the sphere
///////////////////////////////////////////////////////////////////////
REAL SPHERE::distance(const VECTOR3& point) const 
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);
  REAL radius = transformed.norm();

  return fabs(radius - 1.0) * _scale(0,0);
}

///////////////////////////////////////////////////////////////////////
// signed distance to the sphere
// remember that "inside" is negative with signed distance
///////////////////////////////////////////////////////////////////////
REAL SPHERE::signedDistance(const VECTOR3& point) const
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);
  REAL radius = transformed.norm();

  return (radius - 1.0) * _scale(0,0);
}

//////////////////////////////////////////////////////////////////////
// get the closest point on the object, as well as the normal at 
// the point
//////////////////////////////////////////////////////////////////////
void SPHERE::getClosestPoint(const VECTOR3& query, 
                             VECTOR3& closestPointLocal, 
                             VECTOR3& normalLocal) const
{
  const VECTOR3 collisionPoint = worldVertexToLocal(query);
  closestPointLocal = collisionPoint.normalized();

  // this is the one instance where both of these are the same
  normalLocal = closestPointLocal;
}

} // HOBAK
