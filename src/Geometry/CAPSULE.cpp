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
#include "CAPSULE.h"
#include <iostream>

using namespace std;

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
// box positions are defined using R * S * x + t
///////////////////////////////////////////////////////////////////////
CAPSULE::CAPSULE(const VECTOR3& center, const REAL& radius, const REAL& height) :
  CYLINDER(center, radius, height)
{
  _scale = MATRIX3::Identity();
  _scale(0,0) = _scale(2,2) = radius;
  _scale(1,1) = height;
  _rotation = MATRIX3::Identity();
  _translation = center;
  _scaleInverse = _scale.inverse();

  _name = string("CAPSULE");
}

CAPSULE::~CAPSULE()
{
}

///////////////////////////////////////////////////////////////////////
// is a point inside the box?
///////////////////////////////////////////////////////////////////////
bool CAPSULE::inside(const VECTOR3& point) const 
{
  // transform back to local coordinates
  VECTOR3 local = worldVertexToLocal(point);

  // check the top end cap 
  if (local[1] > 0.5 * _height)
  {
    const VECTOR3 topCap(0.0, 0.5 * _height, 0.0);
    const VECTOR3 diff = topCap - local;

    if (diff.norm() > _radius) return false;
    return true;
  }

  // check the bottom end cap 
  if (local[1] < -0.5 * _height)
  {
    const VECTOR3 bottomCap(0.0, -0.5 * _height, 0.0);
    const VECTOR3 diff = bottomCap - local;

    if (diff.norm() > _radius) return false;
    return true;
  }

  // if it's inside the top and bottom slabs, check the radius
  const REAL radius = sqrt(local[0] * local[0] + 
                           local[2] * local[2]);

  if (radius > _radius)
    return false;

  return true;
}

///////////////////////////////////////////////////////////////////////
// distance to the box
///////////////////////////////////////////////////////////////////////
REAL CAPSULE::distance(const VECTOR3& point) const 
{
  // transform back to local coordinates, but keep the scaling
  VECTOR3 local = worldVertexToLocal(point);

  // it must be closest to the top sphere
  if (local[1] > 0.5 * _height)
  {
    const VECTOR3 topCap(0.0, 0.5 * _height, 0.0);
    const VECTOR3 diff = topCap - local;

    return fabs(diff.norm() - _radius);
  }

  // check the bottom end sphere
  if (local[1] < -0.5 * _height)
  {
    const VECTOR3 bottomCap(0.0, -0.5 * _height, 0.0);
    const VECTOR3 diff = bottomCap - local;

    return fabs(diff.norm() - _radius);
  }

  // if it's inside the endcap slabs then it's just radius to the 
  // circular wall
  const REAL radius = sqrt(local[0] * local[0] + 
                           local[2] * local[2]);
  const REAL circularDistance = fabs(radius - _radius);

  return circularDistance;
}

///////////////////////////////////////////////////////////////////////
// signed distance to the box
///////////////////////////////////////////////////////////////////////
REAL CAPSULE::signedDistance(const VECTOR3& point) const
{
  const REAL sign = inside(point) ? -1.0 : 1.0;
  return sign * distance(point);
}

//////////////////////////////////////////////////////////////////////
// get the closest point on the cube, as well as the normal at the point
//////////////////////////////////////////////////////////////////////
void CAPSULE::getClosestPoint(const VECTOR3& query, VECTOR3& closestPoint, VECTOR3& normal) const
{
  const VECTOR3 local = worldVertexToLocal(query);
  
  // if it's above the top end sphere
  if (local[1] > 0.5 * _height)
  {
    const VECTOR3 topCap(0.0, 0.5 * _height, 0.0);
    const VECTOR3 direction = (local - topCap).normalized();

    closestPoint = topCap + direction * _radius;
    normal = closestPoint.normalized();
    return;
  }
  
  // if it's below the bottom end sphere
  if (local[1] < -0.5 * _height)
  {
    const VECTOR3 bottomCap(0.0, -0.5 * _height, 0.0);
    const VECTOR3 direction = (local - bottomCap).normalized();

    closestPoint = bottomCap + direction * _radius;
    normal = closestPoint.normalized();
    return;
  }

  // if it's between the end caps, get the nearest point on the y axis;
  const REAL radiusXZ = sqrt(local[0] * local[0] + local[2] * local[2]);
  closestPoint = local;
  closestPoint[0] *= _radius / radiusXZ;
  closestPoint[2] *= _radius / radiusXZ;

  normal = closestPoint;
  normal[1] = 0.0;
  normal.normalize();
  return;
}

} // HOBAK
