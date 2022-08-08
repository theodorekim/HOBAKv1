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
#ifndef KINEMATIC_SHAPE_H
#define KINEMATIC_SHAPE_H

#include "SETTINGS.h"

namespace HOBAK {

// positions are defined using R * S * x + t,
// starting from a primitive centered at (0,0,0), with radius of 1,
//
// Reminder: to make a new rotation matrix in Eigen, do:
// Eigen::AngleAxisd(0.1, VECTOR3::UnitX())
class KINEMATIC_SHAPE
{
public:
  KINEMATIC_SHAPE() { _name = std::string("KINEMATIC_SHAPE"); };
	virtual ~KINEMATIC_SHAPE() {};

  virtual bool inside(const VECTOR3& point) const = 0;
  virtual REAL distance(const VECTOR3& point) const = 0;

  // remember that "inside" is negative with signed distance
  virtual REAL signedDistance(const VECTOR3& point) const = 0;

  const MATRIX3& rotation() const     { return _rotation; }; 
  const MATRIX3& scale() const        { return _scale; }; 
  const MATRIX3& scaleInverse() const { return _scaleInverse; }; 
  const VECTOR3& translation() const  { return _translation; }; 
  MATRIX3& rotation()     { return _rotation; }; 
  MATRIX3& scale()        { return _scale; }; 
  MATRIX3& scaleInverse() { return _scaleInverse; }; 
  VECTOR3& translation()  { return _translation; }; 
  const std::string& name() const { return _name; }; 

  // transform a vertex from the local space to world.
  // this is used to track how a constraint has moved. If the object
  // moves, we can then see where a specific point on the object has 
  // now moved
  virtual VECTOR3 localVertexToWorld(const VECTOR3& local) const
  {
    return _rotation * _scale * local + _translation;
  };

  // transform a vertex from world space to local.
  // this is used to track constraints. If we want to attach a node
  // to a specific place on the kinematic object, we need a local i
  // coordinate in the object's local space. Then later, when the i
  // object moves, we can get the relative world position by calling 
  // localToWorld
  virtual VECTOR3 worldVertexToLocal(const VECTOR3& world) const
  {
    return _scaleInverse * _rotation.transpose() * (world - _translation);
  };

  // transform a normal from the local space to world.
  // The reasoning is identical to localVertexToWorld
  VECTOR3 localNormalToWorld(const VECTOR3& normal) const { return _rotation * normal; };

  // get the closest point on the object, as well as the normal at 
  // the point
  virtual void getClosestPoint(const VECTOR3& query, 
                               VECTOR3& closestPointLocal, 
                               VECTOR3& normalLocal) const = 0;
protected:
  //VECTOR3 _center;
  MATRIX3 _scale;
  MATRIX3 _scaleInverse;
  MATRIX3 _rotation;
  VECTOR3 _translation;

  std::string _name;
};

}

#endif
