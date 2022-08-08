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
#ifndef CYLINDER_H
#define CYLINDER_H

#include "KINEMATIC_SHAPE.h"

namespace HOBAK {

// cylinder positions are defined using R * S * x + t,
// starting from a cylinder centered at (0,0,0), running up and down
// the y-axis with height 1 and radius 1
//
// Reminder: to make a new rotation matrix in Eigen, do:
// Eigen::AngleAxisd(0.1, VECTOR3::UnitX())
class CYLINDER : public KINEMATIC_SHAPE
{
public:
  CYLINDER(const VECTOR3& center, const REAL& radius, const REAL& height);
	~CYLINDER();

  const REAL radius() const { return _radius; };
  const REAL height() const { return _height; };

  virtual bool inside(const VECTOR3& point) const override;
  virtual REAL distance(const VECTOR3& point) const override;

  // for the cylinder it's slightly easier if we don't apply the scaling
  // via a matrix in the local-to-world transform
  virtual VECTOR3 localVertexToWorld(const VECTOR3& local) const override
  {
    return _rotation * local + _translation;
  };
  virtual VECTOR3 worldVertexToLocal(const VECTOR3& world) const override
  {
    return _rotation.transpose() * (world - _translation);
  };

  // remember that "inside" is negative with signed distance
  virtual REAL signedDistance(const VECTOR3& point) const override;

  // get the closest point on the cube, as well as the normal at the point
  virtual void getClosestPoint(const VECTOR3& query, 
                               VECTOR3& closestPointLocal, 
                               VECTOR3& normalLocal) const override;

protected:
  REAL _radius;
  REAL _height;
};

}

#endif
