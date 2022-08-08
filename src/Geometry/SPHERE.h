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
#ifndef SPHERE_H
#define SPHERE_H

#include "KINEMATIC_SHAPE.h"

namespace HOBAK {

// positions are defined using R * S * x + t,
// starting from a sphere centered at (0,0,0), with radius of 1,
//
// Reminder: to make a new rotation matrix in Eigen, do:
// Eigen::AngleAxisd(0.1, VECTOR3::UnitX())
class SPHERE : public KINEMATIC_SHAPE
{
public:
  SPHERE(const VECTOR3& center, const REAL& scale);
	~SPHERE();

  virtual bool inside(const VECTOR3& point) const override;
  virtual REAL distance(const VECTOR3& point) const override;

  // remember that "inside" is negative with signed distance
  virtual REAL signedDistance(const VECTOR3& point) const override;

  // get the closest point on the object, as well as the normal at 
  // the point
  virtual void getClosestPoint(const VECTOR3& query, 
                               VECTOR3& closestPointLocal, 
                               VECTOR3& normalLocal) const override;
};

}

#endif
