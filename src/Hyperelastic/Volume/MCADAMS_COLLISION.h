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
#ifndef VOLUME_MCADAMS_COLLISION_H
#define VOLUME_MCADAMS_COLLISION_H

#include "VERTEX_FACE_COLLISION.h"

namespace HOBAK {
namespace VOLUME {

///////////////////////////////////////////////////////////////////////////////////
// This is the collision energy described in the SIGGRAPH 2011 paper 
//
// "Efficient elasticity for character skinning with contact and collisions"
//
// by Alex McAdams, Yongning Zhu, Andrew Selle, Mark Empey, Rasmus Tamstorf,
// Joseph Teran, and Eftychios Sifakis
//
// This energy is from the last paragraph of Section 7, and the x_s term there
// appears to be fixed; it's a proxy point sprinkled onto the body surface in
// the style of "Hybrid simulation of deformable solids" by Sifakis et al. 2007.
//
// The VERTEX_FACE_COLLISION class implements the 'unpinned' version where
// x_s is allowed to slide, this version implements the version where x_s
// is fixed and specified by barycentric coordinates. With the normal projection,
// it should work out to the (exact?) same thing.
///////////////////////////////////////////////////////////////////////////////////
class MCADAMS_COLLISION : public VERTEX_FACE_COLLISION
{
public:
  MCADAMS_COLLISION(const REAL& mu, const REAL& eps = 0.0);
  ~MCADAMS_COLLISION() {};

  // get the strain energy, the barycentric coordinate is the position
  // of the original projection
  REAL psi(const std::vector<VECTOR3>& v, const VECTOR3& bary) const;
  virtual REAL psi(const std::vector<VECTOR3>& v) const override;

  // This is the *gradient* of psi. The force is the *negative* gradient of psi.
  VECTOR12 gradient(const std::vector<VECTOR3>& v, const VECTOR3& bary) const;
  virtual VECTOR12 gradient(const std::vector<VECTOR3>& v) const override;

  virtual std::string name() const override;

  MATRIX12 hessian(const std::vector<VECTOR3>& v, const VECTOR3& bary) const;
  virtual MATRIX12 hessian(const std::vector<VECTOR3>& v) const override;

  MATRIX12 clampedHessian(const std::vector<VECTOR3>& v, const VECTOR3& bary) const;
  virtual MATRIX12 clampedHessian(const std::vector<VECTOR3>& v) const override;

private:
  virtual REAL psi(const VECTOR12& x, const VECTOR3& bary) const;
  virtual VECTOR12 gradient(const VECTOR12& x, const VECTOR3& bary) const;
  virtual MATRIX12 hessian(const VECTOR12& x, const VECTOR3& bary) const;
  MATRIX12 clampedHessian(const VECTOR12& x, const VECTOR3& bary) const;

  // gradient of spring length, n' * (v[2] - xs)
  static VECTOR12 barySpringLengthGradient(const std::vector<VECTOR3>& v,
                                           const std::vector<VECTOR3>& e,
                                           const VECTOR3& n,
                                           const VECTOR3& bary);

  // hessian of spring length, n' * (v[2] - xs)
  static MATRIX12 barySpringLengthHessian(const std::vector<VECTOR3>& v,
                                          const std::vector<VECTOR3>& e,
                                          const VECTOR3& n,
                                          const VECTOR3& bary);

protected:
  // partial of (v[0] - xs)
  static MATRIX3x12 tDiffPartial(const VECTOR3& bary);
};

} // VOLUME
} // HOBAK

#endif
