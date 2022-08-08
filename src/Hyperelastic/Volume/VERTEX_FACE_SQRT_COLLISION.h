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
#ifndef VOLUME_VERTEX_FACE_SQRT_COLLISION_H
#define VOLUME_VERTEX_FACE_SQRT_COLLISION_H

#include "MCADAMS_COLLISION.h"

namespace HOBAK {
namespace VOLUME {

///////////////////////////////////////////////////////////////////////////////////
// This is the super-basic, difference based vertex-face collision energy
// described in in "Collision Energies" chapter of
//
// "Dynamic Deformables: Implementation and Production Practicalities"
//
// No tangent sliding, no nothing. Super-basic.
///////////////////////////////////////////////////////////////////////////////////
class VERTEX_FACE_SQRT_COLLISION : public MCADAMS_COLLISION
{
public:
  VERTEX_FACE_SQRT_COLLISION(const REAL& mu, const REAL& eps = 0.0);
  ~VERTEX_FACE_SQRT_COLLISION() {};

  virtual REAL psi(const std::vector<VECTOR3>& v) const override;
  virtual VECTOR12 gradient(const std::vector<VECTOR3>& v) const override;
  virtual MATRIX12 hessian(const std::vector<VECTOR3>& v) const override;
  virtual MATRIX12 clampedHessian(const std::vector<VECTOR3>& v) const override;

  virtual std::string name() const override;
protected:
  virtual REAL psi(const VECTOR12& x, const VECTOR3& bary) const override;
  virtual VECTOR12 gradient(const VECTOR12& x, const VECTOR3& bary) const override;
  virtual MATRIX12 hessian(const VECTOR12& x, const VECTOR3& bary) const override;

  // should we reverse the direction of the force?
  bool reverse(const std::vector<VECTOR3>& v, const std::vector<VECTOR3>& e) const;

  // what's the divide-by-zero threshold where we zero out the force?
  REAL _inverseEps;
};

} // VOLUME
} // HOBAK

#endif
