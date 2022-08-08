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
#ifndef VOLUME_EDGE_SQRT_COLLISION_H
#define VOLUME_EDGE_SQRT_COLLISION_H

#include "SETTINGS.h"
#include <vector>

#include "EDGE_COLLISION.h"

namespace HOBAK {
namespace VOLUME {

///////////////////////////////////////////////////////////////////////////////////
// This is (sort of) the edge-edge collision energy described in:
//
// "Dynamic deformables: implementation and production practicalities"
// Kim and Eberle, 2020, Chapter 11.4
//
// However, unlike the treatments in other work:
//
// "Robust treatement of simultaneous contacts"
// Harmon et al., 2008, Section 2
//
// "Collision and self-collision handling in cloth model dedicated to design garments"
// Provot 1997, Section 3.3
//
// we don't treat the normal as the cross product of the two edges, and instead
// do something simpler where we take the difference between the barycentric
// locations along the two edges, and use that as the collision normal.
///////////////////////////////////////////////////////////////////////////////////
class EDGE_SQRT_COLLISION : public EDGE_COLLISION
{
public:
  EDGE_SQRT_COLLISION(const REAL& mu, const REAL& eps = 0.0);
  ~EDGE_SQRT_COLLISION() {};

  // get the strain energy
  virtual REAL psi(const VECTOR12& x,
                   const VECTOR2& a, const VECTOR2& b) const override;
  virtual REAL psiNegated(const VECTOR12& x, 
                          const VECTOR2& a, const VECTOR2& b) const override;

  // This is the *gradient* of psi. The force is the *negative* gradient of psi.
  virtual VECTOR12 gradient(const VECTOR12& x, 
                            const VECTOR2& a, const VECTOR2& b) const override;
  virtual VECTOR12 gradientNegated(const VECTOR12& x, 
                                   const VECTOR2& a, const VECTOR2& b) const override;

  virtual MATRIX12 hessian(const VECTOR12& x, 
                           const VECTOR2& a, const VECTOR2& b) const override;
  virtual MATRIX12 hessianNegated(const VECTOR12& x, 
                                  const VECTOR2& a, const VECTOR2& b) const override;

  virtual std::string name() const override;

private:

  // is the repulsion direction vector too small, and we should give up?
  REAL _tooSmall;
};

} // VOLUME
} // HOBAK

#endif
