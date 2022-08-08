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
#ifndef VOLUME_EDGE_HYBRID_COLLISION_H
#define VOLUME_EDGE_HYBRID_COLLISION_H

#include "SETTINGS.h"
#include <vector>

#include "EDGE_COLLISION.h"
#include "EDGE_SQRT_COLLISION.h"

namespace HOBAK {
namespace VOLUME {

///////////////////////////////////////////////////////////////////////////////////
// This energy tries to get the best of both worlds between
//  EDGE_COLLISION and EDGE_SQRT_COLLISION
//
// The barycentric version seems the more robust overall, but the 
// cross-product-based EDGE_COLLISION version does better when the
// edges are really close together.
//
// So, we'll call EDGE_SQRT_COLLISION in general, unless the direction is really
// small, and then we'll call EDGE_COLLISION
///////////////////////////////////////////////////////////////////////////////////
class EDGE_HYBRID_COLLISION : public EDGE_COLLISION
{
public:
  EDGE_HYBRID_COLLISION(const REAL& mu, const REAL& eps = 0.0);
  ~EDGE_HYBRID_COLLISION() {};

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
  virtual void setEps(const REAL& eps) override;

private:
  // should be punt computation to the cross-product version?
  bool puntToCrossProduct(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;

  // can call EDGE_COLLISION via the parent class, but need this sibling around too
  EDGE_SQRT_COLLISION _sqrt;

  // separation threshold to punt from EDGE_SQRT_COLLISION to EDGE_COLLISION
  REAL _separationEps;
};

} // VOLUME
} // HOBAK

#endif
