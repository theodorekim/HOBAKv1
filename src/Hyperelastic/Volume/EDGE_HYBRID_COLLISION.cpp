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
#include "EDGE_HYBRID_COLLISION.h"
#include "EDGE_SQRT_COLLISION.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

EDGE_HYBRID_COLLISION::EDGE_HYBRID_COLLISION(const REAL& mu, const REAL& eps) :
  EDGE_COLLISION(mu, eps), _sqrt(mu, eps)
{
  //_separationEps = 1e-8;  // UNSTABLE
  //_separationEps = 1e-6;  // UNSTABLE
  
  // can't set this too small -- when the threshold is hit, it can
  // inject a huge force because the normalization went haywire
  _separationEps = 1e-4;
}

std::string EDGE_HYBRID_COLLISION::name() const
{
  return "EDGE_HYBRID_COLLISION";
}

///////////////////////////////////////////////////////////////////////
// should be punt computation to the cross-product version?
///////////////////////////////////////////////////////////////////////
bool EDGE_HYBRID_COLLISION::puntToCrossProduct(const VECTOR12& x, 
                                               const VECTOR2& a, 
                                               const VECTOR2& b) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x,v,e);

  // get the interpolated vertices
  const VECTOR3 va = (a[0] * v[0] + a[1] * v[1]);
  const VECTOR3 vb = (b[0] * v[2] + b[1] * v[3]);
  const VECTOR3 diff = vb - va;

  // get the normal
  VECTOR3 n = diff;

  if (n.norm() < _separationEps)
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL EDGE_HYBRID_COLLISION::psi(const VECTOR12& x,
                                const VECTOR2& a, 
                                const VECTOR2& b) const
{
  if (puntToCrossProduct(x,a,b))
    return EDGE_COLLISION::psi(x,a,b);
  return _sqrt.psi(x,a,b);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL EDGE_HYBRID_COLLISION::psiNegated(const VECTOR12& x,
                                       const VECTOR2& a, 
                                       const VECTOR2& b) const
{
  if (puntToCrossProduct(x,a,b))
    return EDGE_COLLISION::psiNegated(x,a,b);
  return _sqrt.psiNegated(x,a,b);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 EDGE_HYBRID_COLLISION::gradient(const VECTOR12& x, 
                                         const VECTOR2& a, 
                                         const VECTOR2& b) const
{
  if (puntToCrossProduct(x,a,b))
    return EDGE_COLLISION::gradient(x,a,b);
  return _sqrt.gradient(x,a,b);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 EDGE_HYBRID_COLLISION::gradientNegated(const VECTOR12& x, 
                                                const VECTOR2& a, 
                                                const VECTOR2& b) const
{
  if (puntToCrossProduct(x,a,b))
    return EDGE_COLLISION::gradientNegated(x,a,b);
  return _sqrt.gradientNegated(x,a,b);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 EDGE_HYBRID_COLLISION::hessian(const VECTOR12& x,
                                        const VECTOR2& a, 
                                        const VECTOR2& b) const
{
  if (puntToCrossProduct(x,a,b))
    return EDGE_COLLISION::hessian(x,a,b);
  return _sqrt.hessian(x,a,b);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 EDGE_HYBRID_COLLISION::hessianNegated(const VECTOR12& x,
                                               const VECTOR2& a, 
                                               const VECTOR2& b) const
{
  if (puntToCrossProduct(x,a,b))
    return EDGE_COLLISION::hessianNegated(x,a,b);
  return _sqrt.hessianNegated(x,a,b);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void EDGE_HYBRID_COLLISION::setEps(const REAL& eps)
{
  _eps = eps;
  _sqrt.setEps(eps);
}

} // VOLUME
} // HOBAK
