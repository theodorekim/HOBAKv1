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
#include "EDGE_SQRT_COLLISION.h"
#include "MATRIX_UTIL.h"
#include "COLLISION_UTIL.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

EDGE_SQRT_COLLISION::EDGE_SQRT_COLLISION(const REAL& mu, const REAL& eps) :
  EDGE_COLLISION(mu, eps)
{
  _tooSmall = 1e-7;
}

std::string EDGE_SQRT_COLLISION::name() const
{
  return "EDGE_SQRT_COLLISION";
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL EDGE_SQRT_COLLISION::psi(const VECTOR12& x,
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
  if ((vb - va).norm() < _tooSmall)
    return 0.0;

  const REAL springLength = _eps - (vb - va).norm();
  return _mu * springLength * springLength;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL EDGE_SQRT_COLLISION::psiNegated(const VECTOR12& x,
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
  if ((vb - va).norm() < _tooSmall)
    return 0.0;

  const REAL springLength = _eps + (vb - va).norm();
  return _mu * springLength * springLength;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 EDGE_SQRT_COLLISION::gradient(const VECTOR12& x, 
                                       const VECTOR2& a, 
                                       const VECTOR2& b) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x,v,e);

  assert(v.size() == 4);
  assert(e.size() == 2);
 
  // get the interpolated vertices
  const VECTOR3 va = (a[0] * v[0] + a[1] * v[1]);
  const VECTOR3 vb = (b[0] * v[2] + b[1] * v[3]);
  const VECTOR3 diff = vb - va;

  // if the two are co-linear, give up
  // should probably fall back to cross-product formula here
  // (see EDGE_HYBRID_COLLISION)
  if (diff.norm() < _tooSmall)
    return VECTOR12::Zero();

  // get the normal
  VECTOR3 n = diff;
  n = n / n.norm();

  const REAL springLength = _eps - diff.norm();
  return -2.0 * _mu * springLength * (vDiffPartial(a,b).transpose() * n);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 EDGE_SQRT_COLLISION::gradientNegated(const VECTOR12& x, 
                                              const VECTOR2& a, 
                                              const VECTOR2& b) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x,v,e);

  assert(v.size() == 4);
  assert(e.size() == 2);
 
  // get the interpolated vertices
  const VECTOR3 va = (a[0] * v[0] + a[1] * v[1]);
  const VECTOR3 vb = (b[0] * v[2] + b[1] * v[3]);
  const VECTOR3 diff = vb - va;

  // if the two are co-linear, give up
  // should probably fall back to cross-product formula here
  // (see EDGE_HYBRID_COLLISION)
  if (diff.norm() < _tooSmall)
    return VECTOR12::Zero();

  // get the direction
  VECTOR3 d = diff;
  d = d / d.norm();

  const REAL springLength = _eps + diff.norm();
  const MATRIX3x12 vPartial = vDiffPartial(a,b);
 
  return 2.0 * _mu * springLength * (vPartial.transpose() * d);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 EDGE_SQRT_COLLISION::hessian(const VECTOR12& x,
                                      const VECTOR2& a, 
                                      const VECTOR2& b) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  assert(v.size() == 4);
  assert(e.size() == 2);

  // get the interpolated vertices
  const VECTOR3 va = (a[0] * v[0] + a[1] * v[1]);
  const VECTOR3 vb = (b[0] * v[2] + b[1] * v[3]);
  const VECTOR3 diff = vb - va;
  const REAL diffNorm = diff.norm();

  // if the two are co-linear, give up
  // should probably fall back to cross-product formula here
  // (see EDGE_HYBRID_COLLISION)
  if (diffNorm < _tooSmall)
    return MATRIX12::Zero();

  // get the normal
  VECTOR3 d = diff;
  d = d / d.norm();

  const MATRIX3x12 vPartial = vDiffPartial(a,b);
  const REAL invNorm = (diffNorm >= 1e-8) ? 1.0 / diffNorm : 1.0;
  const REAL invNorm3 = invNorm * invNorm * invNorm;

  const VECTOR12 normPartial = -invNorm * (vPartial.transpose() * diff);
  const MATRIX3x12 dGrad = invNorm * vPartial -
                           invNorm3 * diff * (vPartial.transpose() * diff).transpose();

  return -2.0 * _mu * ((_eps - diffNorm) * (vPartial.transpose() * dGrad) +
                       (normPartial) * (vPartial.transpose() * d).transpose());
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 EDGE_SQRT_COLLISION::hessianNegated(const VECTOR12& x,
                                             const VECTOR2& a, 
                                             const VECTOR2& b) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  assert(v.size() == 4);
  assert(e.size() == 2);

  // get the interpolated vertices
  const VECTOR3 va = (a[0] * v[0] + a[1] * v[1]);
  const VECTOR3 vb = (b[0] * v[2] + b[1] * v[3]);
  const VECTOR3 diff = vb - va;
  const REAL diffNorm = diff.norm();
  const REAL diffNorm3 = diffNorm * diffNorm * diffNorm;

  // if the two are co-linear, give up
  // should probably fall back to cross-product formula here
  // (see EDGE_HYBRID_COLLISION)
  if (diffNorm < _tooSmall)
    return MATRIX12::Zero();

  // get the normal
  VECTOR3 n = diff;
  n = n / n.norm();

  const MATRIX3x12 vPartial = vDiffPartial(a,b);
  const VECTOR12 normPartial = (-1.0 / diffNorm) * (vPartial.transpose() * diff);

  const MATRIX3x12 nGrad = (1.0 / diffNorm) * vPartial -
                           (1.0 / diffNorm3) * diff * (vPartial.transpose() * diff).transpose();

  // this is the energetically consistent one
  return 2.0 * _mu * ((_eps + diffNorm) * (vPartial.transpose() * nGrad) -
                       (normPartial) * (vPartial.transpose() * n).transpose());
}

} // VOLUME
} // HOBAK
