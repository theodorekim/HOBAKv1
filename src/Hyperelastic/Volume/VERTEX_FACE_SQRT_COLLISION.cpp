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
#include "VERTEX_FACE_SQRT_COLLISION.h"
#include "MATRIX_UTIL.h"
#include "COLLISION_UTIL.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

VERTEX_FACE_SQRT_COLLISION::VERTEX_FACE_SQRT_COLLISION(const REAL& mu, const REAL& eps) :
  MCADAMS_COLLISION(mu, eps)
{
  _inverseEps = 1e-8;
}

std::string VERTEX_FACE_SQRT_COLLISION::name() const
{
  return "VERTEX_FACE_SQRT_COLLISION";
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL VERTEX_FACE_SQRT_COLLISION::psi(const vector<VECTOR3>& v) const
{
  const VECTOR3 bary = getInsideBarycentricCoordinates(v);
  return psi(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 VERTEX_FACE_SQRT_COLLISION::gradient(const vector<VECTOR3>& v) const
{
  const VECTOR3 bary = getInsideBarycentricCoordinates(v);
  return gradient(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 VERTEX_FACE_SQRT_COLLISION::hessian(const vector<VECTOR3>& v) const
{
  const VECTOR3 bary = getInsideBarycentricCoordinates(v);
  return hessian(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 VERTEX_FACE_SQRT_COLLISION::clampedHessian(const vector<VECTOR3>& v) const
{
  const MATRIX12 H = hessian(v);
  return clampEigenvalues(H);
}

///////////////////////////////////////////////////////////////////////
// should we reverse the direction of the force?
///////////////////////////////////////////////////////////////////////
bool VERTEX_FACE_SQRT_COLLISION::reverse(const vector<VECTOR3>& v,
                                         const vector<VECTOR3>& e) const
{
  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();

  // e[1] is already the collision vertex recentered to the origin
  // (v[0] - v[2])
  const REAL dotted = n.dot(e[1]);

  return (dotted < 0) ? true : false;
}

#if USING_SOFT_MU
///////////////////////////////////////////////////////////////////////
// Catmull-Rom cubic
///////////////////////////////////////////////////////////////////////
REAL smoothstep(const REAL x)
{
  return (3.0 * x - 2.0 * x * x) * x;
}

REAL lerp(const REAL x)
{
  return x;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL softMu(const REAL mu, const REAL tMagnitude, const REAL eps)
{
  // if the collision distance is zero, there's no notion of a
  // soft mu
  if (eps <= 0.0) return mu;
  const REAL diff = (eps - tMagnitude) / eps;

  // if we're imore than halfway past the collision eps, the
  // mu gets applied at full power
  if (diff > 0.5)
    return mu;

  //return mu * smoothstep(2.0 * diff);
  return mu * lerp(2.0 * diff);
}
#endif

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL VERTEX_FACE_SQRT_COLLISION::psi(const VECTOR12& x, const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  const bool reversal = reverse(v,e);

  const VECTOR3 xs = bary[0] * v[1] + bary[1] * v[2] + bary[2] * v[3];
  const VECTOR3 t = v[0] - xs;
  const REAL tMagnitude = sqrt(t.dot(t));
  const REAL springDiff = (reversal) ? tMagnitude + _eps : tMagnitude - _eps;

  return _mu * springDiff * springDiff;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 VERTEX_FACE_SQRT_COLLISION::gradient(const VECTOR12& x, const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  const bool reversal = reverse(v,e);
  
  // remember we had to reorder vertices in a wonky way
  const VECTOR3 xs = bary[0] * v[1] + bary[1] * v[2] + bary[2] * v[3];
  const VECTOR3 t = v[0] - xs;
  const REAL tMagnitude = sqrt(t.dot(t));
  //const REAL springDiff = tMagnitude - _eps;
  const REAL springDiff = (reversal) ? tMagnitude + _eps : tMagnitude - _eps;
  const MATRIX3x12 tDiff = tDiffPartial(bary); 

  // if everything has become undefined, just give up
  const REAL tDott = t.dot(t);
  if (fabs(tMagnitude) <= _inverseEps || fabs(tDott) < _inverseEps)
    return VECTOR12::Zero();

  const VECTOR12 result = 2.0 * _mu * springDiff * (1.0 / tMagnitude) * tDiff.transpose() * t;

  // could instead try to trap all the inverses and hand back something fixed up,
  // but consistency is not guaranteed, so let's just zero it out at the first
  // sign of trouble
  //const REAL tMagnitudeInv = (fabs(tMagnitude) > _inverseEps) ? 1.0 / tMagnitude : 0.0;
  //const VECTOR12 result = 2.0 * _mu * springDiff * tMagnitudeInv * tDiff.transpose() * t;

#if ENABLE_DEBUG_TRAPS
  if (result.hasNaN())
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    cout << " springDiff: " << springDiff << endl;
    cout << " tMagnitude: " << tMagnitude << endl;
    cout << " tDiff: " << endl << tDiff << endl;
    cout << " result: " << result << endl;
  }
#endif

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 VERTEX_FACE_SQRT_COLLISION::hessian(const VECTOR12& x, const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  const bool reversal = reverse(v,e);
  
  // remember we had to reorder vertices in a wonky way
  const VECTOR3 xs = bary[0] * v[1] + bary[1] * v[2] + bary[2] * v[3];
  const VECTOR3 t = v[0] - xs;
  const REAL tDott = t.dot(t);
  const REAL tMagnitude = sqrt(tDott);
  //const REAL springDiff = tMagnitude - _eps;
  const REAL springDiff = (reversal) ? tMagnitude + _eps : tMagnitude - _eps;
  const MATRIX3x12 tDiff = tDiffPartial(bary); 

  // get the spring length, non-zero rest-length
  const VECTOR12 product = tDiff.transpose() * t;

  // if everything has become undefined, just give up
  if (fabs(tMagnitude) <= _inverseEps || fabs(tDott) < _inverseEps)
    return MATRIX12::Zero();

  return 2.0 * _mu * ((1.0 / tDott - springDiff / (tDott * tMagnitude)) * (product * product.transpose()) +
                      (springDiff / tMagnitude) * tDiff.transpose() * tDiff); 

  // could instead try to trap all the inverses and hand back something fixed up,
  // but consistency is not guaranteed, so let's just zero it out at the first
  // sign of trouble
  //const REAL tMagnitudeInv = (fabs(tMagnitude) > _inverseEps) ? 1.0 / tMagnitude : 0.0;
  //const REAL tDottInv = (fabs(tDott) > _inverseEps) ? 1.0 / tDott : 1.0;
  //return 2.0 * _mu * ((tDottInv - springDiff / (tDott * tMagnitude)) * (product * product.transpose()) +
  //                    (springDiff * tMagnitudeInv) * tDiff.transpose() * tDiff); 
}

} // VOLUME
} // HOBAK
