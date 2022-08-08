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
#include "MCADAMS_COLLISION.h"
#include "MATRIX_UTIL.h"
#include "COLLISION_UTIL.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

MCADAMS_COLLISION::MCADAMS_COLLISION(const REAL& mu, const REAL& eps) :
  VERTEX_FACE_COLLISION(mu, eps)
{
}

std::string MCADAMS_COLLISION::name() const
{
  return "MCADAMS_COLLISION";
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL MCADAMS_COLLISION::psi(const VECTOR12& x, const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);

  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();

  const VECTOR3 xs = bary[0] * v[1] + bary[1] * v[2] + bary[2] * v[3];
  const VECTOR3 t = v[0] - xs;
  
  // get the spring length, non-zero rest-length
  REAL springLength = t.dot(n) - _eps;
  return _mu * springLength * springLength;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL MCADAMS_COLLISION::psi(const vector<VECTOR3>& v,
                            const VECTOR3 & bary) const
{
  return psi(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL MCADAMS_COLLISION::psi(const vector<VECTOR3>& v) const
{
  const VECTOR3 bary = getBarycentricCoordinates(v);
  return psi(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 MCADAMS_COLLISION::gradient(const vector<VECTOR3>& v,
                                     const VECTOR3 & bary) const
{
  return gradient(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 MCADAMS_COLLISION::gradient(const vector<VECTOR3>& v) const
{
  const VECTOR3 bary = getBarycentricCoordinates(v);
  return gradient(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 MCADAMS_COLLISION::hessian(const vector<VECTOR3>& v,
                                    const VECTOR3 & bary) const
{
  return hessian(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 MCADAMS_COLLISION::hessian(const vector<VECTOR3>& v) const
{
  const VECTOR3 bary = getBarycentricCoordinates(v);
  return hessian(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 MCADAMS_COLLISION::clampedHessian(const vector<VECTOR3>& v, 
                                           const VECTOR3& bary) const
{
  return clampedHessian(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 MCADAMS_COLLISION::clampedHessian(const vector<VECTOR3>& v) const
{
  const VECTOR3 bary = getBarycentricCoordinates(v);
  return clampedHessian(flattenVertices(v), bary);
}

///////////////////////////////////////////////////////////////////////
// gradient of spring length, n' * (v[2] - v[0])
///////////////////////////////////////////////////////////////////////
VECTOR12 MCADAMS_COLLISION::barySpringLengthGradient(const std::vector<VECTOR3>& v,
                                                     const std::vector<VECTOR3>& e,
                                                     const VECTOR3& n,
                                                     const VECTOR3& bary)
{
  MATRIX3x12 nPartial = normalGradientVF(e);
  MATRIX3x12 tPartial = tDiffPartial(bary);

  // remember we had to reorder vertices in a wonky way
  const VECTOR3 xs = bary[0] * v[1] + bary[1] * v[2] + bary[2] * v[3];
  const VECTOR3 t = v[0] - xs;

  //f = nPartial' * (v2 - v0) + vPartial' * n;
  return nPartial.transpose() * t + tPartial.transpose() * n;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 MCADAMS_COLLISION::gradient(const VECTOR12& x, const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  
  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();
  
  // remember we had to reorder vertices in a wonky way
  const VECTOR3 xs = bary[0] * v[1] + bary[1] * v[2] + bary[2] * v[3];
  const VECTOR3 t = v[0] - xs;

  // get the spring length, non-zero rest-length
  REAL springLength = t.dot(n) - _eps;
  return 2.0 * _mu * springLength * barySpringLengthGradient(v,e,n,bary);
}

///////////////////////////////////////////////////////////////////////
// hessian of spring length, n' * (v[2] - v[0])
///////////////////////////////////////////////////////////////////////
MATRIX12 MCADAMS_COLLISION::barySpringLengthHessian(const vector<VECTOR3>& v,
                                                    const vector<VECTOR3>& e,
                                                    const VECTOR3& n,
                                                    const VECTOR3& bary)
{
  // remember we had to reorder vertices in a wonky way
  const VECTOR3 xs = bary[0] * v[1] + bary[1] * v[2] + bary[2] * v[3];
  const VECTOR3 t = v[0] - xs;

  MATRIX3x12 tPartial = tDiffPartial(bary);

  //% mode-3 contraction
  //[nx ny nz] = normal_hessian(x);
  //final = nx * delta(1) + ny * delta(2) + nz * delta(3);
  vector<MATRIX12> normalH = normalHessianVF(e);

  MATRIX12 contracted = t[0] * normalH[0] + 
                        t[1] * normalH[1] + 
                        t[2] * normalH[2];
  
  //nGrad= normal_gradient(x);
  MATRIX3x12 nGrad = normalGradientVF(e);

  //product = nGrad' * vGrad;
  //final = final + product + product';
  MATRIX12 product = nGrad.transpose() * tPartial;

  return contracted + product + product.transpose();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 MCADAMS_COLLISION::hessian(const VECTOR12& x, const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  
  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();

  // remember we had to reorder vertices in a wonky way
  const VECTOR3 xs = bary[0] * v[1] + bary[1] * v[2] + bary[2] * v[3];
  const VECTOR3 t = v[0] - xs;
  
  // get the spring length, non-zero rest-length
  REAL springLength = t.dot(n) - _eps;

  // ndotGrad    = ndot_gradient(x);
  VECTOR12 springLengthGrad = barySpringLengthGradient(v,e,n,bary);

  // ndotHessian = ndot_hessian(x);
  MATRIX12 springLengthH = barySpringLengthHessian(v,e,n,bary);
  
  // final = 2 * k * (ndotGrad * ndotGrad' + ndot * ndotHessian);
  return 2.0 * _mu * (springLengthGrad * springLengthGrad.transpose() + 
                      springLength * springLengthH);
  // Gauss-Newton approximation
  //return 2.0 * _mu * (springLengthGrad * springLengthGrad.transpose());
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 MCADAMS_COLLISION::clampedHessian(const VECTOR12& x, const VECTOR3& bary) const
{
  return clampEigenvalues(hessian(x, bary));
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX3x12 MCADAMS_COLLISION::tDiffPartial(const VECTOR3& bary)
{
  MATRIX3x12 tPartial;
  tPartial.setZero();
  tPartial(0,0) = tPartial(1,1)  = tPartial(2,2) = 1.0;
  tPartial(0,3) = tPartial(1,4)  = tPartial(2,5) = -bary[0];
  tPartial(0,6) = tPartial(1,7)  = tPartial(2,8) = -bary[1];
  tPartial(0,9) = tPartial(1,10) = tPartial(2,11) = -bary[2];

  return tPartial;
}

} // VOLUME
} // HOBAK
