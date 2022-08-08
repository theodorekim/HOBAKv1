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
#include "ANISOTROPIC_FUNG.h"
#include "MATRIX_UTIL.h"
#include <iostream>

namespace HOBAK {
namespace VOLUME {

using namespace std;

ANISOTROPIC_FUNG::ANISOTROPIC_FUNG(const REAL& mu, const VECTOR3& a) :
    _mu(mu), _a(a)
{
  _a.normalize();
  _A = _a * _a.transpose();
}

std::string ANISOTROPIC_FUNG::name() const
{
  return "ANISOTROPIC_FUNG";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
// Psi = mu * (exp(I5 - 1) - I5);
///////////////////////////////////////////////////////////////////////
REAL ANISOTROPIC_FUNG::psi(const MATRIX3& F) const
{
  const REAL I5 = (F * _a).squaredNorm();
  return _mu * (exp(I5 - 1.0) - I5);
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// PK1 = mu * (exp(I5 - 1) - 1) * 2 * F * A;
///////////////////////////////////////////////////////////////////////
MATRIX3 ANISOTROPIC_FUNG::PK1(const MATRIX3& F) const
{
  const REAL I5 = (F * _a).squaredNorm();
  return 2.0 * _mu * (exp(I5 - 1.0) - 1.0) * F * _A;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
// fa = vec(F * A);
// I3 = eye(3,3);
// H5 = [A(1,1) * I3 A(1,2) * I3 A(1,3) * I3;
//       A(2,1) * I3 A(2,2) * I3 A(2,3) * I3;
//       A(3,1) * I3 A(3,2) * I3 A(3,3) * I3];
// H = 2 * mu * ((exp(I5 - 1) - 1) * H5 + 2 * exp(I5 - 1) * fa * fa');
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_FUNG::hessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  const VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);
  const REAL I5minusOne = I5 - 1.0;
  MATRIX9 pPpF = (exp(I5minusOne) - 1.0) * H5 + 2.0 * exp(I5minusOne) * (fa * fa.transpose());
  pPpF *= 2.0 * _mu;

  return pPpF;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_FUNG::clampedHessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  const VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);
  const REAL lambda1 =           2.0 * _mu * (exp(I5 - 1.0) - 1.0);
  const REAL lambda0 = lambda1 + 4.0 * _mu * I5 * exp(I5 - 1.0);

  // try the case where it's all positive
  if (lambda0 > 0.0 && lambda1 > 0.0)
  {
    const REAL I5minusOne = I5 - 1.0;
    MATRIX9 pPpF = (exp(I5minusOne) - 1.0) * H5 + 2.0 * exp(I5minusOne) * (fa * fa.transpose());
    pPpF *= 2.0 * _mu;
    return pPpF;
  }

  // if this one's negative, you're out of luck, just return all zeros
  if (lambda0 < 0.0) return MATRIX9::Zero();

  // the only possibility left is that a rank-one subspace is positive definite
  const VECTOR9 faNormed = fa.normalized();
  return lambda0 * (faNormed * faNormed.transpose());
}

bool ANISOTROPIC_FUNG::energyNeedsSVD() const
{
  return true;
}

bool ANISOTROPIC_FUNG::PK1NeedsSVD() const
{
  return true;
}

} // VOLUME
} // HOBAK
