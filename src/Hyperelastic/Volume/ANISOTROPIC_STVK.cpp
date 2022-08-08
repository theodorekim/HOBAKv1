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
#include "ANISOTROPIC_STVK.h"
#include "MATRIX_UTIL.h"

// DEBUG
#include <iostream>

namespace HOBAK {
namespace VOLUME {

using namespace std;

ANISOTROPIC_STVK::ANISOTROPIC_STVK(const REAL& mu, const VECTOR3& a) :
    _mu(mu), _a(a)
{
  _a.normalize();
  _A = _a * _a.transpose();
}

std::string ANISOTROPIC_STVK::name() const
{
  return "ANISOTROPIC_STVK";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL ANISOTROPIC_STVK::psi(const MATRIX3& F) const
{
  const REAL I5minusOne = (F * _a).squaredNorm() - 1.0;
  return 0.5 * _mu * I5minusOne * I5minusOne;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 ANISOTROPIC_STVK::PK1(const MATRIX3& F) const
{
  const REAL I5minusOne = (F * _a).squaredNorm() - 1.0;
  return (2.0 * _mu * I5minusOne) * (F * _A);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_STVK::hessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);
  MATRIX9 pPpF = (I5 - 1.0) * H5 + 2.0 * fa * fa.transpose();
  pPpF *= 2.0 * _mu;

  return pPpF;
}


///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_STVK::clampedHessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  const VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);

  const REAL lambda1 =           2.0 * _mu * (I5 - 1.0);
  const REAL lambda0 = lambda1 + 4.0 * _mu * I5;

  // try the case where it's all positive
  if (lambda0 > 0.0 && lambda1 > 0.0)
  {
    MATRIX9 pPpF = (I5 - 1.0) * H5 + 2.0 * fa * fa.transpose();
    pPpF *= 2.0 * _mu;
    return pPpF;
  }

  // if this one's negative, you're out of luck, just return all zeros
  if (lambda0 < 0.0) return MATRIX9::Zero();

  // the only possibility left is that a rank-one subspace is positive definite
  const VECTOR9 faNormed = fa.normalized();
  return lambda0 * (faNormed * faNormed.transpose());
}

bool ANISOTROPIC_STVK::energyNeedsSVD() const
{
  return true;
}

bool ANISOTROPIC_STVK::PK1NeedsSVD() const
{
  return true;
}

} // VOLUME
} // HOBAK
