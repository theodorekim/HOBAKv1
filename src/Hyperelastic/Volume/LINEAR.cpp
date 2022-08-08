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
#include "LINEAR.h"
#include "MATRIX_UTIL.h"

namespace HOBAK {
namespace VOLUME {

LINEAR::LINEAR(const REAL& mu, const REAL& lambda)
: _mu(mu)
, _lambda(lambda)
{
  assert(_mu >= 0.0);
  assert(_lambda >= 0.0);
}

REAL LINEAR::psi(const MATRIX3& F) const
{
  const MATRIX3 eps = 0.5 * (F + F.transpose()) - MATRIX3::Identity();
  const REAL trEps = eps.trace();
  return _mu * eps.squaredNorm() + 0.5 * _lambda * trEps * trEps;
}

MATRIX3 LINEAR::PK1(const MATRIX3& F) const
{
  const MATRIX3 eps = 0.5 * (F + F.transpose()) - MATRIX3::Identity();
  const REAL trEps = eps.trace();
  return 2.0 * _mu * eps + _lambda * trEps * MATRIX3::Identity();
}

MATRIX9 LINEAR::hessian(const MATRIX3& F) const
{
  MATRIX9 H;
  H = _mu * MATRIX9::Identity();
  H(0,0) += _mu + _lambda;
  H(4,4) += _mu + _lambda;
  H(8,8) += _mu + _lambda;

  H(1,3) = H(3,1) = _mu;
  H(2,6) = H(6,2) = _mu;
  H(5,7) = H(7,5) = _mu;

  H(0,4) = H(0,8) = H(4,8) = _lambda;
  H(4,0) = H(8,0) = H(8,4) = _lambda;
  return H;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX9 LINEAR::clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return hessian(U * Sigma.asDiagonal() * V.transpose());
}

std::string LINEAR::name() const
{
  return "LINEAR";
}

bool LINEAR::energyNeedsSVD() const
{
  return false;
}

bool LINEAR::PK1NeedsSVD() const
{
  return false;
}

} // VOLUME
} // HOBAK
