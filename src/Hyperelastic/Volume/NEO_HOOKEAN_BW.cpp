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
#include "NEO_HOOKEAN_BW.h"
#include "MATRIX_UTIL.h"

namespace HOBAK {
namespace VOLUME {

NEO_HOOKEAN_BW::NEO_HOOKEAN_BW(const REAL& mu, const REAL& lambda)
: _mu(mu)
, _lambda(lambda)
, _alpha(1.0 + _mu / _lambda)
{
  assert(_mu > 0.0);
  assert(_lambda > 0.0);
}

REAL NEO_HOOKEAN_BW::psi(const MATRIX3& F) const
{
  const REAL Ic = F.squaredNorm();
  const REAL J = invariant3(F);
  const REAL logJ = log(J);

  return _mu * 0.5 * (Ic - 3.0) - _mu * logJ + _lambda * 0.5 * logJ * logJ;
}

MATRIX3 NEO_HOOKEAN_BW::PK1(const MATRIX3& F) const
{
  const MATRIX3 pJpF = partialJpartialF(F);
  const REAL J = invariant3(F);
  return _mu * (F - (1.0 / J) * pJpF) + _lambda * log(J) * (1.0 / J) * pJpF;
}

MATRIX9 NEO_HOOKEAN_BW::hessian(const MATRIX3& F) const
{
  const VECTOR9 pjpf = flatten(partialJpartialF(F));

  const REAL J = invariant3(F);
  const REAL logJ = log(J);
  const MATRIX3 f0hat = crossProduct(F,0);
  const MATRIX3 f1hat = crossProduct(F,1);
  const MATRIX3 f2hat = crossProduct(F,2);

  // build the fractal cross-product 
  MATRIX9 hessJ;
  hessJ.setZero();
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    {
      hessJ(i, j + 3) = -f2hat(i,j);
      hessJ(i + 3, j) =  f2hat(i,j);

      hessJ(i, j + 6) =  f1hat(i,j);
      hessJ(i + 6, j) = -f1hat(i,j);

      hessJ(i + 3, j + 6) = -f0hat(i,j);
      hessJ(i + 6, j + 3) =  f0hat(i,j);
    }

  return _mu * MATRIX9::Identity() + ((_mu + _lambda * (1.0 - logJ)) / (J * J))* pjpf * pjpf.transpose() + 
         ((_lambda * logJ - _mu) / J) * hessJ;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX9 NEO_HOOKEAN_BW::clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  VECTOR9 eigenvalues;
  MATRIX9 eigenvectors;

  const REAL J = invariant3(Sigma);
  const REAL logJ = log(J);
  const REAL s0s1inv = 1.0 / (Sigma(0) * Sigma(1));
  const REAL s0s2inv = 1.0 / (Sigma(0) * Sigma(2));
  const REAL s1s2inv = 1.0 / (Sigma(1) * Sigma(2));

  // 0-2 are twist
  // 3-5 are flip
  // 6-8 are scaling
  const REAL front = _lambda * logJ - _mu;
  eigenvalues[0] = front * s1s2inv + _mu;
  eigenvalues[1] = front * s0s2inv + _mu;
  eigenvalues[2] = front * s0s1inv + _mu;
  eigenvalues[3] = -front * s1s2inv + _mu;
  eigenvalues[4] = -front * s0s2inv + _mu;
  eigenvalues[5] = -front * s0s1inv + _mu;

  // populate matrix for scaling eigenvalues
  MATRIX3 A;
  const REAL s0s0 = Sigma(0) * Sigma(0);
  const REAL s1s1 = Sigma(1) * Sigma(1);
  const REAL s2s2 = Sigma(2) * Sigma(2);
  const REAL frontDiag = _lambda * (1.0 - logJ) + _mu;
  A(0,0) = frontDiag / s0s0 + _mu;
  A(1,1) = frontDiag / s1s1 + _mu;
  A(2,2) = frontDiag / s2s2 + _mu;

  A(0,1) = _lambda * s0s1inv;
  A(0,2) = _lambda * s0s2inv;
  A(1,2) = _lambda * s1s2inv;
  A(1,0) = A(0,1);
  A(2,0) = A(0,2);
  A(2,1) = A(1,2);

  // get the scaling eigenvalues
  const Eigen::SelfAdjointEigenSolver<MATRIX3> Aeigs(A);
  for (int x = 0; x < 3; x++)
    eigenvalues[x + 6] = Aeigs.eigenvalues()[x];

  // Compute the eigenvectors
  buildTwistAndFlipEigenvectors(U, V, eigenvectors);
  buildScalingEigenvectors(U, Aeigs.eigenvectors(), V, eigenvectors);

  // Clamp the eigenvalues
  for (int i = 0; i < 9; i++)
    if (eigenvalues(i) < 0.0)
      eigenvalues(i) = 0.0;

  return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
}

std::string NEO_HOOKEAN_BW::name() const
{
  return "NEO_HOOKEAN_BW";
}

bool NEO_HOOKEAN_BW::energyNeedsSVD() const
{
  return false;
}

bool NEO_HOOKEAN_BW::PK1NeedsSVD() const
{
  return false;
}

} // VOLUME
} // HOBAK
