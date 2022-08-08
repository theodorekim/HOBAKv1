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
#include "SNH_WITH_BARRIER.h"
#include "MATRIX_UTIL.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

SNH_WITH_BARRIER::SNH_WITH_BARRIER(const REAL& mu, const REAL& lambda)
// reparamaterizing to be consistent with linear, see Section 3.4
// of "Stable Neo-Hookean Flesh Simulation", right before Eqn. 16
: _mu((4.0 / 3.0) * mu)
, _lambda(lambda + (5.0 / 6.0) * mu)
, _alpha(1.0 + _mu / _lambda - _mu / (4.0 * _lambda))
{
  assert(_mu > 0.0);
  assert(_lambda > 0.0);
  assert(_alpha > 0.0);
}

//////////////////////////////////////////////////////////////////////////////
// Eqn. 14 from Section 3.3 in "Stable Neo-Hookean Flesh Simulation"
//////////////////////////////////////////////////////////////////////////////
REAL SNH_WITH_BARRIER::psi(const MATRIX3& F) const
{
  const REAL Ic = F.squaredNorm();
  const REAL Jminus1 = F.determinant() - _alpha;
  return 0.5 * (_mu * (Ic - 3.0) + _lambda * Jminus1 * Jminus1 - _mu * log(Ic + 1.0));
}


//////////////////////////////////////////////////////////////////////////////
// Eqn. 18 from Section 4.2 in "Stable Neo-Hookean Flesh Simulation"
//////////////////////////////////////////////////////////////////////////////
MATRIX3 SNH_WITH_BARRIER::PK1(const MATRIX3& F) const
{
  const REAL Ic = F.squaredNorm();
  const MATRIX3 pJpF = partialJpartialF(F);
  const REAL Jminus1 = F.determinant() - _alpha;
  return _mu * F + _lambda * Jminus1 * pJpF - _mu / (Ic + 1) * F;
}

//////////////////////////////////////////////////////////////////////////////
// Eqn. 41 from 4.6 in "Stable Neo-Hookean Flesh Simulation"
//////////////////////////////////////////////////////////////////////////////
MATRIX9 SNH_WITH_BARRIER::hessian(const MATRIX3& F) const
{
  const REAL I2 = invariant2(F);
  const REAL I3 = invariant3(F);
  const REAL I2PlusOne = I2 + 1.0;
  const VECTOR9 pjpf = flatten(partialJpartialF(F));
  const VECTOR9 f = flatten(F);

  const REAL scale = _lambda * (I3 - _alpha);
  const MATRIX3 f0hat = crossProduct(F,0) * scale;
  const MATRIX3 f1hat = crossProduct(F,1) * scale;
  const MATRIX3 f2hat = crossProduct(F,2) * scale;

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

  return _mu * (1.0 - 1.0 / (I2 + 1.0)) * MATRIX9::Identity() + 
         _lambda * pjpf * pjpf.transpose() + 
         _mu * (2.0 / (I2PlusOne * I2PlusOne)) * f * f.transpose() + hessJ;
}

//////////////////////////////////////////////////////////////////////////////
// This is not what is actually what is used in production at Pixar.
// For more details, see the SNH class.
//////////////////////////////////////////////////////////////////////////////
MATRIX9 SNH_WITH_BARRIER::clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  VECTOR9 eigenvalues;
  MATRIX9 eigenvectors;

  const REAL J = invariant3(Sigma);
  const REAL I2 = invariant2(Sigma);

  // 0-2 are twist
  // 3-5 are flip
  // 6-8 are scaling
  const REAL front = -_mu / (I2 + 1.0) + _mu;
  const REAL back = _lambda * (J - _alpha);
  eigenvalues[0] = front + back * Sigma[0];
  eigenvalues[1] = front + back * Sigma[1];
  eigenvalues[2] = front + back * Sigma[2];
  eigenvalues[3] = front - back * Sigma[0];
  eigenvalues[4] = front - back * Sigma[1];
  eigenvalues[5] = front - back * Sigma[2];

  // Scaling eigenvalues
  MATRIX3 A;
  const REAL s0s0 = Sigma(0) * Sigma(0);
  const REAL s1s1 = Sigma(1) * Sigma(1);
  const REAL s2s2 = Sigma(2) * Sigma(2);
  const REAL I2plusOneSq = (I2 + 1.0) * (I2 + 1.0);
  const REAL I2plusOneSqInv = 1.0 / I2plusOneSq;
  const REAL twoMu = 2.0 * _mu;
  const REAL frontDiag = (- _mu * (I2 + 1.0)) * I2plusOneSqInv + _mu;
  const REAL frontOffDiag = (2.0 * J - _alpha) * _lambda;

  A(0,0) = frontDiag + twoMu * s0s0 * I2plusOneSqInv + _lambda * s1s1 * s2s2;
  A(1,1) = frontDiag + twoMu * s1s1 * I2plusOneSqInv + _lambda * s0s0 * s2s2;
  A(2,2) = frontDiag + twoMu * s2s2 * I2plusOneSqInv + _lambda * s0s0 * s1s1;
  A(0,1) = frontOffDiag * Sigma[2] + (twoMu * Sigma[0] * Sigma[1]) * I2plusOneSqInv;
  A(0,2) = frontOffDiag * Sigma[1] + (twoMu * Sigma[0] * Sigma[2]) * I2plusOneSqInv;
  A(1,2) = frontOffDiag * Sigma[0] + (twoMu * Sigma[1] * Sigma[2]) * I2plusOneSqInv;
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

std::string SNH_WITH_BARRIER::name() const
{
  return "SNH_WITH_BARRIER";
}

bool SNH_WITH_BARRIER::energyNeedsSVD() const
{
  return false;
}

bool SNH_WITH_BARRIER::PK1NeedsSVD() const
{
  return false;
}

} // VOLUME
} // HOBAK
