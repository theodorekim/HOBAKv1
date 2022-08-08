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
#include "STVK.h"
#include "MATRIX_UTIL.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace VOLUME {

using namespace std;

STVK::STVK(const REAL& mu, const REAL& lambda ) :
    _lambda(lambda), _mu(mu)
{
}

std::string STVK::name() const
{
  return "STVK";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL STVK::psi(const MATRIX3& F) const
{
  const MATRIX3 E = 0.5 * (F.transpose() * F - MATRIX3::Identity());
  return _mu * E.squaredNorm() + _lambda * 0.5 * pow(E.trace(), 2.0); 
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 STVK::PK1(const MATRIX3& F) const
{
  // build the second PK
  const MATRIX3 E = 0.5 * (F.transpose() * F - MATRIX3::Identity());
  MATRIX3 S = _lambda * E.trace() * MATRIX3::Identity() + 2.0 * _mu * E;

  // convert to PK1
  return F * S;
}

///////////////////////////////////////////////////////////////////////
// Build out the eigensystem
//
// This is from Appendix F.1 in
// "Analytic Eigensystems for Isotropic Distortion Energies"
///////////////////////////////////////////////////////////////////////
void STVK::buildEigensystem(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V,
                            VECTOR9& eigenvalues, MATRIX9& eigenvectors) const
{
  const REAL I2 = invariant2(Sigma);
  const REAL front = -_mu + _lambda * 0.5 * (I2 - 3.0);
  const REAL s0Sq = Sigma[0] * Sigma[0];
  const REAL s1Sq = Sigma[1] * Sigma[1];
  const REAL s2Sq = Sigma[2] * Sigma[2];
  const REAL s0s1 = Sigma[0] * Sigma[1];
  const REAL s0s2 = Sigma[0] * Sigma[2];
  const REAL s1s2 = Sigma[1] * Sigma[2];

  // twists
  eigenvalues[0] = front + _mu * (s1Sq + s2Sq - s1s2);
  eigenvalues[1] = front + _mu * (s0Sq + s2Sq - s0s2);
  eigenvalues[2] = front + _mu * (s0Sq + s1Sq - s0s1);

  // flips
  eigenvalues[3] = front + _mu * (s1Sq + s2Sq + s1s2);
  eigenvalues[4] = front + _mu * (s0Sq + s2Sq + s0s2);
  eigenvalues[5] = front + _mu * (s0Sq + s1Sq + s0s1);

  // populate the scaling mode matrix
  MATRIX3 A;
  for (int i = 0; i < 3; i++)
    A(i,i) = front + (_lambda + 3.0 * _mu) * Sigma[i] * Sigma[i];

  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
      if (i != j)
        A(i,j) = _lambda * Sigma[i] * Sigma[j];

  // get the scaling modes
  MATRIX3 Q;
  VECTOR3 Lambda;
  eigensystem(A, Q, Lambda);

  eigenvalues[6] = Lambda[0];
  eigenvalues[7] = Lambda[1];
  eigenvalues[8] = Lambda[2];
    
  // Compute the eigenvectors
  buildTwistAndFlipEigenvectors(U, V, eigenvectors);
  buildScalingEigenvectors(U, Q, V, eigenvectors);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 STVK::hessian(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);

  VECTOR9 eigenvalues;
  MATRIX9 eigenvectors;
  buildEigensystem(U, Sigma, V, eigenvalues, eigenvectors);

  return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 STVK::clampedHessian(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);

  VECTOR9 eigenvalues;
  MATRIX9 eigenvectors;
  buildEigensystem(U, Sigma, V, eigenvalues, eigenvectors);

  for (int x = 0; x < 9; x++)
    if (eigenvalues[x] < 0.0)
      eigenvalues[x] = 0.0;

  return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
}

bool STVK::energyNeedsSVD() const
{
    return false;
}

bool STVK::PK1NeedsSVD() const
{
    return false;
}

} // VOLUME
} // HOBAK
