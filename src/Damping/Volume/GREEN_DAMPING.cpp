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
#include "GREEN_DAMPING.h"
#include "MATRIX_UTIL.h"

namespace HOBAK {
namespace VOLUME {

GREEN_DAMPING::GREEN_DAMPING(const REAL& mu)
{
  _mu = mu;
}

std::string GREEN_DAMPING::name() const
{
  return "GREEN_DAMPING";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL GREEN_DAMPING::psi(const MATRIX3& F, const MATRIX3& Fdot) const
{
  const MATRIX3 FdotF = Fdot.transpose() * F;
  const MATRIX3 Edot = 0.5 * (FdotF + FdotF.transpose());
  return _mu * Edot.squaredNorm();
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 GREEN_DAMPING::PK1(const MATRIX3& F, const MATRIX3& Fdot) const
{
  return _mu * F * (F.transpose() * Fdot + Fdot.transpose() * F);
}

///////////////////////////////////////////////////////////////////////
// there's probably a cleaner outer product form to be discovered
// here
///////////////////////////////////////////////////////////////////////
MATRIX9 GREEN_DAMPING::hessian(const MATRIX3& F, const MATRIX3& Fdot) const
{
  MATRIX9 pPpF = MATRIX9::Zero();
  int index = 0; 
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++, index++)
    {
      MATRIX3 pFpF;
      partialFpartialF(i, j, pFpF);
      MATRIX3 column =  F * F.transpose() * pFpF + F * pFpF.transpose() * F;
      pPpF.col(index) = flatten(column);
    }

  return _mu * pPpF;
}

///////////////////////////////////////////////////////////////////////
// filtered energy hessian
///////////////////////////////////////////////////////////////////////
MATRIX9 GREEN_DAMPING::clampedHessian(const MATRIX3& F, const MATRIX3& Fdot) const
{
  // analysis shows that it can never be indefinite, so no need to
  // filter
  //return hessian(F, Fdot);

  // doing some debugging ...
  return clampEigenvalues(hessian(F, Fdot));
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX9 GREEN_DAMPING::positionGradient(const MATRIX3& F, const MATRIX3& Fdot) const
{
  MATRIX9 pPpF = MATRIX9::Zero();
  int index = 0;
  MATRIX3 Fsum = F.transpose() * Fdot + Fdot.transpose() * F; 
  MATRIX3 Fproduct = F * Fdot.transpose();
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++, index++)
    {
      MATRIX3 pFpF;
      partialFpartialF(i, j, pFpF);
      MATRIX3 column =  pFpF * Fsum + F * pFpF.transpose() * Fdot + Fproduct * pFpF;
      pPpF.col(index) = flatten(column);
    }

  return _mu * pPpF;
}

}
}
