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
#include "ANISOTROPIC_DIRICHLET.h"
#include "MATRIX_UTIL.h"

namespace HOBAK {
namespace VOLUME {

ANISOTROPIC_DIRICHLET::ANISOTROPIC_DIRICHLET(const REAL& mu, const VECTOR3& a) :
    _mu(mu), _a(a)
{
  _a.normalize();
  _A = _a * _a.transpose();
}

std::string ANISOTROPIC_DIRICHLET::name() const
{
  return "ANISOTROPIC_DIRICHLET";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL ANISOTROPIC_DIRICHLET::psi(const MATRIX3& F) const
{
  return _mu * invariant5(F, _a);
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 ANISOTROPIC_DIRICHLET::PK1(const MATRIX3& F) const
{
  return 2.0 * _mu * F * _A;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_DIRICHLET::hessian(const MATRIX3& F) const
{
  MATRIX9 H = MATRIX9::Zero();
  
  // iterate over 3x3 blocks
  for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
      {
          const int x3 = 3 * x;
          const int y3 = 3 * y;

          for (int i = 0; i < 3; i++)
              H(x3 + i, y3 + i) = 2.0 * _mu * _A(x,y);
      }

  return H;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_DIRICHLET::clampedHessian(const MATRIX3& F) const
{
  // it can't go indefinite, so just return the Hessian
  return hessian(F);
}

bool ANISOTROPIC_DIRICHLET::energyNeedsSVD() const
{
    return false;
}

bool ANISOTROPIC_DIRICHLET::PK1NeedsSVD() const
{
    return false;
}

}
}
