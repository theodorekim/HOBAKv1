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
///////////////////////////////////////////////////////////////////////////////////////////////////
// This is a file in the HOBAK library
//
// May 26, 2021 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "HYPERELASTIC.h"
#include "util/MATRIX_UTIL.h"

namespace HOBAK {
namespace VOLUME {

HYPERELASTIC::~HYPERELASTIC() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, then just recombine everything into F and call that version of Psi
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL HYPERELASTIC::psi(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return psi(U * Sigma.asDiagonal() * V.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call Psi
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL HYPERELASTIC::psi(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  return psi(U, Sigma, V);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, then just recombine everything into F and call that version of PK1
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3 HYPERELASTIC::PK1(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return PK1(U * Sigma.asDiagonal() * V.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call PK1
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3 HYPERELASTIC::PK1(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  return PK1(U, Sigma, V);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, then just recombine everything into F and call that version of Hessian
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 HYPERELASTIC::hessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return hessian(U * Sigma.asDiagonal() * V.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call Hessian
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 HYPERELASTIC::hessian(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  return hessian(U, Sigma, V);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, then just recombine everything into F and call that version of 
// ClampedHessian
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 HYPERELASTIC::clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return clampedHessian(U * Sigma.asDiagonal() * V.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call ClampedHessian
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 HYPERELASTIC::clampedHessian(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  return clampedHessian(U, Sigma, V);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL HYPERELASTIC::computeMu(const REAL E, const REAL nu)
{
  return E / (2.0 * (1.0 + nu));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL HYPERELASTIC::computeLambda(const REAL E, const REAL nu)
{
  return (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

} // VOLUME
} // HOBAK
