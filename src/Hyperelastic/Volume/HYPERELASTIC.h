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

#ifndef VOLUME_HYPERELASTIC_H
#define VOLUME_HYPERELASTIC_H

#include "SETTINGS.h"

namespace HOBAK {
namespace VOLUME {

class HYPERELASTIC
{
public:
  virtual ~HYPERELASTIC() = 0;

  // Computes the strain energy density
  virtual REAL psi(const MATRIX3& F) const;
  virtual REAL psi(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const;

  // Computes the first Piola-Kirchoff PK1 stress
  virtual MATRIX3 PK1(const MATRIX3& F) const;
  virtual MATRIX3 PK1(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const;

  // Computes the derivative of the PK1 stress
  virtual MATRIX9 hessian(const MATRIX3& F) const;
  virtual MATRIX9 hessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const; 

  // Computes the derivative of the PK1 stress, clamped to semi-positive definiteness
  virtual MATRIX9 clampedHessian(const MATRIX3& F) const;
  virtual MATRIX9 clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const;

  // The name of the material
  virtual std::string name() const = 0;

  // True if the energy computation requires the SVD of F
  virtual bool energyNeedsSVD() const = 0;

  // True if the PK1 computation requires the SVD of F
  virtual bool PK1NeedsSVD() const = 0;

  // convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
  static REAL computeMu(const REAL E, const REAL nu);
  static REAL computeLambda(const REAL E, const REAL nu);
};

}
}

#endif
