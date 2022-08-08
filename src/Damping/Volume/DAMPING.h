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
#ifndef VOLUME_DAMPING_H
#define VOLUME_DAMPING_H

#include "SETTINGS.h"

namespace HOBAK {
namespace VOLUME {

class DAMPING
{
public:
  // Computes the strain energy density
  virtual REAL psi(const MATRIX3& F, const MATRIX3& Fdot) const = 0;

  // Computes the first Piola-Kirchoff PK1 stress
  virtual MATRIX3 PK1(const MATRIX3& F, const MATRIX3& Fdot) const = 0;

  // The name of the material
  virtual std::string name() const = 0;

  // Computes the derivative of the PK1 stress
  virtual MATRIX9 hessian(const MATRIX3& F, const MATRIX3& Fdot) const = 0;

  // Computes the derivative of the PK1 stress, clamped to semi-positive definiteness
  virtual MATRIX9 clampedHessian(const MATRIX3& F, const MATRIX3& Fdot) const = 0;

  // The asymmetric term in damping that can occur because we take
  // a mixed derivative w.r.t. position and velocity. Still not sure how
  // valuable it is, but will include it here so we can do some testing
  virtual MATRIX9 positionGradient(const MATRIX3& F, const MATRIX3& Fdot) const { return MATRIX9::Zero(); };

  REAL& mu() { return _mu; };
  const REAL mu() const { return _mu; };

protected:
  REAL _mu;
};

}
}

#endif
