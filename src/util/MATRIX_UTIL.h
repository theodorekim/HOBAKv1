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
#ifndef MATRIX_UTIL_H
#define MATRIX_UTIL_H

#include "SETTINGS.h"

namespace HOBAK {

  // convert a MATRIX3 to a VECTOR9 in a consistent way
  VECTOR9 flatten(const MATRIX3& A);

  // convert a VECTOR9 to a MATRIX3 in a consistent way
  MATRIX3 unflatten(const VECTOR9& v);

  // rotation variant of the SVD where the reflections are loaded into
  // Sigma and not U and V
  void svd_rv(const MATRIX3& F, MATRIX3& U, VECTOR3& Sigma, MATRIX3& V);

  // get the polar decomposition of matrix A = RS
  void polarDecomposition(const MATRIX3& A, MATRIX3& R, MATRIX3& S);

  // clamp the eigenvalues of a 9x9 to semi-positive-definite
  MATRIX9 clampEigenvalues(const MATRIX9& A);
  MATRIX12 clampEigenvalues(const MATRIX12& A);
  MATRIX12 clampEigenvaluesToSemiNegative(const MATRIX12& A);

  // get the eigensystem of a 3x3 matrix
  void eigensystem(const MATRIX3& A, MATRIX3& Q, VECTOR3& Lambda);
  void eigensystem(const MATRIX9& A, MATRIX9& Q, VECTOR9& Lambda);
  VECTOR9 eigenvalues(const MATRIX9& A);
  VECTOR12 eigenvalues(const MATRIX12& A);
  VECTOR eigenvalues(const MATRIX& A);
  VECTOR eigenvalues(const SPARSE_MATRIX& A);

  // use Spectra to get the largest eig of a big matrix
  REAL largestEigenvalue(const SPARSE_MATRIX& A);
  
  // use Spectra to get the smallest eig of a big matrix
  REAL smallestEigenvalue(const SPARSE_MATRIX& A);

  // Let's make some random deformation gradients
  MATRIX3 randomMatrix3(const REAL scaling = 3.0);

  // Let's make some random positive-definite deformation gradients
  MATRIX3 randomPositiveDefiniteMatrix3(const REAL scaling = 3.0);

  // Let's make some random directions
  VECTOR3 randomVector3(const REAL scaling = 3.0);

  // Let's make some random directions
  VECTOR12 randomVector12(const REAL scaling = 10.0);

  // Let's make some random rotations
  MATRIX3 randomRotation();

  // Let's make a random barycentric coordinate
  VECTOR2 randomBarycentric();

  // Matrix double-contraction
  REAL ddot(const MATRIX3& A, const MATRIX3& B);

  // eigenvectors 0-2 are the twist modes
  // eigenvectors 3-5 are the flip modes
  void buildTwistAndFlipEigenvectors(const MATRIX3& U, const MATRIX3& V, 
                                     MATRIX9& Q);

  // eigenvectors 6-8 are the scaling modes, jackpot version
  void buildScalingEigenvectors(const MATRIX3& U, const MATRIX3& V, 
                                MATRIX9& Q9);

  // eigenvectors 6-8 are the scaling modes, non-jackpot version
  void buildScalingEigenvectors(const MATRIX3& U, const MATRIX3& Q,
                                const MATRIX3& V, MATRIX9& Q9);

  // get the Kronecker product of matrix with respect to 3x3 identity,
  // used a lot in anisotropic materials
  MATRIX9 kronIdentity(const MATRIX3& A);

  // Tensor invariants
  REAL invariant1(const MATRIX3& F);
  REAL invariant2(const MATRIX3& F);
  REAL invariant3(const MATRIX3& F);
  REAL invariant4(const MATRIX3& F, const VECTOR3& a);
  REAL invariant5(const MATRIX3& F, const VECTOR3& a);

  REAL invariant2(const VECTOR3& Sigma);
  REAL invariant3(const VECTOR3& Sigma);

  // rotation gradient, w.r.t. deformation gradient F
  // \frac{\partial R}{\partial F}
  MATRIX9 rotationGradient(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V);

  // time derivative of rotation
  // \frac{\partial R}{\partial t}
  MATRIX3 rotationDot(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V, const MATRIX3& Fdot);

  // Eqn. 19 from Section 4.2 in "Stable Neo-Hookean Flesh Simulation"
  MATRIX3 partialJpartialF(const MATRIX3& F);

  // Eqn. 29 from Section 4.5 in "Stable Neo-Hookean Flesh Simulation"
  MATRIX3 crossProduct(const MATRIX3& F, const int col);

  // 3rd order tensor derivative of deformation gradient F with respect to itself
  void partialFpartialF(const int i, const int j, MATRIX3& pFpF);
}

#endif
