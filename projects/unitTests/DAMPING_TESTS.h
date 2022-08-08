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
//////////////////////////////////////////////////////////////////////////////
// do a damping convergence test on dPdF
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestDampingHessian(const HOBAK::VOLUME::DAMPING* material, const MATRIX3 &F, const MATRIX3& Fdot)
{
  cout << "=============================================================== " << endl;
  std::cout << " VERIFYING Damping dPdF for " << material->name().c_str() << std::endl;
  cout << "=============================================================== " << endl;

  MATRIX9 dPdF = material->hessian(F, Fdot);
  MATRIX3 P = material->PK1(F, Fdot);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-9)
  {
    MATRIX9 finiteDiff;
    int column = 0;

    // for each of the degrees of the freedom
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++, column++)
      {
        MATRIX3 FdotNew = Fdot;
        FdotNew(x,y) += eps;

        // get the new psi
        MATRIX3 Pnew = material->PK1(F, FdotNew); 

        // store the finite difference
        MATRIX3 diff = (Pnew - P) / eps;
        finiteDiff.col(column) = flatten(diff);
      }

    MATRIX9 diff = dPdF - finiteDiff;
    REAL diffNorm = fabs(diff.norm() / P.norm());
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    std::cout << "eps: " << eps << " diff: " << diffNorm << std::endl;

    if (e == 5 && minSeen > 1e-6)
    {
      std::cout << " TEST FAILED!!!!!" << endl;
      std::cout << " dPdF: " << std::endl << dPdF << std::endl;
      std::cout << " finite diff: " << std::endl << finiteDiff << std::endl;
      std::cout << " diff: " << std::endl << diff << std::endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    std::cout << " TEST PASSED. " << endl;
    return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////
// do a damping convergence test on d^2Psi / dFdFdot
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestDampingPositionGradient(const HOBAK::VOLUME::DAMPING* material, const MATRIX3 &F, const MATRIX3& Fdot)
{
  cout << "=============================================================== " << endl;
  std::cout << " VERIFYING Damping position gradient for " << material->name().c_str() << std::endl;
  cout << "=============================================================== " << endl;

  MATRIX9 H = material->positionGradient(F, Fdot);
  MATRIX3 P = material->PK1(F, Fdot);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-9)
  {
    MATRIX9 finiteDiff;
    int column = 0;

    // for each of the degrees of the freedom
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++, column++)
      {
        MATRIX3 Fnew = F;
        Fnew(x,y) += eps;

        // get the new psi
        MATRIX3 Pnew = material->PK1(Fnew, Fdot); 

        // store the finite difference
        MATRIX3 diff = (Pnew - P) / eps;
        finiteDiff.col(column) = flatten(diff);
      }

    MATRIX9 diff = H - finiteDiff;
    REAL diffNorm = fabs(diff.norm() / P.norm());
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    std::cout << "eps: " << eps << " diff: " << diffNorm << std::endl;

    if (e == 5 && minSeen > 1e-6)
    {
      std::cout << " TEST FAILED!!!!!" << endl;
      std::cout << " H: " << std::endl << H << std::endl;
      std::cout << " finite diff: " << std::endl << finiteDiff << std::endl;
      std::cout << " diff: " << std::endl << diff << std::endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    std::cout << " TEST PASSED. " << endl;
    return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on damping energy's PK1
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestDampingPK1(const HOBAK::VOLUME::DAMPING* material, const MATRIX3 &F, const MATRIX3& Fdot)
{
  cout << "=============================================================== " << endl;
  std::cout << " VERIFYING Damping P for " << material->name().c_str() << std::endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psi(F, Fdot);
  MATRIX3 P = material->PK1(F, Fdot);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-9)
  {
    MATRIX3 finiteDiffP;

    // for each of the degrees of the freedom
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
      {
        MATRIX3 FdotNew = Fdot;
        FdotNew(x,y) += eps;

        double psi = material->psi(F, FdotNew);

        // store the finite difference
        finiteDiffP(x, y) = (psi - psi0) / eps;
      }

    MATRIX3 diff = P - finiteDiffP;
    REAL diffNorm = fabs(diff.norm() / P.norm());
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    std::cout << "eps: " << eps << " diff: " << diffNorm << std::endl;

    if (e == 5 && minSeen > 1e-6)
    {
      std::cout << " TEST FAILED!!!!!" << endl;
      std::cout << " P: " << std::endl << P << std::endl;
      std::cout << " finite diff: " << std::endl << finiteDiffP << std::endl;
      std::cout << " diff: " << std::endl << diff << std::endl;

      MATRIX3 div = P;
      for (int y = 0; y < 3; y++)
        for (int x = 0; x < 3; x++)
          div(x,y) = P(x,y) / finiteDiffP(x,y);
      std::cout << " div: " << std::endl << div << std::endl;
      return false;
    }
      
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    std::cout << " TEST PASSED. " << endl;
    return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works
//////////////////////////////////////////////////////////////////////////////
bool testClampedDampingHessian(const HOBAK::VOLUME::DAMPING* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  const MATRIX3 UF = randomRotation();
  const MATRIX3 VF = randomRotation();
  
  const MATRIX3 Udot = randomRotation();
  const MATRIX3 Vdot = randomRotation();
  VECTOR3 sigmasF(2.0, 3.0, 5.0);
  VECTOR3 sigmasDot(2.0, 3.0, 5.0);

  MATRIX3 F, Fdot;
  MATRIX9 original, clamped, diff;
  REAL diffNorm;

  // verify that with all positive, we still get the same Hessian
  F = UF * sigmasF.asDiagonal() * VF.transpose();
  Fdot = Udot * sigmasDot.asDiagonal() * Vdot.transpose();
  original = clampEigenvalues(material->hessian(F, Fdot));
  clamped  = material->clampedHessian(F, Fdot);
  diff = original - clamped;
  diffNorm = fabs(diff.norm() / original.norm());
  cout << " Diff: " << diffNorm << "\t";

  if (diffNorm < 1e-8 || original.norm() < 1e-8)
    cout << " All positive Sigma test: PASSED " << endl;
  else
  {
    cout << " All positive Sigma test: FAILED " << endl;
    cout << " diff: " << endl;
    cout << diff << endl;
    cout << " Numerical Clamp: " << endl;
    cout << original << endl;
    cout << " Analytic Clamp: " << endl;
    cout << clamped << endl;

    cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
    cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
    return false;
  }
  
  // try the single reflection case, we still get the same Hessian
  sigmasF = VECTOR3(-2.0, 3.0, 5.0);
  sigmasDot = VECTOR3(2.0, -3.0, 5.0);
  F = UF * sigmasF.asDiagonal() * VF.transpose();
  Fdot = Udot * sigmasDot.asDiagonal() * Vdot.transpose();
  original = clampEigenvalues(material->hessian(F, Fdot));
  clamped  = material->clampedHessian(F, Fdot);
  diff = original - clamped;
  diffNorm = fabs(diff.norm() / original.norm());
  cout << " Diff: " << diffNorm << "\t";

  if (diffNorm < 1e-8 || original.norm() < 1e-8)
    cout << " Single negative Sigma test: PASSED " << endl;
  else
  {
    cout << " Single negative Sigma test: FAILED " << endl;
    cout << " diff: " << endl;
    cout << diff << endl;
    cout << " Numerical Clamp: " << endl;
    cout << original << endl;
    cout << " Analytic Clamp: " << endl;
    cout << clamped << endl;

    cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
    cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
    return false;
  }

  // note that the double reflection case is a rotation, so don't 
  // need to check it

  // try the triple reflection case, we still get the same Hessian
  sigmasF = VECTOR3(-2.0, -3.0, -5.0);
  sigmasDot = VECTOR3(-5.0, -2.0, -3.0);
  F = UF * sigmasF.asDiagonal() * VF.transpose();
  Fdot = Udot * sigmasDot.asDiagonal() * Vdot.transpose();
  original = clampEigenvalues(material->hessian(F, Fdot));
  clamped  = material->clampedHessian(F, Fdot);
  diff = original - clamped;
  diffNorm = fabs(diff.norm() / original.norm());
  cout << " Diff: " << diffNorm << "\t";

  if (diffNorm < 1e-8 || original.norm() < 1e-8)
    cout << " Triple negative Sigma test: PASSED " << endl;
  else
  {
    cout << " Triple negative Sigma test: FAILED " << endl;
    cout << " diff: " << endl;
    cout << diff << endl;
    cout << " Numerical Clamp: " << endl;
    cout << original << endl;
    cout << " Analytic Clamp: " << endl;
    cout << clamped << endl;

    cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
    cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
    return false;
  }

  // try the all-under-compression case, see if it still holds 
  sigmasF = VECTOR3(0.2, 0.3, 0.5);
  sigmasDot = VECTOR3(0.5, 0.2, 0.3);
  F = UF * sigmasF.asDiagonal() * VF.transpose();
  Fdot = Udot * sigmasDot.asDiagonal() * Vdot.transpose();
  original = clampEigenvalues(material->hessian(F, Fdot));
  clamped  = material->clampedHessian(F, Fdot);
  diff = original - clamped;
  diffNorm = fabs(diff.norm() / original.norm());
  cout << " Diff: " << diffNorm << "\t";

  if (diffNorm < 1e-8 || original.norm() < 1e-8)
    cout << " Compression test: PASSED " << endl;
  else
  {
    MATRIX9 unclamped = material->hessian(F, Fdot);
    cout << " Compression test: FAILED " << endl;
    cout << " Unclamped: " << endl << unclamped << endl;
    cout << " Numerical clamped: " << endl << original << endl;
    cout << " Analytic clamped: " << endl << clamped << endl;
    cout << " Diff: " << endl << diff << endl;

    cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
    cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
    cout << " Unclamped eigs: " << eigenvalues(unclamped).transpose() << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// test out damping energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Damping Energy tests", "[damping]" )
{
  using namespace HOBAK::VOLUME;
  MATRIX3 F = randomMatrix3();
  MATRIX3 Fdot = randomMatrix3();
  SECTION("Testing Green Damping")
  {
    // Fung explodes easily for large deformation
    GREEN_DAMPING green(1.0);
    REQUIRE(convergenceTestDampingPK1(&green, F, Fdot) == true);
    REQUIRE(convergenceTestDampingHessian(&green, F, Fdot) == true);
    REQUIRE(convergenceTestDampingPositionGradient(&green, F, Fdot) == true);
    REQUIRE(testClampedDampingHessian(&green) == true);
  }
}
