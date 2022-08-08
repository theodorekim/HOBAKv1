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
using namespace HOBAK;
using namespace std;

//////////////////////////////////////////////////////////////////////////////
// do a convergence test to see if hessians are working on tet meshes
// this does assume that convergenceTestTetMeshForces passed. It is possible
// that things can become inconsistent with the original energy
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestTetMeshHessian(const HOBAK::VOLUME::HYPERELASTIC& material, TET_MESH& tetMesh)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING TET_MESH Hessian " << endl;
  cout << "=============================================================== " << endl;
  vector<VECTOR3>& vertices = tetMesh.vertices();

  // randomize the vertices
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] += randomVector3(1.0);
  
  // cache the current state
  const vector<VECTOR3> v0 = vertices;

  // get the analytic forces
  tetMesh.computeFs();
  MATRIX hessian = -1.0 * MATRIX(tetMesh.computeHyperelasticHessian(material));

  // get the current force
  VECTOR f0 = tetMesh.computeHyperelasticForces(material);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-9)
  {
    MATRIX finiteDiff(f0.size(), f0.size());
    finiteDiff.setZero();

    int index = 0;
    for (unsigned int y = 0; y < vertices.size(); y++)
    {
      for (unsigned int x = 0; x < 3; x++, index++)
      {
        // reset the other components
        vertices[y] = v0[y];

        // just perturb one component
        vertices[y][x] += eps;
       
        // compute the energy 
        tetMesh.computeFs();
        VECTOR f1 = tetMesh.computeHyperelasticForces(material);

        finiteDiff.col(index) = -(f1 - f0) / eps;
      }

      // last reset
      vertices[y] = v0[y];
    }

    MATRIX diff = hessian - finiteDiff;
    REAL diffNorm = fabs(diff.norm() / hessian.norm());
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 5 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      return false;
    }
    else
      
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test to see if forces are working on tet meshes
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestTetMeshForces(const HOBAK::VOLUME::HYPERELASTIC& material, TET_MESH& tetMesh)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING TET_MESH forces " << endl;
  cout << "=============================================================== " << endl;
  vector<VECTOR3>& vertices = tetMesh.vertices();

  // randomize the vertices
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] += randomVector3(1.0);
  
  // cache the current state
  const vector<VECTOR3> v0 = vertices;

  // get the analytic forces
  tetMesh.computeFs();
  VECTOR forces = tetMesh.computeHyperelasticForces(material);

  // get the current energy
  REAL psi0 = tetMesh.computeHyperelasticEnergy(material);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-9)
  {
    VECTOR finiteDiff(forces.size());
    finiteDiff.setZero();

    int index = 0;
    for (unsigned int y = 0; y < vertices.size(); y++)
    {
      for (unsigned int x = 0; x < 3; x++, index++)
      {
        // reset the other components
        vertices[y] = v0[y];

        // just perturb one component
        vertices[y][x] += eps;
       
        // compute the energy 
        tetMesh.computeFs();
        REAL psi1 = tetMesh.computeHyperelasticEnergy(material);

        finiteDiff[index] = -(psi1 - psi0) / eps;
      }

      // last reset
      vertices[y] = v0[y];
    }

    VECTOR diff = forces - finiteDiff;
    REAL diffNorm = fabs(diff.norm() / forces.norm());
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 5 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      //cout << " forces: " << endl << forces.transpose() << endl;
      //cout << " finite diff: " << endl << finiteDiff.transpose() << endl;
      //cout << " diff: " << endl << diff.transpose() << endl;
      return false;
    }
    else
      
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works
//
// invertible means that the material supports inversion (i.e. poking
// inside out), so we can test with negative singular values
//////////////////////////////////////////////////////////////////////////////
bool testClampedHessian(const HOBAK::VOLUME::HYPERELASTIC* material, bool invertible = true)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  const MATRIX3 U = randomRotation();
  const MATRIX3 V = randomRotation();
  VECTOR3 sigmas(2.0, 3.0, 5.0);

  MATRIX9 original, clamped, diff;
  REAL diffNorm;

  // verify that with all positive, we still get the same Hessian
  original = clampEigenvalues(material->hessian(U, sigmas, V));
  clamped  = material->clampedHessian(U, sigmas, V);
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
  
  if (invertible)
  {
    // try the single reflection case, we still get the same Hessian
    sigmas = VECTOR3(-2.0, 3.0, 5.0);
    original = clampEigenvalues(material->hessian(U, sigmas, V));
    clamped  = material->clampedHessian(U, sigmas, V);
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
    sigmas = VECTOR3(-2.0, -3.0, -5.0);
    original = clampEigenvalues(material->hessian(U, sigmas, V));
    clamped  = material->clampedHessian(U, sigmas, V);
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
  }
  else
  {
    cout << " Material does NOT support inversion, so skipping those cases." << endl;
  }

  // if it's not invertible, don't pass it an inverted F
  const REAL leftBound = (invertible) ? -1.0 : 0.01;

  // try lots of Fs here, as this case tends to be tricky
  std::mt19937 gen(314159);
  std::uniform_real_distribution<REAL> dist(leftBound, 1.0);
  for (int x = 0; x < 20; x++)
  {
    // try the all-under-compression case, see if it still holds 
    sigmas = VECTOR3(dist(gen), dist(gen), dist(gen));
    original = clampEigenvalues(material->hessian(U, sigmas, V));
    clamped  = material->clampedHessian(U, sigmas, V);
    diff = original - clamped;
    diffNorm = fabs(diff.norm() / original.norm());
    cout << " Diff: " << diffNorm << "\t";

    if (diffNorm < 1e-8 || original.norm() < 1e-8)
      cout << " Compression test: PASSED " << endl;
    else
    {
      MATRIX9 unclamped = material->hessian(U, sigmas, V);
      cout << " Compression test: FAILED " << endl;
      cout << " Unclamped: " << endl << unclamped << endl;
      cout << " Numerical clamped: " << endl << original << endl;
      cout << " Analytic clamped: " << endl << clamped << endl;
      cout << " Diff: " << endl << diff << endl;

      cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
      cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
      cout << " Unclamped eigs: " << eigenvalues(unclamped).transpose() << endl;
      cout << " Sigmas:         " << sigmas.transpose() << endl;
      return false;
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessian(const HOBAK::VOLUME::HYPERELASTIC* material, const MATRIX3 &F)
{
  using namespace HOBAK;
  using namespace std;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  MATRIX9 dPdF = material->hessian(U, Sigma, V);
  MATRIX3 P = material->PK1(U, Sigma, V);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
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
        MATRIX3 Unew,Vnew;
        VECTOR3 Snew;
        svd_rv(Fnew, Unew, Snew, Vnew);
        MATRIX3 Pnew = material->PK1(Unew, Snew, Vnew);

        // store the finite difference
        MATRIX3 diff = (Pnew - P) / eps;
        finiteDiff.col(column) = flatten(diff);
      }

    MATRIX9 diff = dPdF - finiteDiff;
    REAL diffNorm = (fabs(diff.norm() / P.norm())) / 81.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    MATRIX9 div = finiteDiff;
    for (int y = 0; y < 9; y++)
      for (int x = 0; x < 9; x++)
        div(x,y) = div(x,y) / dPdF(x,y);
    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " dPdF: " << endl << dPdF << endl;
      cout << " finite diff: " << endl << finiteDiff << endl;
      cout << " diff: " << endl << diff << endl;
      cout << " div: " << endl << div << endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    cout << " TEST PASSED. " << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the PK1
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestPK1(const HOBAK::VOLUME::HYPERELASTIC* material, const MATRIX3 &F)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING P for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  REAL psi0 = material->psi(U, Sigma, V);

  MATRIX3 P = material->PK1(U, Sigma, V);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX3 finiteDiffP;

    // for each of the degrees of the freedom
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
      {
        MATRIX3 Fnew = F;
        Fnew(x,y) += eps;

        // get the new psi
        MATRIX3 Unew,Vnew;
        VECTOR3 Snew;
        svd_rv(Fnew, Unew, Snew, Vnew);
        double psi = material->psi(Unew, Snew, Vnew);

        // store the finite difference
        finiteDiffP(x, y) = (psi - psi0) / eps;
      }

    MATRIX3 diff = P - finiteDiffP;
    REAL diffNorm = (fabs(diff.norm() / P.norm())) / 9.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " P: " << endl << P << endl;
      cout << " finite diff: " << endl << finiteDiffP << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// Make sure that the material reproduces linear elasticity under
// linearization
//////////////////////////////////////////////////////////////////////////////
bool validateAgainstLinear(const HOBAK::VOLUME::HYPERELASTIC* material,
                           const HOBAK::VOLUME::LINEAR* linear, 
                           const MATRIX3& F)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING linearization consistency for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  // build the linearized strain
  MATRIX3 I = MATRIX3::Identity();
  MATRIX3 eps = 0.5 * (F + F.transpose()) - I;
  VECTOR9 epsVector = flatten(eps);

  // get the linearized force
  MATRIX9 H = material->hessian(I);
  VECTOR9 linearized = H * epsVector;

  // get the force from the LINEAR material
  VECTOR9 groundTruth = flatten(linear->PK1(F));

  // are they the same?
  VECTOR9 diff = groundTruth - linearized;

  if (diff.norm() / 9.0 > 1e-7)
  {
    cout << " TEST FAILED! " << endl;
    cout << " Material PK1: " << linearized.transpose() << endl;
    cout << " Linear PK1:   " << groundTruth.transpose() << endl;
  }
  else
  {
    cout << " TEST PASSED. " << endl;
  }

  return diff.norm() / 9.0 <= 1e-7;
}

//////////////////////////////////////////////////////////////////////////////
// test out isotropic energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE( "Isotropic Energy tests", "[isotropic]" ) 
{
  using namespace HOBAK::VOLUME;
  MATRIX3 F = randomMatrix3();

  SECTION("Testing ARAP") 
  {
    ARAP arap(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&arap, F) == true);
    REQUIRE(convergenceTestHessian(&arap, F) == true);
    REQUIRE(testClampedHessian(&arap) == true);
  }

  SECTION("Testing SNH") 
  {
    SNH snh(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&snh, F) == true);
    REQUIRE(convergenceTestHessian(&snh, F) == true);
    REQUIRE(testClampedHessian(&snh) == true);
  }

  SECTION("Testing StVK") 
  {
    STVK stvk(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&stvk, F) == true);
    REQUIRE(convergenceTestHessian(&stvk, F) == true);
    REQUIRE(testClampedHessian(&stvk) == true);
  }
  
  SECTION("Testing SNH With Barrier term")
  {
    SNH_WITH_BARRIER snh(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&snh, F) == true);
    REQUIRE(convergenceTestHessian(&snh, F) == true);
    REQUIRE(testClampedHessian(&snh) == true);
  }
 
  // Bonet-Wood-style Neo-Hookean can't handle inversion 
  MATRIX3 Fpd = randomPositiveDefiniteMatrix3();
  SECTION("Testing Bonet-Wood Neo-Hookean")
  {
    NEO_HOOKEAN_BW nh(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&nh, Fpd) == true);
    REQUIRE(convergenceTestHessian(&nh, Fpd) == true);
    REQUIRE(testClampedHessian(&nh, false) == true);
  }
  
  SECTION("Testing Linear")
  {
    LINEAR linear(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&linear, F) == true);
    REQUIRE(convergenceTestHessian(&linear, F) == true);
    REQUIRE(testClampedHessian(&linear) == true);
  } 
}

//////////////////////////////////////////////////////////////////////////////
// test out anisotropic energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Anisotropic Energy tests", "[anisotropic]" ) 
{
  using namespace HOBAK::VOLUME;
  MATRIX3 F = randomMatrix3();
  VECTOR3 a = randomVector3().normalized();
  SECTION("Testing Anisotropic ARAP")
  {
    ANISOTROPIC_ARAP anisoARAP(1.0, a);
    REQUIRE(convergenceTestPK1(&anisoARAP, F) == true);
    REQUIRE(convergenceTestHessian(&anisoARAP, F) == true);
    REQUIRE(testClampedHessian(&anisoARAP) == true);
  }
  SECTION("Testing Anisotropic StVK")
  {
    ANISOTROPIC_STVK anisoStVK(1.0, a);
    REQUIRE(convergenceTestPK1(&anisoStVK, F) == true);
    REQUIRE(convergenceTestHessian(&anisoStVK, F) == true);
    REQUIRE(testClampedHessian(&anisoStVK) == true);
  }
  SECTION("Testing Anisotropic Dirichlet")
  {
    ANISOTROPIC_DIRICHLET anisoDirichlet(1.0, a);
    REQUIRE(convergenceTestPK1(&anisoDirichlet, F) == true);
    REQUIRE(convergenceTestHessian(&anisoDirichlet, F) == true);
    REQUIRE(testClampedHessian(&anisoDirichlet) == true);
  }
  SECTION("Testing Anisotropic Fung")
  {
    // Fung explodes easily for large deformation
    ANISOTROPIC_FUNG anisoFung(1.0, a);
    REQUIRE(convergenceTestPK1(&anisoFung, F) == true);
    REQUIRE(convergenceTestHessian(&anisoFung, F) == true);
    REQUIRE(testClampedHessian(&anisoFung) == true);
  }
}

//////////////////////////////////////////////////////////////////////////////
// test out tet mesh 
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Tet Mesh tests", "[tetMesh]" )
{
  using namespace HOBAK::VOLUME;

  // read in the test tet mesh
  vector<VECTOR3> vertices;
  vector<VECTOR4I> tets;
  string testFile("../data/cube_6.tobj");
  bool success = TET_MESH::readTobjFile(testFile, vertices, tets);
  if (!success)
    cout << " Test file " << testFile.c_str() << " not found! " << endl;
  REQUIRE(success);

  // run a convergence test to see if the tet mesh forces are working
  TET_MESH tetMesh(vertices, tets);
  SNH snh(1.0, 1.0);
  REQUIRE(convergenceTestTetMeshForces(snh, tetMesh));
  REQUIRE(convergenceTestTetMeshHessian(snh, tetMesh));
}

//////////////////////////////////////////////////////////////////////////////
// see if Lame reparametrization worked, and is now consistent with linear
// elasticity
//////////////////////////////////////////////////////////////////////////////
TEST_CASE( "Lame reparameterization", "[isotropic]" ) 
{
  using namespace HOBAK::VOLUME;

  const REAL mu = 2.0;
  const REAL lambda = 3.0;
  LINEAR linear(mu, lambda);
  MATRIX3 F = randomPositiveDefiniteMatrix3();

  SECTION("Testing SNH")
  {
    SNH snh(mu, lambda);
    REQUIRE(validateAgainstLinear(&snh, &linear, F));
  } 
  SECTION("Testing SNH_WITH_BARRIER")
  {
    SNH_WITH_BARRIER snh(mu, lambda);
    REQUIRE(validateAgainstLinear(&snh, &linear, F));
  }
  SECTION("Testing STVK")
  {
    STVK stvk(mu, lambda);
    REQUIRE(validateAgainstLinear(&stvk, &linear, F));
  }
}
