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
// do a convergence test on an x-based edge-edge collision Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestEdgeHessian(const HOBAK::VOLUME::EDGE_COLLISION* material, 
                                const vector<VECTOR3>& vertices,
                                const VECTOR2& a, const VECTOR2& b)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING vertex-based Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  VECTOR12 gradient0 = material->gradient(vertices,a,b);
  MATRIX12 H = material->hessian(vertices,a,b);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX12 finiteDiffH;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++, entry++)
      {
        vector<VECTOR3> verticesNew = vertices;
        verticesNew[j][i] += eps;

        // get the new psi
        VECTOR12 gradient = material->gradient(verticesNew,a,b);

        // store the finite difference
        finiteDiffH.col(entry) = (gradient - gradient0) / eps;
      }

    MATRIX12 diff = H - finiteDiffH;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " Hessian: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
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
// do a convergence test on an x-based edge-edge collision Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestEdgeHessianNegated(const HOBAK::VOLUME::EDGE_COLLISION* material, 
                                       const vector<VECTOR3>& vertices,
                                       const VECTOR2& a, const VECTOR2& b)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING vertex-based negated Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  VECTOR12 gradient0 = material->gradientNegated(vertices,a,b);
  MATRIX12 H = material->hessianNegated(vertices,a,b);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX12 finiteDiffH;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++, entry++)
      {
        vector<VECTOR3> verticesNew = vertices;
        verticesNew[j][i] += eps;

        // get the new psi
        VECTOR12 gradient = material->gradientNegated(verticesNew,a,b);

        // store the finite difference
        finiteDiffH.col(entry) = (gradient - gradient0) / eps;
      }

    MATRIX12 diff = H - finiteDiffH;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " Hessian: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
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
// do a convergence test on an x-based edge-edge collision gradient
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestEdgeGradient(const HOBAK::VOLUME::EDGE_COLLISION* material, 
                                 const vector<VECTOR3>& vertices,
                                 const VECTOR2& a, const VECTOR2& b)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING vertex-based gradients for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psi(vertices,a,b);
  const VECTOR12 gradient = material->gradient(vertices,a,b);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR12 finiteDiffGradient;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++, entry++)
      {
        vector<VECTOR3> verticesNew = vertices;
        verticesNew[j][i] += eps;

        // get the new psi
        double psi = material->psi(verticesNew,a,b);

        // store the finite difference
        finiteDiffGradient[entry] = (psi - psi0) / eps;
        
        // the force is the negative gradient, so should take into account here
        finiteDiffGradient[entry] *= 1.0;
      }

    VECTOR12 diff = gradient - finiteDiffGradient;
    REAL diffNorm = (fabs(diff.norm() / gradient.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " gradient: " << endl << gradient << endl;
      cout << " finite diff: " << endl << finiteDiffGradient << endl;
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
// do a convergence test on an x-based edge-edge collision gradient
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestEdgeGradientNegated(const HOBAK::VOLUME::EDGE_COLLISION* material, 
                                        const vector<VECTOR3>& vertices,
                                        const VECTOR2& a, const VECTOR2& b)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING vertex-based negated gradients for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psiNegated(vertices,a,b);
  const VECTOR12 gradient = material->gradientNegated(vertices,a,b);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR12 finiteDiffGradient;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++, entry++)
      {
        vector<VECTOR3> verticesNew = vertices;
        verticesNew[j][i] += eps;

        // get the new psi
        double psi = material->psiNegated(verticesNew,a,b);

        // store the finite difference
        finiteDiffGradient[entry] = (psi - psi0) / eps;
        
        // the force is the negative gradient, so should take into account here
        finiteDiffGradient[entry] *= 1.0;
      }

    VECTOR12 diff = gradient - finiteDiffGradient;
    REAL diffNorm = (fabs(diff.norm() / gradient.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " gradient: " << endl << gradient << endl;
      cout << " finite diff: " << endl << finiteDiffGradient << endl;
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
// do a convergence test on an x-based gradient
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCollisionGradient(const HOBAK::VOLUME::VERTEX_FACE_COLLISION* material, const vector<VECTOR3>& vertices)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING vertex-based gradients for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psi(vertices);
  VECTOR12 gradient = material->gradient(vertices);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR12 finiteDiffGradient;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++, entry++)
      {
        vector<VECTOR3> verticesNew = vertices;
        verticesNew[j][i] += eps;

        // get the new psi
        double psi = material->psi(verticesNew);

        // store the finite difference
        finiteDiffGradient[entry] = (psi - psi0) / eps;
      }

    VECTOR12 diff = gradient - finiteDiffGradient;
    REAL diffNorm = (fabs(diff.norm() / gradient.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " gradient: " << endl << gradient << endl;
      cout << " finite diff: " << endl << finiteDiffGradient << endl;
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
// do a convergence test on an x-based gradient
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCollisionGradient(const HOBAK::VOLUME::VERTEX_FACE_COLLISION* material, const VECTOR12& x)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING gradients for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psi(x);

  VECTOR12 gradient = material->gradient(x);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR12 finiteDiffGradient;

    // for each of the degrees of the freedom
    for (int i = 0; i < 12; i++)
    {
      VECTOR12 xNew = x;
      xNew[i] += eps;

      // get the new psi
      double psi = material->psi(xNew);

      // store the finite difference
      finiteDiffGradient[i] = (psi - psi0) / eps;
    }

    VECTOR12 diff = gradient - finiteDiffGradient;
    REAL diffNorm = (fabs(diff.norm() / gradient.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " gradient: " << endl << gradient << endl;
      cout << " finite diff: " << endl << finiteDiffGradient << endl;
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
// do a convergence test on a barycentric x-based gradient
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestPinnedCollisionGradient(const HOBAK::VOLUME::MCADAMS_COLLISION* material, 
                                            const vector<VECTOR3>& vertices)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING barycentric gradients for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  const VECTOR3 bary = HOBAK::getBarycentricCoordinates(vertices);

  const REAL psi0 = material->psi(vertices,bary);
  const VECTOR12 gradient = material->gradient(vertices,bary);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR12 finiteDiffGradient;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++, entry++)
      {
        vector<VECTOR3> vNew = vertices;
        vNew[i][j] += eps;

        // get the new psi
        double psi = material->psi(vNew, bary);

        // store the finite difference
        finiteDiffGradient[entry] = (psi - psi0) / eps;
      }

    VECTOR12 diff = gradient - finiteDiffGradient;
    REAL diffNorm = (fabs(diff.norm() / gradient.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " gradient: " << endl << gradient << endl;
      cout << " finite diff: " << endl << finiteDiffGradient << endl;
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
// do a convergence test on an x-based Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCollisionHessian(const HOBAK::VOLUME::VERTEX_FACE_COLLISION* material, const vector<VECTOR3>& vertices)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING vertex-based Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  VECTOR12 gradient0 = material->gradient(vertices);
  MATRIX12 H = material->hessian(vertices);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX12 finiteDiffH;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++, entry++)
      {
        vector<VECTOR3> verticesNew = vertices;
        verticesNew[j][i] += eps;

        // get the new psi
        VECTOR12 gradient = material->gradient(verticesNew);

        // store the finite difference
        finiteDiffH.col(entry) = (gradient - gradient0) / eps;
      }

    MATRIX12 diff = H - finiteDiffH;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " Hessian: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
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
// do a convergence test on an x-based Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestCollisionHessian(const HOBAK::VOLUME::VERTEX_FACE_COLLISION* material, const VECTOR12& x)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  VECTOR12 gradient0 = material->gradient(x);
  MATRIX12 H = material->hessian(x);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX12 finiteDiffH;

    // for each of the degrees of the freedom
    for (int i = 0; i < 12; i++)
    {
      VECTOR12 xNew = x;
      xNew[i] += eps;

      // get the new psi
      VECTOR12 gradient = material->gradient(xNew);

      // store the finite difference
      finiteDiffH.col(i) = (gradient - gradient0) / eps;
    }

    MATRIX12 diff = H - finiteDiffH;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " Hessian: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
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
// do a convergence test on an x-based Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestPinnedCollisionHessian(const HOBAK::VOLUME::MCADAMS_COLLISION* material, 
                                           const vector<VECTOR3>& vertices)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING barycentric Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  const VECTOR3 bary = HOBAK::getBarycentricCoordinates(vertices);

  const VECTOR12 gradient0 = material->gradient(vertices, bary);
  const MATRIX12 H = material->hessian(vertices, bary);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX12 finiteDiffH;
    finiteDiffH.setZero();

    // for each of the degrees of the freedom
    int entry = 0;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++, entry++)
      {
        vector<VECTOR3> vNew = vertices;
        vNew[i][j] += eps;
        //VECTOR12 xNew = x0;
        //xNew[entry] += eps;

        // get the new psi
        const VECTOR12 gradient = material->gradient(vNew,bary);
        //const VECTOR12 gradient = material->gradient(xNew,bary);

        // store the finite difference
        finiteDiffH.col(entry) = (gradient - gradient0) / eps;
      }

    MATRIX12 diff = H - finiteDiffH;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " Hessian: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
      cout << " diff: " << endl << diff << endl;

      cout << "symmetric: " << endl << finiteDiffH - finiteDiffH.transpose() << endl;
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
// do a convergence test on an x-based Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergencePsiTestPinnedCollisionHessian(const HOBAK::VOLUME::MCADAMS_COLLISION* material, 
                                              const vector<VECTOR3>& vertices)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING barycentric Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  const VECTOR3 bary = HOBAK::getBarycentricCoordinates(vertices);

  const REAL psi0 = material->psi(vertices, bary);
  const MATRIX12 H = material->hessian(vertices, bary);

  VECTOR12 x0;
  x0[0] = vertices[0][0];
  x0[1] = vertices[0][1];
  x0[2] = vertices[0][2];
  x0[3] = vertices[1][0];
  x0[4] = vertices[1][1];
  x0[5] = vertices[1][2];
  x0[6] = vertices[2][0];
  x0[7] = vertices[2][1];
  x0[8] = vertices[2][2];
  x0[9] = vertices[3][0];
  x0[10] = vertices[3][1];
  x0[11] = vertices[3][2];

  //double eps = 1e-4;
  double eps = 1e-1;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-7)
  {
    MATRIX12 finiteDiffH;
    finiteDiffH.setZero();
    vector<VECTOR3> vertices1(4);

    for (int i = 0; i < 12; i++)
    {
      VECTOR12 x1 = x0;
      x1 = x0;
      x1[i] += eps;
      int entry = 0;
      for (int y = 0; y < 4; y++)
        for (int x = 0; x < 3; x++,entry++)
          vertices1[y][x] = x1[entry];
      REAL p1 = material->psi(vertices1, bary);

      x1 = x0;
      x1[i] -= eps;
      entry = 0;
      for (int y = 0; y < 4; y++)
        for (int x = 0; x < 3; x++,entry++)
          vertices1[y][x] = x1[entry];
      REAL m1 = material->psi(vertices1, bary);

      finiteDiffH(i,i) = (-2.0 * psi0 + p1 + m1) / (eps * eps);
    }

    MATRIX12 diff = H - finiteDiffH;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " Hessian: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
      cout << " diff: " << endl << diff << endl;

      cout << "symmetric: " << endl << finiteDiffH - finiteDiffH.transpose() << endl;
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
// test out collision energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Collision Energy tests", "[collision]" )
{
  using namespace HOBAK::VOLUME;
  VECTOR12 x = randomVector12();
  SECTION("Testing vertex-face collisions, flat vector")
  {
    VERTEX_FACE_COLLISION vertexFace(1.0, 0.01);
    REQUIRE(convergenceTestCollisionGradient(&vertexFace, x) == true);
    REQUIRE(convergenceTestCollisionHessian(&vertexFace, x) == true);
  }
  SECTION("Testing vertex-face collisions, vertices")
  {
    vector<VECTOR3> vertices(4);
    for (int i = 0; i < 4; i++)
      vertices[i] = randomVector3();
    VERTEX_FACE_COLLISION vertexFace(1.0, 0.01);
    REQUIRE(convergenceTestCollisionGradient(&vertexFace, vertices) == true);
    REQUIRE(convergenceTestCollisionHessian(&vertexFace, vertices) == true);
    
    // make sure the clamping didn't accidentally zero out the Hessian
    MATRIX12 H = vertexFace.clampedHessian(vertices);
    REQUIRE(H.norm() > 1e-8);
  }

  vector<VECTOR3> vertices(4);
  for (int x = 0; x < 4; x++)
    vertices[x] = randomVector3();
  SECTION("Testing Barycentric computation")
  {
    cout << "========================================================== " << endl;
    cout << " VERIFYING Barycentric projection for collisions" << endl;
    cout << "========================================================== " << endl;
    VECTOR3 bary = HOBAK::getBarycentricCoordinates(vertices);
    VECTOR3 guess = bary[0] * vertices[1] + bary[1] * vertices[2] + bary[2] * vertices[3];

    // get the projected position
    const VECTOR3 v0 = vertices[1];
    const VECTOR3 v1 = vertices[2];
    const VECTOR3 v2 = vertices[3];
    const VECTOR3 e1 = v1 - v0;
    const VECTOR3 e2 = v2 - v0;
    const VECTOR3 n = e1.cross(e2);
    const VECTOR3 nHat = n / n.norm();
    const VECTOR3 v = vertices[0] - (nHat.dot(vertices[0] - v0)) * nHat;
    
    // make sure we get the projected point back
    REQUIRE((guess - v).norm() < 1e-8);
  }

  SECTION("Testing McAdams barycentric collisions")
  {
    vector<VECTOR3> vertices(4);
    for (int x = 0; x < 4; x++)
      vertices[x] = randomVector3();
    MCADAMS_COLLISION mcAdams(1.0, 0.01);
    
    REQUIRE(convergenceTestPinnedCollisionGradient(&mcAdams, vertices) == true);
    REQUIRE(convergenceTestPinnedCollisionHessian(&mcAdams, vertices) == true);

    // make sure the clamping didn't accidentally zero out the Hessian
    const VECTOR3 bary = HOBAK::getBarycentricCoordinates(vertices);
    MATRIX12 H = mcAdams.clampedHessian(vertices,bary);
    REQUIRE(H.norm() > 1e-8);
  }
  SECTION("Testing sqrt-based barycentric collisions")
  {
    vector<VECTOR3> vertices(4);
    for (int x = 0; x < 4; x++)
      vertices[x] = randomVector3();
    VERTEX_FACE_SQRT_COLLISION sqrtEnergy(1.0, 0.01);
    
    REQUIRE(convergenceTestPinnedCollisionGradient(&sqrtEnergy, vertices) == true);
    REQUIRE(convergenceTestPinnedCollisionHessian(&sqrtEnergy, vertices) == true);

    // make sure the clamping didn't accidentally zero out the Hessian
    //const VECTOR3 bary = HOBAK::getBarycentricCoordinates(vertices);
    //MATRIX12 H = sqrtEnergy.clampedHessian(vertices,bary);
    MATRIX12 H = sqrtEnergy.clampedHessian(vertices);
    REQUIRE(H.norm() > 1e-8);
  }
  SECTION("Testing sqrt-based reversed barycentric collisions")
  {
    vector<VECTOR3> vertices(4);
    vertices[0] = VECTOR3(0.5, -1.0, sqrt(1.0 / 12.0));
    vertices[1] = VECTOR3(1,0,0);
    vertices[2] = VECTOR3(0,0,0);
    vertices[3] = VECTOR3(0.5,0,sqrt(3.0) / 2.0);
    VERTEX_FACE_SQRT_COLLISION sqrtEnergy(1.0, 0.01);
    
    REQUIRE(convergenceTestPinnedCollisionGradient(&sqrtEnergy, vertices) == true);
    REQUIRE(convergenceTestPinnedCollisionHessian(&sqrtEnergy, vertices) == true);

    // make sure the clamping didn't accidentally zero out the Hessian
    //const VECTOR3 bary = HOBAK::getBarycentricCoordinates(vertices);
    //MATRIX12 H = sqrtEnergy.clampedHessian(vertices,bary);
    MATRIX12 H = sqrtEnergy.clampedHessian(vertices);
    REQUIRE(H.norm() > 1e-8);
  }

  SECTION("Testing barycentric edge-edge collisions energy")
  {
    vector<VECTOR3> v(4);
    for (int x = 0; x < 4; x++)
      v[x] = randomVector3();

    VECTOR2 a = randomBarycentric();
    VECTOR2 b = randomBarycentric();

    EDGE_SQRT_COLLISION edge(1.0, 0.01);
    REQUIRE(convergenceTestEdgeGradient(&edge, v, a, b) == true);
    REQUIRE(convergenceTestEdgeHessian(&edge, v, a, b) == true);

    // make sure the clamping didn't accidentally zero out the Hessian
    MATRIX12 H = edge.EDGE_COLLISION::clampedHessian(v,a,b);
    REQUIRE(H.norm() > 1e-8);
  }

  SECTION("Testing edge-edge collisions energy")
  {
    vector<VECTOR3> v(4);
    for (int x = 0; x < 4; x++)
      v[x] = randomVector3();

    VECTOR2 a = randomBarycentric();
    VECTOR2 b = randomBarycentric();

    EDGE_COLLISION edge(1.0, 0.01);
    REQUIRE(convergenceTestEdgeGradient(&edge, v, a, b) == true);
    REQUIRE(convergenceTestEdgeHessian(&edge, v, a, b) == true);
    
    // make sure the clamping didn't accidentally zero out the Hessian
    REQUIRE(edge.clampedHessian(v,a,b).norm() > 1e-8);
  }
  SECTION("Testing negated edge-edge collisions energy")
  {
    vector<VECTOR3> v(4);
    for (int x = 0; x < 4; x++)
      v[x] = randomVector3();

    VECTOR2 a = randomBarycentric();
    VECTOR2 b = randomBarycentric();

    EDGE_COLLISION edge(1.0, 0.01);
    REQUIRE(convergenceTestEdgeGradientNegated(&edge, v, a, b) == true);
    REQUIRE(convergenceTestEdgeHessianNegated(&edge, v, a, b) == true);
    
    // make sure the clamping didn't accidentally zero out the Hessian
    REQUIRE(edge.EDGE_COLLISION::clampedHessianNegated(v,a,b).norm() > 1e-8);
  }

  SECTION("Testing barycentric edge-edge collisions energy")
  {
    vector<VECTOR3> v(4);
    for (int x = 0; x < 4; x++)
      v[x] = randomVector3();

    VECTOR2 a = randomBarycentric();
    VECTOR2 b = randomBarycentric();

    EDGE_SQRT_COLLISION edge(1.0, 0.01);
    REQUIRE(convergenceTestEdgeGradient(&edge, v, a, b) == true);
    REQUIRE(convergenceTestEdgeHessian(&edge, v, a, b) == true);
    
    // make sure the clamping didn't accidentally zero out the Hessian
    REQUIRE(edge.EDGE_COLLISION::clampedHessian(v,a,b).norm() > 1e-8);
  }
  
  SECTION("Testing barycentric, negated edge-edge collisions energy")
  {
    vector<VECTOR3> v(4);
    for (int x = 0; x < 4; x++)
      v[x] = randomVector3();

    VECTOR2 a = randomBarycentric();
    VECTOR2 b = randomBarycentric();

    EDGE_SQRT_COLLISION edge(1.0, 0.01);
    REQUIRE(convergenceTestEdgeGradientNegated(&edge, v, a, b) == true);
    REQUIRE(convergenceTestEdgeHessianNegated(&edge, v, a, b) == true);
    
    // make sure the clamping didn't accidentally zero out the Hessian
    REQUIRE(edge.EDGE_COLLISION::clampedHessianNegated(v,a,b).norm() > 1e-8);
  }
}
