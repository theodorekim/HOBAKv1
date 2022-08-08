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
// perform a stretch test with the quasistatic integrator
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Quasistatic stretch test", "[quasistatic stretch]" )
{
  using namespace HOBAK::VOLUME;
  using namespace HOBAK::TIMESTEPPER;

  // read in the tet mesh file
  string testFile("../data/cube_20.tobj");
  vector<VECTOR3> restVertices;
  vector<VECTOR4I> tets;
  bool success = TET_MESH::readTobjFile(testFile, restVertices, tets);
  if (!success)
    cout << " FILE LOAD FAILED. You should be running from the ./bin directory!" << endl;

  // this probably means you didn't run from the right directory
  REQUIRE(success);

  // build the tet mesh object
  TET_MESH tetMesh(restVertices, tets);
  ARAP hyperelastic(1.0, 1.0);

  // build the time integrator
  QUASISTATIC quasistaticSolver(tetMesh, hyperelastic);

  // attach it to a cube in the scene
  VECTOR3 center0(0.56,0.48, -0.24);
  CUBE cube0(center0, 1.0);
  quasistaticSolver.attachKinematicSurfaceConstraints(&cube0);

  // peek at the bounding box
  VECTOR3 mins, maxs;
  tetMesh.getBoundingBox(mins, maxs);

  // translate the cube's left wall to z = 0
  VECTOR3 translate = (mins + maxs) * 0.5;
  translate[2] = mins[2];

  // apply a stretch to the tet mesh
  for (unsigned int x = 0; x < tetMesh.vertices().size(); x++)
  {
    tetMesh.vertices()[x] -= translate;
    tetMesh.vertices()[x][2] *= 10.0;
    tetMesh.vertices()[x] += translate;
  }
  
  // run the quasistatic solver, it should snap everything back to the rest pose,
  // since there are no fores applied, and the boundary conditions correspond
  // to the rest pose
  quasistaticSolver.solve(true, 1.0);

  // see if it snapped back to rest pose
  REAL totalDiff = 0.0;
  for (unsigned int x = 0; x < tetMesh.vertices().size(); x++)
  {
    VECTOR3 diff = tetMesh.vertices()[x] - restVertices[x];
    totalDiff += diff.norm();
  }

  // it should be really small, as it should be snapped back to the rest pose
  REQUIRE(totalDiff <= 1e-8);
}

//////////////////////////////////////////////////////////////////////////////
// see if Baraff-Witkin constraints allow for clean separation with no
// sticking in quasistatics
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Quasistatic separation test", "[quasistatic separation]" )
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING quasistatic constraint separation" << endl;
  cout << "=============================================================== " << endl;
  using namespace HOBAK::VOLUME;
  using namespace HOBAK::TIMESTEPPER;

  // read in the tet mesh file
  string testFile("../data/cube_10.tobj");
  vector<VECTOR3> vertices;
  vector<VECTOR4I> tets;
  bool success = TET_MESH::readTobjFile(testFile, vertices, tets);

  // this probably means you didn't run from the right directory
  REQUIRE(success);

  ///////////////////////////////////////////////////////////////
  // BUILDING THE SCENE
  ///////////////////////////////////////////////////////////////

  // build the tet mesh object
  TET_MESH tetMesh(vertices, tets);
  ARAP hyperelastic(1.0, 1.0);

  // build the time integrator
  QUASISTATIC quasistaticSolver(tetMesh, hyperelastic);

  // build the kinematic cubes
  const VECTOR3 center0(0.56,0.48, -0.25);
  CUBE cube0(center0, 1.0);
  const VECTOR3 center1(0.56,0.48,1.23);
  CUBE cube1(center1, 1.0);

  // attach it to the rest of the scene
  quasistaticSolver.attachKinematicSurfaceConstraints(&cube1);
  quasistaticSolver.addKinematicCollisionObject(&cube0);

  ///////////////////////////////////////////////////////////////
  // RUNNING THE TEST
  ///////////////////////////////////////////////////////////////
  
  // translate the cube back and forth. When it's done it
  // should have separated from the tet mesh again.
  for (int x = 0; x < 4; x ++)
  {
    cube0.translation()[2] += 0.01;
    quasistaticSolver.solve(false, 1.0);
  }
  for (int x = 0; x < 5; x ++)
  {
    cube0.translation()[2] -= 0.01;
    quasistaticSolver.solve(false, 1.0);
  }

  // the mesh should have cleanly separated from the colliding cubes,
  // and should just retain its original positions
  const VECTOR displacement = tetMesh.getDisplacement();
  const REAL displacementNorm = displacement.norm();

  // mesh should go back to rest pose
  cout << " Displacement norm is: " << displacementNorm << endl;
  const REAL threshold = 1e-6;
  if (displacementNorm > threshold)
    cout << " TEST FAILED. Cube did not separate cleanly" << endl;
  REQUIRE(displacementNorm < threshold);

  // there shouldn't be any plane constraints left
  const int totalConstraints = quasistaticSolver.totalPlaneConstraints();
  cout << " Total plane constraints: " << totalConstraints << endl;
  if (totalConstraints > 0)
    cout << " TEST FAILED. All plane constraints should have been deleted. " << endl;
  REQUIRE(totalConstraints == 0);

  cout << " TEST PASSED. " << endl;
}

//////////////////////////////////////////////////////////////////////////////
// material-agnostic test to see if pop-through is fixed
//////////////////////////////////////////////////////////////////////////////
bool runPopThroughTest(TET_MESH_FASTER& tetMesh,
                       VOLUME::HYPERELASTIC* hyperelastic)
{
  ///////////////////////////////////////////////////////////////
  // BUILDING THE SCENE
  ///////////////////////////////////////////////////////////////
  using namespace HOBAK::VOLUME;
  using namespace HOBAK::TIMESTEPPER;

  // build the time integrator
  QUASISTATIC quasistaticSolver(tetMesh, *hyperelastic);

  // build the kinematic cubes
  VECTOR3 center0(0.56,0.48, -0.25);
  CUBE cube0(center0, 1.0);
  VECTOR3 center1(0.56,0.48,1.23);
  CUBE cube1(center1, 1.0);

  // attach it to the rest of the scene
  quasistaticSolver.attachKinematicSurfaceConstraints(&cube1);
  quasistaticSolver.addKinematicCollisionObject(&cube0);

  ///////////////////////////////////////////////////////////////
  // RUNNING THE TEST
  ///////////////////////////////////////////////////////////////
  
  // translate the cube back and forth. When it's done it
  // should have separated from the tet mesh again.
  for (int x = 0; x < 4; x ++)
  {
    cube0.translation()[2] += 0.01;
    quasistaticSolver.solve(false, 1.0);
  }
  for (int x = 0; x < 4; x ++)
  {
    cube0.translation()[2] -= 0.01;
    quasistaticSolver.solve(false, 1.0);
  }

  // now compress again, see if we get pop-through
  for (int x = 0; x < 2; x ++)
  {
    cube0.translation()[2] += 0.01;
    quasistaticSolver.solve(false, 1.0);
  }
  cube0.translation()[2] -= 0.01;
  quasistaticSolver.solve(false, 1.0);

  // the z-component of all the plane-constrained 
  // nodes should be at 0.26
  const vector<PLANE_CONSTRAINT>& constraints = quasistaticSolver.planeConstraints();

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const int vertexID = constraints[x].vertexID;
    const vector<VECTOR3>& vertices = tetMesh.vertices();
    const VECTOR3& vertex = vertices[vertexID];

    const REAL diff = fabs(vertex[2] - 0.26);
    if (diff > 1e-8)
    {
      cout << " Constraint " << x << " is WAY OFF: " << diff << endl;
      cout << " TEST FAILED. " << endl;
      return false;
    }
  }
  cout << " TEST PASSED. " << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if Baraff-Witkin constraints can be preserved without pop-through
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Quasistatic pop-through test", "[quasistatic pop-through]" )
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING quasistatic pop-through fixed" << endl;
  cout << "=============================================================== " << endl;
  using namespace HOBAK::VOLUME;
  using namespace HOBAK::TIMESTEPPER;

  // read in the tet mesh file
  string testFile("../data/cube_10.tobj");
  vector<VECTOR3> inputVertices;
  vector<VECTOR4I> inputTets;
  bool success = TET_MESH::readTobjFile(testFile, inputVertices, inputTets);
  TET_MESH_FASTER tetMesh(inputVertices, inputTets);

  // this probably means you didn't run from the right directory
  REQUIRE(success);

  // try the test using ARAP
  ARAP* arap= new ARAP(1.0, 1.0);
  cout << " Testing ARAP ... " << endl;
  REQUIRE(runPopThroughTest(tetMesh, arap));
  tetMesh.vertices() = tetMesh.restVertices();
  delete arap;
  
  SNH* snh = new SNH(1.0, 10.0);
  cout << " Testing SNH ... " << endl;
  REQUIRE(runPopThroughTest(tetMesh, snh));
  tetMesh.vertices() = tetMesh.restVertices();
  delete snh;
}

//////////////////////////////////////////////////////////////////////////////
// see if velocity filtering of kinematic constraints is working properly
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Newmark kinematic velocity filtering", "[kinematic filtering]" )
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Newmark velocity filtering works" << endl;
  cout << "=============================================================== " << endl;
  using namespace HOBAK::VOLUME;
  using namespace HOBAK::TIMESTEPPER;

  // read in the tet mesh file
  string testFile("../data/cube_6.tobj");
  vector<VECTOR3> inputVertices;
  vector<VECTOR4I> inputTets;
  bool success = TET_MESH::readTobjFile(testFile, inputVertices, inputTets);

  // this probably means you didn't run from the right directory
  REQUIRE(success);

  // build kinematic cubes
  VECTOR3 center0(0.56,0.48, -0.25);
  CUBE cube0(center0, 1.0);
  VECTOR3 center1(0.56,0.48,1.23);
  CUBE cube1(center1, 1.0);

  // build the tet mesh object
  TET_MESH tetMesh(inputVertices, inputTets);
  ARAP arap(1.0, 1.0);

  // build the time integrator
  NEWMARK newmarkSolver(tetMesh, arap);

  // attach it to the rest of the scene
  newmarkSolver.attachKinematicSurfaceConstraints(&cube1);
  newmarkSolver.addKinematicCollisionObject(&cube0);
  newmarkSolver.setRayeligh(0.01, 0.01);

  REQUIRE(fabs(newmarkSolver.rayleighAlpha() - 0.01) < 1e-8);
  REQUIRE(fabs(newmarkSolver.rayleighBeta() - 0.01) < 1e-8);
  REQUIRE(fabs(newmarkSolver.dt() - (1.0 / 60.0)) < 1e-8);

  newmarkSolver.maxNewtonIterations() = 10;
  newmarkSolver.vertexFaceSelfCollisionsOn() = false;
  newmarkSolver.edgeEdgeSelfCollisionsOn() = false;

  bool verbose = false;

  // crush it
  cube0.translation()[2] += 0.01;
  newmarkSolver.solveRayleighDamped(verbose);
  newmarkSolver.solveRayleighDamped(verbose);

  // back off
  for (int x = 0; x < 4; x++)
  {
    cube0.translation()[2] -= 0.01;
    newmarkSolver.solveRayleighDamped(verbose);
  }

  // crush it more
  for (int x = 0; x < 6; x++)
  {
    cube0.translation()[2] += 0.01;
    newmarkSolver.solveRayleighDamped(verbose);
  }

  // let it settle
  for (int x = 0; x < 5; x++)
    newmarkSolver.solveRayleighDamped(verbose);

  // back off again -- does it stick?
  for (int x = 0; x < 3; x++)
  {
    cube0.translation()[2] -= 0.01;
    newmarkSolver.solveRayleighDamped(verbose);
  }

  if (newmarkSolver.velocity().norm() < 2.0)
  {
    cout << " Newmark velocity filter test PASSED " << endl;
  }
  else
  {
    cout << " Newmark velocity filter test FAILED" << endl;
    cout << " Velocity norm is very big, suggesting that filtering did not properly";
    cout << " occur, and simulation is exploding!!!" << endl << endl;
    cout << " Velocity norm is: " << newmarkSolver.velocity().norm() << endl;
    cout << " (should be less than 2.0, around 1.19826, NOT 277.711)" << endl;
  }
  
  REQUIRE(newmarkSolver.velocity().norm() < 2.0);
}

//////////////////////////////////////////////////////////////////////////////
// see if objects will separate after a collision
//////////////////////////////////////////////////////////////////////////////
bool postCollisionSeparationTest(HOBAK::TIMESTEPPER::TIMESTEPPER* stepper)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING " << stepper->name().c_str() << " collision separation works" << endl;
  cout << "=============================================================== " << endl;
  using namespace HOBAK::VOLUME;
  using namespace HOBAK::TIMESTEPPER;

  // the floor
  VECTOR3 center0(0.0,-5.0, 0.0);
  CUBE cube0(center0, 10.0);
  stepper->addKinematicCollisionObject(&cube0);

  const TET_MESH& tetMesh = stepper->tetMesh();

  // give it a big ballistic velocity
  VECTOR velocity(tetMesh.DOFs());
  for (int x = 0; x < velocity.size() / 3; x++)
    velocity[3 * x + 1] = -3.0;

  stepper->velocity() = velocity;
  
  for (int x = 0; x < 42; x++)
  {
    stepper->externalForces().setZero();
    stepper->addGravity(VECTOR3(0.0, -1.0, 0.0));
    stepper->solve(false);
  }
  
  const vector<VECTOR3>& vertices = tetMesh.vertices();
  const vector<int>& surfaceVertices = tetMesh.surfaceVertices();

  // if the cube jumped off the ground after the collision, all the 
  // y-coordinates of the surface vertices should be greater than zero
  bool surfaceIsAboveZero = true;
  for (unsigned int x = 0; x < surfaceVertices.size(); x++)
  {
    if (vertices[surfaceVertices[x]][1] <= 1e-4)
      surfaceIsAboveZero = false;
  }
  
  if (surfaceIsAboveZero)
  {
    cout << " " << stepper->name().c_str() << " collision separation test PASSED " << endl;
  }
  else
  {
    cout << " " << stepper->name().c_str() << " collision separation test FAILED" << endl;
    cout << " Did the tet mesh stick to the ground? Did you mess with the velocity " << endl;
    cout << " separation conditions in findNewSurfaceConstraints or " << endl;
    cout << " findSeparatingSurfaceConstraints?" << endl;
  }

  return surfaceIsAboveZero;
}

//////////////////////////////////////////////////////////////////////////////
// see if objects will separate after a collision
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Collision separation tests", "[Collision separation]" )
{
  using namespace HOBAK::VOLUME;
  //using namespace HOBAK::TIMESTEPPER;

  // read in the tet mesh file
  string testFile("../data/cube_6.tobj");
  vector<VECTOR3> inputVertices;
  vector<VECTOR4I> inputTets;
  bool success = TET_MESH::readTobjFile(testFile, inputVertices, inputTets);

  REQUIRE(success);

  // build the tet mesh object
  TET_MESH tetMesh(inputVertices, inputTets);
  SNH hyperelastic(1.0, 10.0);

  // build the time integrator
  TIMESTEPPER::BACKWARD_EULER_VELOCITY backwardEuler(tetMesh, hyperelastic);
  backwardEuler.vertexFaceSelfCollisionsOn() = false;
  backwardEuler.edgeEdgeSelfCollisionsOn() = false;
  REQUIRE(postCollisionSeparationTest(&backwardEuler));
  
  TIMESTEPPER::NEWMARK newmarkSolver(tetMesh, hyperelastic);
  newmarkSolver.vertexFaceSelfCollisionsOn() = false;
  newmarkSolver.edgeEdgeSelfCollisionsOn() = false;
  REQUIRE(postCollisionSeparationTest(&newmarkSolver));
}
