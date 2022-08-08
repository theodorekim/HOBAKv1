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
#include "BACKWARD_EULER_VELOCITY.h"
#include "TIMER.h"
#include "Damping/Volume/GREEN_DAMPING.h"
#include "Geometry/TET_MESH_FASTER.h"
#include <float.h>
#include <iostream>

using namespace std;

namespace HOBAK {
namespace TIMESTEPPER {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
BACKWARD_EULER_VELOCITY::BACKWARD_EULER_VELOCITY(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic) :
  TIMESTEPPER(tetMesh, hyperelastic)
{
  _rayleighAlpha = 0.01;
  _rayleighBeta = 0.01;

  _dt = 1.0 / 60.0;

  _time = 0;
  _currentTimestep = 0;

  _name = string("BACKWARD_EULER_VELOCITY");
  
  //_damping = new VOLUME::GREEN_DAMPING(0.001);
  //_damping = NULL;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the displacement targets the the Baraff-Witkin-style constraints
// are trying to hit. Assumes that buildConstraintMatrix() has already been called
///////////////////////////////////////////////////////////////////////////////////////////////////////
void BACKWARD_EULER_VELOCITY::updateConstraintTargets()
{
  TIMER functionTimer(__FUNCTION__);
  _constraintTargets.setZero();
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    // should ignore if we've tagged it for deletion
    if (_planeConstraints[x].isSeparating)
      continue;

    // retrieve collision information
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const int vertexID = constraint.vertexID;
    const int index = 3 * vertexID;

    // compute the target displacement
    const VECTOR3& vertex = _tetMesh.vertices()[vertexID];
    const VECTOR3& localClosestPoint = _planeConstraints[x].localClosestPoint;
    const VECTOR3& closestPoint = shape->localVertexToWorld(localClosestPoint);

    const VECTOR3& xDelta = closestPoint - vertex;

    const VECTOR3 vertexVelocity = velocity(vertexID);
    VECTOR3 vDelta = (xDelta) / _dt - vertexVelocity;

    for (int i = 0; i < 3; i++)
      _constraintTargets[index + i] = vDelta[i];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void printEntry(const VECTOR& v, const int i, const string& varname)
{
  VECTOR3 v3;
  v3[0] = v[3 * i];
  v3[1] = v[3 * i + 1];
  v3[2] = v[3 * i + 2];
  cout << varname.c_str() << ": " << v3.transpose() << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool BACKWARD_EULER_VELOCITY::solve(const bool verbose)
{
  if (_damping != NULL)
    return solveEnergyDamped(verbose);
  return solveRayleighDamped(verbose);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// [BW98] refers to "Large Steps in Cloth Simulation", SIGGRAPH 1998
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool BACKWARD_EULER_VELOCITY::solveRayleighDamped(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " BACKWARD_EULER_VELOCITY RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // get the damping matrix
  SPARSE_MATRIX C = buildRayleighDampingMatrix();

  // only caching this for visualization purposes
  _positionOld = _position;

  // should need to call once, but then preserved throughout
  applyKinematicConstraints();

  // store the filtered b for later
  VECTOR unfiltered;

  // build new constraints and see if we should break any
  findNewSurfaceConstraints(verbose);
  buildConstraintMatrix();

  _tetMesh.setDisplacement(_position);
  _tetMesh.computeFs();
  _tetMesh.computeSVDs();

  // do collision detection, including spatial data structure updates 
  computeCollisionDetection();

  // "z is a vector of the desired values for the constrained variables",
  // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
  // _constraintTargets did not project off the null directions
  updateConstraintTargets();
  VECTOR z =_IminusS * _constraintTargets;

  // get the internal forces
  VECTOR R = _tetMesh.computeHyperelasticForces(_hyperelastic);

  // get the reduced stiffness matrix
  SPARSE_MATRIX K = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);

  // compute collision forces and stiffnesses
  computeCollisionResponse(R,K,C);

  // assemble RHS from Eqn. 18 in [BW98], but
  TIMER systemTimer("Forming linear system");
  _b = _dt * (R + _dt * K * _velocity + _externalForces);

  // assemble system matrix A, LHS from Eqn. 18 in [BW98], but
  // we negate the K term because of how hyperelastic energies work
  _A = _M - _dt * C - _dt * _dt * K;

  // from [TJM15], this is c = b - Az (page 8, top of column 2)
  VECTOR c = _b - _A * z;
  systemTimer.stop();
 
  // just for clarity, go ahead and form the RHS and LHS explicitly
  TIMER projectionTimer("PPCG projection");
  VECTOR rhs = _S * c;
  SPARSE_MATRIX LHS = _S * _A * _S + _IminusS;
  projectionTimer.stop();

#if 1
  TIMER pcgTimer("PCG Solve");
  _cgSolver.compute(LHS);
  VECTOR y = _cgSolver.solve(rhs);
  pcgTimer.stop();

  if (verbose)
    printf("  PCG iters: %3i err: %6.4e \n", (int)_cgSolver.iterations(), (float)_cgSolver.error());
#else
  TIMER choleskyTimer("Cholesky Solve");
  Eigen::SimplicialLDLT<SPARSE_MATRIX> solver;
  solver.compute(LHS);
  VECTOR y = solver.solve(rhs);
  choleskyTimer.stop();
#endif

  // aliasing _solution to \Delta v just to make clear what we're doing here
  VECTOR& vDelta = _solution;
  vDelta = y + z;

  // update velocity
  _velocity = _velocity + vDelta;
  _position = _position + _dt * _velocity;

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;

  const bool constraintsChanged = findSeparatingSurfaceConstraints(_b);

  // ONLY change the constraints here. Trying to do it inside the Newton loop
  // is too oscillatory, and the direction that the forces point in can oscilliate
  // if something was separating
  //
  // Delete constraints this AFTER the Newmark update. Otherwise, we don't 
  // know how to filter off the position changes due to purely the constraints, 
  // and we see huge accelerations
  if (constraintsChanged)
  {
    deleteSurfaceConstraints(verbose);
    updateSurfaceConstraints();
    buildConstraintMatrix();
    updateConstraintTargets();
  }
  // otherwise, update the targets, but the constraint matrix should not have changed.
  else
  {
    updateSurfaceConstraints();
    updateConstraintTargets();
  }
  
  // update node positions
  _tetMesh.setDisplacement(_position);

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;


  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// [BW98] refers to "Large Steps in Cloth Simulation", SIGGRAPH 1998
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool BACKWARD_EULER_VELOCITY::solveEnergyDamped(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " BACKWARD_EULER_VELOCITY ENERGY-DAMPED SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // only caching this for visualization purposes
  _positionOld = _position;

  // should need to call once, but then preserved throughout
  applyKinematicConstraints();

  // store the filtered b for later
  VECTOR unfiltered;

  // do collision detection, including spatial data structure updates 
  computeCollisionDetection();

  // build new constraints and see if we should break any
  findNewSurfaceConstraints(verbose);
  buildConstraintMatrix();

  _tetMesh.setDisplacement(_position);
  _tetMesh.computeFs();
  _tetMesh.computeSVDs();

  // "z is a vector of the desired values for the constrained variables",
  // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
  // _constraintTargets did not project off the null directions
  updateConstraintTargets();
  VECTOR z =_IminusS * _constraintTargets;

  // get the internal forces
  VECTOR R = _tetMesh.computeInternalForces(_hyperelastic, *_damping);

  // get the stiffness matrix
  SPARSE_MATRIX K = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);
  SPARSE_MATRIX C = _tetMesh.computeDampingHessian(*_damping);

  // compute collision forces and stiffnesses
  computeCollisionResponse(R,K,C);

  // assemble RHS from Eqn. 18 in [BW98], but
  TIMER systemTimer("Forming linear system");
  _b = _dt * (R + _dt * K * _velocity + _externalForces);

  // assemble system matrix A, LHS from Eqn. 18 in [BW98], but
  // we negate the K term because of how hyperelastic energies work
  _A = _M - _dt * C - _dt * _dt * K;

  // from [TJM15], this is c = b - Az (page 8, top of column 2)
  VECTOR c = _b - _A * z;
  systemTimer.stop();
 
  // just for clarity, go ahead and form the RHS and LHS explicitly
  TIMER projectionTimer("PPCG projection");
  VECTOR rhs = _S * c;
  SPARSE_MATRIX LHS = _S * _A * _S + _IminusS;
  projectionTimer.stop();

#if 1
  TIMER pcgTimer("PCG Solve");
  _cgSolver.compute(LHS);
  VECTOR y = _cgSolver.solve(rhs);
  pcgTimer.stop();

  if (verbose)
    printf("  PCG iters: %3i err: %6.4e \n", (int)_cgSolver.iterations(), (float)_cgSolver.error());
#else
  TIMER choleskyTimer("Cholesky Solve");
  Eigen::SimplicialLDLT<SPARSE_MATRIX> solver;
  solver.compute(LHS);
  VECTOR y = solver.solve(rhs);
  choleskyTimer.stop();
#endif

  // aliasing _solution to \Delta v just to make clear what we're doing here
  VECTOR& vDelta = _solution;
  vDelta = y + z;

  // update velocity
  _velocity = _velocity + vDelta;
  _position = _position + _dt * _velocity;

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;

  const bool constraintsChanged = findSeparatingSurfaceConstraints(_b);

  // ONLY change the constraints here. Trying to do it inside the Newton loop
  // is too oscillatory, and the direction that the forces point in can oscilliate
  // if something was separating
  //
  // Delete constraints this AFTER the Newmark update. Otherwise, we don't 
  // know how to filter off the position changes due to purely the constraints, 
  // and we see huge accelerations
  if (constraintsChanged)
  {
    deleteSurfaceConstraints(verbose);
    updateSurfaceConstraints();
    buildConstraintMatrix();
    updateConstraintTargets();
  }
  // otherwise, update the targets, but the constraint matrix should not have changed.
  else
  {
    updateSurfaceConstraints();
    updateConstraintTargets();
  }
  
  // update node positions
  _tetMesh.setDisplacement(_position);

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;

  return true;
}

} // HOBAK
} // TIMESTEPPER
