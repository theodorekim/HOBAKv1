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
#include "BACKWARD_EULER_POSITION.h"
#include "TIMER.h"
#include "Damping/Volume/GREEN_DAMPING.h"
#include "Geometry/TET_MESH_FASTER.h"
#include "util/MATRIX_UTIL.h"
#include <float.h>
#include <iostream>

using namespace std;

namespace HOBAK {
namespace TIMESTEPPER {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
BACKWARD_EULER_POSITION::BACKWARD_EULER_POSITION(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic) :
  TIMESTEPPER(tetMesh, hyperelastic)
{
  initialize();
}

BACKWARD_EULER_POSITION::BACKWARD_EULER_POSITION(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, VOLUME::DAMPING& damping) :
  TIMESTEPPER(tetMesh, hyperelastic, damping)
{
  initialize();
}

BACKWARD_EULER_POSITION::~BACKWARD_EULER_POSITION()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void BACKWARD_EULER_POSITION::initialize()
{
  cout << " Initializing position-based backward Euler ... " << flush;

  _dt = 1.0 / 60.0;
  _rayleighAlpha = 0.01;
  _rayleighBeta = 0.01;
  //_rayleighAlpha = 0.001;
  //_rayleighBeta = 0.001;

  _time = 0;
  _currentTimestep = 0;
  
  _name = string("BACKWARD_EULER_POSITION");

  _acceleration.resize(_DOFs);
  _velocityOld.resize(_DOFs);
  
  _acceleration.setZero();
  _velocityOld.setZero();

  //_damping = new VOLUME::ARAP_DAMPING(0.01);
  //_damping = new VOLUME::ARAP_DAMPING(0.1);
  //_damping = new VOLUME::GREEN_DAMPING(0.001);
  //_damping = new VOLUME::GREEN_DAMPING(0.01);
  //_damping = new VOLUME::GREEN_DAMPING(0.1);
  _damping = NULL;

  cout << "done. " << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the displacement targets the the Baraff-Witkin-style constraints
// are trying to hit. Assumes that buildConstraintMatrix() has already been called
///////////////////////////////////////////////////////////////////////////////////////////////////////
void BACKWARD_EULER_POSITION::updateConstraintTargets()
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

    const VECTOR3& displacement = closestPoint - vertex;
    for (int i = 0; i < 3; i++)
      _constraintTargets[index + i] = displacement[i];
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
bool BACKWARD_EULER_POSITION::solve(const bool verbose)
{
  // if there's an energy-based damping, use it
  if (_damping != NULL)
    return solveEnergyDamped(verbose);
  return solveRayleighDamped(verbose);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
//
// solves with self collisions, and collision forces are Rayleigh damped
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool BACKWARD_EULER_POSITION::solveRayleighDamped(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " BACKWARD_EULER_POSITION RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // get the damping matrix
  SPARSE_MATRIX C = buildRayleighDampingMatrix();

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;

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
#if VERY_VERBOSE
  SPARSE_MATRIX elasticK = -K;
#endif

  // collision damping only appears on LHS
  const int rank = R.size();
  SPARSE_MATRIX collisionC(rank, rank);

  // compute collision forces and stiffnesses
  computeCollisionResponse(R,K,collisionC, true);

  // compute the RHS of the residual:
  TIMER rhsTimer("Forming Initial RHS");
  const REAL invDt = 1.0 / _dt;
  const REAL invDt2 = invDt * invDt;
  _b = (invDt * _M - C) * _velocity + R + _externalForces;

  // collisionC does *not* appear here, since the damping only depends on
  // v^{t+1}, which only depends on \Delta x, which is the variable
  // we are solving for
  //_b = (invDt * _M - C - collisionC) * _velocity + R + _externalForces;
  rhsTimer.stop();

  // assemble system matrix A
  TIMER lhsTimer("Forming Initial LHS");
  _A = _M * invDt2 - (C + collisionC) * invDt - K; // w/ Rayleigh
  lhsTimer.stop();

  // in [TJM15], this is c = b - Az (page 8, top of column 2)
  TIMER projectionTimer("PPCG projection");
  VECTOR c = _b - _A * z;

#if VERY_VERBOSE
  if (verbose)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    cout << " b:      " << _b.norm() << endl;
    cout << " R:      " << R.norm() << endl;
    cout << " f_ext:  " << _externalForces.norm() << endl;
    cout << " v term: " << ((invDt * _M - C) * _velocity).norm() << endl;
    cout << " K, elastic (2-norm):   " << elasticK.norm() << endl;
    cout << " K, elastic (inf-norm): " << MATRIX(elasticK).lpNorm<Eigen::Infinity>() << endl;
    cout << " K, elastic, max eig:   " << largestEigenvalue(elasticK) << endl;
    // slow, and almost certainly zero
    //cout << " K, elastic, min eig:   " << smallestEigenvalue(elasticK) << endl;
    cout << " K:      " << K.norm() << endl;
    cout << " C:      " << C.norm() << endl;
    cout << " A:      " << _A.norm() << endl;
    cout << " A max eig: " << largestEigenvalue(_A) << endl;
    //cout << " A min eig: " << smallestEigenvalue(_A) << endl;
    cout << " alpha:  " << _rayleighAlpha << endl;
    cout << " beta:   " << _rayleighBeta << endl;
  }
#endif

  // just for clarity, go ahead and form the RHS and LHS explicitly
  //
  // Since _S is sparse, this multipy could be accelerated significantly, 
  // but leaving it as it is for now
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
  
  // aliasing _solution to \Delta x just to make clear what we're doing here
  VECTOR& xDelta = _solution;
  xDelta = y + z;

#if VERY_VERBOSE
  if (verbose)
  {
    cout << " xDelta (2-norm):     " << xDelta.norm() << endl;
    cout << " xDelta (inf-norm):   " << xDelta.lpNorm<Eigen::Infinity>() << endl;
  }
#endif

  // update positions
  _position += xDelta;
  /*
  cout << " x delta: " << xDelta.norm() << endl;
  cout << " b        " << _b.norm() << endl;
  cout << " R        " << R.norm() << endl;
  cout << " ext      " << _externalForces.norm() << endl;
  cout << " Cv:      " << (-C * _velocity).norm() << endl;
  cout << " Ma:      " << (invDt * _M * _velocity).norm() << endl;
  cout << " v:       " << _velocity.norm() << endl;
  cout << " C:       " << C.norm() << endl;
  */

  // when checking against normals, unfiltered should be negated for Newmark
  const bool constraintsChanged = findSeparatingSurfaceConstraints(_b);

  // see if any of the constraints changed. Used to be that this was outside the Newton loop
  // because the behavior was too oscillatory, but causes too many penetrations to slip
  // through when the Poisson's ratio gets high
  if (constraintsChanged)
  {
    deleteSurfaceConstraints(verbose);
    updateSurfaceConstraints();
    buildConstraintMatrix();
    updateConstraintTargets();
  }
  // update the targets, but the constraint matrix should not have changed.
  else
  {
    updateSurfaceConstraints();
    updateConstraintTargets();
  }
  // update node positions
  _tetMesh.setDisplacement(_position);

  // update velocity
  _velocity = invDt * (_position - _positionOld);

	// update acceleration
  _acceleration = invDt * (_velocity - _velocityOld);

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;
  _acceleration = _S * _acceleration;

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool BACKWARD_EULER_POSITION::solveEnergyDamped(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " BACKWARD_EULER_POSITION ENERGY DAMPED SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;

  // should need to call once, but then preserved throughout
  applyKinematicConstraints();

  // store the filtered b for later
  VECTOR unfiltered;

  // store the internal forces for later
  VECTOR R;
  VECTOR z;

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
  z =_IminusS * _constraintTargets;

  // get the internal forces
  R = _tetMesh.computeInternalForces(_hyperelastic, *_damping);

  // get the reduced stiffness matrix
  SPARSE_MATRIX K = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);
  SPARSE_MATRIX C = _tetMesh.computeDampingHessian(*_damping);

  // compute collision forces and stiffnesses
  computeCollisionResponse(R,K,C);

  // compute the RHS of the residual:
  TIMER rhsTimer("Forming Initial RHS");
  const REAL invDt = 1.0 / _dt;
  const REAL invDt2 = invDt * invDt;
  _b = (invDt * _M) * _velocity + R + _externalForces;
  rhsTimer.stop();

  // assemble system matrix A
  TIMER lhsTimer("Forming Initial LHS");
  _A = _M * invDt2 - C * invDt - K; // w/ Rayleigh
  lhsTimer.stop();

  // in [TJM15], this is c = b - Az (page 8, top of column 2)
  TIMER projectionTimer("PPCG projection");
  VECTOR c = _b - _A * z;
 
  // just for clarity, go ahead and form the RHS and LHS explicitly
  //
  // Since _S is sparse, this multipy could be accelerated significantly, 
  // but leaving it as it is for now
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

  // aliasing _solution to \Delta x just to make clear what we're doing here
  VECTOR& xDelta = _solution;
  xDelta = y + z;

  // update positions
  _position += xDelta;

  // when checking against normals, unfiltered should be negated for Newmark
  const bool constraintsChanged = findSeparatingSurfaceConstraints(_b);

  // see if any of the constraints changed. Used to be that this was outside the Newton loop
  // because the behavior was too oscillatory, but causes too many penetrations to slip
  // through when the Poisson's ratio gets high
  if (constraintsChanged)
  {
    deleteSurfaceConstraints(verbose);
    updateSurfaceConstraints();
    buildConstraintMatrix();
    updateConstraintTargets();
  }
  // update the targets, but the constraint matrix should not have changed.
  else
  {
    updateSurfaceConstraints();
    updateConstraintTargets();
  }
  // update node positions
  _tetMesh.setDisplacement(_position);
  
  // update velocity
  _velocity = invDt * (_position - _positionOld);

	// update acceleration
  _acceleration = invDt * (_velocity - _velocityOld);

  // In addition to filtering by _S here, the right thing is to pick up the velocity of the kinematic
  // object in the constraint direction. I.e. we've implemented the _S part, but not the _IminusS part
  // of this update. For now, stomping these components to zero will at least keep things stable,
  // so keeping it for future work when somebody wants to paddle wheel
  _velocity = _S * _velocity;
  _acceleration = _S * _acceleration;

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;

  return true;
}

} // HOBAK
} // TIMESTEPPER
