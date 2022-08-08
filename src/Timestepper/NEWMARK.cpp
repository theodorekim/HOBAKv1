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
#include "NEWMARK.h"
#include "TIMER.h"
#include "Damping/Volume/GREEN_DAMPING.h"
#include "util/MATRIX_UTIL.h"
#include "Geometry/TET_MESH_FASTER.h"
#include <float.h>
#include <iostream>

#define VF_TEST_ENABLED 0
#define EE_TEST_ENABLED 1

using namespace std;

namespace HOBAK {
namespace TIMESTEPPER {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
NEWMARK::NEWMARK(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic) :
  TIMESTEPPER(tetMesh, hyperelastic)
{
  cout << " Initializing Newmark ... " << flush;
  _minNewtonIterations  = 1;
  //_maxNewtonIterations  = 3;
  _maxNewtonIterations  = 10;
  //_maxNewtonIterations  = 20;
  _residualTolerance = 1e-7;

  _dt = 1.0 / 60.0;
  _rayleighAlpha = 0.001;
  _rayleighBeta = 0.001;
  //_rayleighBeta = 0.25;

  _seenNewtonIterations = -1;

  _time = 0;
  _currentTimestep = 0;
  
  _name = string("NEWMARK");

  _acceleration.resize(_DOFs);
  _positionOld.resize(_DOFs);
  _velocityOld.resize(_DOFs);
  _accelerationOld.resize(_DOFs);
  
  _acceleration.setZero();
  _positionOld.setZero();
  _velocityOld.setZero();
  _accelerationOld.setZero();

  // implicit Newmark
	_beta = 0.25;
	_gamma = 0.50;

	_alpha[0] = 1.0 / (_beta * _dt * _dt);
	_alpha[1] = 1.0 / (_beta * _dt);
	_alpha[2] = (1.0 - 2.0 * _beta) / (2.0 * _beta);
	_alpha[3] = _gamma / (_beta * _dt);
	_alpha[4] = 1.0 - _gamma / _beta;
	_alpha[5] = (1.0 - _gamma / (2.0 * _beta)) * _dt;

  /*
  // DEBUG
  _alpha[0] = 1.0 / (_dt * _dt);
  _alpha[1] = 1.0 / _dt;
  _alpha[2] = 0.0;
  _alpha[3] = 1.0 / _dt;
  _alpha[4] = 0.0;
  _alpha[5] = 0.0;
  */

  /*
  // alphas for acceleration-level update
  _accelerationAlpha[0] = _dt * (1.0 - _gamma);
  _accelerationAlpha[1] = _dt * _gamma;
  _accelerationAlpha[2] = _dt;
  _accelerationAlpha[3] = _dt * _dt * (0.5 - _beta);
  _accelerationAlpha[4] = _dt * _dt * _beta;
  */

  _outputNewton = false;
  //_damping = new VOLUME::GREEN_DAMPING(0.001);
  //_damping = new VOLUME::GREEN_DAMPING(0.01);
  _damping = NULL;

  cout << "done. " << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void NEWMARK::setDt(const REAL dt)
{
  _dt = dt;

	_alpha[0] = 1.0 / (_beta * _dt * _dt);
	_alpha[1] = 1.0 / (_beta * _dt);
	_alpha[2] = (1.0 - 2.0 * _beta) / (2.0 * _beta);
	_alpha[3] = _gamma / (_beta * _dt);
	_alpha[4] = 1.0 - _gamma / _beta;
	_alpha[5] = (1.0 - _gamma / (2.0 * _beta)) * _dt;

  /*
  // alphas for acceleration-level update
  _accelerationAlpha[0] = _dt * (1.0 - _gamma);
  _accelerationAlpha[1] = _dt * _gamma;
  _accelerationAlpha[2] = _dt;
  _accelerationAlpha[3] = _dt * _dt * (0.5 - _beta);
  _accelerationAlpha[4] = _dt * _dt * _beta;
  */
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the displacement targets the the Baraff-Witkin-style constraints
// are trying to hit. Assumes that buildConstraintMatrix() has already been called
///////////////////////////////////////////////////////////////////////////////////////////////////////
void NEWMARK::updateConstraintTargets()
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
bool NEWMARK::solve(const bool verbose)
{
  // if there's an energy-based damping, use it
  if (_damping != NULL)
    return solveEnergyDamped(verbose);
  return solveRayleighDamped(verbose);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool NEWMARK::solveRayleighDamped(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " NEWMARK RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // get the damping matrix
  SPARSE_MATRIX C = buildRayleighDampingMatrix();

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;

  // should need to call once, but then preserved throughout
  applyKinematicConstraints();

  // store the filtered b for later
  VECTOR unfiltered;

  // store the internal forces for later
  VECTOR R;
  VECTOR z;

  // do Newton-Raphson
  REAL eps = _residualTolerance;
  REAL maxR = eps * 10;
  int step = 0;
  
  while (step < _maxNewtonIterations && maxR > eps)
  {
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
    R = _tetMesh.computeHyperelasticForces(_hyperelastic);

    // get the reduced stiffness matrix
    SPARSE_MATRIX K = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);

    // compute collision forces and stiffnesses
    computeCollisionResponse(R,K,C);

    // compute the LHS of the residual:
    // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
    TIMER rhsTimer("Forming RHS");
    _temp = _alpha[0] * (_position - _positionOld) - _alpha[2] * _accelerationOld -_alpha[1] * _velocityOld;
    _temp = _M * _temp;

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    _b = _alpha[3] * (_position - _positionOld) + _alpha[4] * _velocityOld + _alpha[5] * _accelerationOld;
    _b = C * (-1.0 * _b);

    // assemble full residual: LHS + RHS + R - F
    _b += _temp - R - _externalForces;
    rhsTimer.stop();

    // store so we can detect contact breaking later
    unfiltered = _b;

    // assemble system matrix A
    TIMER lhsTimer("Forming LHS");
    _A = _alpha[0] * _M - _alpha[3] * C - K;

    // in [TJM15], this is c = b - Az (page 8, top of column 2)
    // but things are negated in Newmark
    VECTOR c = _b + _A * z;
   
    // just for clarity, go ahead and form the RHS and LHS explicitly
    VECTOR rhs = _S * c;
    SPARSE_MATRIX LHS = _S * _A * _S + _IminusS;
    lhsTimer.stop();

    maxR = rhs.squaredNorm();

    TIMER pcgTimer("PCG Solve");
    _cgSolver.compute(LHS);
    VECTOR y = _cgSolver.solve(rhs);
    pcgTimer.stop();

    if (verbose)
    {
      printf("  Step: %2i PCG iters: %3i err: %6.4e residual %g\n", step, (int)_cgSolver.iterations(), 
             (float)_cgSolver.error(), maxR);
    }

    // aliasing _solution to \Delta x just to make clear what we're doing here
    VECTOR& xDelta = _solution;
    xDelta = y - z;
  
    // update positions
    _position -= xDelta;

    // when checking against normals, unfiltered should be negated for Newmark
    unfiltered *= -1.0;
    const bool constraintsChanged = findSeparatingSurfaceConstraints(unfiltered);

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
   
    step++;
  }
	
  // update velocity
  _velocity = _alpha[3] * (_position - _positionOld) + _alpha[4] * _velocityOld + _alpha[5] * _accelerationOld;

	// update acceleration
  _acceleration = _alpha[0] * (_position - _positionOld) - _alpha[1] * _velocityOld - _alpha[2] * _accelerationOld;

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
bool NEWMARK::solveEnergyDamped(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " NEWMARK ENERGY DAMPED SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
	  
  // should need to call once, but then preserved throughout
  applyKinematicConstraints();
 
  // store the filtered b for later
  VECTOR unfiltered;

  // store the internal elastic and damping forces for later
  VECTOR R;

  VECTOR z;

  // do Newton-Raphson
  REAL eps = _residualTolerance;
  REAL maxR = eps * 10;
  int step = 0;

  while (step < _maxNewtonIterations && maxR > eps)
  {
    // do collision detection, including spatial data structure updates 
    computeCollisionDetection();

    // build new constraints and see if we should break any
    findNewSurfaceConstraints(verbose);
    buildConstraintMatrix();

    // update acceleration, in anticipation of f = M * a
    _velocity = _alpha[3] * (_position - _positionOld) + _alpha[4] * _velocityOld + _alpha[5] * _accelerationOld;
    _acceleration = _alpha[0] * (_position - _positionOld) - _alpha[1] * _velocityOld - _alpha[2] * _accelerationOld;
    _velocity = _S * _velocity;
    _acceleration = _S * _acceleration;

    _tetMesh.setDisplacement(_position);
    _tetMesh.computeFs();
    _tetMesh.computeFdots(_velocity);
    _tetMesh.computeSVDs();

    // "z is a vector of the desired values for the constrained variables",
    // from [TJM15], Section 8, paragraph 3. We apply the _IminusS because
    // _constraintTargets did not project off the null directions
    updateConstraintTargets();
    z =_IminusS * _constraintTargets;

    // get the internal forces, both elastic and damping
    R = _tetMesh.computeInternalForces(_hyperelastic, *_damping);
    
    // get the stiffness matrices
    SPARSE_MATRIX K = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);
    SPARSE_MATRIX Kdamping = _tetMesh.computeDampingHessian(*_damping);

    // compute collision forces and stiffnesses
    computeCollisionResponse(R,K,Kdamping);

    // assemble full residual: RHS + R - F
    TIMER rhsTimer("Forming RHS");
    _b = _M * _acceleration - R - _externalForces;
    rhsTimer.stop();

    // store so we can detect contact breaking later
    unfiltered = _b;

    // assemble system matrix A. Since Kdamping is a derivative w.r.t. velocity,
    // so it should get the velocity Newmark constant
    TIMER lhsTimer("Forming LHS");
    _A = _alpha[0] * _M + _alpha[3] * Kdamping - K;

    // in [TJM15], this is c = b - Az (page 8, top of column 2)
    // but things are negated in Newmark
    VECTOR c = _b + _A * z;
   
    // just for clarity, go ahead and form the RHS and LHS explicitly
    VECTOR rhs = _S * c;
    SPARSE_MATRIX LHS = _S * _A * _S + _IminusS;
    lhsTimer.stop();

    maxR = rhs.squaredNorm();

    TIMER pcgTimer("PCG Solve");
    _cgSolver.compute(LHS);
    VECTOR y = _cgSolver.solve(rhs);
    pcgTimer.stop();

    if (verbose)
    {
      printf("  Step: %2i PCG iters: %3i err: %6.4e residual %g\n", step, (int)_cgSolver.iterations(), 
             (float)_cgSolver.error(), maxR);
    }

    // aliasing _solution to \Delta x just to make clear what we're doing here
    VECTOR& xDelta = _solution;
    xDelta = y - z;

    // update everybody
    _position -= xDelta;
    _velocity = _alpha[3] * (_position - _positionOld) + _alpha[4] * _velocityOld + _alpha[5] * _accelerationOld;
    _acceleration = _alpha[0] * (_position - _positionOld) - _alpha[1] * _velocityOld - _alpha[2] * _accelerationOld;

    // when checking against normals, unfiltered should be negated for Newmark
    unfiltered *= -1.0;
    const bool constraintsChanged = findSeparatingSurfaceConstraints(unfiltered);

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
   
    step++;
  }

  // one last time, since the while loop exited without updating
  _velocity = _alpha[3] * (_position - _positionOld) + _alpha[4] * _velocityOld + _alpha[5] * _accelerationOld;
  _acceleration = _alpha[0] * (_position - _positionOld) - _alpha[1] * _velocityOld - _alpha[2] * _accelerationOld;
  _velocity = _S * _velocity;
  _acceleration = _S * _acceleration;

  if (verbose)
  {
    /*
    cout << " Position norm:         " << _position.norm() << endl;
    cout << " Velocity norm:         " << _velocity.norm() << endl;
    cout << " Velocity old norm:     " << _velocityOld.norm() << endl;
    cout << " Acceleration norm:     " << _acceleration.norm() << endl;
    cout << " Acceleration old norm: " << _accelerationOld.norm() << endl;
    printEntry(_velocity, 0, "_velocity");
    printEntry(_velocityOld, 0, "_velocityOld");
    printEntry(_b, 0, "_b");
    printEntry(R, 0, "R");
    //printEntry(_alpha[3] * _S * (_position - _positionOld), 0, "position term");
    */
  }

  // record which timestep we're on
  _time += _dt;
  _currentTimestep++;

  return true;
}

} // HOBAK
} // TIMESTEPPER
