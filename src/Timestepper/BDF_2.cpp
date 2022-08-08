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
#include "BDF_2.h"
#include "TIMER.h"
#include "Damping/Volume/GREEN_DAMPING.h"
#include "util/MATRIX_UTIL.h"
#include "Geometry/TET_MESH_FASTER.h"
#include <float.h>
#include <iostream>

using namespace std;

namespace HOBAK {
namespace TIMESTEPPER {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
BDF_2::BDF_2(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic) :
  NEWMARK(tetMesh, hyperelastic)
{
  cout << " Initializing BDF-2 ... " << flush;
  //_maxNewtonIterations  = 1;
  _maxNewtonIterations  = 3;
  _rayleighAlpha = 0.01;
  _rayleighBeta = 0.01;
  //_rayleighAlpha = 0.001;
  //_rayleighBeta = 0.001;

  _positionOlder.resize(_DOFs);
  _velocityOlder.resize(_DOFs);

  _positionOlder.setZero();
  _velocityOlder.setZero();

  // implicit Newmark
	_beta = 0.25;
	_gamma = 0.50;

  // these are the BDF-2 settings of the Newmark constants
  // this is here mostly for illustrative purposes. Many of the alphas
  // never get called, since we can assume they are zero.
	_alpha[0] = (9.0 / 4.0) / (_dt * _dt);
	_alpha[1] = 0;
	_alpha[2] = 0;
	_alpha[3] = (3.0 / 2.0) / _dt;
	_alpha[4] = 0;
	_alpha[5] = 0;

  cout << "done. " << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void BDF_2::setDt(const REAL dt)
{
  _dt = dt;

  // these are the BDF-2 settings of the Newmark constants
  // this is here mostly for illustrative purposes. Many of the alphas
  // never get called, since we can assume they are zero.
	_alpha[0] = (9.0 / 4.0) / (_dt * _dt);
	_alpha[1] = 0;
	_alpha[2] = 0;
	_alpha[3] = (3.0 / 2.0) / _dt;
	_alpha[4] = 0;
	_alpha[5] = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool BDF_2::solveRayleighDamped(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " BDF-2 RAYLEIGH SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }

  // get the damping matrix
  SPARSE_MATRIX C = buildRayleighDampingMatrix();

  // copy the current values into the old values
  _positionOlder = _positionOld;
  _velocityOlder = _velocityOld;
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
    const int rank = R.size();
    SPARSE_MATRIX collisionC(rank, rank);
    computeCollisionResponse(R,K,collisionC);

    // compute the LHS of the residual:
    TIMER rhsTimer("Forming RHS");
    VECTOR v = (1.0 / _dt) * (1.5 * _position - 2.0 * _positionOld + 0.5 * _positionOlder);
    VECTOR a = (1.0 / _dt) * (1.5 * v         - 2.0 * _velocityOld + 0.5 * _velocityOlder);
    _temp = -_M * a;

    // compute the RHS of the residual:
    //_b = C * v; // Newton way
    _b = C * (v - _velocityOld); // Baraff-Witkin way (livelier)

    // assemble full residual: LHS + RHS + R - F
    _b += _temp + R + _externalForces;
    rhsTimer.stop();

    // store so we can detect contact breaking later
    unfiltered = _b;

    // assemble system matrix A
    TIMER lhsTimer("Forming LHS");
    _A = _alpha[0] * _M - _alpha[3] * (C + collisionC) - K;

    // in [TJM15], this is c = b - Az (page 8, top of column 2)
    VECTOR c = _b - _A * z;
   
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
      printf("  Step: %2i PCG iters: %3i err: %6.4e Newton residual %g\n", step, (int)_cgSolver.iterations(), 
             (float)_cgSolver.error(), maxR);
    }

    // aliasing _solution to \Delta x just to make clear what we're doing here
    VECTOR& xDelta = _solution;
    xDelta = y + z;
  
    // update positions
    _position += xDelta;

    // when checking against normals, unfiltered should be negated for Newmark
    //unfiltered *= -1.0;
    //const bool constraintsChanged = findSeparatingSurfaceConstraints(unfiltered);
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
   
    step++;
  }
	
  // update velocity
  _velocity     = (1.0 / _dt) * (1.5 * _position - 2.0 * _positionOld + 0.5 * _positionOlder);

	// update acceleration
  _acceleration = (1.0 / _dt) * (1.5 * _velocity - 2.0 * _velocityOld + 0.5 * _velocityOlder);

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
bool BDF_2::solveEnergyDamped(const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " BDF-2 ENERGY DAMPED SOLVE " << _currentTimestep << endl;
    cout << "=================================================" << endl;
  }
  cout << " NOT IMPLEMENTED! " << endl;
  exit(0);
}

} // HOBAK
} // TIMESTEPPER
