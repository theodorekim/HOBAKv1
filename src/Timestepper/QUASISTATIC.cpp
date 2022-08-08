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
#include "QUASISTATIC.h"
#include "TIMER.h"
#include <float.h>
#include <iostream>

using namespace std;

namespace HOBAK {
namespace TIMESTEPPER {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
QUASISTATIC::QUASISTATIC(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic) :
  TIMESTEPPER(tetMesh, hyperelastic)
{
  _minIterations  = 1;
  _maxIterations  = 3;
  _residualTolerance = 1e-12;

  _seenNewtonIterations = -1;
  _name = string("QUASISTATIC");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// filter positions to incorporate Baraff-Witkin-style constraints
///////////////////////////////////////////////////////////////////////////////////////////////////////
void QUASISTATIC::applyKinematicConstraints()
{
  vector<VECTOR3>& vertices = _tetMesh.vertices();
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];

    // directly pin the tet mesh position according to the constraint
    const VECTOR3& localPosition = constraint.localPosition;
    VECTOR3 world = constraint.shape->localVertexToWorld(localPosition);
    vertices[constraint.vertexID] = world;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the displacement targets the the Baraff-Witkin-style constraints
// are trying to hit. Assumes that buildConstraintMatrix() has already been called
//
// unlike the dynamics case, the constraint is just set to the goal position,
// not the delta w.r.t. the goal position
///////////////////////////////////////////////////////////////////////////////////////////////////////
void QUASISTATIC::updateConstraintTargets()
{
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
// The Baraff-Witkin plane constraints are being imposed at the beginning and then frozen during
// the Newton solve. At the end of the solve, we check if the forces are separating from the
// constraint and if they are, we delete the constraint and run one more Newton iteration.
//
// This seems to prevent both sticking and pop-through, while nothing else does.
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool QUASISTATIC::solve(const bool verbose, const REAL stepSize)
#if 1
{
  TIMER functionTimer(__FUNCTION__);
  // shouldn't be taking a step size per Newton iteration that is bigger
  // than the solution direction itself. (Not using momentum methods here,
  // so what would it even mean?)
  assert(stepSize <= 1.0);

  _seenPCGIterations    = 0;
  _seenNewtonIterations = 0;
  _residual             = FLT_MAX;
  _solution             = _tetMesh.getDisplacement();

  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " QUASISTATIC SOLVE" << endl;
    cout << "=================================================" << endl;
  }

  /*
  cout << " vertex 0: " << _tetMesh.vertices()[0].transpose() << endl;
  cout << " constraint target for vertex 0: " << _constraintTargets[0] << ", " << _constraintTargets[1] << ", " << _constraintTargets[2] << endl;
  cout << " b for vertex 0: " << _b[0] << ", " << _b[1] << ", " << _b[2] << endl;
  */

  bool success = false;
  int iterations = 0;
 
  // build new constraints and see if we should break any
  findNewSurfaceConstraints(verbose);
  buildConstraintMatrix();
  applyKinematicConstraints();  // should be redundant after first call

  // don't change the constraints inside the Newton loop, it's too oscillatory
  while (((iterations < _minIterations) || (iterations < _maxIterations)) && 
         (_residual > _residualTolerance))
  {
    if (iterations == 0)
    {
      _tetMesh.computeFs();
      _tetMesh.computeSVDs();
      _forces.setZero();
      _forces += _externalForces;
      _forces += _tetMesh.computeHyperelasticForces(_hyperelastic);
    }

    updateConstraintTargets();
    _b = _S * _forces + _IminusS * _constraintTargets;
    _H = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);

    // filter _forces and _H to incorporate Baraff-Witkin-style constraints
    //SPARSE_MATRIX A = _S * _H * _S + _IminusS;
    SPARSE_MATRIX A = _IminusS - _S * _H * _S;

    _cgSolver.compute(A);
    _solution = _cgSolver.solve(_b);
    
    _tetMesh.addDisplacement(stepSize * _solution);

    _tetMesh.computeFs();
    _tetMesh.computeSVDs();
    _forces.setZero();
    _forces += _externalForces;
    _forces += _tetMesh.computeHyperelasticForces(_hyperelastic);

    // apply the Baraff-Witkin-style constraint targets
    updateConstraintTargets();
    _residual = _b.dot(_b);
    //cout << " constraint target for vertex 0: " << _constraintTargets[0] << ", " << _constraintTargets[1] << ", " << _constraintTargets[2] << endl;
    //cout << " b for vertex 0: " << _b[0] << ", " << _b[1] << ", " << _b[2] << endl;

    iterations++;

    if (verbose)
    {
      cout << "Newton iteration " << iterations << endl;
      cout << "\t CG iterations: " << _cgSolver.iterations() << endl;
      cout << "\t Residual: " << _residual << endl << endl;
    }
  }

  // find which constraints are separating and need to be deleted
  VECTOR unfiltered = _forces;
  const bool constraintsChanged = findSeparatingSurfaceConstraints(unfiltered);

  // ONLY change the constraints here. Trying to do it inside the Newton loop
  // is too oscillatory, and the direction that the forces point in can oscilliate
  // if something was separating, delete those constraints and do one more solve
  if (constraintsChanged)
  {
    deleteSurfaceConstraints();
    buildConstraintMatrix();
    updateConstraintTargets();
    _b = _S * _forces + _IminusS * _constraintTargets;
    _H = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);

    // filter _forces and _H to incorporate Baraff-Witkin-style constraints
    //SPARSE_MATRIX A = _S * _H * _S + _IminusS;
    SPARSE_MATRIX A = _IminusS - _S * _H * _S;
    _cgSolver.compute(A);
    _solution = _cgSolver.solve(_b);
    
    _tetMesh.addDisplacement(stepSize * _solution);
  }

  _position = _tetMesh.getDisplacement();

  success = _residual <= _residualTolerance;
  return success;
}
#else
{
  // shouldn't be taking a step size per Newton iteration that is bigger
  // than the solution direction itself. (Not using momentum methods here,
  // so what would it even mean?)
  assert(stepSize <= 1.0);

  // find and build all the surface constraints
  findNewSurfaceConstraints(verbose);

  // build _S and _IminusS
  buildConstraintMatrix();
  applyKinematicConstraints();
  _tetMesh.computeFs();

  _forces.setZero();
  _forces += _externalForces;
  _forces += _tetMesh.computeHyperelasticForces(_hyperelastic);

  updateConstraintTargets();
  _b = _S * _forces + _IminusS * _constraintTargets;

  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " QUASISTATIC SOLVE" << endl;
    cout << "=================================================" << endl;
  }

  _seenPCGIterations    = 0;
  _seenNewtonIterations = 0;
  _residual             = _b.dot(_b);
  _solution             = _tetMesh.getDisplacement();

  cout << " Beginning force residual: " << _residual << endl;

  // DEBUG
  cout << " vertex 0: " << _tetMesh.vertices()[0].transpose() << endl;
  cout << " constraint target for vertex 0: " << _constraintTargets[0] << ", " << _constraintTargets[1] << ", " << _constraintTargets[2] << endl;
  cout << " b for vertex 0: " << _b[0] << ", " << _b[1] << ", " << _b[2] << endl;

  bool success = false;
  int iterations = 0;
  
  // if the constraints changed in the inner loop, the run again, because
  // the residual will be different
  //bool constraintsChanged = false;
  while (((iterations < _minIterations) || (iterations < _maxIterations)) && 
         (_residual > _residualTolerance))
         //((_residual > _residualTolerance) || constraintsChanged))
  {
    cout << "===== NEWTON ITERATION " << iterations << " ====== " << endl;
    // be on the lookout for constraints changing
    //constraintsChanged = false;

    // build new constraints and see if we should break any
    findNewSurfaceConstraints(verbose);
    deleteSurfaceConstraints();
    buildConstraintMatrix();
    applyKinematicConstraints();  // should be redundant
    _H = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);

    // filter _forces and _H to incorporate Baraff-Witkin-style constraints
    SPARSE_MATRIX A = _S * _H * _S + _IminusS;

    _cgSolver.compute(A);
    _solution = _cgSolver.solve(_b);
    
    _tetMesh.addDisplacement(stepSize * _solution);
    _tetMesh.computeFs();

    cout << " vertex 0: " << _tetMesh.vertices()[0].transpose() << endl;

    _forces.setZero();
    _forces += _externalForces;
    _forces += _tetMesh.computeHyperelasticForces(_hyperelastic);

    // find which constraints are separating and need to be deleted
    //VECTOR unfiltered = _H * _solution - _b;
    //VECTOR unfiltered = _H * _solution - _forces;
    VECTOR unfiltered = _forces;
    findSeparatingSurfaceConstraints(unfiltered);

    // go ahead and delete the ones that are separating
    //deleteSurfaceConstraints();

    // apply the Baraff-Witkin-style constraint targets
    updateConstraintTargets();
    _b = _S * _forces + _IminusS * _constraintTargets;
    _residual = _b.dot(_b);
    cout << " constraint target for vertex 0: " << _constraintTargets[0] << ", " << _constraintTargets[1] << ", " << _constraintTargets[2] << endl;
    cout << " b for vertex 0: " << _b[0] << ", " << _b[1] << ", " << _b[2] << endl;

    iterations++;

    if (verbose)
    {
      cout << "Newton iteration " << iterations << endl;
      cout << "\t CG iterations: " << _cgSolver.iterations() << endl;
      cout << "\t Residual: " << _residual << endl;
    }
  }

  success = _residual <= _residualTolerance;
  return success;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Note: The Baraff-Witkin constraints don't work that well with a line search, because it gets
// treated like another energy request. May be better to do a separated line search where the 
// elastic DOFs can be eased in, but the Baraff-Witkin DOFs are always hit dead on.
//
// We already have the target displacements stored, so it should be straightforward to factor.
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool QUASISTATIC::solveWithBacktracking(const bool verbose)
{
#if 1
  // filter _forces and _H to incorporate Baraff-Witkin-style constraints
  buildConstraintMatrix();
  applyKinematicConstraints();
  _tetMesh.computeFs();

  // cache the starting psi and the original positions
  //vector<VECTOR3> startingVertices = _tetMesh.vertices();
  REAL startingPsi = _tetMesh.computeHyperelasticEnergy(_hyperelastic);
  vector<VECTOR3> startingVertices;
  //REAL startingPsi;

  _forces.setZero();
  _forces += _tetMesh.computeHyperelasticForces(_hyperelastic);
  updateConstraintTargets();
  _b = _S * _forces + _IminusS * _constraintTargets;

  if (verbose)
  {
    cout << "=================================================" << endl;
    cout << " QUASISTATIC SOLVE" << endl;
    cout << "=================================================" << endl;
  }

  _seenPCGIterations    = 0;
  _seenNewtonIterations = 0;
  _residual        = _b.dot(_b);
  _solution        = _tetMesh.getDisplacement();

  cout << " Beginning force residual: " << _residual << endl;
  cout << " Beginning psi:            " << startingPsi << endl;

  bool success = false;
  int iterations = 0;
  while (((iterations < _minIterations) || (iterations < _maxIterations)) && (_residual > _residualTolerance))
  {
    _H = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);

    // filter _forces and _H to incorporate Baraff-Witkin-style constraints
    //_H = _S * _H * _S + _IminusS;
    _H = _IminusS - _S * _H * _S;

    _cgSolver.compute(_H);
    _solution = _cgSolver.solve(_b);
  
    // pin the constraint part of the solution, don't line search those
    VECTOR constraints = _IminusS * _solution;
    _tetMesh.addDisplacement(constraints);
    startingVertices = _tetMesh.vertices();
    startingPsi = _tetMesh.computeHyperelasticEnergy(_hyperelastic);

    // do a line search on the leftovers
    _solution = _S * _solution;

    // find a step that decreases the energy
    int backtracks = 0;
    REAL stepSize = 1.0;
    REAL bestPsiFound = FLT_MAX;
    REAL bestStepSizeFound = 1.0;
    while (backtracks < 20 && bestPsiFound > startingPsi)
    {
      // apply the perturbation
      _tetMesh.vertices() = startingVertices;
      _tetMesh.addDisplacement(stepSize * _solution);
      applyKinematicConstraints();
      _tetMesh.computeFs();

      // compute the new energy, see if it's better
      REAL currentPsi = _tetMesh.computeHyperelasticEnergy(_hyperelastic);
      if (currentPsi < bestPsiFound)
      {
        bestPsiFound = currentPsi;
        bestStepSizeFound = stepSize;
      }
      cout << " Step size " << stepSize << "\t = " << currentPsi << " starting = " << startingPsi << endl;
      backtracks++;
      //stepSize *= 0.9;
      stepSize *= 0.5;
    }

    if (bestPsiFound > startingPsi)
    {
      cout << " LINE SEARCH FAILED. " << endl;
      bestStepSizeFound = stepSize;
    }

    // set it to the best energy found
    _tetMesh.vertices() = startingVertices;
    _tetMesh.addDisplacement(bestStepSizeFound * _solution);
    applyKinematicConstraints();
    _tetMesh.computeFs();

    // build forces for next time
    _forces.setZero();
    _forces += _tetMesh.computeHyperelasticForces(_hyperelastic);

    // build the final RHS
    updateConstraintTargets();
    _b = _S * _forces + _IminusS * _constraintTargets;

    // check the residual
    _residual = _b.dot(_b);
    
    // beat this energy in next iteration
    startingPsi = bestPsiFound;
    startingVertices = _tetMesh.vertices();

    iterations++;

    if (verbose)
    {
      cout << "Newton iteration " << iterations << endl;
      cout << "\t CG iterations: " << _cgSolver.iterations() << endl;
      cout << "\t Residual: " << _residual << endl;
    }
  }

  success = _residual <= _residualTolerance;
  return success;
#else
  return true;
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// find all the surface vertices that are in collision and create constraints
//
// unlike the dynamics case, velocity is not taken into account here
///////////////////////////////////////////////////////////////////////////////////////////////////////
void QUASISTATIC::findNewSurfaceConstraints(const bool verbose)
{
  const vector<VECTOR3> vertices = _tetMesh.vertices();
  const vector<int> surfaceVertices = _tetMesh.surfaceVertices();

  // build any new constraints
  int newConstraints = 0;
  for (unsigned int y = 0; y < _collisionObjects.size(); y++)
  {
    const KINEMATIC_SHAPE* shape = _collisionObjects[y];
    for (unsigned int x = 0; x < surfaceVertices.size(); x++)
    {
      // get the vertex
      assert(surfaceVertices[x] < int(vertices.size()));
      int vertexID = surfaceVertices[x];

      // if it's already in collision, skip it
      if (_inCollision[vertexID]) continue;

      // see if it's inside the shape
      const VECTOR3& vertex = vertices[vertexID];
      if (!shape->inside(vertex)) continue;

      VECTOR3 closestPoint;
      VECTOR3 closestNormal;
      shape->getClosestPoint(vertex, closestPoint, closestNormal);

      // store the constraint
      PLANE_CONSTRAINT constraint;
      constraint.shape = shape;
      constraint.vertexID = vertexID;
      constraint.localClosestPoint = closestPoint;
      constraint.localNormal = closestNormal;
      constraint.isSeparating = false;
      addPlaneConstraint(constraint);

      _inCollision[vertexID] = true;
      newConstraints++;
    }
  }
  if (verbose)
    cout << " Found " << newConstraints << " new constraints " << endl;
}

} // HOBAK
} // TIMESTEPPER
