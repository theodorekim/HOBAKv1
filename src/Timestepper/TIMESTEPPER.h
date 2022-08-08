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
#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Damping/Volume/DAMPING.h"

#include <iostream>

namespace HOBAK {
namespace TIMESTEPPER {

////////////////////////////////////////////////////////////////////////////////////////////////////
// This implements Baraff-Witkin-style constraints from
//  "Large Steps in Cloth Simulation", SIGGRAPH 1998
// by building the system described in 
//  "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
////////////////////////////////////////////////////////////////////////////////////////////////////
class TIMESTEPPER
{
public:
  TIMESTEPPER(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic);
  TIMESTEPPER(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, VOLUME::DAMPING& damping);
  virtual ~TIMESTEPPER();

  VECTOR& externalForces()             { return _externalForces; };
  const VECTOR& externalForces() const { return _externalForces; };
  const vector<PLANE_CONSTRAINT>& planeConstraints() const { return _planeConstraints; };
  const VOLUME::HYPERELASTIC& material() const             { return _hyperelastic; };
  const string materialName() const                        { return _hyperelastic.name(); };
  const string& name() const                               { return _name; };
  const REAL dt() const                                    { return _dt; };

  const TET_MESH& tetMesh() const       { return _tetMesh; };
  const VECTOR position() const         { return _position; };
  const VECTOR positionOld() const      { return _positionOld; };
  const VECTOR velocity() const         { return _velocity; };
  VECTOR& position()                    { return _position; };
  VECTOR& positionOld()                 { return _positionOld; };
  VECTOR& velocity()                    { return _velocity; };
  // Newmark needs to recompute things if this is set differently
  //REAL& dt()                            { return _dt; };
  //
  const bool& vertexFaceSelfCollisionsOn() const { return _vertexFaceSelfCollisionsOn; };
  const bool& edgeEdgeSelfCollisionsOn() const   { return _edgeEdgeSelfCollisionsOn; };
  const REAL& collisionStiffness() const         { return _collisionStiffness; };
  const REAL& collisionDampingBeta() const       { return _collisionDampingBeta; };
  bool& vertexFaceSelfCollisionsOn()    { return _vertexFaceSelfCollisionsOn; };
  bool& edgeEdgeSelfCollisionsOn()      { return _edgeEdgeSelfCollisionsOn; };
  REAL& collisionStiffness()            { return _collisionStiffness; };
  REAL& collisionDampingBeta()          { return _collisionDampingBeta; };
  virtual void setDt(const REAL dt)     { _dt = dt; };
  void setRayeligh(const REAL alpha, const REAL beta);

  // velocity at a specific vertex
  const VECTOR3 velocity(unsigned int index) const;

  // take a timestep
  virtual bool solve(const bool verbose) = 0;

  // add a gravity body force to the simulation
  void addGravity(const VECTOR3& bodyForce);

  // add a plane constraint
  void addPlaneConstraint(const PLANE_CONSTRAINT& constraint) { _planeConstraints.push_back(constraint); };
  void clearPlaneConstraints()                                { _planeConstraints.clear(); };
  int totalPlaneConstraints()                                 { return _planeConstraints.size(); };

  // constrain surface nodes inside a kinematic body to move along with that body
  void attachKinematicSurfaceConstraints(const KINEMATIC_SHAPE* shape);

  // constrain all nodes inside a kinematic body to move along with that body
  void attachKinematicConstraints(const KINEMATIC_SHAPE* shape);

  // which nodes are the constrained ones?
  vector<int> constrainedNodes() const;

  // add kinematic collision object to system
  void addKinematicCollisionObject(const KINEMATIC_SHAPE* shape);

  // make all objects lighter or heavier
  void scaleMass(const REAL& scalar)  { _M *= scalar; };

protected:
  // shared initialization across constructors
  void initialize();

  // build the constraint matrix to incorporate Baraff-Witkin-style constraints, 
  // but using the [TJM15] projection matrix
  void buildConstraintMatrix();

  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  //
  // the target formulations are different in dynamics vs. quasistatics, and
  // position vs. velocity updates, so it is pure virtual here
  virtual void updateConstraintTargets() = 0;

  // filter positions to incorporate Baraff-Witkin-style constraints,
  // this one is slightly different in QUASISTATIC, so this functions has been
  // made virtual
  virtual void applyKinematicConstraints();

  // find all the surface vertices that are in collision and create constraints
  // this one is slightly different in QUASISTATIC, i.e. that one can't take
  // velocity into account as a separation condition, so this function has
  // been made virtual
  virtual void findNewSurfaceConstraints(const bool verbose = false);

  // update the closest point positions on surface constraints
  void updateSurfaceConstraints();

  // find all the constraints tagged for deletion and delete them
  void deleteSurfaceConstraints(const bool verbose = false);

  // find the surface constraints that are separating
  bool findSeparatingSurfaceConstraints(const VECTOR& unfiltered);

  // build the mass matrix based on the one-ring volumes
  SPARSE_MATRIX buildMassMatrix();
  
  // build the damping matrix based on the rest pose stiffness
  SPARSE_MATRIX buildRayleighDampingMatrix();

  // do the collision detection, in anticipation of collision response
  void computeCollisionDetection();

  // compute collision forces, add them to the forces and stiffness matrix
  // R = forces, K = stiffness, C = damping
  void computeCollisionResponse(VECTOR& R, SPARSE_MATRIX& K, SPARSE_MATRIX& C, const bool verbose = false);

  REAL _residual;
  int _seenPCGIterations;

  TET_MESH& _tetMesh;
  VOLUME::HYPERELASTIC& _hyperelastic;
  VOLUME::DAMPING* _damping;

  int _DOFs;
  VECTOR _forces;
  VECTOR _externalForces;

  // RHS of the solve
  VECTOR _b;

  // result of the solve
  VECTOR _solution;

  // everybody needs a scratchpad sometimes
  VECTOR _temp;

  // constraint matrix
  SPARSE_MATRIX _S;
  SPARSE_MATRIX _IminusS;

  // constraint targets
  VECTOR _constraintTargets;

  // is the vertex already experiencing a kinematic collision?
  vector<bool> _inCollision;

  // constraints to have vertices move with a kinematic body
  vector<KINEMATIC_CONSTRAINT> _kinematicConstraints;

  // constraints to have vertex slide along a plane
  vector<PLANE_CONSTRAINT> _planeConstraints;

  // kinematic collision objects
  vector<const KINEMATIC_SHAPE*> _collisionObjects;

  // variables to solve for
  VECTOR _position;
  VECTOR _velocity;

  // in case the user wants to rewind to the previous positions 
  VECTOR _positionOld;

  // timestep
  REAL _dt;
  REAL _rayleighAlpha;
  REAL _rayleighBeta;
  
  // solver vars
  SPARSE_MATRIX _A;
  SPARSE_MATRIX _M;
  
  // global Hessian matrix
  SPARSE_MATRIX _H;
  Eigen::ConjugateGradient<SPARSE_MATRIX, Eigen::Lower|Eigen::Upper> _cgSolver;

  // what's this timestepper called?
  string _name;

  // are self-collisions activated?
  bool _vertexFaceSelfCollisionsOn;
  bool _edgeEdgeSelfCollisionsOn;

  // collision spring and damping constants
  REAL _collisionStiffness;
  REAL _collisionDampingBeta;
};

} // HOBAK
} // TIMESTEPPER

#endif
