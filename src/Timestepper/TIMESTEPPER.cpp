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
#include "TIMESTEPPER.h"
#include "TIMER.h"
#include "Geometry/TET_MESH_FASTER.h"
#include <float.h>
#include <iostream>

using namespace std;

namespace HOBAK {
namespace TIMESTEPPER {

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
TIMESTEPPER::TIMESTEPPER(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic) :
  _tetMesh(tetMesh), _hyperelastic(hyperelastic), _damping(NULL)
{
  initialize();
}

TIMESTEPPER::TIMESTEPPER(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, VOLUME::DAMPING& damping) :
  _tetMesh(tetMesh), _hyperelastic(hyperelastic), _damping(&damping)
{
  initialize();
}

TIMESTEPPER::~TIMESTEPPER()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::initialize()
{
  _residual = FLT_MAX;
  _seenPCGIterations    = -1;

  _DOFs = _tetMesh.DOFs();
  _b.resize(_DOFs);
  _forces.resize(_DOFs);
  _externalForces.resize(_DOFs);
  _constraintTargets.resize(_DOFs);

  _b.setZero();
  _forces.setZero();
  _externalForces.setZero();
  _constraintTargets.setZero();

  int totalVertices = _tetMesh.vertices().size();
  _inCollision.resize(totalVertices);
  for (int x = 0; x < totalVertices; x++)
    _inCollision[x] = false;

  _position.resize(_DOFs);
  _positionOld.resize(_DOFs);
  _velocity.resize(_DOFs);
  _temp.resize(_DOFs);
  _solution.resize(_DOFs);
  
  _position.setZero();
  _positionOld.setZero();
  _velocity.setZero();
  _temp.setZero();
  _solution.setZero();

  _name = string("UNKNOWN");

  _vertexFaceSelfCollisionsOn = true;
  _edgeEdgeSelfCollisionsOn   = true;
  _collisionStiffness = 1.0;
  _collisionDampingBeta = 0.001;

  _dt = 1.0 / 30.0;

  // build the mass matrix once and for all
  _M = buildMassMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// filter positions to incorporate Baraff-Witkin-style constraints
//
// NOTE: this is different from the quasistatic case, in that it stores the kinematic
// corrections in _position, and not directly the vertices of _tetMesh. When the
// displacement of _tetMesh is pinned during the Newton solve, then the _teMesh
// will see the constraints.
//
// The QUASISTATIC class overrides this one with its own implementation, but the
// dynamics solve all share this version
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::applyKinematicConstraints()
{
  const vector<VECTOR3>& restVertices = _tetMesh.restVertices();
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];

    // pin the tet mesh position according to the constraint
    const VECTOR3& localPosition = constraint.localPosition;
    VECTOR3 world = constraint.shape->localVertexToWorld(localPosition);

    const int vertexID = constraint.vertexID;
    const VECTOR3 diff = world - restVertices[vertexID];
    _position[3 * vertexID] = diff[0];
    _position[3 * vertexID + 1] = diff[1];
    _position[3 * vertexID + 2] = diff[2];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the constraint matrix to incorporate Baraff-Witkin-style constraints, but using 
// the "Smoothed aggregation multigrid for cloth simulation" projection matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::buildConstraintMatrix()
{
  SPARSE_MATRIX I(_DOFs, _DOFs);
  I.setIdentity();
  _S = I;
 
  // build the plane constraints for the LHS
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // if this one is tagged for deletion, ignore it
    if (constraint.isSeparating) continue;

    // get the normal direction
    const KINEMATIC_SHAPE* shape = _planeConstraints[x].shape;
    const VECTOR3& localNormal = _planeConstraints[x].localNormal;
    const VECTOR3 normal = shape->localNormalToWorld(localNormal).normalized();

    // build the filter matrix
    const MATRIX3 Sblock = MATRIX3::Identity() - normal * normal.transpose();
    const int vertexID = constraint.vertexID;
    const int index = 3 * vertexID;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        _S.coeffRef(index + i, index + j) = Sblock(i,j);
  }

  // apply the kinematic constraints LAST. These override any prior plane
  // constraints
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = _kinematicConstraints[x];

    // set the filter matrix entries
    const int index = 3 * constraint.vertexID;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        _S.coeffRef(index + i, index + j) = 0.0;
  }

  // store the complement
  _IminusS = I - _S;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the mass matrix based on the one-ring volumes
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TIMESTEPPER::buildMassMatrix()
{
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  const vector<REAL>& volumes = _tetMesh.restOneRingVolumes(); 

  // set diagonal the one-ring volumes
  // TODO: take into account density
  for (int x = 0; x < _tetMesh.totalVertices(); x++)
  {
    const REAL entry = volumes[x];
    for (int y = 0; y < 3; y++)
    {
      TRIPLET triplet(3 * x + y, 3 * x + y, entry);
      triplets.push_back(triplet);
    }
  }
  
  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

#if 0
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " ZEROING MASS MATRIX " << endl;
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  A *= 0.0;
#endif

  return A;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// build the damping matrix based on the rest pose stiffness
///////////////////////////////////////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TIMESTEPPER::buildRayleighDampingMatrix()
{
  // back up current state
  _temp = _tetMesh.getDisplacement();

  // set to zero displacement
  VECTOR zero(_DOFs);
  zero.setZero();
  _tetMesh.setDisplacement(zero);

  // get stiffness matrix at that state
  _tetMesh.computeFs();
  _tetMesh.computeSVDs();
  SPARSE_MATRIX K = _tetMesh.computeHyperelasticClampedHessian(_hyperelastic);

  // restore old state
  _tetMesh.setDisplacement(_temp);
  
  // build out the Rayleigh damping
  SPARSE_MATRIX C = _rayleighAlpha * _M;
  C += _rayleighBeta * K;

  return C;
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
// find the surface constraints that are separating
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool TIMESTEPPER::findSeparatingSurfaceConstraints(const VECTOR& unfiltered)
{
  bool changed = false;

  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    PLANE_CONSTRAINT& constraint = _planeConstraints[x];

    // is the vertex still inside the object? In that case, keep the constraint.
    const KINEMATIC_SHAPE* shape = constraint.shape;
    const int vertexID = constraint.vertexID;
    const REAL signedDistance = shape->signedDistance(_tetMesh.vertices()[vertexID]);

    //bool debug = constraint.vertexID == 0;
    bool debug = false;
    if (debug)
    {
      cout << " Constraint for vertex: " << constraint.vertexID << endl;
      cout << "   Signed distance:     " << signedDistance << endl;
    }

    // if the distance is outside and large, move on
    if (signedDistance > 1e-6) 
    {
      constraint.isSeparating = true;
      changed = true;
      if (debug) cout << " CONSTRAINT IS OUTSIDE " << endl;
      continue;
    }

    // what direction is the solution pointing in?
    const int vectorID = 3 * vertexID;
    VECTOR3 xDirection;
    xDirection[0] = unfiltered[vectorID];
    xDirection[1] = unfiltered[vectorID + 1];
    xDirection[2] = unfiltered[vectorID + 2];

    // make the test material agnostic and only look at the direction; then if it's a big force, 
    // the testing threshold won't get messed up later 
    if (xDirection.norm() > 1.0)
      xDirection.normalize();

    // what direction is the kinematic object's surface normal pointing in?
    VECTOR3 normal = shape->localNormalToWorld(constraint.localNormal);

    // what is the magnitude in the separation direction?
    const REAL separationMagnitude = xDirection.dot(normal);

    if (debug) cout << "  separation magnitude: " << separationMagnitude << endl;

    /*
    // is the velocity pulling away?
    VECTOR3 localVelocity = velocity(vertexID);
    const REAL velocitySeparation = localVelocity.dot(normal);
    */

    // are they the same? in that case, they're separating
    //
    // using a small epsilon threshold, especially for velocitySeparation
    // because otherwise there is jittering during resting contact
    if (separationMagnitude > 1e-6)   // velocity condition doesn't seem to ever get triggered
    //if (separationMagnitude > 1e-6 || velocitySeparation > 1e-6)
    //if (separationMagnitude > 1e-5 || velocitySeparation > 1e-5)
    //if (separationMagnitude > 0.0 || velocitySeparation > 0.0)
    //if (separationMagnitude > 1e-4 || velocitySeparation > 1e-4)  // jitter goes away for SNH lambda = 1000, but gets sticky
    //if (separationMagnitude > 9e-5 || velocitySeparation > 9e-5)
    {
      //cout << " Separation: " << separationMagnitude << " velocity: " << velocitySeparation << endl;
      constraint.isSeparating = true;
      changed = true;

      if (debug)
        cout << " SEPARATING" << endl;
    }
  }

  return changed;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// add a gravity body force to the simulation
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::addGravity(const VECTOR3& bodyForce)
{
  const vector<REAL>& oneRingVolumes = _tetMesh.restOneRingVolumes();

  for (int x = 0; x < _DOFs / 3; x++)
  {
    const VECTOR3 scaledForce = oneRingVolumes[x] * bodyForce;
    _externalForces[3 * x]     += scaledForce[0];
    _externalForces[3 * x + 1] += scaledForce[1];
    _externalForces[3 * x + 2] += scaledForce[2];
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// constrain surface nodes inside a kinematic body to move along with that body
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::attachKinematicSurfaceConstraints(const KINEMATIC_SHAPE* shape)
{
  // get all the nodes inside the shape
  const vector<VECTOR3>& vertices = _tetMesh.vertices();
  const vector<int>& surfaceVertices = _tetMesh.surfaceVertices();
  for (unsigned int x = 0; x < surfaceVertices.size(); x++)
  {
    int whichVertex = surfaceVertices[x];
    VECTOR3 v = vertices[whichVertex];

    // if it's not inside, move on
    if (!shape->inside(v)) continue;

    // if it's inside, get its local coordinates
    VECTOR3 local = shape->worldVertexToLocal(v);

    // record everything the solver will need later
    KINEMATIC_CONSTRAINT constraint;
    constraint.shape = shape;
    constraint.vertexID = whichVertex;
    constraint.localPosition = local;

    // rememeber the constraint for later
    _kinematicConstraints.push_back(constraint);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// constrain all nodes inside a kinematic body to move along with that body
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::attachKinematicConstraints(const KINEMATIC_SHAPE* shape)
{
  // get all the nodes inside the shape
  const vector<VECTOR3>& vertices = _tetMesh.vertices();
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    VECTOR3 v = vertices[x];

    // if it's not inside, move on
    if (!shape->inside(v)) continue;

    // if it's inside, get its local coordinates
    VECTOR3 local = shape->worldVertexToLocal(v);

    // record everything the solver will need later
    KINEMATIC_CONSTRAINT constraint;
    constraint.shape = shape;
    constraint.vertexID = x;
    constraint.localPosition = local;

    // rememeber the constraint for later
    _kinematicConstraints.push_back(constraint);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// which nodes are the constrained ones?
///////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> TIMESTEPPER::constrainedNodes() const
{
  // find the (unique) constrained nodes
  map<int, bool> isConstrained;
  for (unsigned int x = 0; x < _kinematicConstraints.size(); x++)
  {
    const int vertexID = _kinematicConstraints[x].vertexID;
    isConstrained[vertexID] = true;
  }

  // tape out the unique IDs
  vector<int> nodes;
  for (auto iter = isConstrained.begin(); iter != isConstrained.end(); iter++)
    nodes.push_back(iter->first);

  return nodes;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// add kinematic collision object to system
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::addKinematicCollisionObject(const KINEMATIC_SHAPE* shape)
{
  // make sure we didn't already add it
  for (unsigned int x = 0; x < _collisionObjects.size(); x++)
    if (shape == _collisionObjects[x])
      return;

  _collisionObjects.push_back(shape);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// find all the surface vertices that are in collision and create constraints
//
// this one is slightly different in QUASISTATIC, i.e. that one can't take
// velocity into account as a separation condition, so this function has
// been made virtual
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::findNewSurfaceConstraints(const bool verbose)
{
  const vector<VECTOR3> vertices = _tetMesh.vertices();
  const vector<int> surfaceVertices = _tetMesh.surfaceVertices();

  if (verbose)
    cout << " Currently tracking " << _planeConstraints.size() << " constraints " << endl;

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

      //bool debug = (vertexID == 0);
      bool debug = false;

      // if it's already in collision, skip it
      if (_inCollision[vertexID]) 
      {
        if (debug) cout << " vertex is already collision, moving on" << endl;
        continue;
      }

      // see if it's inside the shape
      const VECTOR3& vertex = vertices[vertexID];
      if (!shape->inside(vertex)) 
      {
        if (debug) cout << " vertex is not inside the shape, moving on" << endl;
        continue;
      }

      VECTOR3 closestPoint;
      VECTOR3 closestNormal;
      shape->getClosestPoint(vertex, closestPoint, closestNormal);
     
      // if the velocity is pulling away from the surface, don't constrain it
      VECTOR3 vertexVelocity = velocity(vertexID);
      VECTOR3 normal = shape->localNormalToWorld(closestNormal);
      const REAL velocitySeparation = vertexVelocity.dot(normal);

      if (debug)
      {
        cout << " velocity:   " << vertexVelocity.transpose() << endl;
        cout << " normal:     " << normal.transpose() << endl;
        cout << " separation: " << velocitySeparation << endl;
      }

      // comparing directly against zero here. Trying for a small
      // epsilon just induces sticking.
      //
      // Without this, objects will always stick to a surface after initially
      // sliding
      //
      // BDF-2 sticks unless -FLT_EPSILON is used, but other integrators seem okay
      //if (velocitySeparation >= 0)
      //if (velocitySeparation >= -1e-9)
      if (velocitySeparation >= -FLT_EPSILON)
        continue;

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

///////////////////////////////////////////////////////////////////////////////////////////////////////
// update the closest point positions on the surface constraints
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::updateSurfaceConstraints()
{
  const vector<VECTOR3> vertices = _tetMesh.vertices();
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    const KINEMATIC_SHAPE& shape = *_planeConstraints[x].shape;

    // get the new vertex position
    const int vertexID = _planeConstraints[x].vertexID;
    const VECTOR3& vertex = vertices[vertexID];

    // recompute closes point
    VECTOR3 closestPointLocal, normalLocal;
    shape.getClosestPoint(vertex, closestPointLocal, normalLocal);

    // store result
    _planeConstraints[x].localClosestPoint = closestPointLocal;
    _planeConstraints[x].localNormal = normalLocal;
  }
  // we're not checking whether it's still inside or separating here.
  // That will be handled by findSeparatingSurfaceConstraints
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// find all the constraints tagged for deletion and delete them
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::deleteSurfaceConstraints(const bool verbose)
{
  int totalDeleted = 0;

  // If any constraints were tagged for deletion last time, delete them now.
  //
  // That's right, I'm just building a whole new vector instead of deleting nodes 
  // from a linked list. If it's all too ugly for you, look away.
  // Like I said in the README.md, this library is not optimized yet.
  vector<PLANE_CONSTRAINT> constraints;
  for (unsigned int x = 0; x < _planeConstraints.size(); x++)
  {
    // if it's not separating, just store it
    if (!_planeConstraints[x].isSeparating)
      constraints.push_back(_planeConstraints[x]);
    // if we're deleting this, make sure this surface vertex isn't 
    // tagged as in collision anymore
    else
    {
      _inCollision[_planeConstraints[x].vertexID] = false;
      totalDeleted++;
    }
  }

  if (verbose)
    cout << " Total deleted: " << totalDeleted << endl;

  _planeConstraints = constraints;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// velocity at a specific vertex
///////////////////////////////////////////////////////////////////////////////////////////////////////
const VECTOR3 TIMESTEPPER::velocity(unsigned int index) const
{
  assert(index >= 0);
  assert(index < _velocity.size());
  VECTOR3 vertexVelocity;
  vertexVelocity[0] = _velocity[3 * index];
  vertexVelocity[1] = _velocity[3 * index + 1];
  vertexVelocity[2] = _velocity[3 * index + 2];

  return vertexVelocity;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// reset the Rayleigh damping constants
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::setRayeligh(const REAL alpha, const REAL beta)
{
  _rayleighAlpha = alpha;
  _rayleighBeta = beta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// do the collision detection, in anticipation of collision response
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::computeCollisionDetection()
{
  TIMER functionTimer(__FUNCTION__);

  // if the tet mesh has an AABB accelerator, refit it
  TET_MESH_FASTER* fast = dynamic_cast<TET_MESH_FASTER*>(&_tetMesh);
  if (fast != NULL)
    fast->refitAABB();

  // do the collision processing
  const REAL invDt = 1.0 / _dt;
  if (_vertexFaceSelfCollisionsOn)
  {
    // vertex-face collision detection
    _tetMesh.computeVertexFaceCollisions();

    // build out the vertex-face "collision tets"
    // TODO: this need to get cut down
    _tetMesh.buildVertexFaceCollisionTets(_velocity);
  }
  if (_edgeEdgeSelfCollisionsOn)
  {
    // edge-edge collision detection
    _tetMesh.computeEdgeEdgeCollisions();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute collision forces, add them to the forces and stiffness matrix
// R = forces, K = stiffness, C = damping
///////////////////////////////////////////////////////////////////////////////////////////////////////
void TIMESTEPPER::computeCollisionResponse(VECTOR& R, SPARSE_MATRIX& K, SPARSE_MATRIX& collisionC, const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  // build the collision forces and Hessians
  const int rank = R.size();
  VECTOR collisionForces(rank);
  SPARSE_MATRIX collisionK(rank, rank);
  //SPARSE_MATRIX collisionC(rank, rank);

  collisionForces.setZero();
  collisionK.setZero();
  collisionC.setZero();

  const REAL dampingBeta = _collisionDampingBeta;
  _tetMesh.setCollisionStiffness(_collisionStiffness);

  // vertex-face case
  VECTOR forcesVF;
  SPARSE_MATRIX hessianVF;
  if (_vertexFaceSelfCollisionsOn)
  {
    // get vertex-face collision forces and gradient
    forcesVF = _tetMesh.computeVertexFaceCollisionForces();
    hessianVF = _tetMesh.computeVertexFaceCollisionClampedHessian();
    
    collisionForces += forcesVF;
    collisionK += hessianVF;
    collisionC += dampingBeta * hessianVF;
  }

  // edge-edge case
  VECTOR forcesEE;
  SPARSE_MATRIX hessianEE;
  if (_edgeEdgeSelfCollisionsOn)
  {
    forcesEE = _tetMesh.computeEdgeEdgeCollisionForces();
    hessianEE = _tetMesh.computeEdgeEdgeCollisionClampedHessian();
    
    collisionForces += forcesEE;
    collisionK += hessianEE;
    collisionC += dampingBeta * hessianEE;
  }

#if VERY_VERBOSE
  if (verbose)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    if (_vertexFaceSelfCollisionsOn)
    {
      cout << " VF f: " << forcesVF.norm() << endl;
      cout << " VF K: " << hessianVF.norm() << endl;
      cout << " VF C: " << (dampingBeta * hessianVF).norm() << endl;
    }

    if (_edgeEdgeSelfCollisionsOn)
    {
      cout << " EE f: " << forcesEE.norm() << endl;
      cout << " EE K: " << hessianEE.norm() << endl;
      cout << " EE C: " << (dampingBeta * hessianEE).norm() << endl;
    }

    cout << " collision forces: " << collisionForces.norm() << endl;
    cout << " collision K: " << collisionK.norm() << endl;
    cout << " collision C: " << collisionC.norm() << endl;
  }
#endif

  // add self-collisions to both LHS and RHS
  if (_vertexFaceSelfCollisionsOn || _edgeEdgeSelfCollisionsOn)
  {
    R += collisionForces;
    K += collisionK;
    //C += collisionC;
  }
}

} // HOBAK
} // TIMESTEPPER
