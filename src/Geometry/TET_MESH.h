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
#ifndef TET_MESH_H
#define TET_MESH_H

#include "SETTINGS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Hyperelastic/Volume/VERTEX_FACE_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_COLLISION.h"
#include "Damping/Volume/DAMPING.h"

#include <map>
#include <vector>

namespace HOBAK {

using namespace std;

class TET_MESH
{
public:
  TET_MESH() = default;
  TET_MESH(const vector<VECTOR3>& restVertices, 
           const vector<VECTOR4I>& tets);
  virtual ~TET_MESH();

  // accessors
  const vector<VECTOR3>& vertices() const     { return _vertices; };
  vector<VECTOR3>& vertices()                 { return _vertices; };
  const vector<VECTOR3>& restVertices() const { return _restVertices; };
  vector<VECTOR3>& restVertices()             { return _restVertices; };
  const vector<VECTOR4I>& tets() const        { return _tets; };
  vector<VECTOR4I>& tets()                    { return _tets; };
  const vector<VECTOR4I>& vertexFaceCollisionTets() const    { return _vertexFaceCollisionTets; };
  vector<VECTOR4I>& vertexFaceCollisionTets()                { return _vertexFaceCollisionTets; };
  const vector<REAL>& restOneRingVolumes() const   { return _restOneRingVolumes; };
  vector<REAL>& restOneRingVolumes()               { return _restOneRingVolumes; };
  const VECTOR3& vertex(const int index) const     { return _vertices[index]; };
  VECTOR3& vertex(const int index)                 { return _vertices[index]; };
  const REAL & collisionEps() const                { return _collisionEps; };
  const VECTOR3& restVertex(const int index) const { return _restVertices[index]; };
  const VECTOR4I& tet(const int tetIndex) const    { return _tets[tetIndex]; };
  const vector<int>& surfaceTets() const           { return _surfaceTets; };
  const vector<int>& surfaceVertices() const       { return _surfaceVertices; };
  const vector<VECTOR3I>& surfaceTriangles() const { return _surfaceTriangles; };
  const vector<VECTOR2I>& surfaceEdges() const     { return _surfaceEdges; };
  const vector<pair<int,int> >& vertexFaceCollisions() const { return _vertexFaceCollisions; };
  const vector<pair<int,int> >& edgeEdgeCollisions() const   { return _edgeEdgeCollisions; };
  const vector<bool>& edgeEdgeIntersections() const          { return _edgeEdgeIntersections; };
  const vector<pair<VECTOR2,VECTOR2> >& edgeEdgeCoordinates() const   { return _edgeEdgeCoordinates; };
  const vector<REAL>& surfaceTriangleAreas() const { return _surfaceTriangleAreas; };
  const vector<VECTOR3I>& surfaceTriangleNeighbors() const { return _surfaceTriangleNeighbors; };

  int totalVertices() const { return _vertices.size(); };
  const int DOFs() const    { return _vertices.size() * 3; };

  // get deformation gradient, and its SVD
  MATRIX3 computeF(const int tetIndex) const;
  void computeFs();
  void computeFdots(const VECTOR& velocity);
  void computeSVDs();

  // get volume-weighted global translation
  VECTOR3 getTranslation() const;
  
  // get volume-weighted global translation, for the rest state
  VECTOR3 getRestTranslation() const;
  
  // get Procrustes-style global rotation
  MATRIX3 getRotation() const;

  // get the current displacement in vector form
  VECTOR getDisplacement() const;
  
  // set the vertex displacements to these values exactly
  void setDisplacement(const VECTOR& delta);

  // set the vertex positions directly exactly
  void setPositions(const VECTOR& positions);

  // add the following deltas to the positions
  void addDisplacement(const VECTOR& delta);

  // set collision eps to something new
  void setCollisionEps(const REAL& eps);
  void setCollisionStiffness(const REAL& stiffness);
 
  // set collision pairs (for replays)
  void setCollisionPairs(const vector<pair<int,int> >& vertexFace, 
                         const vector<pair<int,int> >& edgeEdge);

  // compute hyperelastic quantities
  REAL computeHyperelasticEnergy(const VOLUME::HYPERELASTIC& hyperelastic) const;
  VECTOR computeHyperelasticForces(const VOLUME::HYPERELASTIC& hyperelastic) const;
  virtual SPARSE_MATRIX computeHyperelasticClampedHessian(const VOLUME::HYPERELASTIC& hyperelastic) const;
  virtual SPARSE_MATRIX computeHyperelasticHessian(const VOLUME::HYPERELASTIC& hyperelastic) const;

  // compute damping quantities
  VECTOR computeDampingForces(const VOLUME::DAMPING& damping) const;
  virtual SPARSE_MATRIX computeDampingHessian(const VOLUME::DAMPING& damping) const;

  // compute x-based collision quantities
  VECTOR computeVertexFaceCollisionForces() const;
  SPARSE_MATRIX computeVertexFaceCollisionClampedHessian() const;
  REAL computeEdgeEdgeCollisionEnergy() const;
  VECTOR computeEdgeEdgeCollisionForces() const;
  SPARSE_MATRIX computeEdgeEdgeCollisionClampedHessian() const;

  // compuate elastic and damping forces at the same time
  virtual VECTOR computeInternalForces(const VOLUME::HYPERELASTIC& hyperelastic,
                                       const VOLUME::DAMPING& damping) const;

  // get the bounding box for the current mesh
  void getBoundingBox(VECTOR3& mins, VECTOR3& maxs) const;

  // find all the vertex-face collision pairs, using the InFaceRegion test
  virtual void computeVertexFaceCollisions();

  // find all the edge-edge collision pairs
  virtual void computeEdgeEdgeCollisions();

  // debug edge-edge collisions, load up some specific pairs
  //void computeEdgeEdgeCollisionsDebug();

  // based on vertex-face collision pairs, build "collision tets"
  void buildVertexFaceCollisionTets(const VECTOR& velocity);
  
  // write out the surface to OBJ triangle mesh
  static bool writeSurfaceToObj(const string& filename, const TET_MESH& tetMesh);

  // read in OBJ-style tet mesh file
  static bool readTobjFile(const string& filename, 
                           vector<VECTOR3>& vertices, vector<VECTOR4I>& tets);

  // write out OBJ-style tet mesh file, the bool tells whether to write the
  // rest or deformed vertices
  static bool writeTobjFile(const string& filename, const TET_MESH& tetMesh, const bool restVertices);

  // normalize vertices so that they're in a unit box, centered at (0.5, 0.5, 0.5)
  static vector<VECTOR3> normalizeVertices(const vector<VECTOR3>& vertices);

  // compute distance between a point and triangle
  static REAL pointTriangleDistance(const VECTOR3& v0, const VECTOR3& v1, 
                                    const VECTOR3& v2, const VECTOR3& v);

  // see if the projection of v onto the plane of v0,v1,v2 is inside the triangle
  // formed by v0,v1,v2
  static bool pointProjectsInsideTriangle(const VECTOR3& v0, const VECTOR3& v1, 
                                          const VECTOR3& v2, const VECTOR3& v);

  // compute the dihedral angle between two surface faces
  REAL surfaceFaceDihedralAngle(const int surfaceID0, const int surfaceID1) const;

protected:
  // compute the volume of a tet
  static REAL computeTetVolume(const vector<VECTOR3>& tetVertices);

  // compute volumes for tets -- works for rest and deformed, just
  // pass it _restVertices or _vertices
  vector<REAL> computeTetVolumes(const vector<VECTOR3>& vertices);

  // compute volumes for each vertex one-ring -- works for rest
  // and deformed, just pass it _restVertices or _vertices
  vector<REAL> computeOneRingVolumes(const vector<VECTOR3>& vertices);

  // compute material inverses for deformation gradient
  vector<MATRIX3> computeDmInvs();

  // compute change-of-basis from deformation gradient F to positions, x
  vector<MATRIX9x12> computePFpxs();

  // find what's on the surface
  void computeSurfaceVertices();
  void computeSurfaceTriangles();
  void computeSurfaceEdges();
  void computeSurfaceAreas();
  void computeSurfaceTriangleNeighbors();
  void computeSurfaceEdgeTriangleNeighbors();

  // (DEACTIVATED) find out how close all the edges are initially
  void computeEdgeEdgeRestDistance();

  // compute a triangle area
  static REAL triangleArea(const vector<VECTOR3>& triangle);

  // get the normal to a plane, specified by three points
  static VECTOR3 planeNormal(const vector<VECTOR3>& plane);

  // project point onto plane, specific by three points
  static VECTOR3 pointPlaneProjection(const vector<VECTOR3>& plane, const VECTOR3& point);

  // see if the vertex is inside the collision cell described in 
  // Chapter 11: Collision Processing, [Kim and Eberle 2020]
  bool insideCollisionCell(const int surfaceTriangleID, const VECTOR3& vertex);

  // compute distance to collision cell wall, where positive means inside
  // and negative means outside
  REAL distanceToCollisionCellWall(const int surfaceTriangleID, const VECTOR3& vertex);

  // compute whether one vertex is inside the vertex one right of another
  void computeSurfaceVertexOneRings();

  // are these two surface triangles neighbors?
  bool areSurfaceTriangleNeighbors(const int id0, const int id1) const;

  // build a consistent tet/flap ordering from two surface triangles
  VECTOR4I buildSurfaceFlap(const int surfaceID0, const int surfaceID1) const;

  // compute the normal of the surface triangle at _surfaceTriangles[triangleID];
  VECTOR3 surfaceTriangleNormal(const int triangleID) const;

  // see if a current surface triangle has been crushed to degeneracy
  bool surfaceTriangleIsDegenerate(const int surfaceTriangleID);

  // compute which vertices are attached to inverted tets
  void computeInvertedVertices();

  // the core geometry
  vector<VECTOR3>    _vertices;
  vector<VECTOR3>    _restVertices;
  vector<VECTOR4I>   _tets;

  // volumes, computed by computeTetVolumes and computeOneRingVolumes
  //
  // TODO: why aren't these just VECTORs?
  vector<REAL>  _restTetVolumes;
  vector<REAL>  _restOneRingVolumes;
  vector<REAL>  _restOneRingAreas;
  VECTOR        _restEdgeAreas;

  // support for computing deformation gradient F
  vector<MATRIX3>    _DmInvs;

  // change-of-basis to go from deformation gradient (F) to positions (x)
  vector<MATRIX9x12> _pFpxs;

  // deformation gradients, and their SVDs
  vector<MATRIX3> _Fs;
  vector<MATRIX3> _Us;
  vector<VECTOR3> _Sigmas;
  vector<MATRIX3> _Vs;

  // velocity gradients
  vector<MATRIX3> _Fdots;

  // list of tets that are on the surface
  vector<int> _surfaceTets;

  // list of triangles that are on the surface
  // each triplet is ordered counter-clockwise, facing outwards
  // the VECTOR3I indexes into _vertices
  vector<VECTOR3I> _surfaceTriangles;
  vector<REAL> _surfaceTriangleAreas;

  // for each surface triangle, what's the index of the neighboring triangles?
  vector<VECTOR3I> _surfaceTriangleNeighbors;

  // list of edges on the surface
  // each pair is in sorted order, and index into _vertices
  vector<VECTOR2I> _surfaceEdges;

  // list of vertices that are on the surface
  // indexes into _vertices
  vector<int> _surfaceVertices;

  // for each _surfaceEdges, what are the one or two neighboring triangles
  // in _surfaceTriangles?
  vector<VECTOR2I> _surfaceEdgeTriangleNeighbors;

  // for each pair of _surfaceEdges, what _collisionEps should we use? If they started
  // out closer than _collisionEps, then we need to set a smaller tolerance.
  //
  // the entries in the pair<int,int> are:
  //   unsigned int flat = edge[0] + edge[1] * _surfaceEdges.size();
  map<pair<unsigned int, unsigned int>, REAL> _edgeEdgeRestDistance;

  // how close is considered to be in collision?
  REAL _collisionEps;

  // list of vertex-face collisions
  // first indexes into _vertices
  // second indexes into _surfaceTriangles
  vector<pair<int, int> > _vertexFaceCollisions;

  // list of edge-edge collision indices
  // first indexes into _surfaceEdges
  // second indexes into _surfaceEdges
  vector<pair<int, int> > _edgeEdgeCollisions;

  // interpolation coordinates for edge-edge collisions
  vector<pair<VECTOR2, VECTOR2> > _edgeEdgeCoordinates;
  
  // are the edge-edge collisions still separate, or is there already a face-edge intersection?
  vector<bool> _edgeEdgeIntersections;

  // list of collisionEps for each _edgeEdgeIntersections. Usually this will be the default,
  // but if it started out closer than that, we use that instead
  //vector<REAL> _edgeEdgeCollisionEps;

  // list of "collision tets" formed by vertex-face pairs
  vector<VECTOR4I> _vertexFaceCollisionTets;

  // list of "collision tets" formed by edge-edge pairs
  //vector<VECTOR4I> _edgeEdgeCollisionsTets;

  // DEBUG: see if the collision tet exists already
  //map<pair<int, int>, int> _vertexFaceCollisionTetsHash;
  
  // scaling term for vertex-face collision forces
  vector<REAL> _vertexFaceCollisionAreas;

  // scaling term for edge-edge collision forces
  vector<REAL> _edgeEdgeCollisionAreas;

  // convert tet mesh vertexID into a surface mesh vertexID
  // convert into into _vertices into index into _surfaceVertices
  map<int, int> _volumeToSurfaceID;

  // constitutive model for collisions
  VOLUME::HYPERELASTIC* _collisionMaterial;

  // have your computed the SVDs since the last time you computed F?
  bool _svdsComputed;

  // see if two indices in _vertices (in sorted order)
  // are within the one ring of each other
  map<pair<int,int>, bool> _insideSurfaceVertexOneRing;

  // which vertex-face collision force are we using?
  VOLUME::VERTEX_FACE_COLLISION* _vertexFaceEnergy;
  
  // which edge-edge collision force are we using?
  VOLUME::EDGE_COLLISION* _edgeEdgeEnergy;

  // which vertices are inverted?
  vector<bool> _invertedVertices;
};

}

#endif
