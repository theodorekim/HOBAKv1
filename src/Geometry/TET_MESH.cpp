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
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include <iostream>
#include <float.h>
#include "util/TIMER.h"
#include "util/COLLISION_UTIL.h"
#include "util/MATRIX_UTIL.h"
#include "LINE_INTERSECT.h"
#include "TET_MESH.h"

// DEBUG: only here specifically to debug collisions
#include "Hyperelastic/Volume/ARAP.h"
#include "Hyperelastic/Volume/VERTEX_FACE_COLLISION.h"
#include "Hyperelastic/Volume/VERTEX_FACE_SQRT_COLLISION.h"
#include "Hyperelastic/Volume/MCADAMS_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_SQRT_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_HYBRID_COLLISION.h"

// add back the bug where the force points in the wrong direction when two
// colliding edges are on triangles that are intersecting
#define ADD_EDGE_EDGE_PENETRATION_BUG 0

namespace HOBAK {

using namespace std;

TET_MESH::TET_MESH(const vector<VECTOR3>& restVertices, 
                   const vector<VECTOR4I>& tets) :
    _vertices(restVertices),
    _restVertices(restVertices),
    _tets(tets)
{
  _restTetVolumes = computeTetVolumes(_restVertices);
  _restOneRingVolumes = computeOneRingVolumes(_restVertices);
  _DmInvs = computeDmInvs();
  _pFpxs = computePFpxs();

  const int totalTets = _tets.size();
  _Fs.resize(totalTets);
  _Us.resize(totalTets);
  _Sigmas.resize(totalTets);
  _Vs.resize(totalTets);
  _Fdots.resize(totalTets);

  computeSurfaceTriangles();
  computeSurfaceVertices();
  computeSurfaceEdges();

  computeSurfaceAreas();
  computeSurfaceTriangleNeighbors();
  computeSurfaceEdgeTriangleNeighbors();

  // one centimeter
  _collisionEps = 0.01;
  // two centimeters, one seems to get into trouble without CCD
  //_collisionEps = 0.02;
  _collisionMaterial = NULL;  // this is still experimental
  _svdsComputed = false;

  // this gets overwritten by TIMESTEPPER every step, so a dummy is fine
  REAL stiffness = 1000.0;
  
  // must be called after _collisionEps has been set
  //computeEdgeEdgeRestDistance();
  
  // store which surface vertices are within the one rings of
  // each other
  computeSurfaceVertexOneRings();

  // if you want to try out the McAdams energy, 
  // here's the place to swap it in
  //_vertexFaceEnergy = new VOLUME::VERTEX_FACE_COLLISION(stiffness, _collisionEps); // default
  //_vertexFaceEnergy = new VOLUME::MCADAMS_COLLISION(stiffness, _collisionEps);
 
  // if you want to try out a different edge-edge energy, 
  // here's the place to swap it in
  //_edgeEdgeEnergy = new VOLUME::EDGE_COLLISION(stiffness, _collisionEps);
  //_edgeEdgeEnergy = new VOLUME::EDGE_HYBRID_COLLISION(stiffness, _collisionEps);
 
  // preferred, verified on the bunny drop scene 
  _vertexFaceEnergy = new VOLUME::VERTEX_FACE_SQRT_COLLISION(stiffness, _collisionEps); // default
  _edgeEdgeEnergy = new VOLUME::EDGE_SQRT_COLLISION(stiffness, _collisionEps); // default
  
  // verified on the bunny drop scene 
  //_vertexFaceEnergy = new VOLUME::MCADAMS_COLLISION(stiffness, _collisionEps);
  //_edgeEdgeEnergy = new VOLUME::EDGE_SQRT_COLLISION(stiffness, _collisionEps);
  
  // verified on the bunny drop scene 
  //_vertexFaceEnergy = new VOLUME::MCADAMS_COLLISION(stiffness, _collisionEps);
  //_edgeEdgeEnergy = new VOLUME::EDGE_HYBRID_COLLISION(stiffness, _collisionEps);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TET_MESH::~TET_MESH()
{
  delete _vertexFaceEnergy;
  delete _edgeEdgeEnergy;
}

///////////////////////////////////////////////////////////////////////
// compute volumes for tets -- works for rest and deformed, just
// pass it _restVertices or _vertices
///////////////////////////////////////////////////////////////////////
vector<REAL> TET_MESH::computeTetVolumes(const vector<VECTOR3>& vertices)
{
  vector<REAL> tetVolumes(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    vector<VECTOR3> tetVertices;
    for (unsigned int j = 0; j < 4; j++)
      tetVertices.push_back(_vertices[tet[j]]);
    tetVolumes[i] = computeTetVolume(tetVertices);

    if (tetVolumes[i] < 0.0)
    {
      cout << " ERROR: Bad rest volume found: " << tetVolumes[i] << endl;
      // tetVolumes[i] = 0.0;
    }
    assert(tetVolumes[i] >= 0.0);
  }
  return tetVolumes;
}

///////////////////////////////////////////////////////////////////////
// compute the volume of a tet
///////////////////////////////////////////////////////////////////////
REAL TET_MESH::computeTetVolume(const vector<VECTOR3>& tetVertices)
{
  const VECTOR3 diff1 = tetVertices[1] - tetVertices[0];
  const VECTOR3 diff2 = tetVertices[2] - tetVertices[0];
  const VECTOR3 diff3 = tetVertices[3] - tetVertices[0];
  return diff3.dot((diff1).cross(diff2)) / 6.0;
}

///////////////////////////////////////////////////////////////////////
// compute rest volumes for each vertex one-ring
///////////////////////////////////////////////////////////////////////
vector<REAL> TET_MESH::computeOneRingVolumes(const vector<VECTOR3>& vertices)
{
  const vector<REAL> tetVolumes = computeTetVolumes(vertices);
  unsigned int size = vertices.size();

  vector<REAL> oneRingVolumes(size);
  for (unsigned int x = 0; x < size; x++)
    oneRingVolumes[x] = 0.0;

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    const REAL quarter = 0.25 * tetVolumes[x];
    for (int y = 0; y < 4; y++)
      oneRingVolumes[_tets[x][y]] += quarter;
  }

  return oneRingVolumes;
}

///////////////////////////////////////////////////////////////////////
// compute material inverses for deformation gradient
///////////////////////////////////////////////////////////////////////
vector<MATRIX3> TET_MESH::computeDmInvs()
{
  vector<MATRIX3> DmInvs(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    MATRIX3 Dm;
    Dm.col(0) = _vertices[tet[1]] - _vertices[tet[0]];
    Dm.col(1) = _vertices[tet[2]] - _vertices[tet[0]];
    Dm.col(2) = _vertices[tet[3]] - _vertices[tet[0]];
    DmInvs[i] = Dm.inverse();
  }
  return DmInvs;
}

///////////////////////////////////////////////////////////////////////
// compute change-of-basis from deformation gradient F to positions x
// for a single DmInv
///////////////////////////////////////////////////////////////////////
static MATRIX9x12 computePFpx(const MATRIX3 &DmInv)
{
  const REAL m = DmInv(0, 0);
  const REAL n = DmInv(0, 1);
  const REAL o = DmInv(0, 2);
  const REAL p = DmInv(1, 0);
  const REAL q = DmInv(1, 1);
  const REAL r = DmInv(1, 2);
  const REAL s = DmInv(2, 0);
  const REAL t = DmInv(2, 1);
  const REAL u = DmInv(2, 2);

  const REAL t1 = -m - p - s;
  const REAL t2 = -n - q - t;
  const REAL t3 = -o - r - u;

  MATRIX9x12 PFPu = MATRIX9x12::Zero();
  PFPu(0, 0)  = t1;
  PFPu(0, 3)  = m;
  PFPu(0, 6)  = p;
  PFPu(0, 9)  = s;
  PFPu(1, 1)  = t1;
  PFPu(1, 4)  = m;
  PFPu(1, 7)  = p;
  PFPu(1, 10) = s;
  PFPu(2, 2)  = t1;
  PFPu(2, 5)  = m;
  PFPu(2, 8)  = p;
  PFPu(2, 11) = s;
  PFPu(3, 0)  = t2;
  PFPu(3, 3)  = n;
  PFPu(3, 6)  = q;
  PFPu(3, 9)  = t;
  PFPu(4, 1)  = t2;
  PFPu(4, 4)  = n;
  PFPu(4, 7)  = q;
  PFPu(4, 10) = t;
  PFPu(5, 2)  = t2;
  PFPu(5, 5)  = n;
  PFPu(5, 8)  = q;
  PFPu(5, 11) = t;
  PFPu(6, 0)  = t3;
  PFPu(6, 3)  = o;
  PFPu(6, 6)  = r;
  PFPu(6, 9)  = u;
  PFPu(7, 1)  = t3;
  PFPu(7, 4)  = o;
  PFPu(7, 7)  = r;
  PFPu(7, 10) = u;
  PFPu(8, 2)  = t3;
  PFPu(8, 5)  = o;
  PFPu(8, 8)  = r;
  PFPu(8, 11) = u;

  return PFPu;
}

///////////////////////////////////////////////////////////////////////
// compute change-of-basis from deformation gradient F to positions, x
///////////////////////////////////////////////////////////////////////
vector<MATRIX9x12> TET_MESH::computePFpxs()
{
  vector<MATRIX9x12> pFpxs(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
    pFpxs[i] = computePFpx(_DmInvs[i]);
  return pFpxs;
}

///////////////////////////////////////////////////////////////////////
// used by computeSurfaceTriangles as a comparator between two
// triangles
///////////////////////////////////////////////////////////////////////
struct triangleCompare
{
  bool operator()(const VECTOR3I &a, const VECTOR3I &b) const
  {
    if (a[0] < b[0])
      return true;
    if (a[0] > b[0])
      return false;

    if (a[1] < b[1])
      return true;
    if (a[1] > b[1])
      return false;

    if (a[2] < b[2])
      return true;
    if (a[2] > b[2])
      return false;
    return false;
  }
};

///////////////////////////////////////////////////////////////////////
// find which triangles are on the surface
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceTriangles()
{
  map<VECTOR3I, int, triangleCompare> faceCounts;

  // for each tet, add its faces to the face count
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    VECTOR4I t = _tets[x];

    VECTOR3I faces[4];
    faces[0] << t[0], t[1], t[3];
    faces[1] << t[0], t[2], t[1];
    faces[2] << t[0], t[3], t[2];
    faces[3] << t[1], t[2], t[3];

    for (int y = 0; y < 4; y++)
      std::sort(faces[y].data(), faces[y].data() + faces[y].size());

    for (int y = 0; y < 4; y++)
      faceCounts[faces[y]]++;
  }

  // go back through the tets, if any of its faces have a count less than
  // 2, then it must be because it faces outside
  _surfaceTriangles.clear();
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    VECTOR4I t = _tets[x];

    VECTOR3I faces[4];

    // these are consistently ordered counter-clockwise
    faces[0] << t[0], t[1], t[3];
    faces[1] << t[0], t[2], t[1];
    faces[2] << t[0], t[3], t[2];
    faces[3] << t[1], t[2], t[3];

    VECTOR3I facesSorted[4];

    // make a sorted copy, but keep the original around for rendering
    for (int y = 0; y < 4; y++)
    {
      facesSorted[y] = faces[y];
      std::sort(facesSorted[y].data(), facesSorted[y].data() + facesSorted[y].size());
    }

    // see which faces don't have a dual
    for (int y = 0; y < 4; y++)
    {
      if (faceCounts[facesSorted[y]] < 2)
        _surfaceTriangles.push_back(faces[y]);
    }
  }

  cout << " Found " << _surfaceTriangles.size() << " surface triangles out of " 
       << _tets.size() * 4 << " possible. " << endl;
}

///////////////////////////////////////////////////////////////////////
// for each edge on the surface, what's the index of the neighboring 
// triangles?
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceEdgeTriangleNeighbors()
{
  // translate the VEC2I into a index
  map<pair<int,int>, int> edgeToIndex;
  for (unsigned int x = 0; x < _surfaceEdges.size(); x++)
  {
    pair<int,int> toHash;
    toHash.first  = _surfaceEdges[x][0];
    toHash.second = _surfaceEdges[x][1];

    if (toHash.first > toHash.second)
    {
      int temp = toHash.first;
      toHash.first = toHash.second;
      toHash.second = temp;
    }
    edgeToIndex[toHash] = x;
  }

  // look up the edges of each surface face, tabulate the adjacent
  // triangles
  vector<vector<int> > faceHash(_surfaceEdges.size());
  for (unsigned int i = 0; i < _surfaceTriangles.size(); i++)
  {
    const VECTOR3I t = _surfaceTriangles[i];

    // store each edge as pair
    pair<int, int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // store the edge
      assert(edgeToIndex.find(edge) != edgeToIndex.end());
      int edgeIndex = edgeToIndex[edge];
      faceHash[edgeIndex].push_back(i);
    }
  }

  // store the final results
  _surfaceEdgeTriangleNeighbors.resize(_surfaceEdges.size());
  for (unsigned int i = 0; i < _surfaceEdges.size(); i++)
  {
    _surfaceEdgeTriangleNeighbors[i][0] = -1;
    _surfaceEdgeTriangleNeighbors[i][1] = -1;

    assert(faceHash[i].size() > 0);
    _surfaceEdgeTriangleNeighbors[i][0] = faceHash[i][0];

    if (faceHash[i].size() == 2)
      _surfaceEdgeTriangleNeighbors[i][1] = faceHash[i][1];
  }
}

///////////////////////////////////////////////////////////////////////
// for each surface triangle, what's the index of the neighboring 
// triangles?
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceTriangleNeighbors()
{
  multimap<pair<int, int>, unsigned int> edgeNeighboringTriangles;

  // hash all the edges from each surface triangle
  for (unsigned int i = 0; i < _surfaceTriangles.size(); i++)
  {
    const VECTOR3I t = _surfaceTriangles[i];

    // store each edge as pair
    pair<int, int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // hash it
      pair<pair<int, int>, unsigned int> hash(edge, i);
      edgeNeighboringTriangles.insert(hash);
    }
  }
 
  // get the other edge that wasn't the current one 
  _surfaceTriangleNeighbors.clear();
  for (unsigned int i = 0; i < _surfaceTriangles.size(); i++)
  {
    const VECTOR3I t = _surfaceTriangles[i];

    // store results here
    VECTOR3I neighbors(-1,-1,-1);

    // reconstruct the edge again
    pair<int,int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // find the matching triangles
      auto range = edgeNeighboringTriangles.equal_range(edge);
      for (auto it = range.first; it != range.second; it++)
      {
        if (it->second != i)
          neighbors[j] = it->second;
      }
    }

    // store the neighbors
    _surfaceTriangleNeighbors.push_back(neighbors);
  }
}

///////////////////////////////////////////////////////////////////////
// compute surface areas for collision weights
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceAreas()
{
  // compute the areas
  _surfaceTriangleAreas.clear();
  for (unsigned int x = 0; x < _surfaceTriangles.size(); x++)
  {
    vector<VECTOR3> vertices(3);
    vertices[0] = _restVertices[_surfaceTriangles[x][0]];
    vertices[1] = _restVertices[_surfaceTriangles[x][1]];
    vertices[2] = _restVertices[_surfaceTriangles[x][2]];

    _surfaceTriangleAreas.push_back(triangleArea(vertices));
  }
  
  // compute the one-ring areas
  assert(_surfaceVertices.size() != 0);
  _restOneRingAreas.resize(_surfaceVertices.size());
  for (unsigned int x = 0; x < _restOneRingAreas.size(); x++)
    _restOneRingAreas[x] = 0;
  for (unsigned int x = 0; x < _surfaceTriangles.size(); x++)
  {
    assert(x < _surfaceTriangleAreas.size());
    assert(x < _surfaceTriangles.size());
    const REAL& area = _surfaceTriangleAreas[x];
    const VECTOR3I& triangle = _surfaceTriangles[x];

    for (int y = 0; y < 3; y++)
    {
      const int surfaceID = _volumeToSurfaceID[triangle[y]];
      assert(surfaceID < (int)_restOneRingAreas.size());
      _restOneRingAreas[surfaceID] += (1.0 / 3.0) * area;
    }
  }

  // build a mapping from edge index pairs to _surfaceEdges
  map<pair<int, int>, int> edgeHash;
  for (unsigned int x = 0; x < _surfaceEdges.size(); x++)
  {
    pair<int,int> edge(_surfaceEdges[x][0], _surfaceEdges[x][1]);
    edgeHash[edge] = x;
  }

  // compute the edge areas
  assert(_surfaceEdges.size() != 0);
  _restEdgeAreas.resize(_surfaceEdges.size());
  _restEdgeAreas.setZero();
  for (unsigned int x = 0; x < _surfaceTriangles.size(); x++)
  {
    // build each edge
    for (int y = 0; y < 3; y++)
    {
      pair<int,int> edge(_surfaceTriangles[x][y],
                         _surfaceTriangles[x][(y + 1) % 3]);

      // swap them to the order the hash expects
      if (edge.first > edge.second)
      {
        const REAL temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      const int edgeIndex = edgeHash[edge];
      assert(edgeIndex >= 0);
      assert(edgeIndex < _restEdgeAreas.size());
      _restEdgeAreas[edgeIndex] += _surfaceTriangleAreas[x] / 3.0;
    }
  }
}

///////////////////////////////////////////////////////////////////////
// find which vertices are on the surface
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceVertices()
{
  if (_surfaceTriangles.size() == 0)
    computeSurfaceTriangles();

  // hash them all out
  std::map<int, bool> foundVertices;
  for (unsigned int x = 0; x < _surfaceTriangles.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++)
      foundVertices[_surfaceTriangles[x][y]] = true;
  }

  // serialize
  _surfaceVertices.clear();
  for (auto iter = foundVertices.begin(); iter != foundVertices.end(); iter++)
    _surfaceVertices.push_back(iter->first);

  // compute the reverse lookup
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
    _volumeToSurfaceID[_surfaceVertices[x]] = x;

  cout << " Found " << _surfaceVertices.size() << " vertices on the surface " << endl;
}

///////////////////////////////////////////////////////////////////////
// find which edges are on the surface
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceEdges()
{
  // you did call this first, right?
  if (_surfaceTriangles.size() == 0)
    computeSurfaceTriangles();

  // hash all the edges, so we don't store any repeats
  map<pair<int, int>, bool> edgeHash;
  for (unsigned int x = 0; x < _surfaceTriangles.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++)
    {
      const int v0 = _surfaceTriangles[x][y];
      const int v1 = _surfaceTriangles[x][(y + 1) % 3];

      // store them in sorted order
      pair<int, int> edge;
      if (v0 > v1)
      {
        edge.first = v1;
        edge.second = v0;
      }
      else
      {
        edge.first = v0;
        edge.second = v1;
      }

      // hash it out
      edgeHash[edge] = true;
    }
  }

  // store all the unique hashes
  _surfaceEdges.clear();
  for (auto iter = edgeHash.begin(); iter != edgeHash.end(); iter++)
  {
    const pair<int,int> e = iter->first;
    const VECTOR2I edge(e.first, e.second);
    _surfaceEdges.push_back(edge);
  }

  cout << " Found " << _surfaceEdges.size() << " edges on the surface " << endl;
}

///////////////////////////////////////////////////////////////////////
// get all the deformation gradients
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeFs()
{
  TIMER functionTimer(__FUNCTION__);
  assert(_Fs.size() == _tets.size());

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _tets.size(); x++)
    _Fs[x] = computeF(x);

  // SVDs are now stale
  _svdsComputed = false;
}

///////////////////////////////////////////////////////////////////////
// get all the velocity gradients
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeFdots(const VECTOR& velocity)
{
  TIMER functionTimer(__FUNCTION__);
  assert(_Fs.size() == _tets.size());

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    const VECTOR4I& tet = _tets[x];
    VECTOR3 v[4];
    for (int y = 0; y < 4; y++)
    {
      v[y][0] = velocity[3 * tet[y]];
      v[y][1] = velocity[3 * tet[y] + 1];
      v[y][2] = velocity[3 * tet[y] + 2];
    }

    MATRIX3 V;
    V.col(0) = v[1] - v[0];
    V.col(1) = v[2] - v[0];
    V.col(2) = v[3] - v[0];
    _Fdots[x] = V * _DmInvs[x];
  }
}

///////////////////////////////////////////////////////////////////////
// get the SVD of all the deformation gradients
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeSVDs()
{
  TIMER functionTimer(__FUNCTION__);
  assert(_Us.size() == _tets.size());
  assert(_Sigmas.size() == _tets.size());
  assert(_Vs.size() == _tets.size());
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _tets.size(); x++)
    svd_rv(_Fs[x], _Us[x], _Sigmas[x], _Vs[x]);

  _svdsComputed = true;
}

///////////////////////////////////////////////////////////////////////
// get deformation gradient
///////////////////////////////////////////////////////////////////////
MATRIX3 TET_MESH::computeF(const int tetIndex) const
{
  const VECTOR4I& tet = _tets[tetIndex];
  MATRIX3 Ds;
  Ds.col(0) = _vertices[tet[1]] - _vertices[tet[0]];
  Ds.col(1) = _vertices[tet[2]] - _vertices[tet[0]];
  Ds.col(2) = _vertices[tet[3]] - _vertices[tet[0]];
  return Ds * _DmInvs[tetIndex];
}

///////////////////////////////////////////////////////////////////////
// get hyperelastic energy over entire mesh
///////////////////////////////////////////////////////////////////////
REAL TET_MESH::computeHyperelasticEnergy(const VOLUME::HYPERELASTIC& hyperelastic) const
{
  assert(_tets.size() == _restTetVolumes.size());

  VECTOR tetEnergies(_tets.size());
  for (int tetIndex = 0; tetIndex < int(_tets.size()); tetIndex++)
  {
    const MATRIX3 F = computeF(tetIndex);
    tetEnergies[tetIndex] = _restTetVolumes[tetIndex] * hyperelastic.psi(F);
  }

  return tetEnergies.sum();
}

///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the hyperelastic force
///////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::computeHyperelasticForces(const VOLUME::HYPERELASTIC& hyperelastic) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR12> perElementForces(_tets.size());
  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const MATRIX3& F = _Fs[tetIndex];
    const MATRIX3 PK1 = hyperelastic.PK1(F);
    const VECTOR12 forceDensity = _pFpxs[tetIndex].transpose() * flatten(PK1);
    const VECTOR12 force = -_restTetVolumes[tetIndex] * forceDensity;
    perElementForces[tetIndex] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const VECTOR4I& tet = _tets[tetIndex];
    const VECTOR12& tetForce = perElementForces[tetIndex];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the damping force
///////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::computeDampingForces(const VOLUME::DAMPING& damping) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR12> perElementForces(_tets.size());
  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const MATRIX3& F = _Fs[tetIndex];
    const MATRIX3& Fdot = _Fdots[tetIndex];
    const MATRIX3 PK1 = damping.PK1(F, Fdot);
    const VECTOR12 forceDensity = _pFpxs[tetIndex].transpose() * flatten(PK1);
    const VECTOR12 force = -_restTetVolumes[tetIndex] * forceDensity;
    perElementForces[tetIndex] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const VECTOR4I& tet = _tets[tetIndex];
    const VECTOR12& tetForce = perElementForces[tetIndex];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the elastic and damping forces
// at the same time
///////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::computeInternalForces(const VOLUME::HYPERELASTIC& hyperelastic,
                                       const VOLUME::DAMPING& damping) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR12> perElementForces(_tets.size());
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const MATRIX3& U       = _Us[tetIndex];
    const MATRIX3& V       = _Vs[tetIndex];
    const VECTOR3& Sigma   = _Sigmas[tetIndex];
    const MATRIX3& F = _Fs[tetIndex];
    const MATRIX3& Fdot = _Fdots[tetIndex];

    const MATRIX3 elasticPK1 = hyperelastic.PK1(U, Sigma, V);
    const MATRIX3 dampingPK1 = damping.PK1(F, Fdot);
    const VECTOR12 forceDensity = _pFpxs[tetIndex].transpose() * flatten(elasticPK1 + dampingPK1);
    const VECTOR12 force = -_restTetVolumes[tetIndex] * forceDensity;
    perElementForces[tetIndex] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const VECTOR4I& tet = _tets[tetIndex];
    const VECTOR12& tetForce = perElementForces[tetIndex];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the damping gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_MESH::computeDampingHessian(const VOLUME::DAMPING& damping) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perElementHessians(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const MATRIX3& F       = _Fs[i];
    const MATRIX3& Fdot    = _Fdots[i];
    const MATRIX9x12& pFpx = _pFpxs[i];
    const MATRIX9 hessian  = -_restTetVolumes[i] * damping.hessian(F, Fdot);
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_MESH::computeHyperelasticHessian(const VOLUME::HYPERELASTIC& hyperelastic) const
{
  vector<MATRIX12> perElementHessians(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const MATRIX3& F       = _Fs[i];
    const MATRIX9x12& pFpx = _pFpxs[i];
    const MATRIX9 hessian  = -_restTetVolumes[i] * hyperelastic.hessian(F);
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_MESH::computeHyperelasticClampedHessian(const VOLUME::HYPERELASTIC& hyperelastic) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perElementHessians(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const MATRIX3& F       = _Fs[i];
    const MATRIX9x12& pFpx = _pFpxs[i];
    const MATRIX9 hessian  = -_restTetVolumes[i] * hyperelastic.clampedHessian(F);
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

/*
///////////////////////////////////////////////////////////////////////
// Use the material and damping Hessian to compute the force gradient
// this is all super-slow, should be optimized
//
// assumes both computeFs() AND computeFdots() have been called
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_MESH::computeClampedHessian(const VOLUME::HYPERELASTIC& hyperelastic,
                                              const VOLUME::DAMPING& damping) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perElementHessians(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const MATRIX3& F       = _Fs[i];
    const MATRIX3& Fdot    = _Fdots[i];
    const MATRIX9x12& pFpx = _pFpxs[i];
    const MATRIX9 hessian  = _restTetVolumes[i] * (hyperelastic.clampedHessian(F) +
                                                   damping.clampedHessian(F, Fdot));
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}
*/

///////////////////////////////////////////////////////////////////////
// get the current displacement in vector form
///////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::getDisplacement() const
{
  VECTOR delta(_vertices.size() * 3);
  delta.setZero();

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    const VECTOR3 diff = _vertices[x] - _restVertices[x];
    const int x3 = 3 * x;
    delta[x3]     = diff[0];
    delta[x3 + 1] = diff[1];
    delta[x3 + 2] = diff[2];
  }

  return delta;
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void TET_MESH::setPositions(const VECTOR& positions)
{
  assert(positions.size() == int(_vertices.size()) * 3);

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] = positions[3 * x];
    _vertices[x][1] = positions[3 * x + 1];
    _vertices[x][2] = positions[3 * x + 2];
  }
}

///////////////////////////////////////////////////////////////////////
// set the vertex displacements to these values exactly
///////////////////////////////////////////////////////////////////////
void TET_MESH::setDisplacement(const VECTOR& delta)
{
  assert(delta.size() == int(_vertices.size()) * 3);

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] = _restVertices[x][0] + delta[3 * x];
    _vertices[x][1] = _restVertices[x][1] + delta[3 * x + 1];
    _vertices[x][2] = _restVertices[x][2] + delta[3 * x + 2];
  }
}

///////////////////////////////////////////////////////////////////////
// add the following deltas to the positions
///////////////////////////////////////////////////////////////////////
void TET_MESH::addDisplacement(const VECTOR& delta)
{
  assert(delta.size() == int(_vertices.size()) * 3);

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] += delta[3 * x];
    _vertices[x][1] += delta[3 * x + 1];
    _vertices[x][2] += delta[3 * x + 2];
  }
}

///////////////////////////////////////////////////////////////////////
// get the bounding box for the current mesh
///////////////////////////////////////////////////////////////////////
void TET_MESH::getBoundingBox(VECTOR3& mins, VECTOR3& maxs) const
{
  assert(_vertices.size() > 0);
  mins = _vertices[0];
  maxs = _vertices[0];

  for (unsigned int x = 1; x < _vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      mins[y] = (mins[y] < _vertices[x][y]) ? mins[y] : _vertices[x][y];
      maxs[y] = (maxs[y] > _vertices[x][y]) ? maxs[y] : _vertices[x][y];
    }
}

///////////////////////////////////////////////////////////////////////
// write out the surface to OBJ triangle mesh
///////////////////////////////////////////////////////////////////////
bool TET_MESH::writeSurfaceToObj(const string& filename, const TET_MESH& tetMesh)
{
  FILE* file = fopen(filename.c_str(), "w");
  if (file == NULL)
  {
    cout << " Failed to open file " << filename.c_str() << endl;
    return false;
  }
  cout << " Writing out file " << filename.c_str() << " ... " << flush;

  const vector<VECTOR3>& vertices = tetMesh.vertices();
  const vector<VECTOR3I>& surfaceTriangles = tetMesh.surfaceTriangles();

  // do the ugly thing and just write out all the vertices, even
  // the internal ones
  for (unsigned int x = 0; x < vertices.size(); x++)
    fprintf(file, "v %f %f %f\n", vertices[x][0], vertices[x][1], vertices[x][2]);

  // write out the indices for the surface triangles, but remember
  // that Objs are 1-indexed
  for (unsigned int x = 0; x < surfaceTriangles.size(); x++)
    fprintf(file, "f %i %i %i\n", surfaceTriangles[x][0] + 1, 
                                  surfaceTriangles[x][1] + 1,
                                  surfaceTriangles[x][2] + 1);

  fclose(file);
  cout << " done. " << endl;
  return true;
}

///////////////////////////////////////////////////////////////////////
// read in OBJ-style tet mesh file
///////////////////////////////////////////////////////////////////////
bool TET_MESH::readTobjFile(const string& filename, 
                            vector<VECTOR3>& vertices,
                            vector<VECTOR4I>& tets)
{
  // erase whatever was in the vectors before
  vertices.clear();
  tets.clear();

  FILE* file = fopen(filename.c_str(), "r");
  cout << " Reading in *.tobj file " << filename.c_str() << endl;

  if (file == NULL)
  {
    cout << " Failed to open file!" << endl;
    return false;
  }

  char nextChar = getc(file);

  // get the vertices
  while (nextChar == 'v' && nextChar != EOF)
  {
    ungetc(nextChar, file);

    double v[3];
    fscanf(file, "v %lf %lf %lf\n", &v[0], &v[1], &v[2]);
    vertices.push_back(VECTOR3(v[0], v[1], v[2]));

    nextChar = getc(file);
  }
  if (nextChar == EOF)
  {
    cout << " ERROR: File contains only vertices and no tets! " << endl;
    return false;
  }
  cout << " Found " << vertices.size() << " vertices " << endl;

  // get the tets
  while (nextChar == 't' && nextChar != EOF)
  {
    ungetc(nextChar, file);

    VECTOR4I tet;
    fscanf(file, "t %i %i %i %i\n", &tet[0], &tet[1], &tet[2], &tet[3]);
    tets.push_back(tet);

    nextChar = getc(file);
  }
  cout << " Found " << tets.size() << " tets " << endl;
  fclose(file);

  return true;
}

///////////////////////////////////////////////////////////////////////
// write out OBJ-style tet mesh file
///////////////////////////////////////////////////////////////////////
bool TET_MESH::writeTobjFile(const string& filename, const TET_MESH& tetMesh, 
                             const bool restVertices)
{
  FILE* file = fopen(filename.c_str(), "w");
  cout << " Writing out tet mesh file: " << filename.c_str() << endl;
 
  if (file == NULL)
  {
    cout << " Failed to open file!" << endl;
    return false;
  }

  const vector<VECTOR3>& vertices = (restVertices) ? tetMesh.restVertices() : tetMesh.vertices(); 
  const vector<VECTOR4I>& tets = tetMesh.tets(); 

  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    const VECTOR3& v = vertices[x];
    fprintf(file, "v %.17g %.17g %.17g\n", v[0], v[1], v[2]);
  }
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    const VECTOR4I& tet = tets[x];
    fprintf(file, "t %i %i %i %i\n", tet[0], tet[1], tet[2], tet[3]);
  }

  fclose(file);
  return true;
}

///////////////////////////////////////////////////////////////////////
// normalize vertices so that they're in a unit box, 
// centered at (0.5, 0.5, 0.5)
///////////////////////////////////////////////////////////////////////
vector<VECTOR3> TET_MESH::normalizeVertices(const vector<VECTOR3>& vertices)
{
  assert(vertices.size() > 0);
  VECTOR3 mins = vertices[0];
  VECTOR3 maxs = vertices[0];
  for (unsigned int x = 1; x < vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      mins[y] = (mins[y] < vertices[x][y]) ? mins[y] : vertices[x][y];
      maxs[y] = (maxs[y] > vertices[x][y]) ? maxs[y] : vertices[x][y];
    }

  const VECTOR3 lengths = maxs - mins;
  const REAL maxLengthInv = 1.0 / lengths.maxCoeff();

  vector<VECTOR3> normalized = vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    normalized[x] -= mins;
    normalized[x] *= maxLengthInv;
    
    normalized[x] += VECTOR3(0.5, 0.5, 0.5);
  }

  return normalized;
}

///////////////////////////////////////////////////////////////////////
// see if the projection of v onto the plane of v0,v1,v2 is inside 
// the triangle formed by v0,v1,v2
///////////////////////////////////////////////////////////////////////
bool TET_MESH::pointProjectsInsideTriangle(const VECTOR3& v0, const VECTOR3& v1, 
                                           const VECTOR3& v2, const VECTOR3& v)
{
  // get the barycentric coordinates
  const VECTOR3 e1 = v1 - v0;
  const VECTOR3 e2 = v2 - v0;
  const VECTOR3 n = e1.cross(e2);
  const VECTOR3 na = (v2 - v1).cross(v - v1);
  const VECTOR3 nb = (v0 - v2).cross(v - v2);
  const VECTOR3 nc = (v1 - v0).cross(v - v0);
  const VECTOR3 barycentric(n.dot(na) / n.squaredNorm(),
                            n.dot(nb) / n.squaredNorm(),
                            n.dot(nc) / n.squaredNorm());
                            
  const REAL barySum = fabs(barycentric[0]) + fabs(barycentric[1]) + fabs(barycentric[2]);

  // if the point projects to inside the triangle, it should sum to 1
  if (barySum - 1.0 < 1e-8)
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute distance between a point and triangle
///////////////////////////////////////////////////////////////////////
REAL TET_MESH::pointTriangleDistance(const VECTOR3& v0, const VECTOR3& v1, 
                                     const VECTOR3& v2, const VECTOR3& v)
{
  // get the barycentric coordinates
  const VECTOR3 e1 = v1 - v0;
  const VECTOR3 e2 = v2 - v0;
  const VECTOR3 n = e1.cross(e2);
  const VECTOR3 na = (v2 - v1).cross(v - v1);
  const VECTOR3 nb = (v0 - v2).cross(v - v2);
  const VECTOR3 nc = (v1 - v0).cross(v - v0);
  const VECTOR3 barycentric(n.dot(na) / n.squaredNorm(),
                            n.dot(nb) / n.squaredNorm(),
                            n.dot(nc) / n.squaredNorm());
                            
  const REAL barySum = fabs(barycentric[0]) + fabs(barycentric[1]) + fabs(barycentric[2]);

  // if the point projects to inside the triangle, it should sum to 1
  if (barySum - 1.0 < 1e-8)
  {
    const VECTOR3 nHat = n / n.norm();
    const REAL normalDistance = (nHat.dot(v - v0));
    return fabs(normalDistance);
  }

  // project onto each edge, find the distance to each edge
  const VECTOR3 e3 = v2 - v1;
  const VECTOR3 ev = v - v0;
  const VECTOR3 ev3 = v - v1;
  const VECTOR3 e1Hat = e1 / e1.norm();
  const VECTOR3 e2Hat = e2 / e2.norm();
  const VECTOR3 e3Hat = e3 / e3.norm();
  VECTOR3 edgeDistances(FLT_MAX, FLT_MAX, FLT_MAX);

  // see if it projects onto the interval of the edge
  // if it doesn't, then the vertex distance will be smaller,
  // so we can skip computing anything
  const REAL e1dot = e1Hat.dot(ev);
  if (e1dot > 0.0 && e1dot < e1.norm())
  {
    const VECTOR3 projected = v0 + e1Hat * e1dot;
    edgeDistances[0] = (v - projected).norm();
  }
  const REAL e2dot = e2Hat.dot(ev);
  if (e2dot > 0.0 && e2dot < e2.norm())
  {
    const VECTOR3 projected = v0 + e2Hat * e2dot;
    edgeDistances[1] = (v - projected).norm();
  }
  const REAL e3dot = e3Hat.dot(ev3);
  if (e3dot > 0.0 && e3dot < e3.norm())
  {
    const VECTOR3 projected = v1 + e3Hat * e3dot;
    edgeDistances[2] = (v - projected).norm();
  }

  // get the distance to each vertex
  const VECTOR3 vertexDistances((v - v0).norm(), 
                                (v - v1).norm(), 
                                (v - v2).norm());

  // get the smallest of both the edge and vertex distances
  const REAL vertexMin = vertexDistances.minCoeff();
  const REAL edgeMin = edgeDistances.minCoeff();

  // return the smallest of those
  return (vertexMin < edgeMin) ? vertexMin : edgeMin;
}

///////////////////////////////////////////////////////////////////////
// see if the vertex is inside the collision cell described in 
// Chapter 11: Collision Processing, [Kim and Eberle 2020]
//
// this seems to work okay, but need stronger debugging checks
///////////////////////////////////////////////////////////////////////
bool TET_MESH::insideCollisionCell(const int surfaceTriangleID, const VECTOR3& vertex)
{
  const VECTOR3I& t = _surfaceTriangles[surfaceTriangleID];
  vector<VECTOR3> v;
  v.push_back(_vertices[t[0]]);
  v.push_back(_vertices[t[1]]);
  v.push_back(_vertices[t[2]]);
  VECTOR3 n = planeNormal(v);

  // get the normals of the three adjacent faces
  vector<VECTOR3> nNeighbors;
  const VECTOR3I& neighbors = _surfaceTriangleNeighbors[surfaceTriangleID];
  for (int x = 0; x < 3; x++)
  {
    assert(neighbors[x] != -1);
    const VECTOR3I& tNeighbor = _surfaceTriangles[neighbors[x]];
    vector<VECTOR3> vNeighbor;
    vNeighbor.push_back(_vertices[tNeighbor[0]]);
    vNeighbor.push_back(_vertices[tNeighbor[1]]);
    vNeighbor.push_back(_vertices[tNeighbor[2]]);
    VECTOR3 nNeighbor = planeNormal(vNeighbor);

    nNeighbors.push_back(nNeighbor);
  }

  // do the inside check
  for (int x = 0; x < 3; x++)
  {
    const VECTOR3 ne = (nNeighbors[x] + n).normalized();
    const VECTOR3 eij = v[(x+1) % 3] - v[x];
    const VECTOR3 neb = ne.cross(eij);
    const VECTOR3 nebHat = neb.normalized();
    const REAL deplane = nebHat.dot(vertex - v[x]);

    if (deplane < 0.0)
      return false;
  }

  // TODO: some face plane compares, and returning if the vertex is above
  // or below the face
  return true;
}

///////////////////////////////////////////////////////////////////////
// compute distance to collision cell wall, where positive means inside
// and negative means outside
///////////////////////////////////////////////////////////////////////
REAL TET_MESH::distanceToCollisionCellWall(const int surfaceTriangleID, const VECTOR3& vertex)
{
  const VECTOR3I& t = _surfaceTriangles[surfaceTriangleID];
  vector<VECTOR3> v;
  v.push_back(_vertices[t[0]]);
  v.push_back(_vertices[t[1]]);
  v.push_back(_vertices[t[2]]);
  VECTOR3 n = planeNormal(v);

  // get the normals of the three adjacent faces
  vector<VECTOR3> nNeighbors;
  const VECTOR3I& neighbors = _surfaceTriangleNeighbors[surfaceTriangleID];
  for (int x = 0; x < 3; x++)
  {
    assert(neighbors[x] != -1);
    const VECTOR3I& tNeighbor = _surfaceTriangles[neighbors[x]];
    vector<VECTOR3> vNeighbor;
    vNeighbor.push_back(_vertices[tNeighbor[0]]);
    vNeighbor.push_back(_vertices[tNeighbor[1]]);
    vNeighbor.push_back(_vertices[tNeighbor[2]]);
    VECTOR3 nNeighbor = planeNormal(vNeighbor);

    nNeighbors.push_back(nNeighbor);
  }

  // do the inside check
  //int whichNeighbor = -1;
  REAL smallestDistance = FLT_MAX;
  for (int x = 0; x < 3; x++)
  {
    // averaged normal along edge
    const VECTOR3 ne = (nNeighbors[x] + n).normalized();

    // the edge itself
    const VECTOR3 eij = v[(x+1) % 3] - v[x];

    // inward-facing normal into cell
    const VECTOR3 neb = ne.cross(eij);
    const VECTOR3 nebHat = neb.normalized();

    // dot of the current vertex against the inward-facing normal
    const REAL deplane = nebHat.dot(vertex - v[x]);

    if (fabs(deplane) < smallestDistance)
    {
      //whichNeighbor = x;
      smallestDistance = fabs(deplane);
    }
  }

  return smallestDistance;
}

///////////////////////////////////////////////////////////////////////
// find all the vertex-face collision pairs, using the 
// InFaceRegion test from "Collision Processing" chapter of
// "Dynamic Deformables"
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeVertexFaceCollisions()
{
  TIMER functionTimer(__FUNCTION__);

  // if a vertex is part of an inverted tet, don't have it participate 
  // in a self-collision. That tet needs to get its house in order 
  // before it starts bossing around a surface face. Not checking for 
  // this causes faces to get horribly tangled in inverted configurations.
  computeInvertedVertices();

  _vertexFaceCollisions.clear();
  const REAL collisionEps = _collisionEps;

  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
  {
    const int currentID = _surfaceVertices[x];

    // if the vertex is involved in an inverted tet, give up
    if (_invertedVertices[currentID]) 
      continue;

    const VECTOR3& surfaceVertex = _vertices[currentID];

    // find the close triangles
    for (unsigned int y = 0; y < _surfaceTriangles.size(); y++)
    {
      // if the surface triangle is so small the normal could be degenerate, skip it
      if (surfaceTriangleIsDegenerate(y))
        continue;

      const VECTOR3I& t = _surfaceTriangles[y];

      // if it's an inverted face, move on
      if (_invertedVertices[t[0]] && _invertedVertices[t[1]] && _invertedVertices[t[2]])
        continue;

      // if this triangle is in the one-ring of the current vertex, skip it
      if (t[0] == currentID || t[1] == currentID || t[2] == currentID) continue;
      
      const REAL distance = pointTriangleDistance(_vertices[t[0]], _vertices[t[1]],
                                                  _vertices[t[2]], surfaceVertex);

      if (distance < collisionEps)
      {
        // if the point, projected onto the face's plane, is inside the face,
        // then record the collision now
        if (pointProjectsInsideTriangle(_vertices[t[0]], _vertices[t[1]], 
                                        _vertices[t[2]], surfaceVertex))
        {
          pair<int,int> collision(currentID, y);
          _vertexFaceCollisions.push_back(collision);
          continue;
        }
        if (insideCollisionCell(y, surfaceVertex))
        {
          pair<int,int> collision(currentID, y);
          _vertexFaceCollisions.push_back(collision);
        }
      }
    }
  }

#if VERY_VERBOSE
  if (_vertexFaceCollisions.size() > 0)
    cout << " Found " << _vertexFaceCollisions.size() << " vertex-face collisions " << endl;

  cout << " pairs: " << endl;
  for (unsigned int x = 0; x < _vertexFaceCollisions.size(); x++)
  {
    const pair<int,int> collision = _vertexFaceCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")" << endl;
  }
#endif
}

///////////////////////////////////////////////////////////////////////
// find all the edge-edge self collision pairs, using the 
// brute-force tests
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeEdgeEdgeCollisions()
{
  TIMER functionTimer(__FUNCTION__);
  _edgeEdgeCollisions.clear();
  _edgeEdgeIntersections.clear();
  _edgeEdgeCoordinates.clear();
  _edgeEdgeCollisionAreas.clear();

  // build a mapping from edge index pairs to _surfaceEdges
  map<pair<int, int>, int> edgeHash;
  for (unsigned int x = 0; x < _surfaceEdges.size(); x++)
  {
    pair<int,int> edge(_surfaceEdges[x][0], _surfaceEdges[x][1]);
    edgeHash[edge] = x;
  }

  // get the nearest edge to each edge, not including itself
  // and ones where it shares a vertex
  for (unsigned int x = 0; x < _surfaceEdges.size(); x++)
  {
    int closestEdge = -1;
    REAL closestDistance = FLT_MAX;
    VECTOR2 aClosest(-1,-1);
    VECTOR2 bClosest(-1,-1);
    const VECTOR2I outerEdge = _surfaceEdges[x];
    const VECTOR3& v0 = _vertices[outerEdge[0]];
    const VECTOR3& v1 = _vertices[outerEdge[1]];

    // find the closest other edge
    for (unsigned int y = x + 1; y < _surfaceEdges.size(); y++)
    {
      const VECTOR2I innerEdge = _surfaceEdges[y];
      // if they share a vertex, skip it
      if ((outerEdge[0] == innerEdge[0]) || (outerEdge[0] == innerEdge[1]) ||
          (outerEdge[1] == innerEdge[0]) || (outerEdge[1] == innerEdge[1]))
        continue;

      const VECTOR3& v2 = _vertices[innerEdge[0]];
      const VECTOR3& v3 = _vertices[innerEdge[1]];

      VECTOR3 innerPoint, outerPoint;
      IntersectLineSegments(v0, v1, v2, v3,
                            outerPoint, innerPoint);  

      const REAL distance = (innerPoint - outerPoint).norm();

      //if (distance < separationDistance)
      //  separationDistance = distance;

      if (distance > closestDistance) continue;

      // get the line interpolation coordinates
      VECTOR2 a,b;
      const VECTOR3 e0 = v1 - v0;
      const VECTOR3 e1 = v3 - v2;

      // this is a little dicey in general, but if the intersection test isn't
      // total garbage, it should still be robust
      a[1] = (outerPoint - v0).norm() / e0.norm();
      a[0] = 1.0 - a[1];
      b[1] = (innerPoint - v2).norm() / e1.norm();
      b[0] = 1.0 - b[1];

      // if it's really close to an end vertex, skip it
      const REAL skipEps = 1e-4;
      if ((a[0] < skipEps) || (a[0] > 1.0 - skipEps)) continue;
      if ((a[1] < skipEps) || (a[1] > 1.0 - skipEps)) continue;
      if ((b[0] < skipEps) || (b[0] > 1.0 - skipEps)) continue;
      if ((b[1] < skipEps) || (b[1] > 1.0 - skipEps)) continue;

      // it's mid-segment, and closest, so remember it
      closestDistance = distance;
      closestEdge = y;

      aClosest = a;
      bClosest = b;
    }

    // if nothing was close, move on
    if (closestEdge == -1) continue;

    // are they within each other's one rings?
    const VECTOR2I innerEdge = _surfaceEdges[closestEdge];
    bool insideOneRing = false;

    for (int j = 0; j < 2; j++)
    {
      pair<int, int> lookup;
      lookup.first = outerEdge[j];
      for (int i = 0; i < 2; i++)
      {
        lookup.second = innerEdge[i];
        if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end())
          insideOneRing = true;
      }
    }
    if (insideOneRing) continue;

    // if it's within the positive threshold, it's in collision
    if (closestDistance < _collisionEps)
    {
      pair<int,int> collision(x, closestEdge);
      _edgeEdgeCollisions.push_back(collision);
      
      // this was actually set, right?
      assert(aClosest[0] > 0.0 && aClosest[1] > 0.0);
      assert(bClosest[0] > 0.0 && bClosest[1] > 0.0);

      pair<VECTOR2,VECTOR2> coordinate(aClosest, bClosest);
      _edgeEdgeCoordinates.push_back(coordinate);

      // get the areas too
      const VECTOR2I innerEdge = _surfaceEdges[closestEdge];
      const pair<int,int> outerPair(outerEdge[0], outerEdge[1]);
      const pair<int,int> innerPair(innerEdge[0], innerEdge[1]);
      const REAL xArea = _restEdgeAreas[edgeHash[outerPair]];
      const REAL closestArea = _restEdgeAreas[edgeHash[innerPair]];
      _edgeEdgeCollisionAreas.push_back(xArea + closestArea);

      // find out if they are penetrating
      vector<VECTOR3> edge(2);
      edge[0] = v0;
      edge[1] = v1;

      // get the adjacent triangles of the *other* edge
      VECTOR2I adjacentTriangles = _surfaceEdgeTriangleNeighbors[edgeHash[innerPair]];

      // build triangle 0
      const VECTOR3I surfaceTriangle0 = _surfaceTriangles[adjacentTriangles[0]];
      vector<VECTOR3> triangle0;
      triangle0.push_back(_vertices[surfaceTriangle0[0]]);
      triangle0.push_back(_vertices[surfaceTriangle0[1]]);
      triangle0.push_back(_vertices[surfaceTriangle0[2]]);

      // build triangle 1
      vector<VECTOR3> triangle1;
      // if there's another triangle on the other side (this is in case we're looking at cloth)
      // then store that one too
      if (adjacentTriangles[1] != -1)
      {
        const VECTOR3I surfaceTriangle1 = _surfaceTriangles[adjacentTriangles[1]];
        triangle1.push_back(_vertices[surfaceTriangle1[0]]);
        triangle1.push_back(_vertices[surfaceTriangle1[1]]);
        triangle1.push_back(_vertices[surfaceTriangle1[2]]);
      }

      // see if the edges are already penetrating the opposing faces
      bool penetrating = false;
      if (triangle0.size() > 0) penetrating = faceEdgeIntersection(triangle0, edge);
      if (triangle1.size() > 0) penetrating = penetrating || faceEdgeIntersection(triangle1, edge);

      _edgeEdgeIntersections.push_back(penetrating);
      //_edgeEdgeCollisionEps.push_back(_collisionEps);

      // TODO: for completeness, should probably test the other edges against the other
      // pair, just in case we're looking at a degenerate case. In general, seems redundant.
    }
  }
  assert(_edgeEdgeCollisions.size() == _edgeEdgeCoordinates.size());

#if VERY_VERBOSE
  if (_edgeEdgeCollisions.size() > 0)
    cout << " Found " << _edgeEdgeCollisions.size() << " edge-edge collisions " << endl;

  cout << " pairs: " << endl;
  for (unsigned int x = 0; x < _edgeEdgeCollisions.size(); x++)
  {
    const pair<int,int> collision = _edgeEdgeCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")" << endl;
  }
#endif
}

///////////////////////////////////////////////////////////////////////
// get the triangle area
///////////////////////////////////////////////////////////////////////
REAL TET_MESH::triangleArea(const vector<VECTOR3>& triangle)
{
  const VECTOR3 edge1 = triangle[1] - triangle[0];
  const VECTOR3 edge2 = triangle[2] - triangle[0];
  return 0.5 * edge1.cross(edge2).norm();
}

///////////////////////////////////////////////////////////////////////
// get the normal to a plane, specified by three points
///////////////////////////////////////////////////////////////////////
VECTOR3 TET_MESH::planeNormal(const vector<VECTOR3>& plane)
{
  const VECTOR3 edge1 = plane[1] - plane[0];
  const VECTOR3 edge2 = plane[2] - plane[0];
  return edge1.cross(edge2).normalized();
}

///////////////////////////////////////////////////////////////////////
// project point onto plane, specific by three points
///////////////////////////////////////////////////////////////////////
VECTOR3 TET_MESH::pointPlaneProjection(const vector<VECTOR3>& plane, const VECTOR3& point)
{
  const VECTOR3 normal = planeNormal(plane);
  return point - (normal.dot(point - plane[0])) * normal;
}

///////////////////////////////////////////////////////////////////////
// based on vertex-face collision pairs, build "collision tets"
///////////////////////////////////////////////////////////////////////
void TET_MESH::buildVertexFaceCollisionTets(const VECTOR& velocity)
{
  // clear out old ones (not clear if this is smart or dumb)
  _vertexFaceCollisionTets.clear();
  _vertexFaceCollisionAreas.clear();

  // make a tet for each vertex-face pair
  for (unsigned int x = 0; x < _vertexFaceCollisions.size(); x++)
  {
    const int vertexID = _vertexFaceCollisions[x].first;
    const int faceID = _vertexFaceCollisions[x].second;
    const VECTOR3I& face = _surfaceTriangles[faceID];

    // build a tet with the correct vertex ordering
    VECTOR4I tet;
    tet[0] = vertexID;

    // reverse the ordering here because we want the normal to face
    // the opposite direction compared to a surface triangle
    tet[1] = face[2];
    tet[2] = face[1];
    tet[3] = face[0];

    // get the rest area of the triangle
    vector<VECTOR3> restFace(3);
    restFace[0] = _restVertices[face[0]]; 
    restFace[1] = _restVertices[face[1]];
    restFace[2] = _restVertices[face[2]];
    const REAL restFaceArea = triangleArea(restFace);

    // reverse lookup exists, right?
    assert(_volumeToSurfaceID.find(vertexID) != _volumeToSurfaceID.end());
    const REAL surfaceID = _volumeToSurfaceID[vertexID];
    const REAL restVertexArea = _restOneRingAreas[surfaceID];

    // store everything
    _vertexFaceCollisionTets.push_back(tet);
    assert(restFaceArea >= 0.0);
    assert(restVertexArea >= 0.0);
    _vertexFaceCollisionAreas.push_back(restFaceArea + restVertexArea);
  }
}

///////////////////////////////////////////////////////////////////////
// compute collision forces using collision tets
///////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::computeVertexFaceCollisionForces() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perElementForces(_vertexFaceCollisionTets.size());
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    vector<VECTOR3> vs(4);
    for (unsigned int j = 0; j < 4; j++)
      vs[j] = _vertices[_vertexFaceCollisionTets[i][j]];
    const VECTOR12 force = -_vertexFaceCollisionAreas[i] * _vertexFaceEnergy->gradient(vs);
    perElementForces[i] = force;

#if ENABLE_DEBUG_TRAPS
    if (force.hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " NaN in collision tet: " << i << endl;
      for (int j = 0; j < 4; j++)
        cout << " v" << j << ": " << vs[j].transpose() << endl;
      cout << " gradient: " << endl << _vertexFaceEnergy->gradient(vs) << endl;
    }
#endif
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    const VECTOR4I& tet = _vertexFaceCollisionTets[i];
    const VECTOR12& tetForce = perElementForces[i];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision energy using x-based formulation
///////////////////////////////////////////////////////////////////////
REAL TET_MESH::computeEdgeEdgeCollisionEnergy() const
{
  TIMER functionTimer(__FUNCTION__);

  REAL finalEnergy = 0.0;
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _surfaceEdges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _surfaceEdges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];

    const VECTOR2& a = _edgeEdgeCoordinates[i].first;
    const VECTOR2& b = _edgeEdgeCoordinates[i].second;

    const REAL psi = _edgeEdgeEnergy->psi(vs,a,b);
    finalEnergy += _edgeEdgeCollisionAreas[i] * psi;
  }

  return finalEnergy;
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision forces using x-based formulation
///////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::computeEdgeEdgeCollisionForces() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perElementForces(_edgeEdgeCollisions.size());
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _surfaceEdges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _surfaceEdges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];

    const VECTOR2& a = _edgeEdgeCoordinates[i].first;
    const VECTOR2& b = _edgeEdgeCoordinates[i].second;

#if ADD_EDGE_EDGE_PENETRATION_BUG
    const VECTOR12 force = -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->gradient(vs,a,b);
#else
    const VECTOR12 force = (!_edgeEdgeIntersections[i]) ? -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->gradient(vs,a,b)
                                                        : -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->gradientNegated(vs,a,b);
#endif

    perElementForces[i] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _surfaceEdges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _surfaceEdges[_edgeEdgeCollisions[i].second];
    const VECTOR12& edgeForce = perElementForces[i];

    vector<int> vertexIndices(4);
    vertexIndices[0] = edge0[0];
    vertexIndices[1] = edge0[1];
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];

    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * vertexIndices[x];
      assert((int)index < DOFs);
      forces[index]     += edgeForce[3 * x];
      forces[index + 1] += edgeForce[3 * x + 1];
      forces[index + 2] += edgeForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision Hessians using x-based formulation
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_MESH::computeEdgeEdgeCollisionClampedHessian() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<MATRIX12> perElementHessians(_edgeEdgeCollisions.size());
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _surfaceEdges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _surfaceEdges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];

    const VECTOR2& a = _edgeEdgeCoordinates[i].first;
    const VECTOR2& b = _edgeEdgeCoordinates[i].second;
#if ADD_EDGE_EDGE_PENETRATION_BUG
    const MATRIX12 H = -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->clampedHessian(vs,a,b);
#else
    const MATRIX12 H = (!_edgeEdgeIntersections[i]) ? -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->clampedHessian(vs,a,b)
                                                    : -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->clampedHessianNegated(vs,a,b);
#endif
    perElementHessians[i] = H;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const MATRIX12& H = perElementHessians[i];
    const VECTOR2I& edge0 = _surfaceEdges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _surfaceEdges[_edgeEdgeCollisions[i].second];

    vector<int> vertexIndices(4);
    vertexIndices[0] = edge0[0];
    vertexIndices[1] = edge0[1];
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];

    for (int y = 0; y < 4; y++)
    {
      int yVertex = vertexIndices[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = vertexIndices[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }
 
  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the collision force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_MESH::computeVertexFaceCollisionClampedHessian() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<MATRIX12> perElementHessians(_vertexFaceCollisionTets.size());
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    vector<VECTOR3> vs(4);
    for (unsigned int j = 0; j < 4; j++)
      vs[j] = _vertices[_vertexFaceCollisionTets[i][j]];
    const MATRIX12 H = -_vertexFaceCollisionAreas[i] * _vertexFaceEnergy->clampedHessian(vs);
    perElementHessians[i] = H;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    const VECTOR4I& tet = _vertexFaceCollisionTets[i];
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// find out how close all the edges are initially
//
// for each pair of _surfaceEdges, what _collisionEps should we use? 
// If they started out closer than _collisionEps, then we need to set 
// a smaller tolerance.
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeEdgeEdgeRestDistance()
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Computing edge-edge rest distances ... " << flush;
  _edgeEdgeRestDistance.clear();

  // hit all the pairs
  for (unsigned int y = 0; y < _surfaceEdges.size(); y++)
  {
    // cache the vertex positions
    const VECTOR2I yEdge = _surfaceEdges[y];
    const VECTOR3& v0 = _vertices[yEdge[0]];
    const VECTOR3& v1 = _vertices[yEdge[1]];

    // Wonder if meshes will get big enough that this needs to
    // be upgraded to 64-bit? Would be a good problem to have ...
    unsigned int yFlat = yEdge[0] + yEdge[1] * _surfaceEdges.size();

    for (unsigned int x = 0; x < _surfaceEdges.size(); x++)
    {
      // don't compare against yourself
      if (x == y) continue;

      const VECTOR2I xEdge = _surfaceEdges[x];
      unsigned int xFlat = xEdge[0] + xEdge[1] * _surfaceEdges.size();

      const VECTOR3& v2 = _vertices[xEdge[0]];
      const VECTOR3& v3 = _vertices[xEdge[1]];
    
      // do the actual distance test
      VECTOR3 innerPoint, outerPoint;
      IntersectLineSegments(v0, v1, v2, v3,
                            outerPoint, innerPoint);  
      const REAL distance = (innerPoint - outerPoint).norm();

      // store the pair
      pair<unsigned int,unsigned int> edgeEdge(xFlat, yFlat);

      // hash the distance
      _edgeEdgeRestDistance[edgeEdge] = distance;
    }
  }
  cout << "done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// compute whether one vertex is inside the vertex one right of another
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceVertexOneRings()
{
  _insideSurfaceVertexOneRing.clear();
  for (unsigned int x = 0; x < _surfaceEdges.size(); x++)
  {
    // gonna be lazy for a moment here
    const VECTOR2I edge = _surfaceEdges[x];
    _insideSurfaceVertexOneRing[pair<int, int>(edge[0], edge[1])] = true;
    _insideSurfaceVertexOneRing[pair<int, int>(edge[1], edge[0])] = true;
  }
}

///////////////////////////////////////////////////////////////////////
// set collision eps to something new
///////////////////////////////////////////////////////////////////////
void TET_MESH::setCollisionEps(const REAL& eps)
{
  _collisionEps = eps;
  _vertexFaceEnergy->eps() = eps;
  _edgeEdgeEnergy->setEps(eps);
}

///////////////////////////////////////////////////////////////////////
// set collision stiffness to something new
///////////////////////////////////////////////////////////////////////
void TET_MESH::setCollisionStiffness(const REAL& mu)
{
  _vertexFaceEnergy->mu() = mu;
  _edgeEdgeEnergy->mu() = mu;
}

///////////////////////////////////////////////////////////////////////
// are these two surface triangles neighbors?
///////////////////////////////////////////////////////////////////////
bool TET_MESH::areSurfaceTriangleNeighbors(const int id0, const int id1) const
{
  assert(_surfaceTriangleNeighbors.size() > 0);
  assert(id0 < (int)_surfaceTriangleNeighbors.size());
  assert(id1 < (int)_surfaceTriangleNeighbors.size());

  const VECTOR3I neighbors0 = _surfaceTriangleNeighbors[id0];

  for (int x = 0; x < 3; x++)
    if (neighbors0[x] == id1)
      return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute the normal of the surface triangle at 
// _surfaceTriangles[triangleID];
///////////////////////////////////////////////////////////////////////
VECTOR3 TET_MESH::surfaceTriangleNormal(const int triangleID) const
{
  assert(triangleID < (int)_surfaceTriangles.size());

  const VECTOR3I& vertexIDs = _surfaceTriangles[triangleID];
  const VECTOR3& v0 = _vertices[vertexIDs[0]];
  const VECTOR3& v1 = _vertices[vertexIDs[1]];
  const VECTOR3& v2 = _vertices[vertexIDs[2]];

  const VECTOR3& e0 = v1 - v0;
  const VECTOR3& e1 = v2 - v0;

  return e0.cross(e1).normalized();
}

///////////////////////////////////////////////////////////////////////
// set collision pairs (for replays)
///////////////////////////////////////////////////////////////////////
void TET_MESH::setCollisionPairs(const vector<pair<int,int> >& vertexFace, 
                                 const vector<pair<int,int> >& edgeEdge)
{
  _vertexFaceCollisions = vertexFace;
  _edgeEdgeCollisions = edgeEdge;
}

///////////////////////////////////////////////////////////////////////
// compute the dihedral angle between two surface faces
///////////////////////////////////////////////////////////////////////
REAL TET_MESH::surfaceFaceDihedralAngle(const int surfaceID0, const int surfaceID1) const
{
  const VECTOR4I tet = buildSurfaceFlap(surfaceID0, surfaceID1);

  // let's do some cross products ...
  //
  //         1
  //
  //         o
  //        /|\
  //       / | \
  //      /  |  \
  //  0  o   |   o  3
  //      \  |  /
  //       \ | /
  //        \|/
  //         o
  //
  //         2
  //

  const VECTOR3& v0 = _vertices[tet[0]];
  const VECTOR3& v1 = _vertices[tet[1]];
  const VECTOR3& v2 = _vertices[tet[2]];
  const VECTOR3& v3 = _vertices[tet[3]];

  const VECTOR3 e20 = v2 - v0;
  const VECTOR3 e10 = v1 - v0;
  const VECTOR3 n0 = e20.cross(e10) / (e20 - e10).norm();

  const VECTOR3 e13 = v1 - v3;
  const VECTOR3 e23 = v2 - v3;
  const VECTOR3 n1 = e13.cross(e23) / (e13 - e23).norm();

  const VECTOR3 e12 = (v1 - v2) / (v1 - v2).norm();

  const REAL sinTheta = (n0.cross(n1)).dot(e12);
  const REAL cosTheta = n0.dot(n1);

  return atan2(sinTheta, cosTheta);
}

///////////////////////////////////////////////////////////////////////
// build a consistent tet/flap ordering from two surface triangles
///////////////////////////////////////////////////////////////////////
VECTOR4I TET_MESH::buildSurfaceFlap(const int surfaceID0, const int surfaceID1) const
{
  // they're legal indices, right?
  assert(surfaceID0 >= 0);
  assert(surfaceID1 >= 0);
  assert(surfaceID0 < (int)_surfaceTriangles.size());
  assert(surfaceID1 < (int)_surfaceTriangles.size());

  // they are neighbors, right?
  assert(areSurfaceTriangleNeighbors(surfaceID0, surfaceID1));

  const VECTOR3I f0 = _surfaceTriangles[surfaceID0];
  const VECTOR3I f1 = _surfaceTriangles[surfaceID1];

  int firstMatch = -1;
  int secondMatch = -1;
  int unmatched0 = -1;
  int unmatched1 = -1;

  // find the tet indices for the first face
  for (int x = 0; x < 3; x++)
  {
    // let's search for this index
    int i0 = f0[x];

    bool matchFound = false;
    for (int y = 0; y < 3; y++)
    {
      // see if it matches
      if (i0 == f1[y])
      {
        // which match is it?
        if (firstMatch == -1)
          firstMatch = i0;
        else
          secondMatch = i0;
        
        matchFound = true;
      }
    }

    if (!matchFound)
      unmatched0 = i0;
  }

  // find the unmatched vertex from the second face
  for (int x = 0; x < 3; x++)
  {
    // let's search for this index
    int i1 = f1[x];
    if (i1 != firstMatch && i1 != secondMatch)
      unmatched1 = i1;
  }

  // we did find one, right?
  assert(unmatched1 != -1);

  // build a tet/flap
  //
  //         1
  //
  //         o
  //        /|\
  //       / | \
  //      /  |  \
  //  0  o   |   o  3
  //      \  |  /
  //       \ | /
  //        \|/
  //         o
  //
  //         2
  //
  VECTOR4I tet;
  tet[0] = unmatched0;
  tet[1] = secondMatch;
  tet[2] = firstMatch;
  tet[3] = unmatched1;

  return tet;
}

///////////////////////////////////////////////////////////////////////
// get mass-weighted global translation
///////////////////////////////////////////////////////////////////////
VECTOR3 TET_MESH::getTranslation() const
{
  VECTOR3 vertexSum;
  vertexSum.setZero();

  REAL volumeSum = 0;
  assert(_vertices.size() == _restOneRingVolumes.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    volumeSum += _restOneRingVolumes[x];
    vertexSum += _vertices[x] * _restOneRingVolumes[x];
  }

  return vertexSum * (1.0 / volumeSum);
}

///////////////////////////////////////////////////////////////////////
// get volume-weighted global translation, for the rest state
///////////////////////////////////////////////////////////////////////
VECTOR3 TET_MESH::getRestTranslation() const
{
  VECTOR3 vertexSum;
  vertexSum.setZero();

  REAL volumeSum = 0;
  assert(_vertices.size() == _restOneRingVolumes.size());
  for (unsigned int x = 0; x < _restVertices.size(); x++)
  {
    volumeSum += _restOneRingVolumes[x];
    vertexSum += _restVertices[x] * _restOneRingVolumes[x];
  }

  return vertexSum * (1.0 / volumeSum);
}

///////////////////////////////////////////////////////////////////////
// get Procrustes-style global rotation using Eqn. 7 and surrounding
// text from Muller et al's "Meshless Deformations Based on Shape 
// Matching" from SIGGRAPH 2005
///////////////////////////////////////////////////////////////////////
MATRIX3 TET_MESH::getRotation() const
{
  // trying to follow Muller's notation here
  const VECTOR3 x_cm0 = getRestTranslation();
  const VECTOR3 x_cm  = getTranslation();

  // left matrix in Eqn. 7 of Muller's paper
  MATRIX3 Apq;
  Apq.setZero();
  for (unsigned int x = 0; x < _restVertices.size(); x++)
  {
    const VECTOR3 p = _vertices[x] - x_cm;
    const VECTOR3 q = _restVertices[x] - x_cm0;
    Apq += _restOneRingVolumes[x] * (p * q.transpose());
  }

  // get the rotation
  MATRIX3 R,S;
  polarDecomposition(Apq, R, S);
  return R;
}

///////////////////////////////////////////////////////////////////////
// see if a current surface triangle has been crushed to degeneracy
///////////////////////////////////////////////////////////////////////
bool TET_MESH::surfaceTriangleIsDegenerate(const int surfaceTriangleID)
{
  assert(surfaceTriangleID >= 0);
  assert(surfaceTriangleID < (int)_surfaceTriangles.size());

  // get the rest area
  vector<VECTOR3> vertices(3);
  vertices[0] = _restVertices[_surfaceTriangles[surfaceTriangleID][0]];
  vertices[1] = _restVertices[_surfaceTriangles[surfaceTriangleID][1]];
  vertices[2] = _restVertices[_surfaceTriangles[surfaceTriangleID][2]];
  const REAL restArea = triangleArea(vertices);

  // get the deformed area
  vertices[0] = _vertices[_surfaceTriangles[surfaceTriangleID][0]];
  vertices[1] = _vertices[_surfaceTriangles[surfaceTriangleID][1]];
  vertices[2] = _vertices[_surfaceTriangles[surfaceTriangleID][2]];
  const REAL deformedArea = triangleArea(vertices);

  const REAL relativeArea = deformedArea / restArea;

  const REAL degeneracyEps = 1e-4;
  if (relativeArea < degeneracyEps) return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute which vertices are attached to inverted tets
///////////////////////////////////////////////////////////////////////
void TET_MESH::computeInvertedVertices()
{
  TIMER functionTimer(__FUNCTION__);
  // first set them all to false
  _invertedVertices.resize(_vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _invertedVertices[x] = false;

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    // if the tet is not inverted, move on
    if (_Fs[x].determinant() > 0.0)
      continue;

    // if tet is inverted, tags all its vertices
    for (int y = 0; y < 4; y++)
      _invertedVertices[_tets[x][y]] = true;
  }

  int totalInverted = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    if (_invertedVertices[x])
      totalInverted++;

  //cout << " Total inverted vertices: " << totalInverted << endl;
}

}
