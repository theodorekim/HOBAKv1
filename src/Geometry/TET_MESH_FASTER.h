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
#ifndef TET_MESH_FASTER_H
#define TET_MESH_FASTER_H

#include "TET_MESH.h"
#include "AABB_TREE.h"

namespace HOBAK {

class TET_MESH_FASTER : public TET_MESH
{
public:
  TET_MESH_FASTER(const vector<VECTOR3>& restVertices, 
                  const vector<VECTOR4I>& tets);

  // do something so that this doesn't run so slow
  virtual SPARSE_MATRIX computeHyperelasticClampedHessian(const VOLUME::HYPERELASTIC& hyperelastic) const override;
  virtual SPARSE_MATRIX computeDampingHessian(const VOLUME::DAMPING& damping) const override;

  // find all the vertex-face collision pairs, using the InFaceRegion test
  virtual void computeVertexFaceCollisions() override;

  // find all the edge-edge collision pairs
  virtual void computeEdgeEdgeCollisions() override;

  const AABB_TREE& aabbTreeTriangles() const { return _aabbTreeTriangles; };
  void refitAABB() { _aabbTreeTriangles.refit(); _aabbTreeEdges.refit(); };

private:
  // find the compressed index mapping
  void computeCompressedIndices();

  mutable bool _sparsityCached;
  mutable SPARSE_MATRIX _sparseA;

  // for sparse matrix entry (x,y), find the compressed index
  map<pair<int,int>, int> _compressedIndex;

  // cache the hessian for each tet
  mutable vector<MATRIX12> _perElementHessians;

  // mapping from edge index pairs to _surfaceEdges
  map<pair<int, int>, int> _edgeHash;

  // for each entry in the global stiffness matrix, the
  // tet indices to gather entries from
  vector<vector<VECTOR3I> > _hessianGathers;

  // collision detection acceleration structure for triangles
  AABB_TREE _aabbTreeTriangles;
  
  // collision detection acceleration structure for edges
  AABB_TREE _aabbTreeEdges;
};

}

#endif
