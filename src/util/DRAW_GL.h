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
#ifndef DRAW_GL_H
#define DRAW_GL_H

#include "Geometry/CAPSULE.h"
#include "Geometry/CYLINDER.h"
#include "Geometry/CUBE.h"
#include "Geometry/SPHERE.h"
#include "Geometry/AABB_TREE.h"
#include "Timestepper/TIMESTEPPER.h"
#include "util/COLLISION_UTIL.h"
#include <glvu.h>

namespace HOBAK {

// Print a string to the GL window
void printGlString(string output);

void drawPlaneConstraints(const TET_MESH* tetMesh, const TIMESTEPPER::TIMESTEPPER* stepper);
void drawTet(const TET_MESH& mesh, const int tetID);
VECTOR3 planeNormal(const vector<VECTOR3>& plane);

// draw the collision tets that have been built for this mesh
void drawCollisionTets(const TET_MESH& mesh);

// just draw a single vertex of a tet mesh
void drawVertex(const TET_MESH& mesh, const int index);

// just draw the vertices of a tet mesh
void drawVertices(const TET_MESH& mesh, const vector<int>& toDraw);

// draw the collision cell around a specific surface triangle
void drawSurfaceFaceCollisionCell(const TET_MESH& mesh, const int triangleIndex);

// draw a specific surface triangles of a tet mesh
void drawSurfaceFace(const TET_MESH& mesh, const int index);

// draw only the surface triangles of a tet mesh
void drawSurfaceTriangles(const TET_MESH& mesh, bool drawOutlines);
void drawSurfaceTriangles(const TET_MESH& mesh, bool drawOutlines, const VECTOR3& tetColor, const VECTOR3& outlineColor);

// draw lines between differences between the two meshes
void drawDiff(const TET_MESH& mesh0, const TET_MESH& mesh1);

// draw kinematic objects
void drawKinematicShape(const KINEMATIC_SHAPE& shape);
void drawCapsule(const CAPSULE& capsule);
void drawCylinder(const CYLINDER& cylinder);
void drawSphere(const SPHERE& sphere);
void drawCube(const CUBE& cube);

// draw an AABB
void drawAABB(const VECTOR3& minCorner, const VECTOR3& maxCorner);
void drawAABB(const AABB_NODE& node);

// draw a AABB tree at a specific depth
void drawAABBTree(const AABB_NODE* node, const int drawDepth, const int currentDepth);
void drawAABBTree(const AABB_TREE& tree, const int drawDepth);

// draw coordinate axes, xyz = rgb
void drawAxes();

// draw collision information
void drawVertexFacePairs(const TET_MESH& tetMesh, const int highlighted = -1);
void drawEdgeEdgePairs(const TET_MESH& tetMesh, const int pairID = -1);

} // HOBAK

#endif
