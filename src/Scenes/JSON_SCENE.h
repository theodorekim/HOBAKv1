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
#ifndef JSON_SCENE_H
#define JSON_SCENE_H

#include "SIMULATION_SCENE.h"

namespace HOBAK {

// Replaying simulation output from JSON file.
class JSON_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Replaying JSON output from simulateScene" << endl;
    cout << "=====================================================================" << endl;
  }

  // don't do anything for building the scene - FILE_IO should be doing
  // the heavy lifting here
  virtual bool buildScene() override
  {
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    if (_normalizedVertices)
      vertices = TET_MESH::normalizeVertices(vertices);

    // apply any initial transforms the scene requested
    for (unsigned int x = 0; x < vertices.size(); x++)
      vertices[x] = _initialA * vertices[x] + _initialTranslation;

    _tetMesh = new TET_MESH(vertices, tets);

    _drawFrame = 0;
    _drawFeature = true;
    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) {
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    
    _tetMesh->setDisplacement(_positions[_drawFrame]);
    _tetMesh->setCollisionPairs(_vertexFaceCollisions[_drawFrame],
                                _edgeEdgeCollisions[_drawFrame]);

    if (verbose)
      cout << " Setting to frame " << _drawFrame << endl;
    _drawFrame++;
  };

  // jump to a specific frame
  void jumpToFrame(const int frame) {
    _drawFrame = frame;
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    
    _tetMesh->setDisplacement(_positions[_drawFrame]);
    _tetMesh->setCollisionPairs(_vertexFaceCollisions[_drawFrame],
                                _edgeEdgeCollisions[_drawFrame]);

    cout << " Jumping to frame " << _drawFrame << endl;
  };

  // in case _drawFrame changed, you can update the positions
  void updatePositions()
  {
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    _tetMesh->setDisplacement(_positions[_drawFrame]);
    _tetMesh->setCollisionPairs(_vertexFaceCollisions[_drawFrame],
                                _edgeEdgeCollisions[_drawFrame]);
  }

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene() override
  {
    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
   
    glEnable(GL_DEPTH_TEST);
    if (_drawFeature)
      drawVertexFacePairs(*_tetMesh, _arrowCounter);
  };

  // draw only the tet mesh
  virtual void drawTetMesh()
  {
    const VECTOR3 tetColor(0, 0, 0.5);
    const VECTOR3 outlineColor(1, 1, 1);
    drawSurfaceTriangles(*_tetMesh, true, tetColor, outlineColor);
  };

  // subtract off the translation and rotation from the main body
  // and draw that
  void drawBodyCenteredScene(const VECTOR3& tetColor = VECTOR3(0.5, 0.5, 0.5), 
                             const VECTOR3& outlineColor = VECTOR3(0,0,0))
  {
    glEnable(GL_DEPTH_TEST);

    // let's orient ourselves
    drawAxes();

    glPushMatrix();

    // undo the global rotation
    const MATRIX3 RT = _tetMesh->getRotation().transpose();
    const Eigen::AngleAxis<GLfloat> rotation{ RT.cast<GLfloat>() };
    glRotatef((180.0 / M_PI) * rotation.angle(), rotation.axis().x(), rotation.axis().y(), rotation.axis().z());
    
    // undo the global translation
    const VECTOR3 translation = _tetMesh->getTranslation();
    glTranslatef(-translation[0], -translation[1], -translation[2]);

    // now draw
    drawSurfaceTriangles(*_tetMesh, true, tetColor, outlineColor);
    if (_drawFeature)
      drawVertexFacePairs(*_tetMesh, _arrowCounter);

    glPopMatrix();
  }
#endif

  const vector<VECTOR>& positions() const  { return _positions; };
  const vector<VECTOR>& velocities() const { return _velocities; };
  const vector<vector<pair<int,int> > >& vertexFaceCollisions() const { return _vertexFaceCollisions; };
  const vector<vector<pair<int,int> > >& edgeEdgeCollisions() const   { return _edgeEdgeCollisions; };
  vector<VECTOR>& velocities()             { return _velocities; };
  vector<VECTOR>& positions()              { return _positions; };
  vector<vector<pair<int,int> > >& vertexFaceCollisions() { return _vertexFaceCollisions; };
  vector<vector<pair<int,int> > >& edgeEdgeCollisions()   { return _edgeEdgeCollisions; };
  int& drawFrame()                         { return _drawFrame; };
  const int totalFrames() const            { return _positions.size(); };

  // have a setter here for data hiding
  void setInitialA(const MATRIX3& A)           { _initialA = A; };
  void setInitialTranslation(const VECTOR3& t) { _initialTranslation = t; };

private:
  vector<VECTOR> _positions;
  vector<VECTOR> _velocities;
  vector<vector<pair<int,int> > > _vertexFaceCollisions;
  vector<vector<pair<int,int> > > _edgeEdgeCollisions;
  int _drawFrame;
};

}

#endif
