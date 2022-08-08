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
#ifndef DROPPED_E_H
#define DROPPED_E_H

#include "SIMULATION_SCENE.h"

namespace HOBAK {

class DROPPED_E : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " An E-shaped mesh is dropped on the ground, resulting in a variety of " << endl;
    cout << " self-collisions. Both VF and EE collisions are enabled." << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "dropped_e";

    // read in the tet mesh file
    //string testFile("../data/e/e_25.tobj");
    //string testFile("../data/e/e_30.tobj");
    //string testFile("../data/e/e_40.tobj");
    //string testFile("../data/e/e_50.tobj");
    _tetMeshFilename = string("../data/e/e_30.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);
    vertices = TET_MESH::normalizeVertices(vertices);
    _normalizedVertices = true;

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    _gravity = VECTOR3(0, -1.0, 0);

    REAL E = 10.0;
    REAL nu = 0.45; // lambda \approx 10

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::ARAP(mu, 1.0);

    const vector<REAL>& areas = _tetMesh->surfaceTriangleAreas();
    REAL smallest = areas[0];
    REAL largest = areas[0];
    for (unsigned int x = 1; x < areas.size(); x++)
    {
      if (areas[x] > largest) largest  = areas[x];
      if (areas[x] < largest) smallest = areas[x];
    }
    cout << "Largest triangle area:  "  << largest << endl;
    cout << "Smallest triangle area: "  << smallest << endl;

    // build the time integrator
    _solver = new TIMESTEPPER::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic);

    VECTOR3 center0(0.0,-4.77, 0.0);
    VECTOR3 center1(0.0,6.9, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 10.0));
    _kinematicShapes.push_back(new CUBE(center1, 10.0));
    _solver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _solver->addKinematicCollisionObject(_kinematicShapes[0]);

    // collision constants
    const REAL collisionMu = 1000.0;  // default
    _solver->collisionStiffness() = collisionMu;
    _solver->collisionDampingBeta() = 0.001;

    _solver->vertexFaceSelfCollisionsOn() = true;
    _solver->edgeEdgeSelfCollisionsOn() = true;

    _solver->setDt(1.0 / 60.0);  // default

    _tetMesh->setCollisionEps(0.02);

    // for the crumpling E
    _eye    = VECTOR3(0.929159, 0.720205, 3.59228);
    _lookAt = VECTOR3(0.870294, 0.670356, 2.59526);
    _up     = VECTOR3(0.00648923, 0.99871, -0.0503162);

    // stay within github large file limits
    _pauseFrame = 200;

    return true;
  }

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    setToPreviousTimestep();

    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
    
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
   
    glEnable(GL_DEPTH_TEST);
    if (_drawFeature)
      drawVertexFacePairs(*_tetMesh, _arrowCounter);

    restoreToCurrentTimestep();
  };
#endif

};

}

#endif
