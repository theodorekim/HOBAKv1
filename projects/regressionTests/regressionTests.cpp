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
#include <cmath>
#include <iostream>
#include <float.h>

#define CATCH_CONFIG_MAIN

// deactivate all the GL stuff for SIMULATION_SCENE
#undef USING_GL
#include "util/FILE_IO.h"
#include "Scenes/JSON_SCENE.h"
#include "Scenes/TWO_TET_DROP_VF.h"
#include "Scenes/TWO_TET_DROP_EE.h"
#include "Scenes/TWO_TET_KISS_VF.h"
#include "Scenes/TWO_TET_KISS_EE.h"
#include "Scenes/PINNED_C.h"
#include "Scenes/DROPPED_E.h"
#include "Scenes/ARM_BEND.h"
#include "Scenes/BUNNY_DROP.h"
#include "Scenes/DROP_TEST.h"
#include "Scenes/NEWMARK_STRETCH.h"
#include "Scenes/QUASISTATIC_COLLISIONS.h"
#include "Scenes/QUASISTATIC_STRETCH.h"

#include "catch_amalgamated.hpp"

using namespace HOBAK;
using namespace std;

bool sceneRegression(SIMULATION_SCENE& scene, const string& groundTruthFilename)
{
  // read in the ground truth
  JSON_SCENE groundTruth;
  bool success = readSceneJSON(groundTruthFilename.c_str(), groundTruth);
  if (!success) return false;

  const int totalFrames = groundTruth.totalFrames();

  // build the current scene
  scene.buildScene();
  groundTruth.buildScene();

  cout << " =========================================================================== " << endl;
  cout << " Regression testing scene: " << scene.sceneName().c_str() << endl;
  cout << " =========================================================================== " << endl;

  // get the tetMeshes
  const TET_MESH* sceneMesh = scene.tetMesh();
  const TET_MESH* groundTruthMesh = groundTruth.tetMesh();

  // compare, frame-by-frame
  for (int x = 0; x < totalFrames; x++)
  {
    // step both
    scene.stepSimulation(false);
    //scene.stepSimulation(true);
    groundTruth.stepSimulation(false);

    // get the positions for both
    const VECTOR sceneDisplacement = sceneMesh->getDisplacement();
    const VECTOR groundTruthDisplacement = groundTruthMesh->getDisplacement();

    const VECTOR diff = sceneDisplacement - groundTruthDisplacement;
    const REAL relativeDiff = diff.norm() / groundTruthDisplacement.norm();

    if (relativeDiff > 1e-8)
    {
      cout << "Scene " << scene.sceneName().c_str() << " regression FAILED on frame " << x << endl;
      cout << "RELATIVE diff " << x << ": " << relativeDiff << "\t Absolute diff: " << diff.norm() << "\t simulation: " << sceneDisplacement.norm() << "\t ground truth: " << groundTruthDisplacement.norm() << endl;

      cout << " diff:   " << diff.transpose() << endl;
      cout << " scene:  " << sceneDisplacement.transpose() << endl;
      cout << " ground: " << groundTruthDisplacement.transpose() << endl;
      return false;
    }
    else
    {
      //cout << " Relative diff " << x << ": " << relativeDiff << endl;
      cout << " Relative diff " << x << ": " << relativeDiff << "\t Absolute diff: " << diff.norm() << "\t simulation: " << sceneDisplacement.norm() << "\t ground truth: " << groundTruthDisplacement.norm() << endl;
      //cout << " diff:   " << diff.transpose() << endl;
      //cout << " scene:  " << sceneDisplacement.transpose() << endl;
      //cout << " ground: " << groundTruthDisplacement.transpose() << endl;
      //exit(0);
    }
  }
  return true;
}

TEST_CASE("Scene regression tests", "[Scene regression testing]")
{
  DROP_TEST dropTest;
  REQUIRE(sceneRegression(dropTest, "../data/regression/drop_test.json.gz"));

  QUASISTATIC_STRETCH quasistaticStretch;
  REQUIRE(sceneRegression(quasistaticStretch, "../data/regression/quasistatic_stretch.json.gz"));

  QUASISTATIC_COLLISIONS quasistaticCollisions;
  REQUIRE(sceneRegression(quasistaticCollisions, "../data/regression/quasistatic_collisions.json.gz"));

  NEWMARK_STRETCH newmarkStretch;
  REQUIRE(sceneRegression(newmarkStretch, "../data/regression/newmark_stretch.json.gz"));

  TWO_TET_DROP_VF twoTetDropVF;
  REQUIRE(sceneRegression(twoTetDropVF, "../data/regression/two_tet_drop_vf.json"));
  
  TWO_TET_DROP_EE twoTetDropEE;
  REQUIRE(sceneRegression(twoTetDropEE, "../data/regression/two_tet_drop_ee.json"));
  
  TWO_TET_KISS_VF twoTetKissVF;
  REQUIRE(sceneRegression(twoTetKissVF, "../data/regression/two_tet_kiss_vf.json"));

  TWO_TET_KISS_EE twoTetKissEE;
  REQUIRE(sceneRegression(twoTetKissEE, "../data/regression/two_tet_kiss_ee.json"));
  
  PINNED_C pinnedC;
  REQUIRE(sceneRegression(pinnedC, "../data/regression/pinned_c.json.gz"));
  
  DROPPED_E droppedE;
  REQUIRE(sceneRegression(droppedE, "../data/regression/dropped_e.json.gz"));
  
  ARM_BEND armBend;
  REQUIRE(sceneRegression(armBend, "../data/regression/arm_bend.json.gz"));

  BUNNY_DROP bunnyDrop;
  REQUIRE(sceneRegression(bunnyDrop, "../data/regression/bunny_drop.json.gz"));
}
