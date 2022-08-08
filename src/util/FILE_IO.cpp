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
#include "FILE_IO.h"
#include <stringbuffer.h>
#include <prettywriter.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <istreamwrapper.h>

// for flattening and unflattening matrices
#include "MATRIX_UTIL.h"

namespace HOBAK {

using namespace std;
using namespace rapidjson;

// the one and only static scene JSON Document
static rapidjson::Document _sceneDocument;
static rapidjson::Value _allFrames(kArrayType);

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single string
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const std::string& i, Value& container)
{
  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();
  Value key;
  key.SetString(i.c_str(), allocator);
  container.AddMember(StringRef(memberName), key, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single REAL
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const REAL& i, Value& container)
{
  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();
  container.AddMember(StringRef(memberName), i, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single bool
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const bool& i, Value& container)
{
  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();
  container.AddMember(StringRef(memberName), i, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single int
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const int& i, Value& container)
{
  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();
  container.AddMember(StringRef(memberName), i, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single VECTOR
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const VECTOR& v, Value& container)
{
  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();

  Value vArray(kArrayType);
  for (unsigned int x = 0; x < v.size(); x++)
    vArray.PushBack(v[x], allocator);
  
  container.AddMember(StringRef(memberName), vArray, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a vector of pairs 
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const vector<pair<int,int> >& v, Value& container)
{
  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();
  
  Value vArray(kArrayType);
  for (unsigned int x = 0; x < v.size(); x++)
  {
    vArray.PushBack(v[x].first, allocator);
    vArray.PushBack(v[x].second, allocator);
  }
  
  container.AddMember(StringRef(memberName), vArray, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single VECTOR9
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const VECTOR9& v, Value& container)
{
  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();

  Value vArray(kArrayType);
  for (unsigned int x = 0; x < 9; x++)
    vArray.PushBack(v[x], allocator);
  
  container.AddMember(StringRef(memberName), vArray, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// read a single bool
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void read(const char* memberName, bool& s, const Value& container)
{
  assert(container.HasMember(memberName));
  const Value& member = container[memberName];
  s = member.GetBool();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// read a single REAL
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void read(const char* memberName, REAL& s, const Value& container)
{
  assert(container.HasMember(memberName));
  const Value& member = container[memberName];
  s = member.GetFloat();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// read a single VECTOR9
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void read(const char* memberName, VECTOR9& v, const Value& container)
{
  assert(container.HasMember(memberName));
  const Value& member = container[memberName];
  assert(member.IsArray());
  for (unsigned int x = 0; x < 9; x++)
    v[x] = member[x].GetDouble();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// read a vector of integer pairs
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void read(const char* memberName, vector<pair<int,int> >& v, Value& container)
{
  v.clear();

  const Value& member = container[memberName];
  assert(member.IsArray());

  // its for pairs, so there needs to be an even number of integers
  assert(member.Size() % 2 == 0);

  v.resize(member.Size() / 2);
  for (unsigned int x = 0; x < member.Size() / 2; x++)
  {
    v[x].first  = member[2 * x].GetDouble();
    v[x].second = member[2 * x + 1].GetDouble();
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// read a single VECTOR
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void read(const char* memberName, VECTOR& v, Value& container)
{
  assert(container.HasMember(memberName));
  const Value& member = container[memberName];
  assert(member.IsArray());
  v.resize(member.Size());

  for (unsigned int x = 0; x < member.Size(); x++)
    v[x] = member[x].GetDouble();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// read a single VECTOR3
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void read(const char* memberName, VECTOR3& v, const Value& container)
{
  assert(container.HasMember(memberName));
  const Value& member = container[memberName];
  assert(member.IsArray());
  v[0] = member[0].GetDouble();
  v[1] = member[1].GetDouble();
  v[2] = member[2].GetDouble();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// read a single string 
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void read(const char* memberName, string& s, const Value& container)
{
  assert(container.HasMember(memberName));
  const Value& member = container[memberName];
  s = member.GetString();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single MATRIX3
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void read(const char* memberName, MATRIX3& A, const Value& container)
{
  VECTOR9 a;
  read(memberName, a, container);
  A = unflatten(a);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single MATRIX3
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const MATRIX3& A, Value& container)
{
  VECTOR9 a = flatten(A);
  write(memberName, a, container);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a single VECTOR3
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const VECTOR3& v, Value& container)
{
  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();

  Value vArray(kArrayType);
  vArray.PushBack(v[0], allocator);
  vArray.PushBack(v[1], allocator);
  vArray.PushBack(v[2], allocator);
  
  container.AddMember(StringRef(memberName), vArray, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// read in a shape
///////////////////////////////////////////////////////////////////////////////////////////////////////
KINEMATIC_SHAPE* read(const Value& shapeContainer)
{
  string name;
  VECTOR3 translation;
  MATRIX3 rotation, scale, scaleInverse;

  read("name", name, shapeContainer);
  read("translation", translation, shapeContainer);
  read("rotation", rotation, shapeContainer);
  read("scale", scale, shapeContainer);
  read("scaleInverse", scaleInverse, shapeContainer);

  KINEMATIC_SHAPE* shape = NULL;
  if (name.compare("CUBE") == 0)
    shape = new CUBE(translation, scale(0,0));
  else if (name.compare("SPHERE") == 0)
    shape = new SPHERE(translation, scale(0,0));
  else if (name.compare("CYLINDER") == 0 || name.compare("CAPSULE") == 0)
  {
    REAL radius, height;
    read("radius", radius, shapeContainer);
    read("height", height, shapeContainer);

    if (name.compare("CYLINDER") == 0)
      shape = new CYLINDER(translation, radius, height);
    else if (name.compare("CAPSULE") == 0)
      shape = new CAPSULE(translation, radius, height);
  }

  assert(shape != NULL);
  shape->rotation() = rotation;
  return shape;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a sphere
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const KINEMATIC_SHAPE* shape, Value& kinematicsContainer)
{
  Value shapeContainer;
  shapeContainer.SetObject();
 
  write("name", shape->name(), shapeContainer);
  write("translation", shape->translation(), shapeContainer);
  write("scale", shape->scale(), shapeContainer);
  write("scaleInverse", shape->scaleInverse(), shapeContainer);
  write("rotation", shape->rotation(), shapeContainer);

  if ((shape->name().compare("CYLINDER") == 0) ||
      (shape->name().compare("CAPSULE") == 0))
  {
    // CAPSULE is a subclass of CYLINDER, so that's covered too
    const CYLINDER* cylinder = (CYLINDER*)shape;
  
    write("radius", cylinder->radius(), shapeContainer);
    write("height", cylinder->height(), shapeContainer);
  }

  // add the shape to the kinematics list
  kinematicsContainer.PushBack(shapeContainer, _sceneDocument.GetAllocator());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write a kinematic shape
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const vector<KINEMATIC_SHAPE*>& shapes, Value& sceneContainer)
{
  Value kinematicsContainer(kArrayType);

  // add each shape to the kinematics container
  for (unsigned int x = 0; x < shapes.size(); x++)
    write(shapes[x], kinematicsContainer);
  
  // add the kinematics container to the scene
  sceneContainer.AddMember("kinematicShapes", kinematicsContainer, _sceneDocument.GetAllocator());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write an array of VECTOR4I
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const vector<VECTOR4I>& tets, Document& document, Value& container)
{
  Document::AllocatorType& allocator = document.GetAllocator();
  Value tetsArray(kArrayType);
  for (auto t = tets.begin(); t != tets.end(); t++) {
    Value tetArray(kArrayType);
    const VECTOR4I& tet = *t;
    tetArray.PushBack(tet[0], allocator);
    tetArray.PushBack(tet[1], allocator);
    tetArray.PushBack(tet[2], allocator);
    tetArray.PushBack(tet[3], allocator);
    tetsArray.PushBack(tetArray, allocator);
  }
  container.AddMember(StringRef(memberName), tetsArray, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write an array of VECTOR3
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void write(const char* memberName, const vector<VECTOR3>& vertices, Document& document, Value& container)
{
  Document::AllocatorType& allocator = document.GetAllocator();
  Value verticesArray(kArrayType);
  for (auto v = vertices.begin(); v != vertices.end(); v++) {
    Value v3Array(kArrayType);
    const VECTOR3& vertex = *v;
    v3Array.PushBack(vertex[0], allocator);
    v3Array.PushBack(vertex[1], allocator);
    v3Array.PushBack(vertex[2], allocator);
    verticesArray.PushBack(v3Array, allocator);
  }
  container.AddMember(StringRef(memberName), verticesArray, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Write out TET_MESH state to a JSON file
///////////////////////////////////////////////////////////////////////////////////////////////////////
void writeJSON(const std::string& filename, const TET_MESH& tetMesh)
{
  using namespace rapidjson;

  Document d;
  d.SetObject();
  Document::AllocatorType &allocator = d.GetAllocator();

  Value tetMeshContainer;
  tetMeshContainer.SetObject();
  write("vertices", tetMesh.vertices(), d, tetMeshContainer);
  write("restVertices", tetMesh.restVertices(), d, tetMeshContainer);
  write("tets", tetMesh.tets(), d, tetMeshContainer);
  d.AddMember("TetMesh", tetMeshContainer, allocator);

  cout << " Writing out JSON file " << filename.c_str() << endl;
  StringBuffer s;
  PrettyWriter<StringBuffer> writer(s);
  //Writer<StringBuffer> writer(s);
  d.Accept(writer);

  ofstream file(filename.c_str());
  file << s.GetString() << endl;
  file.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize the JSON file for a simulation scene
///////////////////////////////////////////////////////////////////////////////////////////////////////
void initializeSceneJSON(const SIMULATION_SCENE& scene)
{
  _sceneDocument.SetObject();

  // write out the frame-independent data
  Value sceneContainer;
  sceneContainer.SetObject();
  write("eye",          scene.eye(),    sceneContainer);
  write("lookAt",       scene.lookAt(), sceneContainer);
  write("up",           scene.up(),     sceneContainer);
  write("worldCenter",  scene.worldCenter(),  sceneContainer);
  write("gravity",      scene.gravity(),      sceneContainer);
  write("pauseFrame",   scene.pauseFrame(),   sceneContainer);
  write("arrowCounter", scene.arrowCounter(), sceneContainer);
  write("drawFeature",  scene.drawFeature(),  sceneContainer);
  write("tetMeshFilename", scene.tetMeshFilename(), sceneContainer);
  write("normalizedVertices", scene.normalizedVertices(), sceneContainer);
  write("initialA", scene.initialA(), sceneContainer);
  write("initialTranslation",  scene.initialTranslation(),  sceneContainer);
 
  // time stepped frame-independent data
  write("dt",  scene.solver()->dt(),  sceneContainer);
  write("vertexFaceSelfCollisionsOn",  scene.solver()->vertexFaceSelfCollisionsOn(),  sceneContainer);
  write("edgeEdgeSelfCollisionsOn",  scene.solver()->edgeEdgeSelfCollisionsOn(),  sceneContainer);
  write("collisionStiffness",  scene.solver()->collisionStiffness(),  sceneContainer);
  write("collisionDampingBeta",  scene.solver()->collisionDampingBeta(),  sceneContainer);
  write("materialName",  scene.solver()->materialName(),  sceneContainer);

  // write out all the kinematic shapes
  write(scene.kinematicShapes(), sceneContainer);

  Document::AllocatorType& allocator = _sceneDocument.GetAllocator();
  _sceneDocument.AddMember("scene", sceneContainer, allocator);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Record a frame of simulation state to a JSON file
///////////////////////////////////////////////////////////////////////////////////////////////////////
void recordFrameToJSON(const SIMULATION_SCENE& scene)
{
  // make a container for just this frame
  Value singleFrame;
  singleFrame.SetObject();

  const TIMESTEPPER::TIMESTEPPER* solver = scene.solver();
  write("frameNumber", scene.frameNumber(), singleFrame);
  write("position", solver->position(), singleFrame);
  write("positionOld", solver->positionOld(), singleFrame);
  write("velocity", solver->velocity(), singleFrame);
  write("externalForces", solver->externalForces(), singleFrame);
  
  const TET_MESH* tetMesh = scene.tetMesh();
  write("vertexFaceCollisions", tetMesh->vertexFaceCollisions(), singleFrame);
  write("edgeEdgeCollisions", tetMesh->edgeEdgeCollisions(), singleFrame);
 
  _allFrames.PushBack(singleFrame, _sceneDocument.GetAllocator()); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write the final state to a JSON file
///////////////////////////////////////////////////////////////////////////////////////////////////////
void writeSceneJSON(const char* filename)
{
  // add all the frames to the scene
  _sceneDocument.AddMember("frames", _allFrames, _sceneDocument.GetAllocator());

  cout << " Writing out JSON file " << filename << endl;
  StringBuffer s;

  // human-readable
  //PrettyWriter<StringBuffer> writer(s);

  // less space-hogging
  Writer<StringBuffer> writer(s);

  _sceneDocument.Accept(writer);

  ofstream file(filename);
  file << s.GetString() << endl;
  file.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// is it a gzip file?
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool isGzipped(const char* name)
{
  string filename(name);
  const int length = filename.length();

  if (length < 3) return false;

  if (filename[length - 1] != 'z')
    return false;
  if (filename[length - 2] != 'g')
    return false;
  if (filename[length - 3] != '.')
    return false;

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// unzip the JSON file, then read it
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool readGzippedJSON(const char* filename, JSON_SCENE& scene)
{
  // make a copy of the original gzip
  string cp("cp ");
  cp = cp + string(filename) + string(" temp.json.gz");
  system(cp.c_str());

  // unzip the copy
  cout << " Unzipping ... " << flush;
  system("gunzip temp.json.gz");
  cout << "done." << endl;

  // read in the scene
  bool success = readSceneJSON("temp.json", scene);

  // delete the unzipped version
  system("rm temp.json");

  return success;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// write the final state to a JSON file
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool readSceneJSON(const char* filename, JSON_SCENE& scene)
{
  // if it's compressed, punt to decompressor
  if (isGzipped(filename))
    return readGzippedJSON(filename, scene);

  // the file does exist, right?
  ifstream ifs(filename);
  if (ifs.fail())
  {
    cout << " Failed to open file " << filename << "!!!" << endl;
    return false;
  }

  IStreamWrapper isw(ifs);
  Document document;
  cout << " Parsing JSON file ... " << flush;
  document.ParseStream(isw);
  cout << "done." << endl;
  assert(document.IsObject());

  assert(document.HasMember("scene"));
  const Value& sceneContainer = document["scene"];

  // read camera parameters
  read("eye", scene.eye(), sceneContainer);
  read("lookAt", scene.lookAt(), sceneContainer);
  read("up", scene.up(), sceneContainer);
  read("worldCenter", scene.worldCenter(), sceneContainer);

  // read simulation parameters
  read("gravity", scene.gravity(), sceneContainer);
  read("tetMeshFilename", scene.tetMeshFilename(), sceneContainer);
  read("normalizedVertices", scene.normalizedVertices(), sceneContainer);

  MATRIX3 A;
  VECTOR3 t; 
  read("initialA", A, sceneContainer);
  read("initialTranslation", t, sceneContainer);
  scene.setInitialA(A);
  scene.setInitialTranslation(t);

  cout << " eye:         " << scene.eye().transpose() << endl;
  cout << " lookAt:      " << scene.lookAt().transpose() << endl;
  cout << " up:          " << scene.up().transpose() << endl;
  cout << " worldCenter: " << scene.worldCenter().transpose() << endl;
  cout << " gravity:     " << scene.gravity().transpose() << endl;
  cout << " tetMeshFilename: " << scene.tetMeshFilename().c_str() << endl;
  cout.flush();

  // read in kinematic objects
  assert(sceneContainer.HasMember("kinematicShapes"));
  const Value& kinematicsContainer = sceneContainer["kinematicShapes"];
  assert(kinematicsContainer.IsArray());
  vector<KINEMATIC_SHAPE*>& shapes = scene.kinematicShapes();
  for (unsigned int x = 0; x < kinematicsContainer.Size(); x++)
  {
    KINEMATIC_SHAPE* shape = read(kinematicsContainer[x]);
    shapes.push_back(shape);
  }
  cout << " Found " << shapes.size() << " kinematic shapes " << endl;
  cout.flush();

  // read in simulation frames
  assert(document.HasMember("frames"));
  Value& framesContainer = document["frames"];

  // read in position frames
  vector<VECTOR>& positions = scene.positions();
  positions.resize(framesContainer.Size());
  for (unsigned int x = 0; x < framesContainer.Size(); x++)
    read("position", positions[x], framesContainer[x]);

  // they're all the same size, right?
  for (unsigned int x = 1; x < positions.size(); x++)
    assert(positions[x].size() == positions[0].size());
  cout << " Read in " << positions.size() << " position frames " << endl;
  cout.flush();
  
  // read in velocity frames
  vector<VECTOR>& velocities = scene.velocities();
  velocities.resize(framesContainer.Size());
  for (unsigned int x = 0; x < framesContainer.Size(); x++)
    read("velocity", velocities[x], framesContainer[x]);

  // read in collision histories
  vector<vector<pair<int, int> > >& vertexFaceCollisions = scene.vertexFaceCollisions();
  vector<vector<pair<int, int> > >& edgeEdgeCollisions = scene.edgeEdgeCollisions();
  vertexFaceCollisions.resize(framesContainer.Size());
  edgeEdgeCollisions.resize(framesContainer.Size());
  for (unsigned int x = 0; x < framesContainer.Size(); x++)
  {
    read("vertexFaceCollisions", vertexFaceCollisions[x], framesContainer[x]);
    read("edgeEdgeCollisions", edgeEdgeCollisions[x], framesContainer[x]);
  }

  // they're all the same size, right?
  for (unsigned int x = 1; x < velocities.size(); x++)
    assert(velocities[x].size() == velocities[0].size());
  cout << " Read in " << velocities.size() << " velocity frames " << endl;
  cout.flush();

  return true;
}


} // HOBAK
