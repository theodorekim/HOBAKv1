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
#include "SETTINGS.h"
#include "MshLoader.h"
#include <vector>
#include <iostream>

using namespace std;

string meshname;
vector<VECTOR3> vertices;
vector<VECTOR4I> tets;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void writeTetMesh(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "w");
  cout << " Writing out tet mesh file: " << filename.c_str() << endl;
  cout.flush();
 
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    const VECTOR3 v = vertices[x];
    fprintf(file, "v %.17g %.17g %.17g\n", v[0], v[1], v[2]);
  }
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    const VECTOR4I tet = tets[x];
    fprintf(file, "t %i %i %i %i\n", tet[0], tet[1], tet[2], tet[3]);
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void ReadGmshFile(const string& filename)
{
  PyMesh::MshLoader loader(filename);

  const PyMesh::VectorF& nodes = loader.get_nodes();
  const PyMesh::VectorI& elements = loader.get_elements();

  cout << " nodes size: " << nodes.size() << endl;
  cout << " nodes size / 3: " << nodes.size() / 3 << endl;
  cout << " elements size: " << elements.size() << endl;
  cout << " elements size / 4: " << elements.size() / 4 << endl;

  vertices.clear();
  for (unsigned int x = 0; x < nodes.size() / 3; x++)
  {
      int x3 = 3 * x;
      VECTOR3 vertex;
      vertex[0] = nodes[x3];
      vertex[1] = nodes[x3 + 1];
      vertex[2] = nodes[x3 + 2];
      vertices.push_back(vertex);
  }
  
  tets.clear();
  for (unsigned int x = 0; x < elements.size() / 4; x++)
  {
      int x4 = 4 * x;
      VECTOR4I tet;
      tet[0] = elements[x4];
      tet[1] = elements[x4 + 1];
      tet[2] = elements[x4 + 2];
      tet[3] = elements[x4 + 3];
      tets.push_back(tet);
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cout << " USAGE: " << argv[0] << " <mesh name>" << endl;
    return 0;
  }

  meshname = string(argv[1]);
  cout << " Reading in file " << meshname.c_str() << endl;

  cout << " Gmsh file detected " << endl;
  ReadGmshFile(meshname.c_str());
  
  string tobjFilename = meshname + string(".tobj");
  writeTetMesh(tobjFilename.c_str());

  return 0;
}
