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
// MIN_HEAP.cpp: implementation of the MIN_HEAP class.
//
//////////////////////////////////////////////////////////////////////

#include "MIN_HEAP.h"

namespace HOBAK {

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MIN_HEAP::MIN_HEAP() :
  _size(0)
{

}

MIN_HEAP::~MIN_HEAP()
{
}

//////////////////////////////////////////////////////////////////////
// insert an element
//////////////////////////////////////////////////////////////////////
void MIN_HEAP::insert(HEAP_ENTRY& cell)
{
  _size++;
  if ((unsigned int)_size > _heap.size())
  {
    cell.heapIndex = _size - 1;
    _heap.push_back(cell);
    _heapIndex[cell.index] = _size - 1;
  }
  else
  {
    cell.heapIndex = _size - 1;
    _heap[_size - 1] = cell;
    _heapIndex[cell.index] = _size - 1;
  }

  decreaseKey(cell.index, _heap[_size - 1].distance);
}
  
//////////////////////////////////////////////////////////////////////
// pop the max off
//////////////////////////////////////////////////////////////////////
HEAP_ENTRY MIN_HEAP::popMin()
{
  if (_size == 0)
    assert(0);

  HEAP_ENTRY currentMin = heapMin();
  swap(0, _size - 1);
  _size--;
  heapify(0);

  return currentMin;
}
  
//////////////////////////////////////////////////////////////////////
// increase key of entry
//////////////////////////////////////////////////////////////////////
void MIN_HEAP::decreaseKey(int toChange, REAL newKey)
{
  int index = _heapIndex[toChange];
  if (fabs(newKey) > fabs(_heap[index].distance))
    return;

  _heap[index].distance = newKey;

  while (index > 0 && 
         (fabs(_heap[parent(index)].distance) > fabs(_heap[index].distance)) &&
         _heap[parent(index)].index < _heap[index].index) // enforce lexicographic?
  {
    assert(index < int(_heap.size()));
    assert(index > 0);
    assert(parent(index) >= 0);
    assert(parent(index) < int(_heap.size()));

    swap(index, parent(index));
    index = parent(index);
  }
}
  
//////////////////////////////////////////////////////////////////////
// get min 
//////////////////////////////////////////////////////////////////////
HEAP_ENTRY MIN_HEAP::heapMin()
{
  return _heap[0];
}

//////////////////////////////////////////////////////////////////////
// swap two vector elements
//////////////////////////////////////////////////////////////////////
void MIN_HEAP::swap(int index1, int index2)
{
  // update heap
  HEAP_ENTRY temp = _heap[index1];
  _heap[index1] = _heap[index2];
  _heap[index2] = temp;

  // update hash table
  int tempIndex = _heapIndex[index1];
  _heapIndex[index1] = _heapIndex[index2];
  _heapIndex[index2] = tempIndex;
}

//////////////////////////////////////////////////////////////////////
// enforce heap property
//////////////////////////////////////////////////////////////////////
void MIN_HEAP::heapify(int index)
{
  int leftIndex = left(index);
  int rightIndex = right(index);

  int largest = index;
  if (leftIndex < _size)
    if (fabs(_heap[leftIndex].distance) < fabs(_heap[index].distance))
      largest = leftIndex;
        
  if (rightIndex < _size)
    if (fabs(_heap[rightIndex].distance) < fabs(_heap[largest].distance))
      largest = rightIndex;

  if (largest != index)
  {
    swap(largest, index);
    heapify(largest);
  }
}

//////////////////////////////////////////////////////////////////////
// print out heap
//////////////////////////////////////////////////////////////////////
void MIN_HEAP::print()
{
  for (int x = 0; x < _size; x++)
    cout << __FILE__ << " " << __LINE__ << " element " << x << " : " << _heap[x].distance << endl;
}

} // HOBAK
