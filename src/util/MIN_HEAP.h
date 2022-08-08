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
// MIN_HEAP.h: interface for the MIN_HEAP class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MIN_HEAP_H
#define MIN_HEAP_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#include "SETTINGS.h"

namespace HOBAK {

class HEAP_ENTRY {
public:
  HEAP_ENTRY() {
    distance = 0;
    index = 0;
    heapIndex = 0;
  };

  REAL distance;
  int index;
  int heapIndex;
};

// the actual heap
class MIN_HEAP  
{
public:
  
	MIN_HEAP();
	virtual ~MIN_HEAP();

  // heap ops
  void insert(HEAP_ENTRY& cell);
  void decreaseKey(int toChange, REAL newKey);
  HEAP_ENTRY popMin();
  HEAP_ENTRY heapMin();
  int size() { return _size; };

  // debugging
  void print();
 
  void clear() {
    _heap.clear();
    _heapIndex.clear();
    _size = 0;
  };

  bool empty() { return _size == 0; };
  
private:
  // tree traversal
  int parent(int i) { return i / 2; };
  int left(int i) { return i * 2; };
  int right(int i) { return i * 2 + 1; };
 
  // the heap
  std::vector<HEAP_ENTRY> _heap;
  
  // hash table mapping grid index to heap index
  std::map<int,int> _heapIndex;

  // size of current heap
  int _size;
  
  // enforce heap property 
  void heapify(int index);

  // swap two entries
  void swap(int index1, int index2);
};

} // HOBAK

#endif
