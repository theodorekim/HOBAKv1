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
#ifndef BDF_2_H
#define BDF_2_H

#include "Timestepper/NEWMARK.h"

namespace HOBAK {
namespace TIMESTEPPER {

////////////////////////////////////////////////////////////////////////////////////////////////////
// Newton-style solver for second-order Backward Differentiation Formula (BDF-2).
//
// This is from Choi and Ko, "Stable but Responsive Cloth", SIGGRAPH 2002.
// Section 4, Equation 16.
////////////////////////////////////////////////////////////////////////////////////////////////////
class BDF_2 : public NEWMARK
{
public:
  BDF_2(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic);
  virtual void setDt(const REAL dt) override;

  // take a timestep
  virtual bool solveRayleighDamped(const bool verbose) override;
  virtual bool solveEnergyDamped(const bool verbose) override;

  VECTOR& velocityOlder() { return _velocityOlder; };
  const VECTOR velocityOlder() const { return _velocityOlder; };

private:

  // quantities from *two* timesteps ago
  VECTOR _positionOlder;
  VECTOR _velocityOlder;
};

} // HOBAK
} // TIMESTEPPER

#endif
