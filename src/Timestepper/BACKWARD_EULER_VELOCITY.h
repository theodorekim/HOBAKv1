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
#ifndef BACKWARD_EULER_VELOCITY_H
#define BACKWARD_EULER_VELOCITY_H

#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Timestepper/TIMESTEPPER.h"

namespace HOBAK {
namespace TIMESTEPPER {

////////////////////////////////////////////////////////////////////////////////////////////////////
// This is an implementation of the Baraff-Witkin-style velocity-level solver from
// "Large Steps in Cloth Simulation", SIGGRAPH 1998
//
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
////////////////////////////////////////////////////////////////////////////////////////////////////
class BACKWARD_EULER_VELOCITY : public TIMESTEPPER
{
public:
  BACKWARD_EULER_VELOCITY(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic);

  // take a timestep
  virtual bool solve(const bool verbose) override;
  bool solveRayleighDamped(const bool verbose);
  bool solveEnergyDamped(const bool verbose);

private:
  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  virtual void updateConstraintTargets() override;

  // Baraff-Witkin solves for change in velocity
  VECTOR _vDelta;

  // current simulation time
  REAL _time;
 
  // current simulation step
  int _currentTimestep;
};

} // HOBAK
} // TIMESTEPPER

#endif
