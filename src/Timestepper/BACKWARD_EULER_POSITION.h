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
#ifndef BACKWARD_EULER_POSITION_H
#define BACKWARD_EULER_POSITION_H

#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Damping/Volume/DAMPING.h"
#include "Timestepper/TIMESTEPPER.h"

namespace HOBAK {
namespace TIMESTEPPER {

////////////////////////////////////////////////////////////////////////////////////////////////////
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
////////////////////////////////////////////////////////////////////////////////////////////////////
class BACKWARD_EULER_POSITION : public TIMESTEPPER
{
public:
  BACKWARD_EULER_POSITION(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic);
  BACKWARD_EULER_POSITION(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, VOLUME::DAMPING& damping);
  virtual ~BACKWARD_EULER_POSITION();

  const REAL rayleighAlpha() const  { return _rayleighAlpha; };
  const REAL rayleighBeta() const   { return _rayleighBeta; };

  const VECTOR velocityOld() const      { return _velocityOld; };
  VECTOR& velocityOld()      { return _velocityOld; };
  virtual void setDt(const REAL dt) override { _dt = dt; };

  // take a timestep
  virtual bool solve(const bool verbose) override;
  bool solveEnergyDamped(const bool verbose);

  // solves with self collisions, and collision forces are Rayleigh damped
  bool solveRayleighDamped(const bool verbose);

private:
  // shared initialization across constructors
  void initialize();

  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  virtual void updateConstraintTargets() override;

  // current simulation time
  REAL _time;
 
  // current simulation step
  int _currentTimestep;

  // variables to solve for
  VECTOR _acceleration;
  VECTOR _velocityOld;
};

} // HOBAK
} // TIMESTEPPER

#endif
