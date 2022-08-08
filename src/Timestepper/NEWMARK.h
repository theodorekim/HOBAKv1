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
#ifndef NEWMARK_H
#define NEWMARK_H

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
class NEWMARK : public TIMESTEPPER
{
public:
  NEWMARK(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic);
  const REAL rayleighAlpha() const  { return _rayleighAlpha; };
  const REAL rayleighBeta() const   { return _rayleighBeta; };

  const int maxNewtonIterations() const { return _maxNewtonIterations; };
  const VECTOR velocityOld() const      { return _velocityOld; };
  VECTOR& velocityOld()      { return _velocityOld; };
  int& maxNewtonIterations() { return _maxNewtonIterations; };
  virtual void setDt(const REAL dt) override;

  // take a timestep
  virtual bool solve(const bool verbose) override;
  virtual bool solveRayleighDamped(const bool verbose);
  virtual bool solveEnergyDamped(const bool verbose);

  // output all Newton steps to OBJ files?
  // this is for just-for-rfun visualizations
  bool& outputNewton() { return _outputNewton; };

protected:
  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  virtual void updateConstraintTargets() override;

  int _minNewtonIterations;
  int _maxNewtonIterations;
  int _seenNewtonIterations;
  REAL _residualTolerance;
  
  // Newmark constants
  REAL _alpha[6];
  //REAL _accelerationAlpha[5];
  REAL _beta;
  REAL _gamma;
  
  // current simulation time
  REAL _time;
 
  // current simulation step
  int _currentTimestep;

  // variables to solve for
  VECTOR _acceleration;
  VECTOR _velocityOld;
  VECTOR _accelerationOld;

  // output all Newton steps to OBJ files?
  bool _outputNewton;
};

} // HOBAK
} // TIMESTEPPER

#endif
