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
#ifndef QUASISTATIC_H
#define QUASISTATIC_H

#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Timestepper/TIMESTEPPER.h"

namespace HOBAK {
namespace TIMESTEPPER {

class QUASISTATIC : public TIMESTEPPER
{
public:
  QUASISTATIC(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic);

  // take a timestep
  //virtual bool solve(const bool verbose) override { return solve(verbose, 0.01); };
  virtual bool solve(const bool verbose) override { return solve(verbose, 1.0); };

  // take a step with a prescribed stepsize
  bool solve(const bool verbose, const REAL stepSize);

  // take a timestep, but use a backtracking line search
  bool solveWithBacktracking(const bool verbose);

private:
  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  virtual void updateConstraintTargets() override;

  // filter positions to incorporate Baraff-Witkin-style constraints
  virtual void applyKinematicConstraints() override;

  // find all the surface vertices that are in collision and create constraints
  virtual void findNewSurfaceConstraints(const bool verbose) override;

  int _minIterations;
  int _maxIterations;
  REAL _residualTolerance;
  
  int _seenNewtonIterations;
};

} // HOBAK
} // TIMESTEPPER

#endif
