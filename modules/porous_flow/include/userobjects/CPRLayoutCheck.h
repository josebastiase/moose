//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"

/**
 * Checks, once at startup, that the dofs of the nonlinear system are interleaved per node (finite
 * elements) or per element (finite volumes), which is what '-pc_type cpr' assumes, and prints the
 * PETSc options to use.  Errors out with an explanation if the layout does not fit, so that CPR
 * never runs silently on the wrong dofs.
 */
class CPRLayoutCheck : public GeneralUserObject
{
public:
  static InputParameters validParams();
  CPRLayoutCheck(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void finalize() override {}

protected:
  /// The nonlinear variable that CPR treats as the pressure
  const NonlinearVariableName & _pressure_var;
};
