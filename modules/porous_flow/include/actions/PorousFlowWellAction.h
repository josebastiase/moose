//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObjectAction.h"

class PorousFlowDictator;

/**
 * Action that expands a [Wells] block into the DiracKernels, UserObjects and Postprocessors that
 * model the well.  Wells are operational forcing on the PorousFlow equations rather than PDE terms,
 * so this is a plain Action reading only the PorousFlowDictator: it works identically whether the
 * equations came from a PorousFlow Physics block or were built by hand.
 */
class PorousFlowWellAction : public MooseObjectAction
{
public:
  static InputParameters validParams();

  PorousFlowWellAction(const InputParameters & parameters);

  /**
   * MOOSE builds a separate Action instance for each task, so nothing may be carried between
   * add_porous_flow_well_gate and add_porous_flow_well in a member.  The names of the objects the
   * gate task builds are derived from name(), and are recomputed rather than stored.
   */
  virtual void act() override;

protected:
  /// Error on any parameter that requires the implicit scalar-BHP well of phase 1b
  void checkNotYetImplemented();

  /// Whether this well's bottom-hole pressure is solved for rather than implied by a prescribed rate
  bool solvesForBottomHolePressure() const;

  /// The name of the scalar variable holding this well's bottom-hole pressure
  VariableName bottomHolePressureVariable() const;

  /**
   * Add the scalar variable holding the bottom-hole pressure, along with the NullScalarKernel that
   * lets the kernel coverage check see it.  The coverage check builds its list of covered variables
   * from Kernels, ScalarKernels, finite volume objects and IntegratedBCs, never from DiracKernels,
   * so the variable would otherwise look uncovered even though a DiracKernel assembles its
   * equation.
   */
  void addBottomHolePressureVariable();

  /**
   * The variable holding the mass balance of fluid component `component`.  The PorousFlow
   * Physics/Action convention is that component i < N is on mass_fraction_vars[i] and the final
   * component N is on porepressure; the PorousFlowDictator does not record this mapping.
   */
  VariableName componentVariable(unsigned int component) const;

  /// Build the objects modelling a PorousFlowInjectionWell
  void addInjectionWell();

  /**
   * Build an injection well whose bottom-hole pressure is solved for: one borehole owning the rate
   * constraint, and one per fluid component contributing flux, all sharing the scalar variable.
   * The owner carries no mass_fraction_component, so its flux is the total the rate refers to,
   * and since the injected fractions sum to one the component fluxes sum to the same thing.
   * @param rate The well's total rate function, already gated
   * @param mass_fractions The composition of the injected stream
   */
  void addConstrainedInjectionWell(const FunctionName & rate,
                                   const std::vector<Real> & mass_fractions);

  /// Build the objects modelling a PorousFlowProductionWell
  void addProductionWell();

  /// Whether the fluid-component variables are finite-volume variables
  bool finiteVolume() const;

  /**
   * A ParsedFunction that is 1 while the well is open and 0 while it is shut, combining the
   * start_time/end_time window with the shut_if condition.  Returns an empty name if the well is
   * never shut.
   *
   * A well is shut by zeroing its flux, never by disabling its objects with a Control.  A disabled
   * DiracKernel stops contributing, but nothing zeroes the objects that report on it: a disabled
   * Postprocessor is simply not executed and holds its last value, and a disabled
   * PorousFlowLineSink never re-runs the addPoints() that zeroes its PorousFlowSumQuantity.  Either
   * way the rates and cumulatives would go stale rather than to zero.
   */
  FunctionName gateFunction();

  /// The rate function to use: `rate` multiplied by `gate`, unless `gate` is empty
  FunctionName gatedRate(const FunctionName & gate, const FunctionName & rate);

  /// Insert the well's schedule breakpoints into the simulation's sync times
  void addSyncTimes();

  /// Copy the geometry parameters of a production well into a PorousFlowLineGeometry-derived object
  void setLineGeometry(InputParameters & pars) const;

  /**
   * An AuxVariable carrying the gate, for use as a line sink's `multiplying_var`.  A producer's
   * rate comes from its deliverability curve rather than from a rate function, so there is no rate
   * to zero: the flux itself must be multiplied by the gate.  Returns an empty name if the well is
   * never shut.
   */
  VariableName gateAuxVariable(const FunctionName & gate);

  /// Add a FunctionValuePostprocessor evaluating `func`, scaled by `scale`, at the start of the step
  void
  addRatePostprocessor(const PostprocessorName & pp_name, const FunctionName & func, Real scale);

  /// Whether this well is ever shut in, and so needs a gate
  bool isGated() const;

  /// Number of fluid components in the simulation
  unsigned int _num_components;
};
