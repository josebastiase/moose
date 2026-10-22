//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowLineSink.h"

class SinglePhaseFluidProperties;
class MooseVariableFieldBase;

namespace libMesh
{
class System;
}

/**
 * Approximates a borehole by a sequence of Dirac Points, using Peaceman's well index.
 *
 * The borehole is controlled either by its bottom-hole pressure (control = bhp, the default), which
 * is then the function bottom_p_or_t, or by its rate (control = rate), in which case the
 * bottom-hole pressure is a scalar unknown, solved for together with the reservoir so that the flux
 * summed over the completions equals the requested rate, unless that would need more than bhp_max.
 *
 * The variable may be a finite-element or a finite-volume variable.  For finite elements the flux
 * is computed from the nodal PorousFlow materials, with a hand-coded Jacobian.  For finite volumes
 * it is computed from the AD quadpoint (cell) materials of the finite-volume PorousFlow objects,
 * and the residual and Jacobian are assembled from AD.
 */
class PorousFlowPeacemanBorehole : public PorousFlowLineSink
{
public:
  /**
   * Creates a new PorousFlowPeacemanBorehole
   * This reads the file containing the lines of the form
   * radius x y z
   * that defines the borehole geometry.
   * It also calculates segment-lengths and rotation matrices
   * needed for computing the borehole well constant
   */
  static InputParameters validParams();

  /**
   * Peaceman's form of the well constant of one borehole half segment, as used by
   * PorousFlowPeacemanBorehole, for finite elements and finite volumes
   * @param perm permeability (or thermal conductivity) tensor
   * @param rot rotation matrix that rotates the half-segment direction to the z axis
   * @param half_len length of the half segment
   * @param ele element containing the borehole point
   * @param rad borehole radius
   * @param re_constant dimensionless constant used in evaluating the borehole effective radius
   */
  static Real peacemanWellConstant(const RealTensorValue & perm,
                                   const RealTensorValue & rot,
                                   Real half_len,
                                   const Elem * ele,
                                   Real rad,
                                   Real re_constant);

  PorousFlowPeacemanBorehole(const InputParameters & parameters);

  void virtual initialSetup() override;
  void virtual timestepSetup() override;
  void virtual residualSetup() override;
  void virtual jacobianSetup() override;

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar_num) override;
  virtual void computeResidualAndJacobian() override;

  /// Whether the borehole is rate controlled
  bool rateControlled() const { return _control == Control::rate; }

  /// Whether this borehole contributes its flux to the bottom-hole-pressure constraint
  bool contributesToConstraint() const { return rateControlled() && _compute_scalar_residuals; }

  /// Number of the bottom-hole-pressure scalar variable (rate control only)
  unsigned int bottomHolePressureVariable() const { return _kappa_var; }

  /// Whether the borehole is holding bhp_max instead of delivering its rate
  bool clamped() const { return _clamped; }

protected:
  /// Whether the bottom-hole pressure is prescribed or solved for to deliver a rate
  const enum class Control { bhp, rate } _control;

  /**
   * If positive then the borehole acts as a sink (producion well) for porepressure > borehole
   * pressure, and does nothing otherwise
   * If negative then the borehole acts as a source (injection well) for porepressure < borehole
   * pressure, and does nothing otherwise
   * The flow rate to/from the borehole is multiplied by |character|, so usually character = +/- 1
   *
   * A derived class whose borehole pressure is an unknown suppresses this parameter, since the
   * direction of flow is then settled by the constraint rather than prescribed, so this may be
   * nullptr
   */
  const Function * const _character;

  /// Bottomhole pressure of borehole, or nullptr if a derived class solves for it instead
  const Function * const _p_bot;

  /// Unit weight of fluid in borehole (for calculating bottomhole pressure at each Dirac Point).
  /// Not used if _use_density_from_temperature is true.
  const RealVectorValue _unit_weight;

  /// Whether the wellbore pressure profile is built from a temperature-dependent fluid density
  /// (true if the 'unit_weight_fp' parameter was supplied) rather than from the constant _unit_weight
  const bool _use_density_from_temperature;

  /// Fluid properties used to evaluate the in-well fluid density.  nullptr unless
  /// _use_density_from_temperature
  const SinglePhaseFluidProperties * const _fp;

  /// The coupled unit_weight_temperature variable.  nullptr unless
  /// _use_density_from_temperature
  const MooseVariableFieldBase * const _temperature_var;

  /// The libMesh system holding _temperature_var, used to sample that variable at borehole
  /// points that may lie outside this processor's mesh partition.  nullptr unless
  /// _use_density_from_temperature
  const libMesh::System * const _temperature_system;

  /// Variable number of unit_weight_temperature within _temperature_system
  const unsigned int _temperature_var_number;

  /// Gravitational acceleration (in the units used elsewhere in the input file), pointing
  /// downwards.  Only used if _use_density_from_temperature
  const RealVectorValue _gravity;

  /// Fixed pressure (Pa) at which the in-well fluid density is evaluated.  Only used if
  /// _use_density_from_temperature
  const Real _density_reference_pressure;

  /// Conversion of unit_weight_temperature's values to Kelvin (0 for Kelvin, 273.15 for Celsius)
  const Real _t_c2k;

  /**
   * Wellbore pressure at each well point, indexed by Dirac point ID.  Only used if
   * _use_density_from_temperature.  Recomputed once per residual evaluation and once per
   * Jacobian evaluation (see residualSetup()/jacobianSetup()) so that the residual always uses
   * the fluid density implied by the current nonlinear iterate's temperature.
   */
  std::vector<Real> _bh_pressure;

  /**
   * (Re)computes _bh_pressure from the temperature at each well point, when
   * _use_density_from_temperature is true.  Does nothing otherwise.
   */
  void computeWellborePressures();

  /**
   * The wellbore pressure (or temperature, for function_of=temperature) at the given Dirac
   * point, ie P_bot + unit_weight.(x_i - x_bottom) in the constant-unit_weight case, or the
   * temperature-dependent-density-integrated profile otherwise.
   */
  Real wellborePressure(unsigned current_dirac_ptid) const;

  /// Borehole constant
  const Real _re_constant;

  /// Well constant
  const Real _well_constant;

  /// Whether there is a quadpoint permeability material (for error checking)
  const bool _has_permeability;

  /// Whether there is a quadpoint thermal conductivity material (for error checking)
  const bool _has_thermal_conductivity;

  /// Permeability or conductivity of porous material (finite elements)
  const MaterialProperty<RealTensorValue> * const _perm_or_cond;

  /// d(Permeability)/d(PorousFlow variable) (finite elements)
  const MaterialProperty<std::vector<RealTensorValue>> * const _dperm_or_cond_dvar;

  /// Rotation matrix used in well_constant calculation
  std::vector<RealTensorValue> _rot_matrix;

  /**
   * Calculates Peaceman's form of the borehole well constant
   * Z Chen, Y Zhang, Well flow models for various numerical methods, Int J Num Analysis and
   * Modeling, 3 (2008) 375-388
   */
  Real wellConstant(const RealTensorValue & perm,
                    const RealTensorValue & rot,
                    const Real & half_len,
                    const Elem * ele,
                    const Real & rad) const;

  /**
   * The sum of the Peaceman well constants of the half-segments either side of a Dirac point, which
   * is the conductance of that completion per unit pressure difference and per unit weighting.
   * Unlike computeQpBaseOutflow this carries no directional gate, so a derived class whose borehole
   * pressure is an unknown can use it to form a flux that is differentiable in that unknown.
   * @param current_dirac_ptid The ID assigned to the point by addPoints()
   */
  Real wellConstantSum(unsigned current_dirac_ptid) const;

  /// Permeability (or thermal conductivity) at the current quadpoint, from the FE or the FV material
  RealTensorValue permOrCond() const;

  Real computeQpBaseOutflow(unsigned current_dirac_ptid) const override;
  void computeQpBaseOutflowJacobian(unsigned jvar,
                                    unsigned current_dirac_ptid,
                                    Real & outflow,
                                    Real & outflowp) const override;

  ///@{ Rate control (control = rate)

  /**
   * The flux out of the porespace at the current Dirac point per unit pressure difference and per
   * unit weighting, for the test function of node _i (finite elements)
   */
  Real conductance(unsigned current_dirac_ptid) const;

  /// The pressure in the borehole at the current Dirac point, from the scalar unknown
  Real boreholePressure() const;

  virtual Real computeQpScalarResidual() override;
  virtual Real computeQpScalarJacobian() override;
  virtual Real computeQpScalarOffDiagJacobian(unsigned int jvar_num) override;
  virtual Real computeQpJacobianScalar() override;

  /**
   * Choose between delivering the rate and holding bhp_max, using what the last residual evaluation
   * of the previous timestep recorded.  Taken once per timestep, so the equation being solved never
   * changes part way through a solve
   */
  void updateConstraint();

  /// Records a contribution to the rate-constraint residual, and the bottom-hole pressure, in the owner
  void recordRateMode(Real rate_mode);

  /// The rate this borehole delivers, positive for production
  const Function * const _rate;

  /// The bottom-hole pressure this borehole may not exceed, or nullptr if it is unlimited
  const Function * const _bhp_max;

  /// How many times per timestep the borehole may switch between its two constraints
  const unsigned int _max_swaps;

  /**
   * The borehole object that owns the constraint: of the objects sharing the bottom-hole pressure
   * that contribute to the constraint, the first.  It adds the rate, decides whether bhp_max holds,
   * and gathers the rate-constraint residual of all of them
   */
  PorousFlowPeacemanBorehole * _owner = nullptr;

  /// Whether this object adds the rate to the constraint (ie it is the owner)
  bool _adds_rate = false;

  /// Whether the borehole is currently holding _bhp_max rather than delivering _rate (owner)
  bool _clamped = false;

  /// Switches made so far this timestep (owner)
  unsigned int _n_swaps = 0;

  /**
   * The residual the rate constraint would have (summed flux - rate), accumulated over the
   * contributing objects and completions during the latest residual evaluation (owner).  It is
   * recorded even while the borehole holds bhp_max, because it says when the limit stops binding:
   * the flux falls with rising bottom-hole pressure, so a negative value means the rate would be
   * delivered below the limit
   */
  Real _rate_mode_residual = 0.0;

  /// The bottom-hole pressure seen during the latest residual evaluation (owner)
  Real _pbh_last = 0.0;
  ///@}

  ///@{ Finite volumes: AD quadpoint (cell) materials and the AD assembly

  /// Pressure (function_of = pressure) or temperature (function_of = temperature) in the cell
  ADReal fvPorT() const;

  /// The flux out of the porespace at borehole point i, which lies in the current cell
  ADReal fvPointOutflow(unsigned int i);

  /**
   * Sums over the borehole points in the current cell: the flux into this object's variable, its
   * contribution to the bottom-hole-pressure constraint, and to the rate-constraint residual
   */
  void fvCellResiduals(ADReal & outflow,
                       ADReal & scalar_residual,
                       Real & rate_mode,
                       bool record_point_fluxes);

  void fvComputeResidual();
  void fvComputeJacobian();
  void fvComputeResidualAndJacobian();

  const ADMaterialProperty<std::vector<Real>> * const _ad_pp;
  const ADMaterialProperty<Real> * const _ad_temperature;
  const ADMaterialProperty<RealTensorValue> * const _ad_perm_or_cond;
  const ADMaterialProperty<std::vector<Real>> * const _ad_density;
  const ADMaterialProperty<std::vector<Real>> * const _ad_viscosity;
  const ADMaterialProperty<std::vector<Real>> * const _ad_relperm;
  const ADMaterialProperty<std::vector<std::vector<Real>>> * const _ad_mass_fractions;
  const ADMaterialProperty<std::vector<Real>> * const _ad_enthalpy;
  const ADMaterialProperty<std::vector<Real>> * const _ad_internal_energy;
  const ADVariableValue * const _ad_multiplying_var;

  /// The bottom-hole pressure as an AD value (finite volumes, control = rate)
  const ADVariableValue * const _ad_bottom_hole_pressure;
  ///@}
};
