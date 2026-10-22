//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPeacemanBorehole.h"
#include "MooseVariableScalar.h"
#include "NonlinearSystemBase.h"
#include "Assembly.h"
#include "RotationMatrix.h"
#include "Function.h"
#include "SinglePhaseFluidProperties.h"
#include "MooseVariableFieldBase.h"
#include "SystemBase.h"
#include "libmesh/system.h"

registerMooseObject("PorousFlowApp", PorousFlowPeacemanBorehole);

InputParameters
PorousFlowPeacemanBorehole::validParams()
{
  InputParameters params = PorousFlowLineSink::validParams();
  // A borehole point that the point locator cannot place in the mesh must not be silently
  // dropped: computeWellborePressures() samples every borehole point's temperature directly
  // (not just the ones that ended up with a Dirac source), and a missing element there corrupts
  // the trapezoidal pressure integral for every other point up the well, not just the missing
  // one.  DiracKernelBase's own default of IGNORE would let that happen silently, so raise it to
  // ERROR here.
  params.set<MooseEnum>("point_not_found_behavior") = "ERROR";
  MooseEnum control("bhp rate", "bhp");
  params.addParam<MooseEnum>(
      "control",
      control,
      "bhp: the bottom-hole pressure is prescribed by bottom_p_or_t.  rate: the bottom-hole "
      "pressure is the scalar variable bottom_hole_pressure, solved for so that the flux summed "
      "over the completions equals 'rate' (unless that would need more than bhp_max)");
  params.addParam<FunctionName>(
      "character",
      "Needed with control = bhp, not used with control = rate.  If zero then borehole does "
      "nothing.  If positive the borehole acts as a sink "
      "(production well) for porepressure > borehole pressure, and does nothing "
      "otherwise.  If negative the borehole acts as a source (injection well) for "
      "porepressure < borehole pressure, and does nothing otherwise.  The flow rate "
      "to/from the borehole is multiplied by |character|, so usually character = +/- "
      "1, but you can specify other quantities to provide an overall scaling to the "
      "flow if you like.");
  params.addParam<FunctionName>("bottom_p_or_t",
                                "Needed with control = bhp.  For function_of=pressure, this "
                                "function is the pressure at the bottom of the borehole, "
                                "otherwise it is the temperature at the bottom of the borehole.");
  params.renameCoupledVar("scalar_variable",
                          "bottom_hole_pressure",
                          "The scalar variable holding the pressure at the bottom of the borehole "
                          "(control = rate).  Its value is solved for, not prescribed: it is "
                          "whatever delivers 'rate'.  It needs a NullScalarKernel, because the "
                          "kernel coverage check does not see that DiracKernels assemble its "
                          "equation");
  params.addParam<FunctionName>(
      "rate",
      "The rate the borehole delivers (control = rate), summed over all of its completions, "
      "positive for production and negative for injection.  The units are those of the equation "
      "the borehole contributes to, so kg/s for a fluid mass equation.  Where a borehole is made "
      "of "
      "several objects sharing bottom_hole_pressure (eg one per fluid phase and component), each "
      "object that contributes to the constraint (compute_scalar_residuals = true) adds its flux "
      "to "
      "it, and the first adds the rate, so they must all give the same rate");
  params.addParam<FunctionName>(
      "bhp_max",
      "The bottom-hole pressure the borehole may not exceed (control = rate).  Where delivering "
      "'rate' would need more than this, the borehole holds this pressure instead and delivers "
      "whatever rate follows.  The constraint is decided once per timestep, so the bottom-hole "
      "pressure may exceed bhp_max for the timestep in which the limit is first reached");
  params.addParam<unsigned int>(
      "max_constraint_swaps",
      5,
      "The number of times in one timestep the borehole may switch between delivering 'rate' and "
      "holding 'bhp_max'.  Once this many switches have been made the borehole keeps whichever "
      "constraint it holds until the next timestep");
  params.addParam<RealVectorValue>(
      "unit_weight",
      "(fluid_density*gravitational_acceleration) as a vector pointing downwards.  "
      "Note that the borehole pressure at a given z position is bottom_p_or_t + "
      "unit_weight*(q - q_bottom), where q=(x,y,z) and q_bottom=(x,y,z) of the "
      "bottom point of the borehole.  The analogous formula holds for "
      "function_of=temperature.  If you don't want bottomhole pressure (or "
      "temperature) to vary in the borehole just set unit_weight=0.  Typical value "
      "is = (0,0,-1E4), for water.  Exactly one of 'unit_weight' or 'unit_weight_fp' must be "
      "given.  Use 'unit_weight_fp' instead of 'unit_weight' if you want the fluid unit weight "
      "to vary along the borehole according to a temperature-dependent fluid "
      "density, rather than being a single constant value.");
  params.addParam<UserObjectName>(
      "unit_weight_fp",
      "SinglePhaseFluidProperties UserObject used to evaluate the in-well fluid density from "
      "'unit_weight_temperature' at each borehole point, in order to build a wellbore pressure "
      "profile that accounts for a thermal gradient along the borehole.  Providing this "
      "parameter activates this temperature-dependent unit-weight mode instead of the constant "
      "'unit_weight'.  Not compatible with function_of=temperature, since in that mode "
      "bottom_p_or_t is a temperature, not a pressure, so there is no pressure profile to build.  "
      "If given, 'unit_weight_temperature' and 'unit_weight_gravity' are also required.  "
      "(Deliberately not named 'fp'/'gravity'/'temperature_variable', even though that mirrors "
      "convention elsewhere in PorousFlow, because those are common GlobalParams names: an "
      "input file that sets 'gravity' or 'fp' in [GlobalParams] for unrelated Darcy kernels or "
      "fluid materials would otherwise silently activate, or fail to validate, this mode on "
      "every PorousFlowPeacemanBorehole in the input.)");
  params.addCoupledVar(
      "unit_weight_temperature",
      "The (nonlinear or auxiliary) variable holding temperature, sampled at each borehole point "
      "to compute the in-well fluid density used to build the wellbore pressure profile.  Must "
      "be a variable, not a constant value.  This is unrelated to function_of=temperature "
      "(which instead selects whether the *outflow* driving this DiracKernel is a function of "
      "porepressure or temperature).  Only used, and required, if 'unit_weight_fp' is given.");
  params.addParam<RealVectorValue>(
      "unit_weight_gravity",
      "Gravitational acceleration, pointing downwards, in the units used elsewhere in this "
      "input file (eg (0,0,-9.81) for SI units and lengths in metres).  Only used, and "
      "required, if 'unit_weight_fp' is given.  Densities computed from 'unit_weight_fp' are "
      "always in kg/m^3, so if lengths in this input file are not metres, scale "
      "'unit_weight_gravity' accordingly (eg (0,0,-9.81E-6) if pressures are in MPa and lengths "
      "in metres, matching the 'gravity' convention used by PorousFlow Darcy kernels).");
  MooseEnum temperature_unit_choice("Kelvin=0 Celsius=1", "Kelvin");
  params.addParam<MooseEnum>(
      "unit_weight_temperature_unit",
      temperature_unit_choice,
      "The unit of 'unit_weight_temperature'.  Only used if 'unit_weight_fp' is given.");
  MooseEnum pressure_unit_choice("Pa MPa", "Pa");
  params.addParam<MooseEnum>(
      "unit_weight_pressure_unit",
      pressure_unit_choice,
      "The unit of 'unit_weight_reference_pressure'.  Only used if 'unit_weight_fp' is given.");
  params.addRangeCheckedParam<Real>(
      "unit_weight_reference_pressure",
      101325.0, // standard atmosphere: liquid-water density is only weakly pressure-dependent
                // (~0.9% per 20MPa), so a fixed reference pressure is used instead of the local
                // (coupled, and hence more expensive to sample) porepressure
      "unit_weight_reference_pressure > 0",
      "The fixed pressure (in the units given by 'unit_weight_pressure_unit') at which the "
      "in-well fluid density is evaluated by 'unit_weight_fp'.  Only used if 'unit_weight_fp' "
      "is given.  Choose a value close to the expected wellbore pressure for the most accurate "
      "density.");
  params.addParam<Real>("re_constant",
                        0.28,
                        "The dimensionless constant used in evaluating the borehole effective "
                        "radius.  This depends on the meshing scheme.  Peacemann "
                        "finite-difference calculations give 0.28, while for rectangular finite "
                        "elements the result is closer to 0.1594.  (See  Eqn(4.13) of Z Chen, Y "
                        "Zhang, Well flow models for various numerical methods, Int J Num "
                        "Analysis and Modeling, 3 (2008) 375-388.)");
  params.addParam<Real>("well_constant",
                        -1.0,
                        "Usually this is calculated internally from the element geometry, the "
                        "local borehole direction and segment length, and the permeability.  "
                        "However, if this parameter is given as a positive number then this "
                        "number is used instead of the internal calculation.  This speeds up "
                        "computation marginally.  re_constant becomes irrelevant");
  params.addClassDescription(
      "Approximates a borehole in the mesh using the Peaceman approach, ie "
      "using a number of point sinks with given radii whose positions are "
      "read from a file.  NOTE: if you are using PorousFlowPorosity that depends on volumetric "
      "strain, you should set strain_at_nearest_qp=true in your GlobalParams, to ensure the nodal "
      "Porosity Material uses the volumetric strain at the Dirac quadpoints, and can therefore be "
      "computed.  The wellbore pressure profile is built either from a constant fluid unit "
      "weight ('unit_weight') or, if a thermal gradient along the borehole makes a single "
      "constant unit weight a poor approximation, from a fluid density computed at each "
      "borehole point from a temperature-dependent fluid-properties UserObject "
      "('unit_weight_fp')");
  return params;
}

PorousFlowPeacemanBorehole::PorousFlowPeacemanBorehole(const InputParameters & parameters)
  : PorousFlowLineSink(parameters),
    _control(getParam<MooseEnum>("control").getEnum<Control>()),
    _character(isParamValid("character") ? &getFunction("character") : nullptr),
    _p_bot(isParamValid("bottom_p_or_t") ? &getFunction("bottom_p_or_t") : nullptr),
    _unit_weight(isParamValid("unit_weight") ? getParam<RealVectorValue>("unit_weight")
                                             : RealVectorValue()),
    _use_density_from_temperature(isParamValid("unit_weight_fp")),
    _fp(_use_density_from_temperature ? &getUserObject<SinglePhaseFluidProperties>("unit_weight_fp")
                                      : nullptr),
    _temperature_var(_use_density_from_temperature && isCoupled("unit_weight_temperature")
                         ? getFieldVar("unit_weight_temperature", 0)
                         : nullptr),
    _temperature_system(_temperature_var ? &_temperature_var->sys().system() : nullptr),
    _temperature_var_number(_temperature_var ? _temperature_var->number() : libMesh::invalid_uint),
    _gravity(isParamValid("unit_weight_gravity") ? getParam<RealVectorValue>("unit_weight_gravity")
                                                 : RealVectorValue()),
    _density_reference_pressure(
        getParam<Real>("unit_weight_reference_pressure") *
        (getParam<MooseEnum>("unit_weight_pressure_unit") == 0 ? 1.0 : 1.0E6)),
    _t_c2k(getParam<MooseEnum>("unit_weight_temperature_unit") == 0 ? 0.0 : 273.15),
    _re_constant(getParam<Real>("re_constant")),
    _well_constant(getParam<Real>("well_constant")),
    _has_permeability(
        !_fv && hasMaterialProperty<RealTensorValue>("PorousFlow_permeability_qp") &&
        hasMaterialProperty<std::vector<RealTensorValue>>("dPorousFlow_permeability_qp_dvar")),
    _has_thermal_conductivity(
        !_fv && hasMaterialProperty<RealTensorValue>("PorousFlow_thermal_conductivity_qp") &&
        hasMaterialProperty<std::vector<RealTensorValue>>(
            "dPorousFlow_thermal_conductivity_qp_dvar")),
    _perm_or_cond(
        _fv ? nullptr
            : (_p_or_t == PorTchoice::pressure
                   ? &getMaterialProperty<RealTensorValue>("PorousFlow_permeability_qp")
                   : &getMaterialProperty<RealTensorValue>("PorousFlow_thermal_conductivity_qp"))),
    _dperm_or_cond_dvar(_fv ? nullptr
                            : (_p_or_t == PorTchoice::pressure
                                   ? &getMaterialProperty<std::vector<RealTensorValue>>(
                                         "dPorousFlow_permeability_qp_dvar")
                                   : &getMaterialProperty<std::vector<RealTensorValue>>(
                                         "dPorousFlow_thermal_conductivity_qp_dvar"))),
    _rate(isParamValid("rate") ? &getFunction("rate") : nullptr),
    _bhp_max(isParamValid("bhp_max") ? &getFunction("bhp_max") : nullptr),
    _max_swaps(getParam<unsigned int>("max_constraint_swaps")),
    _ad_pp(_fv && _p_or_t == PorTchoice::pressure
               ? &getADMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")
               : nullptr),
    _ad_temperature(_fv && _p_or_t == PorTchoice::temperature
                        ? &getADMaterialProperty<Real>("PorousFlow_temperature_qp")
                        : nullptr),
    _ad_perm_or_cond(
        _fv ? (_p_or_t == PorTchoice::pressure
                   ? &getADMaterialProperty<RealTensorValue>("PorousFlow_permeability_qp")
                   : &getADMaterialProperty<RealTensorValue>("PorousFlow_thermal_conductivity_qp"))
            : nullptr),
    _ad_density(_fv && _use_mobility
                    ? &getADMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")
                    : nullptr),
    _ad_viscosity(_fv && _use_mobility
                      ? &getADMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_qp")
                      : nullptr),
    _ad_relperm(
        _fv && (_use_mobility || _use_relative_permeability)
            ? &getADMaterialProperty<std::vector<Real>>("PorousFlow_relative_permeability_qp")
            : nullptr),
    _ad_mass_fractions(
        _fv && _use_mass_fraction
            ? &getADMaterialProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_qp")
            : nullptr),
    _ad_enthalpy(_fv && _use_enthalpy ? &getADMaterialProperty<std::vector<Real>>(
                                            "PorousFlow_fluid_phase_enthalpy_qp")
                                      : nullptr),
    _ad_internal_energy(
        _fv && _use_internal_energy
            ? &getADMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_internal_energy_qp")
            : nullptr),
    _ad_multiplying_var(_fv ? &adCoupledValue("multiplying_var") : nullptr),
    _ad_bottom_hole_pressure(_fv && _control == Control::rate && _use_scalar
                                 ? &adCoupledScalarValue("scalar_variable")
                                 : nullptr)
{
  if (!_fv && _p_or_t == PorTchoice::pressure && !_has_permeability)
    mooseError("PorousFlowPeacemanBorehole: You have specified function_of=porepressure, but you "
               "do not have a quadpoint permeability material");
  if (!_fv && _p_or_t == PorTchoice::temperature && !_has_thermal_conductivity)
    mooseError("PorousFlowPeacemanBorehole: You have specified function_of=temperature, but you do "
               "not have a quadpoint thermal_conductivity material");

  // The wellbore pressure profile is built from either a single constant unit_weight, or a
  // fluid density computed from temperature at each point (via 'unit_weight_fp') - never both,
  // and never neither.  Note this deliberately does not error if 'unit_weight_temperature' or
  // 'unit_weight_gravity' are valid while 'unit_weight_fp' is not: unlike 'unit_weight' and
  // 'unit_weight_fp' themselves, those two are not required to detect the user's intended mode,
  // so treating their unexpected presence as an error would only serve to reject input files
  // that harmlessly set an identically-named GlobalParam for some unrelated object.
  const int checkWellborePressureFormat =
      int(isParamValid("unit_weight")) + int(isParamValid("unit_weight_fp"));
  if (checkWellborePressureFormat > 1)
    paramError("unit_weight",
               "PorousFlowPeacemanBorehole: must specify only one of 'unit_weight' (a constant "
               "fluid unit weight) or 'unit_weight_fp' (a fluid-properties UserObject, so that "
               "the fluid unit weight is instead computed from the temperature at each borehole "
               "point)");
  else if (checkWellborePressureFormat == 0)
    paramError("unit_weight",
               "PorousFlowPeacemanBorehole: must specify at least one of 'unit_weight' or "
               "'unit_weight_fp'");

  if (_control == Control::bhp)
  {
    if (!_p_bot)
      paramError("bottom_p_or_t", "bottom_p_or_t is needed with control = bhp");
    if (!_character)
      paramError("character", "character is needed with control = bhp");
    if (_use_scalar || _rate || _bhp_max)
      paramError("control",
                 "bottom_hole_pressure, rate and bhp_max are only used with control = rate");
  }
  else
  {
    if (!_use_scalar)
      paramError(
          "bottom_hole_pressure",
          "with control = rate the borehole solves for its bottom-hole pressure, so it needs "
          "a scalar variable to hold it");
    if (_k_order != 1)
      paramError("bottom_hole_pressure",
                 "a borehole has one bottom-hole pressure, so this must be a scalar variable "
                 "carrying a single degree of freedom, but the one supplied carries ",
                 _k_order);
    if (_compute_scalar_residuals && !_rate)
      paramError("rate",
                 "this object contributes to the borehole's constraint, so it needs the rate.  "
                 "Objects that contribute only flux set compute_scalar_residuals = false");
    if (_p_bot)
      paramError("bottom_p_or_t",
                 "bottom_p_or_t is not used with control = rate: the bottom-hole pressure is the "
                 "bottom_hole_pressure variable");
    if (_character)
      paramError("character",
                 "character is not used with control = rate: fluid flows in either direction at "
                 "every completion, as the solved bottom-hole pressure dictates");
    if (_use_density_from_temperature)
      paramError("unit_weight_fp", "unit_weight_fp is not available with control = rate");
  }

  if (_fv && _use_density_from_temperature)
    paramError("unit_weight_fp", "unit_weight_fp is not available for finite-volume variables");

  if (_use_density_from_temperature)
  {
    if (_p_or_t == PorTchoice::temperature)
      paramError("unit_weight_fp",
                 "PorousFlowPeacemanBorehole: 'unit_weight_fp' computes a fluid density to "
                 "build a hydrostatic *pressure* profile along the borehole, which is "
                 "meaningless when function_of=temperature (bottom_p_or_t is then a "
                 "temperature, not a pressure)");
    if (!isParamValid("unit_weight_temperature"))
      paramError("unit_weight_temperature",
                 "PorousFlowPeacemanBorehole: 'unit_weight_temperature' must be supplied when "
                 "'unit_weight_fp' is supplied");
    if (!isCoupled("unit_weight_temperature"))
      paramError("unit_weight_temperature",
                 "PorousFlowPeacemanBorehole: 'unit_weight_temperature' must be a nonlinear or "
                 "auxiliary variable, not a constant value.  The wellbore pressure profile is "
                 "built by sampling this variable at each borehole point, so a spatially "
                 "constant temperature has no profile to sample: use 'unit_weight' instead if "
                 "the in-well fluid density really is constant");
    if (!isParamValid("unit_weight_gravity"))
      paramError("unit_weight_gravity",
                 "PorousFlowPeacemanBorehole: 'unit_weight_gravity' must be supplied when "
                 "'unit_weight_fp' is supplied");
  }
}

void
PorousFlowPeacemanBorehole::initialSetup()
{
  PorousFlowLineGeometry::initialSetup();

  if (!_point_file.empty() && _zs[0] < _zs.back())
    mooseError("PorousFlowPeacemanBorehole: The last entry in the point_file needs to be at the "
               "bottom of the well_bore because this is the point where the function bottom_p_or_t "
               "is evaluated.  The depth of the first point is z=",
               _zs[0],
               " and the last point is z=",
               _zs.back());

  // construct the rotation matrix needed to rotate the permeability
  const unsigned int num_pts = _zs.size();
  _rot_matrix.resize(std::max(num_pts - 1, (unsigned)1));
  for (unsigned int i = 0; i + 1 < num_pts; ++i)
  {
    const RealVectorValue v2(_xs[i + 1] - _xs[i], _ys[i + 1] - _ys[i], _zs[i + 1] - _zs[i]);
    _rot_matrix[i] = RotationMatrix::rotVecToZ(v2);
  }
  if (num_pts == (unsigned)1)
    _rot_matrix[0] = RotationMatrix::rotVecToZ(_line_direction);

  if (contributesToConstraint())
  {
    // The objects sharing the bottom-hole pressure that contribute to the constraint each add
    // their flux to it, and the first of them adds the rate and decides whether bhp_max holds
    const auto rate_name = getParam<FunctionName>("rate");
    const auto bhp_max_name = isParamValid("bhp_max") ? getParam<FunctionName>("bhp_max") : "";
    const auto & dirac_kernels = _fe_problem.getNonlinearSystemBase(_sys.number())
                                     .getDiracKernelWarehouse()
                                     .getActiveObjects();
    for (const auto & dk : dirac_kernels)
    {
      auto * bh = dynamic_cast<PorousFlowPeacemanBorehole *>(dk.get());
      if (!bh || !bh->contributesToConstraint() ||
          bh->bottomHolePressureVariable() != bottomHolePressureVariable())
        continue;
      if (bh->getParam<FunctionName>("rate") != rate_name)
        paramError("rate",
                   "the objects ",
                   name(),
                   " and ",
                   bh->name(),
                   " share bottom_hole_pressure and both contribute to its constraint, so they "
                   "are one borehole and must use the same rate");
      if ((bh->isParamValid("bhp_max") ? bh->getParam<FunctionName>("bhp_max") : "") !=
          bhp_max_name)
        paramError("bhp_max",
                   "the objects ",
                   name(),
                   " and ",
                   bh->name(),
                   " share bottom_hole_pressure and both contribute to its constraint, so they "
                   "must use the same bhp_max");
      if (!_owner)
        _owner = bh;
    }
    _adds_rate = (_owner == this);
  }
}

void
PorousFlowPeacemanBorehole::timestepSetup()
{
  PorousFlowLineSink::timestepSetup();
  if (!_adds_rate || !_bhp_max)
    return;

  // The constraint is decided here, once per timestep, so the equation being solved never changes
  // part way through a solve.  Deciding it between Newton iterations instead was measured to cost a
  // failed timestep wherever the limit engages, while never engaging the limit within the step
  // that first breaches it
  _n_swaps = 0;
  _communicator.sum(_rate_mode_residual);
  updateConstraint();
}

void
PorousFlowPeacemanBorehole::updateConstraint()
{
  // Once the borehole has switched often enough in one timestep, stop: an operating point sitting
  // on the limit can otherwise chatter between the two constraints and never converge
  if (_n_swaps >= _max_swaps)
    return;

  const Real bhp_max = _bhp_max->value(_t, _bottom_point);
  if (!_clamped && _pbh_last > bhp_max)
  {
    // Delivering the rate needs a pressure the borehole is not allowed to reach
    _clamped = true;
    ++_n_swaps;
  }
  else if (_clamped && _rate_mode_residual < 0.0)
  {
    // The flux falls as the bottom-hole pressure rises, so a negative rate-constraint residual at
    // the limit means the rate would be met below it, and the limit no longer binds
    _clamped = false;
    ++_n_swaps;
  }
}

void
PorousFlowPeacemanBorehole::recordRateMode(const Real rate_mode)
{
  if (!_owner || !_owner->_bhp_max)
    return;
  _owner->_rate_mode_residual += rate_mode;
  _owner->_pbh_last = _kappa[0];
}

void
PorousFlowPeacemanBorehole::residualSetup()
{
  PorousFlowLineSink::residualSetup();
  computeWellborePressures();
  // The value left by the last residual evaluation of a timestep decides the next
  if (_adds_rate)
    _rate_mode_residual = 0.0;
}

void
PorousFlowPeacemanBorehole::jacobianSetup()
{
  PorousFlowLineSink::jacobianSetup();
  computeWellborePressures();
}

void
PorousFlowPeacemanBorehole::computeWellborePressures()
{
  if (!_use_density_from_temperature)
    return;

  const std::size_t num_pts = _z_coord->size();
  _bh_pressure.assign(num_pts, 0.0);
  if (num_pts == 0)
    return;

  // Sample the temperature, and hence the in-well fluid density, at every well point (not just
  // the points owned by this processor), because the wellbore pressure at any point is a
  // cumulative integral over all the points between it and the bottom point.
  // System::point_value performs the parallel point-location and broadcast internally, so this
  // is safe to call even for points this processor's mesh partition does not own.
  std::vector<Real> density(num_pts);
  for (const auto i : make_range(num_pts))
  {
    const Point p(_x_coord->at(i), _y_coord->at(i), _z_coord->at(i));
    const Real temperature = _temperature_system->point_value(_temperature_var_number, p, false);
    density[i] = _fp->rho_from_p_T(_density_reference_pressure, temperature + _t_c2k);
  }

  // Integrate the fluid unit weight (density*gravity) up the wellbore from the bottom point,
  // where the pressure is prescribed by bottom_p_or_t, using the trapezoidal rule on each
  // polyline segment.  This is exact for a density varying linearly along the well (the target
  // use case: a linear thermal gradient with a locally-linear density(temperature)), and it
  // degenerates exactly to the constant-unit_weight formula when the density is constant, since
  // the sum then telescopes to density*gravity.(x_i - x_bottom).
  _bh_pressure[num_pts - 1] = _p_bot->value(_t, _bottom_point);
  for (std::size_t i = num_pts - 1; i > 0; --i)
  {
    const RealVectorValue segment(_x_coord->at(i - 1) - _x_coord->at(i),
                                  _y_coord->at(i - 1) - _y_coord->at(i),
                                  _z_coord->at(i - 1) - _z_coord->at(i));
    _bh_pressure[i - 1] =
        _bh_pressure[i] + 0.5 * (density[i - 1] + density[i]) * (_gravity * segment);
  }
}

Real
PorousFlowPeacemanBorehole::wellborePressure(unsigned current_dirac_ptid) const
{
  if (!_use_density_from_temperature)
    return _p_bot->value(_t, _bottom_point) + _unit_weight * (_q_point[_qp] - _bottom_point);

  mooseAssert(current_dirac_ptid < _bh_pressure.size(),
              "PorousFlowPeacemanBorehole: the wellbore pressure profile has not been computed "
              "for this Dirac point");
  return _bh_pressure[current_dirac_ptid];
}

Real
PorousFlowPeacemanBorehole::wellConstant(const RealTensorValue & perm,
                                         const RealTensorValue & rot,
                                         const Real & half_len,
                                         const Elem * ele,
                                         const Real & rad) const
{
  if (_well_constant > 0)
    return _well_constant;
  return peacemanWellConstant(perm, rot, half_len, ele, rad, _re_constant);
}

Real
PorousFlowPeacemanBorehole::peacemanWellConstant(const RealTensorValue & perm,
                                                 const RealTensorValue & rot,
                                                 Real half_len,
                                                 const Elem * ele,
                                                 Real rad,
                                                 Real re_constant)
// Peaceman's form for the borehole well constant
{
  // rot_perm has its "2" component lying along the half segment.
  // We want to determine the eigenvectors of rot(0:1, 0:1), since, when
  // rotated back to the original frame we will determine the element
  // lengths along these directions
  const RealTensorValue rot_perm = (rot * perm) * rot.transpose();
  const Real trace2D = rot_perm(0, 0) + rot_perm(1, 1);
  const Real det2D = rot_perm(0, 0) * rot_perm(1, 1) - rot_perm(0, 1) * rot_perm(1, 0);
  const Real sq = std::sqrt(std::max(0.25 * trace2D * trace2D - det2D,
                                     0.0)); // the std::max accounts for wierdo precision loss
  const Real eig_val1 = 0.5 * trace2D + sq;
  const Real eig_val2 = 0.5 * trace2D - sq;
  RealVectorValue eig_vec1, eig_vec2;
  if (sq > std::abs(trace2D) * 1E-7) // matrix is not a multiple of the identity (1E-7 accounts for
                                     // precision in a crude way)
  {
    if (rot_perm(1, 0) != 0)
    {
      eig_vec1(0) = eig_val1 - rot_perm(1, 1);
      eig_vec1(1) = rot_perm(1, 0);
      eig_vec2(0) = eig_val2 - rot_perm(1, 1);
      eig_vec2(1) = rot_perm(1, 0);
    }
    else if (rot_perm(0, 1) != 0)
    {
      eig_vec1(0) = rot_perm(0, 1);
      eig_vec1(1) = eig_val1 - rot_perm(0, 0);
      eig_vec2(0) = rot_perm(0, 1);
      eig_vec2(1) = eig_val2 - rot_perm(0, 0);
    }
    else // off diagonal terms are both zero
    {
      eig_vec1(0) = 1.0;
      eig_vec2(1) = 1.0;
    }
  }
  else // matrix is basically a multiple of the identity
  {
    eig_vec1(0) = 1.0;
    eig_vec2(1) = 1.0;
  }

  // finally, rotate these to original frame and normalise
  eig_vec1 = rot.transpose() * eig_vec1;
  eig_vec1 /= std::sqrt(eig_vec1 * eig_vec1);
  eig_vec2 = rot.transpose() * eig_vec2;
  eig_vec2 /= std::sqrt(eig_vec2 * eig_vec2);

  // find the "length" of the element in these directions
  // TODO - maybe better to use variance than max&min
  Real max1 = eig_vec1 * ele->point(0);
  Real max2 = eig_vec2 * ele->point(0);
  Real min1 = max1;
  Real min2 = max2;
  Real proj;
  for (unsigned int i = 1; i < ele->n_nodes(); i++)
  {
    proj = eig_vec1 * ele->point(i);
    max1 = (max1 < proj) ? proj : max1;
    min1 = (min1 < proj) ? min1 : proj;

    proj = eig_vec2 * ele->point(i);
    max2 = (max2 < proj) ? proj : max2;
    min2 = (min2 < proj) ? min2 : proj;
  }
  const Real ll1 = max1 - min1;
  const Real ll2 = max2 - min2;

  Real r0;
  if (eig_val1 <= 0.0)
    r0 = re_constant * ll1;
  else if (eig_val2 <= 0.0)
    r0 = re_constant * ll2;
  else
    r0 = re_constant *
         std::sqrt(std::sqrt(eig_val1 / eig_val2) * std::pow(ll2, 2) +
                   std::sqrt(eig_val2 / eig_val1) * std::pow(ll1, 2)) /
         (std::pow(eig_val1 / eig_val2, 0.25) + std::pow(eig_val2 / eig_val1, 0.25));

  const Real effective_perm = (det2D >= 0.0 ? std::sqrt(det2D) : 0.0);

  const Real halfPi = acos(0.0);

  if (r0 <= rad)
    ::mooseError(
        "The effective element size (about 0.2-times-true-ele-size) for an element "
        "containing a Peaceman-type borehole must be (much) larger than the borehole radius "
        "for the Peaceman formulation to be correct.  Your element has effective size ",
        r0,
        " and the borehole radius is ",
        rad,
        "\n");

  return 4 * halfPi * effective_perm * half_len / std::log(r0 / rad);
}

Real
PorousFlowPeacemanBorehole::wellConstantSum(unsigned current_dirac_ptid) const
{
  Real wc_sum = 0.0;

  if (current_dirac_ptid > 0)
    // contribution from half-segment "behind" this point (must have >1 point for
    // current_dirac_ptid>0)
    wc_sum += wellConstant(permOrCond(),
                           _rot_matrix[current_dirac_ptid - 1],
                           _half_seg_len[current_dirac_ptid - 1],
                           _current_elem,
                           _weight->at(current_dirac_ptid));

  if (current_dirac_ptid + 1 < _zs.size() || _zs.size() == 1)
    // contribution from half-segment "ahead of" this point, or we only have one point
    wc_sum += wellConstant(permOrCond(),
                           _rot_matrix[current_dirac_ptid],
                           _half_seg_len[current_dirac_ptid],
                           _current_elem,
                           _weight->at(current_dirac_ptid));

  return wc_sum;
}

RealTensorValue
PorousFlowPeacemanBorehole::permOrCond() const
{
  return _fv ? MetaPhysicL::raw_value((*_ad_perm_or_cond)[_qp]) : (*_perm_or_cond)[_qp];
}

Real
PorousFlowPeacemanBorehole::computeQpBaseOutflow(unsigned current_dirac_ptid) const
{
  if (_control == Control::rate)
    // No directional gate: fluid flows in either direction at every completion, and the flux is
    // differentiable in the bottom-hole pressure at the operating point the constraint settles on
    return conductance(current_dirac_ptid) * (ptqp() - boreholePressure());

  const Real character = _character->value(_t, _q_point[_qp]);
  if (character == 0.0)
    return 0.0;

  const Real bh_pressure = wellborePressure(current_dirac_ptid);
  const Real pp = ptqp();

  // this is the flow rate from porespace out of the system.  The borehole only acts in the
  // direction given by the sign of character: a production well does nothing where the borehole
  // pressure exceeds the porepressure, and an injection well does nothing where it does not
  if (!((character < 0.0 && pp < bh_pressure) || (character > 0.0 && pp > bh_pressure)))
    return 0.0;

  const Real outflow = wellConstantSum(current_dirac_ptid) * (pp - bh_pressure);

  return outflow * _test[_i][_qp] * std::abs(character);
}

void
PorousFlowPeacemanBorehole::computeQpBaseOutflowJacobian(unsigned jvar,
                                                         unsigned current_dirac_ptid,
                                                         Real & outflow,
                                                         Real & outflowp) const
{
  outflow = 0.0;
  outflowp = 0.0;

  if (_control == Control::rate)
  {
    if (_dictator.notPorousFlowVariable(jvar))
      return;
    const unsigned pvar = _dictator.porousFlowVariableNum(jvar);
    const Real c = conductance(current_dirac_ptid);
    outflow = c * (ptqp() - boreholePressure());
    outflowp = c * dptqp(pvar) * _phi[_j][_qp];
    return;
  }

  const Real character = _character->value(_t, _q_point[_qp]);
  if (character == 0.0)
    return;

  if (_dictator.notPorousFlowVariable(jvar))
    return;
  const unsigned pvar = _dictator.porousFlowVariableNum(jvar);

  // When _use_density_from_temperature is true, bh_pressure also depends on the temperature at
  // every well point between here and the bottom point (via computeWellborePressures()), not
  // just on the porous flow variables at this quadpoint.  That dependence is deliberately not
  // differentiated here: a DiracKernel can only assemble into the (test, phi) dofs of the
  // element containing its own quadpoint, so the coupling to temperature dofs in other elements
  // along the well cannot be represented in this Jacobian.  residualSetup()/jacobianSetup()
  // still recompute bh_pressure from the current nonlinear iterate before every evaluation, so
  // the converged solution is unaffected; only the Newton convergence rate may be mildly slower.
  // A proper fix (a relationship-managed line/segment kernel base class) is tracked separately
  // in https://github.com/idaholab/moose/issues/33757, to be picked up after this PR merges.
  const Real bh_pressure = wellborePressure(current_dirac_ptid);
  const Real pp = ptqp();
  const Real pp_prime = dptqp(pvar) * _phi[_j][_qp];

  if (current_dirac_ptid > 0)
  // contribution from half-segment "behind" this point
  {
    if ((character < 0.0 && pp < bh_pressure) || (character > 0.0 && pp > bh_pressure))
    {
      // injection, so outflow<0 || // production, so outflow>0
      const Real wc = wellConstant(permOrCond(),
                                   _rot_matrix[current_dirac_ptid - 1],
                                   _half_seg_len[current_dirac_ptid - 1],
                                   _current_elem,
                                   _weight->at(current_dirac_ptid));
      outflowp += wc * pp_prime;
      outflow += wc * (pp - bh_pressure);
    }
  }

  if (current_dirac_ptid < _zs.size() - 1 || _zs.size() == 1)
  // contribution from half-segment "ahead of" this point
  {
    if ((character < 0.0 && pp < bh_pressure) || (character > 0.0 && pp > bh_pressure))
    {
      // injection, so outflow<0 || // production, so outflow>0
      const Real wc = wellConstant(permOrCond(),
                                   _rot_matrix[current_dirac_ptid],
                                   _half_seg_len[current_dirac_ptid],
                                   _current_elem,
                                   _weight->at(current_dirac_ptid));
      outflowp += wc * pp_prime;
      outflow += wc * (pp - bh_pressure);
    }
  }

  outflowp *= _test[_i][_qp] * std::abs(character);
  outflow *= _test[_i][_qp] * std::abs(character);
}

// ------------------------------------------------------------------------------------------------
// Rate control, finite elements

Real
PorousFlowPeacemanBorehole::conductance(unsigned current_dirac_ptid) const
{
  return wellConstantSum(current_dirac_ptid) * _test[_i][_qp];
}

Real
PorousFlowPeacemanBorehole::boreholePressure() const
{
  return _kappa[0] + _unit_weight * (_q_point[_qp] - _bottom_point);
}

Real
PorousFlowPeacemanBorehole::computeQpScalarResidual()
{
  // The flux of the whole completion is the sum over the element's nodes, since the residual of
  // each node carries only its own share via the test function
  const unsigned current_dirac_ptid = currentPointCachedID();
  const Real fraction = completionFraction(current_dirac_ptid);

  const auto i_save = _i;
  Real flux = 0.0;
  for (_i = 0; _i < _test.size(); ++_i)
    flux += weightedOutflow();
  _i = i_save;

  // The rate is a property of the borehole, not of any one completion or object, so it is added
  // by the owner alone, shared between the completions so that it enters exactly once however
  // they are spread over elements and processors
  const Real rate_mode = flux - (_adds_rate ? _rate->value(_t, _q_point[_qp]) * fraction : 0.0);
  recordRateMode(rate_mode);

  if (_owner && _owner->clamped())
    // The limit has taken over: the bottom-hole pressure is pinned to bhp_max
    return _adds_rate ? (_kappa[0] - _bhp_max->value(_t, _bottom_point)) * fraction : 0.0;

  return rate_mode;
}

Real
PorousFlowPeacemanBorehole::computeQpJacobianScalar()
{
  // outflow = weighting * conductance * (porepressure - boreholePressure()), and neither the
  // weighting nor the conductance depends on the borehole pressure
  const unsigned current_dirac_ptid = currentPointCachedID();
  const Real base_outflow = computeQpBaseOutflow(current_dirac_ptid);
  return -weightingFactor(base_outflow) * conductance(current_dirac_ptid);
}

Real
PorousFlowPeacemanBorehole::computeQpScalarJacobian()
{
  if (_owner && _owner->clamped())
    return _adds_rate ? completionFraction(currentPointCachedID()) : 0.0;

  const auto i_save = _i;
  Real jacobian = 0.0;
  for (_i = 0; _i < _test.size(); ++_i)
    jacobian += computeQpJacobianScalar();
  _i = i_save;
  return jacobian;
}

Real
PorousFlowPeacemanBorehole::computeQpScalarOffDiagJacobian(const unsigned int jvar_num)
{
  if (_owner && _owner->clamped())
    // A pinned bottom-hole pressure does not depend on the state of the reservoir
    return 0.0;

  const auto i_save = _i;
  Real jacobian = 0.0;
  for (_i = 0; _i < _test.size(); ++_i)
    jacobian += jac(jvar_num);
  _i = i_save;
  return jacobian;
}

// ------------------------------------------------------------------------------------------------
// Finite volumes: AD cell materials and AD assembly

ADReal
PorousFlowPeacemanBorehole::fvPorT() const
{
  return _p_or_t == PorTchoice::pressure ? (*_ad_pp)[_qp][_ph] : (*_ad_temperature)[_qp];
}

ADReal
PorousFlowPeacemanBorehole::fvPointOutflow(unsigned int i)
{
  Real character = 1.0;
  ADReal borehole_p = 0.0;
  if (_control == Control::rate)
    borehole_p = (*_ad_bottom_hole_pressure)[0] + _unit_weight * (_q_point[_qp] - _bottom_point);
  else
  {
    character = _character->value(_t, _q_point[_qp]);
    if (character == 0.0)
      return 0.0;
    borehole_p = wellborePressure(i);
  }

  const ADReal p_or_t = fvPorT();
  const Real drawdown = MetaPhysicL::raw_value(p_or_t - borehole_p);
  // With control = bhp the borehole only acts in the direction given by the sign of character, as
  // for finite elements
  if (_control == Control::bhp &&
      !((character < 0.0 && drawdown < 0.0) || (character > 0.0 && drawdown > 0.0)))
    return 0.0;
  const bool injecting = drawdown < 0.0;

  ADReal multiplier = std::abs(character) * wellConstantSum(i) * (*_ad_multiplying_var)[_qp];
  if (_use_relative_permeability)
    multiplier *= (*_ad_relperm)[_qp][_ph];
  if (_use_mobility)
  {
    if (_mobility_type == MobilityTypeChoice::phase)
      multiplier *=
          (*_ad_relperm)[_qp][_ph] * (*_ad_density)[_qp][_ph] / (*_ad_viscosity)[_qp][_ph];
    else
    {
      // the density of the nominated (injected) phase times the total mobility of the cell
      ADReal total_mobility = 0.0;
      for (const auto p : make_range(_dictator.numPhases()))
        total_mobility += (*_ad_relperm)[_qp][p] / (*_ad_viscosity)[_qp][p];
      multiplier *= (*_ad_density)[_qp][_ph] * total_mobility;
    }
  }
  if (_use_mass_fraction)
  {
    // Injected fluid carries the composition of the injected stream, outflow that of the cell
    if (injecting && _use_injected_mass_fraction)
      multiplier *= _injected_mass_fraction->value(_t, _q_point[_qp]);
    else
      multiplier *= (*_ad_mass_fractions)[_qp][_ph][_sp];
  }
  if (_use_enthalpy)
    multiplier *= (*_ad_enthalpy)[_qp][_ph];
  if (_use_internal_energy)
    multiplier *= (*_ad_internal_energy)[_qp][_ph];

  return multiplier * (p_or_t - borehole_p);
}

void
PorousFlowPeacemanBorehole::fvCellResiduals(ADReal & outflow,
                                            ADReal & scalar_residual,
                                            Real & rate_mode,
                                            const bool record_point_fluxes)
{
  outflow = 0.0;
  scalar_residual = 0.0;
  rate_mode = 0.0;
  const bool contributes = contributesToConstraint();
  const bool clamped = contributes && _owner && _owner->clamped();

  forEachPoint(
      [&](const Real point_value)
      {
        const unsigned int i = currentPointCachedID();
        const ADReal q = point_value * fvPointOutflow(i);
        outflow += q;
        if (record_point_fluxes && _point_fluxes)
          _point_fluxes->add(i, MetaPhysicL::raw_value(q));
        if (!contributes)
          return;
        const Real fraction = completionFraction(i);
        const Real rate_part = _adds_rate ? _rate->value(_t, _q_point[_qp]) * fraction : 0.0;
        rate_mode += MetaPhysicL::raw_value(q) - rate_part;
        if (clamped)
        {
          if (_adds_rate)
            scalar_residual +=
                ((*_ad_bottom_hole_pressure)[0] - _bhp_max->value(_t, _bottom_point)) * fraction;
        }
        else
          scalar_residual += q - rate_part;
      });
}

void
PorousFlowPeacemanBorehole::fvComputeResidual()
{
  ADReal outflow, scalar_residual;
  Real rate_mode;
  fvCellResiduals(outflow, scalar_residual, rate_mode, true);

  // A finite-volume variable has one dof per cell, whose test function is 1 in the cell
  if (_compute_field_residuals)
    addResiduals(_assembly,
                 std::array<Real, 1>{{MetaPhysicL::raw_value(outflow)}},
                 _var.dofIndices(),
                 _var.scalingFactor());
  _total_outflow_mass.add(MetaPhysicL::raw_value(outflow) * _dt);

  if (contributesToConstraint())
  {
    addResiduals(_assembly,
                 std::array<Real, 1>{{MetaPhysicL::raw_value(scalar_residual)}},
                 _kappa_var_ptr->dofIndices(),
                 _kappa_var_ptr->scalingFactor());
    recordRateMode(rate_mode);
  }

  // The Dirac residual loop does not add the residuals cached with explicit dof indices, so add
  // them here.  This does not double count: Assembly clears the cache as it adds
  _fe_problem.addCachedResidual(_tid);
}

void
PorousFlowPeacemanBorehole::fvComputeJacobian()
{
  // AD gives the derivatives with respect to all the variables at once (PorousFlow variables and
  // the bottom-hole pressure)
  ADReal outflow, scalar_residual;
  Real rate_mode;
  fvCellResiduals(outflow, scalar_residual, rate_mode, false);

  if (_compute_field_residuals)
    addJacobian(
        _assembly, std::array<ADReal, 1>{{outflow}}, _var.dofIndices(), _var.scalingFactor());
  if (contributesToConstraint())
    addJacobian(_assembly,
                std::array<ADReal, 1>{{scalar_residual}},
                _kappa_var_ptr->dofIndices(),
                _kappa_var_ptr->scalingFactor());
}

void
PorousFlowPeacemanBorehole::fvComputeResidualAndJacobian()
{
  ADReal outflow, scalar_residual;
  Real rate_mode;
  fvCellResiduals(outflow, scalar_residual, rate_mode, true);

  if (_compute_field_residuals)
    addResidualsAndJacobian(
        _assembly, std::array<ADReal, 1>{{outflow}}, _var.dofIndices(), _var.scalingFactor());
  _total_outflow_mass.add(MetaPhysicL::raw_value(outflow) * _dt);

  if (contributesToConstraint())
  {
    addResidualsAndJacobian(_assembly,
                            std::array<ADReal, 1>{{scalar_residual}},
                            _kappa_var_ptr->dofIndices(),
                            _kappa_var_ptr->scalingFactor());
    recordRateMode(rate_mode);
  }
  _fe_problem.addCachedResidual(_tid);
}

// ------------------------------------------------------------------------------------------------
// Dispatch between finite elements and finite volumes

void
PorousFlowPeacemanBorehole::computeResidual()
{
  if (_fv)
    fvComputeResidual();
  else
    PorousFlowLineSink::computeResidual();
}

void
PorousFlowPeacemanBorehole::computeJacobian()
{
  if (_fv)
    fvComputeJacobian();
  else
    PorousFlowLineSink::computeJacobian();
}

void
PorousFlowPeacemanBorehole::computeOffDiagJacobian(const unsigned int jvar_num)
{
  if (!_fv)
    PorousFlowLineSink::computeOffDiagJacobian(jvar_num);
  else if (jvar_num == _var.number())
    // ComputeDiracThread calls this once per coupled variable.  The AD Jacobian holds every
    // column, so it is assembled once, for the diagonal block
    fvComputeJacobian();
}

void
PorousFlowPeacemanBorehole::computeResidualAndJacobian()
{
  if (_fv)
    fvComputeResidualAndJacobian();
  else
    PorousFlowLineSink::computeResidualAndJacobian();
}
