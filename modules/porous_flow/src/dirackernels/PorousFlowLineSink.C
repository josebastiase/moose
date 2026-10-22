//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowLineSink.h"
#include "Function.h"
#include "libmesh/utility.h"

InputParameters
PorousFlowLineSink::validParams()
{
  InputParameters params = PorousFlowLineGeometry::validParams();
  MooseEnum p_or_t_choice("pressure=0 temperature=1", "pressure");
  params.addParam<MooseEnum>("function_of",
                             p_or_t_choice,
                             "Modifying functions will be a function of either pressure and "
                             "permeability (eg, for boreholes that pump fluids) or "
                             "temperature and thermal conductivity (eg, for boreholes that "
                             "pump pure heat with no fluid flow)");
  params.addRequiredParam<UserObjectName>(
      "SumQuantityUO",
      "User Object of type=PorousFlowSumQuantity in which to place the total "
      "outflow from the line sink for each time step.");
  params.addParam<UserObjectName>(
      "PointFluxUO",
      "Optional UserObject of type=PorousFlowPointFluxQuantity in which to record the "
      "instantaneous flux (eg kg.s^-1 for fluid, J.s^-1 for heat) at each individual Dirac "
      "point of this line sink, as computed during the most recent residual evaluation.  Use a "
      "PorousFlowPlotPointFluxQuantity VectorPostprocessor to output the recorded values.  "
      "Unlike SumQuantityUO, this is not multiplied by the timestep size.  Use a separate "
      "UserObject for each line sink.");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addParam<unsigned int>(
      "fluid_phase",
      0,
      "The fluid phase whose pressure (and potentially mobility, enthalpy, etc) "
      "controls the flux to the line sink.  For p_or_t=temperature, and without "
      "any use_*, this parameter is irrelevant");
  params.addParam<unsigned int>("mass_fraction_component",
                                "The index corresponding to a fluid "
                                "component.  If supplied, the flux will "
                                "be multiplied by the nodal mass "
                                "fraction for the component");
  params.addParam<bool>(
      "use_relative_permeability", false, "Multiply the flux by the fluid relative permeability");
  params.addParam<bool>("use_mobility", false, "Multiply the flux by the fluid mobility");
  MooseEnum mobility_type_choice("phase=0 total=1", "phase");
  params.addParam<MooseEnum>(
      "mobility_type",
      mobility_type_choice,
      "The mobility used when use_mobility=true.  'phase' is the mobility of 'fluid_phase' alone, "
      "k_r*rho/mu.  'total' is the density of 'fluid_phase' multiplied by the total mobility of "
      "all phases, rho*sum_p(k_r_p/mu_p).  'phase' is correct for a production well, whose stream "
      "carries what is mobile at the node.  An injection well needs 'total': the phase mobility of "
      "a phase that is absent from the block is zero, so with 'phase' such a well has no "
      "injectivity at all and can never begin to flow.");
  params.addParam<FunctionName>(
      "injected_mass_fraction",
      "The mass fraction of 'mass_fraction_component' in the fluid this object injects.  Fluid "
      "leaving the porespace carries the composition of the porespace, but fluid entering it "
      "carries the composition of whatever is being injected, which the reservoir has no knowledge "
      "of.  Without this parameter inflow is weighted by the resident mass fraction, which is "
      "correct for a producer but not for an injector: injecting a gas into brine that holds "
      "almost none of it, the resident weighting injects almost nothing.  May be a function of "
      "time, for an injectate whose composition changes.  Only meaningful with "
      "'mass_fraction_component'.");
  params.addParam<bool>("use_enthalpy", false, "Multiply the flux by the fluid enthalpy");
  params.addParam<bool>(
      "use_internal_energy", false, "Multiply the flux by the fluid internal energy");
  params.addCoupledVar("multiplying_var", 1.0, "Fluxes will be moultiplied by this variable");
  params.addClassDescription("Approximates a line sink in the mesh by a sequence of weighted Dirac "
                             "points whose positions are read from a file");
  return params;
}

PorousFlowLineSink::PorousFlowLineSink(const InputParameters & parameters)
  : PorousFlowLineGeometry(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _total_outflow_mass(
        const_cast<PorousFlowSumQuantity &>(getUserObject<PorousFlowSumQuantity>("SumQuantityUO"))),
    _point_fluxes(isParamValid("PointFluxUO")
                      ? &const_cast<PorousFlowPointFluxQuantity &>(
                            getUserObject<PorousFlowPointFluxQuantity>("PointFluxUO"))
                      : nullptr),
    _fv(_var.isFV()),

    _has_porepressure(
        !_fv && hasMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp") &&
        hasMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_porepressure_qp_dvar")),
    _has_temperature(!_fv && hasMaterialProperty<Real>("PorousFlow_temperature_qp") &&
                     hasMaterialProperty<std::vector<Real>>("dPorousFlow_temperature_qp_dvar")),
    _has_mass_fraction(
        !_fv && hasMaterialProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_nodal") &&
        hasMaterialProperty<std::vector<std::vector<std::vector<Real>>>>(
            "dPorousFlow_mass_frac_nodal_dvar")),
    _has_relative_permeability(
        !_fv && hasMaterialProperty<std::vector<Real>>("PorousFlow_relative_permeability_nodal") &&
        hasMaterialProperty<std::vector<std::vector<Real>>>(
            "dPorousFlow_relative_permeability_nodal_dvar")),
    _has_mobility(
        !_fv && hasMaterialProperty<std::vector<Real>>("PorousFlow_relative_permeability_nodal") &&
        hasMaterialProperty<std::vector<std::vector<Real>>>(
            "dPorousFlow_relative_permeability_nodal_dvar") &&
        hasMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal") &&
        hasMaterialProperty<std::vector<std::vector<Real>>>(
            "dPorousFlow_fluid_phase_density_nodal_dvar") &&
        hasMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_nodal") &&
        hasMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_viscosity_nodal_dvar")),
    _has_enthalpy(!_fv &&
                  hasMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_enthalpy_nodal") &&
                  hasMaterialProperty<std::vector<std::vector<Real>>>(
                      "dPorousFlow_fluid_phase_enthalpy_nodal_dvar")),
    _has_internal_energy(
        !_fv &&
        hasMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_internal_energy_nodal") &&
        hasMaterialProperty<std::vector<std::vector<Real>>>(
            "dPorousFlow_fluid_phase_internal_energy_nodal_dvar")),

    _p_or_t(getParam<MooseEnum>("function_of").getEnum<PorTchoice>()),
    _use_mass_fraction(isParamValid("mass_fraction_component")),
    _use_relative_permeability(getParam<bool>("use_relative_permeability")),
    _use_mobility(getParam<bool>("use_mobility")),
    _use_injected_mass_fraction(isParamValid("injected_mass_fraction")),
    _injected_mass_fraction(_use_injected_mass_fraction ? &getFunction("injected_mass_fraction")
                                                        : nullptr),
    _mobility_type(getParam<MooseEnum>("mobility_type").getEnum<MobilityTypeChoice>()),
    _use_enthalpy(getParam<bool>("use_enthalpy")),
    _use_internal_energy(getParam<bool>("use_internal_energy")),

    _ph(getParam<unsigned int>("fluid_phase")),
    _sp(_use_mass_fraction ? getParam<unsigned int>("mass_fraction_component") : 0),

    _pp((_p_or_t == PorTchoice::pressure && _has_porepressure)
            ? &getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")
            : nullptr),
    _dpp_dvar((_p_or_t == PorTchoice::pressure && _has_porepressure)
                  ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                        "dPorousFlow_porepressure_qp_dvar")
                  : nullptr),
    _temperature((_p_or_t == PorTchoice::temperature && _has_temperature)
                     ? &getMaterialProperty<Real>("PorousFlow_temperature_qp")
                     : nullptr),
    _dtemperature_dvar(
        (_p_or_t == PorTchoice::temperature && _has_temperature)
            ? &getMaterialProperty<std::vector<Real>>("dPorousFlow_temperature_qp_dvar")
            : nullptr),
    _fluid_density_node(
        (_use_mobility && _has_mobility)
            ? &getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")
            : nullptr),
    _dfluid_density_node_dvar((_use_mobility && _has_mobility)
                                  ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                        "dPorousFlow_fluid_phase_density_nodal_dvar")
                                  : nullptr),
    _fluid_viscosity((_use_mobility && _has_mobility)
                         ? &getMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_nodal")
                         : nullptr),
    _dfluid_viscosity_dvar((_use_mobility && _has_mobility)
                               ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                     "dPorousFlow_viscosity_nodal_dvar")
                               : nullptr),
    _relative_permeability(
        ((_use_mobility && _has_mobility) ||
         (_use_relative_permeability && _has_relative_permeability))
            ? &getMaterialProperty<std::vector<Real>>("PorousFlow_relative_permeability_nodal")
            : nullptr),
    _drelative_permeability_dvar(((_use_mobility && _has_mobility) ||
                                  (_use_relative_permeability && _has_relative_permeability))
                                     ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                           "dPorousFlow_relative_permeability_nodal_dvar")
                                     : nullptr),
    _mass_fractions(
        (_use_mass_fraction && _has_mass_fraction)
            ? &getMaterialProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_nodal")
            : nullptr),
    _dmass_fractions_dvar((_use_mass_fraction && _has_mass_fraction)
                              ? &getMaterialProperty<std::vector<std::vector<std::vector<Real>>>>(
                                    "dPorousFlow_mass_frac_nodal_dvar")
                              : nullptr),
    _enthalpy(_has_enthalpy ? &getMaterialPropertyByName<std::vector<Real>>(
                                  "PorousFlow_fluid_phase_enthalpy_nodal")
                            : nullptr),
    _denthalpy_dvar(_has_enthalpy ? &getMaterialPropertyByName<std::vector<std::vector<Real>>>(
                                        "dPorousFlow_fluid_phase_enthalpy_nodal_dvar")
                                  : nullptr),
    _internal_energy(_has_internal_energy ? &getMaterialPropertyByName<std::vector<Real>>(
                                                "PorousFlow_fluid_phase_internal_energy_nodal")
                                          : nullptr),
    _dinternal_energy_dvar(_has_internal_energy
                               ? &getMaterialPropertyByName<std::vector<std::vector<Real>>>(
                                     "dPorousFlow_fluid_phase_internal_energy_nodal_dvar")
                               : nullptr),
    _multiplying_var(coupledValue("multiplying_var"))
{
  // zero the outflow mass
  _total_outflow_mass.zero();

  if (_ph >= _dictator.numPhases())
    paramError("fluid_phase",
               "The Dictator proclaims that the maximum phase index in this simulation is ",
               _dictator.numPhases() - 1,
               " whereas you have used ",
               _ph,
               ". Remember that indexing starts at 0. You must try harder.");

  if (_use_mass_fraction && _sp >= _dictator.numComponents())
    paramError(
        "mass_fraction_component",
        "The Dictator proclaims that the maximum fluid component index in this simulation is ",
        _dictator.numComponents() - 1,
        " whereas you have used ",
        _sp,
        ". Remember that indexing starts at 0. Please be assured that the Dictator has noted your "
        "error.");

  if (!_fv && _p_or_t == PorTchoice::pressure && !_has_porepressure)
    mooseError("PorousFlowLineSink: You have specified function_of=porepressure, but you do not "
               "have a quadpoint porepressure material");

  if (!_fv && _p_or_t == PorTchoice::temperature && !_has_temperature)
    mooseError("PorousFlowLineSink: You have specified function_of=temperature, but you do not "
               "have a quadpoint temperature material");

  if (!_fv && _use_mass_fraction && !_has_mass_fraction)
    mooseError("PorousFlowLineSink: You have specified a fluid component, but do not have a nodal "
               "mass-fraction material");

  if (!_fv && _use_relative_permeability && !_has_relative_permeability)
    mooseError("PorousFlowLineSink: You have set use_relative_permeability=true, but do not have a "
               "nodal relative permeability material");

  if (!_fv && _use_mobility && !_has_mobility)
    mooseError("PorousFlowLineSink: You have set use_mobility=true, but do not have nodal density, "
               "relative permeability or viscosity material");

  if (_mobility_type == MobilityTypeChoice::total && !_use_mobility)
    paramError("mobility_type",
               "You have set mobility_type=total, but use_mobility=false, so the flux is not "
               "multiplied by a mobility of any type");

  if (!_fv && _use_enthalpy && !_has_enthalpy)
    mooseError("PorousFlowLineSink: You have set use_enthalpy=true, but do not have a nodal "
               "enthalpy material");

  if (!_fv && _use_internal_energy && !_has_internal_energy)
    mooseError("PorousFlowLineSink: You have set use_internal_energy=true, but do not have a nodal "
               "internal-energy material");

  // To correctly compute the Jacobian terms,
  // tell MOOSE that this DiracKernel depends on all the PorousFlow Variables
  const std::vector<MooseVariableFEBase *> & coupled_vars = _dictator.getCoupledMooseVars();
  for (unsigned int i = 0; i < coupled_vars.size(); i++)
    addMooseVariableDependency(coupled_vars[i]);
}

void
PorousFlowLineSink::addPoints()
{
  // This function gets called just before the DiracKernel is evaluated
  // so this is a handy place to zero this out.
  _total_outflow_mass.zero();
  if (_point_fluxes)
    _point_fluxes->zero(*_x_coord, *_y_coord, *_z_coord);

  PorousFlowLineGeometry::addPoints();
}

Real
PorousFlowLineSink::massFraction(const Real base_outflow) const
{
  // A negative outflow is fluid entering the porespace, which carries the injected composition
  if (_use_injected_mass_fraction && base_outflow < 0.0)
    return _injected_mass_fraction->value(_t, _q_point[_qp]);

  return (*_mass_fractions)[_i][_ph][_sp];
}

Real
PorousFlowLineSink::weightingFactor(const Real base_outflow) const
{
  Real weight = _multiplying_var[_qp];

  if (_use_relative_permeability)
    weight *= (*_relative_permeability)[_i][_ph];

  if (_use_mobility)
    weight *= mobility();

  if (_use_mass_fraction)
    weight *= massFraction(base_outflow);

  if (_use_enthalpy)
    weight *= (*_enthalpy)[_i][_ph];

  if (_use_internal_energy)
    weight *= (*_internal_energy)[_i][_ph];

  return weight;
}

Real
PorousFlowLineSink::weightedOutflow()
{
  // Get the ID we initially assigned to this point
  const unsigned current_dirac_ptid = currentPointCachedID();
  const Real outflow = computeQpBaseOutflow(current_dirac_ptid);
  if (outflow == 0.0)
    return 0.0;

  return outflow * weightingFactor(outflow);
}

Real
PorousFlowLineSink::computeQpResidual()
{
  const Real outflow = weightedOutflow();

  _total_outflow_mass.add(
      outflow * _dt); // this is not thread safe, but DiracKernel's aren't currently threaded
  if (_point_fluxes)
    _point_fluxes->add(currentPointCachedID(), outflow);

  return outflow;
}

Real
PorousFlowLineSink::computeQpJacobian()
{
  return jac(_var.number());
}

Real
PorousFlowLineSink::computeQpOffDiagJacobian(unsigned int jvar)
{
  return jac(jvar);
}

Real
PorousFlowLineSink::jac(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  const unsigned pvar = _dictator.porousFlowVariableNum(jvar);

  Real outflow;
  Real outflowp;
  const unsigned current_dirac_ptid = currentPointCachedID();
  computeQpBaseOutflowJacobian(jvar, current_dirac_ptid, outflow, outflowp);
  if (outflow == 0.0 && outflowp == 0.0)
    return 0.0;

  // Which way the fluid moves, and so which composition it carries, is settled by the unweighted
  // flux, before any of the weightings below have multiplied into it
  const Real base_outflow = outflow;

  outflow *= _multiplying_var[_qp];
  outflowp *= _multiplying_var[_qp];

  if (_use_relative_permeability)
  {
    const Real relperm_prime = (_i != _j ? 0.0 : (*_drelative_permeability_dvar)[_i][_ph][pvar]);
    outflowp = (*_relative_permeability)[_i][_ph] * outflowp + relperm_prime * outflow;
    outflow *= (*_relative_permeability)[_i][_ph];
  }

  if (_use_mobility)
  {
    const Real mob = mobility();
    const Real mob_prime = (_i != _j ? 0.0 : dmobility(pvar));
    outflowp = mob * outflowp + mob_prime * outflow;
    outflow *= mob;
  }

  if (_use_mass_fraction)
  {
    // The injected composition is prescribed rather than a property of the reservoir, so where it
    // is used the derivative of the weighting vanishes
    const bool injecting = (base_outflow < 0.0);
    const Real mass_frac = massFraction(base_outflow);
    const Real mass_fractions_prime = (injecting && _use_injected_mass_fraction) || _i != _j
                                          ? 0.0
                                          : (*_dmass_fractions_dvar)[_i][_ph][_sp][pvar];
    outflowp = mass_frac * outflowp + mass_fractions_prime * outflow;
    outflow *= mass_frac;
  }

  if (_use_enthalpy)
  {
    const Real enthalpy_prime = (_i != _j ? 0.0 : (*_denthalpy_dvar)[_i][_ph][pvar]);
    outflowp = (*_enthalpy)[_i][_ph] * outflowp + enthalpy_prime * outflow;
    outflow *= (*_enthalpy)[_i][_ph];
  }

  if (_use_internal_energy)
  {
    const Real internal_energy_prime = (_i != _j ? 0.0 : (*_dinternal_energy_dvar)[_i][_ph][pvar]);
    outflowp = (*_internal_energy)[_i][_ph] * outflowp + internal_energy_prime * outflow;
    // this multiplication was performed, but the code does not need to know: outflow *=
    // (*_internal_energy)[_i][_ph];
  }

  return outflowp;
}

Real
PorousFlowLineSink::mobility() const
{
  if (_mobility_type == MobilityTypeChoice::phase)
    return (*_relative_permeability)[_i][_ph] * (*_fluid_density_node)[_i][_ph] /
           (*_fluid_viscosity)[_i][_ph];

  // Injected fluid displaces whatever is resident, so the resistance it meets is that of all the
  // phases together.  Unlike the phase mobility, this is nonzero in a block containing none of the
  // injected phase, which is what allows an injection well to begin flowing.  The mass rate is
  // carried at the density of the injected phase, so that factor remains the nominated phase's.
  Real total_mobility = 0.0;
  for (const auto p : make_range(_dictator.numPhases()))
    total_mobility += (*_relative_permeability)[_i][p] / (*_fluid_viscosity)[_i][p];
  return (*_fluid_density_node)[_i][_ph] * total_mobility;
}

Real
PorousFlowLineSink::dmobility(unsigned pvar) const
{
  if (_mobility_type == MobilityTypeChoice::phase)
    return (*_drelative_permeability_dvar)[_i][_ph][pvar] * (*_fluid_density_node)[_i][_ph] /
               (*_fluid_viscosity)[_i][_ph] +
           (*_relative_permeability)[_i][_ph] * (*_dfluid_density_node_dvar)[_i][_ph][pvar] /
               (*_fluid_viscosity)[_i][_ph] -
           (*_relative_permeability)[_i][_ph] * (*_fluid_density_node)[_i][_ph] *
               (*_dfluid_viscosity_dvar)[_i][_ph][pvar] /
               Utility::pow<2>((*_fluid_viscosity)[_i][_ph]);

  Real total_mobility = 0.0;
  Real dtotal_mobility = 0.0;
  for (const auto p : make_range(_dictator.numPhases()))
  {
    total_mobility += (*_relative_permeability)[_i][p] / (*_fluid_viscosity)[_i][p];
    dtotal_mobility += (*_drelative_permeability_dvar)[_i][p][pvar] / (*_fluid_viscosity)[_i][p] -
                       (*_relative_permeability)[_i][p] * (*_dfluid_viscosity_dvar)[_i][p][pvar] /
                           Utility::pow<2>((*_fluid_viscosity)[_i][p]);
  }
  return (*_dfluid_density_node_dvar)[_i][_ph][pvar] * total_mobility +
         (*_fluid_density_node)[_i][_ph] * dtotal_mobility;
}

Real
PorousFlowLineSink::ptqp() const
{
  return (_p_or_t == PorTchoice::pressure ? (*_pp)[_qp][_ph] : (*_temperature)[_qp]);
}

Real
PorousFlowLineSink::dptqp(unsigned pvar) const
{
  return (_p_or_t == PorTchoice::pressure ? (*_dpp_dvar)[_qp][_ph][pvar]
                                          : (*_dtemperature_dvar)[_qp][pvar]);
}
