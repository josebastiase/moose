//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowWellAction.h"
#include "PorousFlowDictator.h"
#include "FEProblem.h"
#include "Factory.h"
#include "OutputWarehouse.h"
#include "MooseVariableFieldBase.h"

registerMooseAction("PorousFlowApp", PorousFlowWellAction, "add_porous_flow_well_gate");
registerMooseAction("PorousFlowApp", PorousFlowWellAction, "add_porous_flow_well");

InputParameters
PorousFlowWellAction::validParams()
{
  InputParameters params = MooseObjectAction::validParams();
  params.addClassDescription("Action that expands a [Wells] block into the DiracKernels, "
                             "UserObjects and Postprocessors that model the well.");
  return params;
}

PorousFlowWellAction::PorousFlowWellAction(const InputParameters & parameters)
  : MooseObjectAction(parameters), _num_components(0)
{
}

void
PorousFlowWellAction::act()
{
  const InputParameters & obj_pars = getObjectParams();

  if (_current_task == "add_porous_flow_well_variable")
  {
    if (solvesForBottomHolePressure())
      addBottomHolePressureVariable();
    return;
  }

  if (_current_task == "add_porous_flow_well_gate")
  {
    checkNotYetImplemented();
    addSyncTimes();

    // The PorousFlowDictator does not exist yet at this task, so nothing here may consult it
    const FunctionName gate = gateFunction();
    if (_type == "PorousFlowInjectionWell")
      gatedRate(gate, obj_pars.get<FunctionName>("rate"));
    else
      gateAuxVariable(gate);
    return;
  }

  const PorousFlowDictator & dictator = _problem->getUserObject<PorousFlowDictator>(
      obj_pars.get<UserObjectName>("PorousFlowDictator"));
  _num_components = dictator.numComponents();

  const auto & mass_fraction_vars = obj_pars.get<std::vector<VariableName>>("mass_fraction_vars");
  if (mass_fraction_vars.size() + 1 != _num_components)
    mooseError(name(),
               ": the simulation has ",
               _num_components,
               " fluid components, so 'mass_fraction_vars' must have ",
               _num_components - 1,
               " entries (the final component's equation lives on 'porepressure'), but ",
               mass_fraction_vars.size(),
               " were given");

  if (_type == "PorousFlowInjectionWell")
    addInjectionWell();
  else if (_type == "PorousFlowProductionWell")
  {
    if (finiteVolume())
      mooseError(name(),
                 ": a PorousFlowProductionWell is a deliverability-curve line sink, which is only "
                 "available for finite-element variables.  With finite-volume variables use a "
                 "PorousFlowPeacemanBorehole DiracKernel instead");
    addProductionWell();
  }
  // NOTE: no reason to catch unknown types here - MOOSE will catch that much earlier
}

bool
PorousFlowWellAction::solvesForBottomHolePressure() const
{
  // The mechanism follows from the constraints the well declares, and is never named by the user: a
  // well that must respect a pressure limit cannot be a prescribed-rate source, because the rate it
  // achieves is whatever the limit allows
  return getObjectParams().isParamValid("bhp_max");
}

VariableName
PorousFlowWellAction::bottomHolePressureVariable() const
{
  return name() + "_bhp";
}

void
PorousFlowWellAction::addBottomHolePressureVariable()
{
  const VariableName var_name = bottomHolePressureVariable();

  auto var_params = _factory.getValidParams("MooseVariableScalar");
  var_params.set<MooseEnum>("family") = "SCALAR";
  var_params.set<MooseEnum>("order") = "FIRST";
  _problem->addVariable("MooseVariableScalar", var_name, var_params);

  auto pars = _factory.getValidParams("NullScalarKernel");
  pars.set<NonlinearVariableName>("variable") = var_name;
  _problem->addScalarKernel("NullScalarKernel", var_name + "_null", pars);
}

void
PorousFlowWellAction::checkNotYetImplemented()
{
  const InputParameters & obj_pars = getObjectParams();
  for (const std::string & p : {"bhp", "bhp_min"})
    if (obj_pars.isParamValid(p))
      mooseError(name(),
                 ": '",
                 p,
                 "' requires the implicit scalar-BHP well, which is not yet implemented (Task 14 "
                 "phase 1b).  A phase-1a well is driven by a prescribed rate or deliverability "
                 "curve only.");

  for (const std::string & p : {"group", "reopen_if"})
    if (obj_pars.isParamValid(p))
      mooseError(name(), ": '", p, "' is not yet implemented.");
}

bool
PorousFlowWellAction::isGated() const
{
  const InputParameters & obj_pars = getObjectParams();
  return obj_pars.isParamValid("start_time") || obj_pars.isParamValid("end_time") ||
         obj_pars.isParamValid("shut_if");
}

VariableName
PorousFlowWellAction::componentVariable(unsigned int component) const
{
  const InputParameters & obj_pars = getObjectParams();
  if (component + 1 == _num_components)
    return obj_pars.get<VariableName>("porepressure");
  return obj_pars.get<std::vector<VariableName>>("mass_fraction_vars")[component];
}

void
PorousFlowWellAction::addSyncTimes()
{
  const InputParameters & obj_pars = getObjectParams();
  std::set<Real> & sync_times = _app.getOutputWarehouse().getSyncTimes();

  // So that a rate change never lands part-way through a timestep
  const auto & times = obj_pars.get<std::vector<Real>>("sync_times");
  sync_times.insert(times.begin(), times.end());

  // So that the solve lands exactly on the times this well switches on and off
  if (obj_pars.isParamValid("start_time"))
    sync_times.insert(obj_pars.get<Real>("start_time"));
  if (obj_pars.isParamValid("end_time"))
    sync_times.insert(obj_pars.get<Real>("end_time"));
}

void
PorousFlowWellAction::addRatePostprocessor(const PostprocessorName & pp_name,
                                           const FunctionName & func,
                                           Real scale)
{
  const std::string class_name = "FunctionValuePostprocessor";
  InputParameters pars = _factory.getValidParams(class_name);
  pars.set<FunctionName>("function") = func;
  pars.set<Real>("scale_factor") = scale;
  // The DiracKernel reads this value during the residual evaluation, so it must be current and
  // frozen across the nonlinear iterations of the step
  pars.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_BEGIN};
  _problem->addPostprocessor(class_name, pp_name, pars);
}

void
PorousFlowWellAction::addInjectionWell()
{
  const InputParameters & obj_pars = getObjectParams();
  const std::string nm = name();

  // A well completed at a single point is modelled by the point sources, and one completed over an
  // interval by the line sources, which divide the prescribed rate between the completions.  The
  // user declares the geometry and never chooses between the two.
  const bool single_completion = obj_pars.isParamValid("point");
  const bool completion_interval =
      obj_pars.isParamValid("point_file") || obj_pars.isParamValid("line_base");

  if (single_completion && completion_interval)
    mooseError(nm,
               ": the completions are given both by 'point' and by 'point_file' or 'line_base'.  "
               "Use 'point' for a well completed at a single point, and 'point_file' or "
               "'line_base' for one completed over an interval.");
  if (!single_completion && !completion_interval)
    mooseError(nm,
               ": a well's completions must be given by 'point', 'point_file' or 'line_base'.  "
               "[Wells] models a well as a point/line source in the model; areal or radial "
               "injection across a boundary is a boundary condition, so use PorousFlowSink "
               "instead.");

  // The gate task built this function; its name is derived from name(), not carried in a member
  const FunctionName rate =
      isGated() ? FunctionName(nm + "_rate_sched") : obj_pars.get<FunctionName>("rate");
  const auto & mass_fractions = obj_pars.get<std::vector<Real>>("injected_mass_fractions");

  if (mass_fractions.size() != _num_components)
    mooseError(nm,
               ": 'injected_mass_fractions' must have one entry per fluid component (",
               _num_components,
               "), but ",
               mass_fractions.size(),
               " were given");

  Real sum = 0.0;
  for (const auto w : mass_fractions)
  {
    if (w < 0.0)
      mooseError(nm, ": entries of 'injected_mass_fractions' cannot be negative");
    sum += w;
  }
  if (std::abs(sum - 1.0) > 1.0E-8)
    mooseError(nm, ": 'injected_mass_fractions' must sum to 1, but they sum to ", sum);

  if (solvesForBottomHolePressure())
  {
    if (finiteVolume() && obj_pars.isParamValid("temperature"))
      mooseError(nm,
                 ": a nonisothermal well with 'bhp_max' is not available for finite-volume "
                 "variables");
    addConstrainedInjectionWell(rate, mass_fractions);
    return;
  }

  // One point source per component, carrying q_i = w_i * q_total.  A positive mass_flux is already
  // a source in PorousFlowPointSourceFromPostprocessor, so the natural-sense sign passes straight
  // through
  for (const auto i : make_range(_num_components))
  {
    if (mass_fractions[i] == 0.0)
      continue;

    const PostprocessorName rate_pp = nm + "_q_" + Moose::stringify(i);
    addRatePostprocessor(rate_pp, rate, mass_fractions[i]);

    const std::string class_name = single_completion ? "PorousFlowPointSourceFromPostprocessor"
                                                     : "PorousFlowLineSourceFromPostprocessor";
    InputParameters pars = _factory.getValidParams(class_name);
    pars.set<NonlinearVariableName>("variable") = componentVariable(i);
    if (single_completion)
      pars.set<Point>("point") = obj_pars.get<Point>("point");
    else
      setLineGeometry(pars);
    pars.set<PostprocessorName>("mass_flux") = rate_pp;
    const std::string kernel_name = nm + "_source_" + Moose::stringify(i);
    _problem->addDiracKernel(class_name, kernel_name, pars);

    const std::string cum_class = "TimeIntegratedPostprocessor";
    InputParameters cum_pars = _factory.getValidParams(cum_class);
    cum_pars.set<PostprocessorName>("value") = rate_pp;
    _problem->addPostprocessor(cum_class, nm + "_cum_" + Moose::stringify(i), cum_pars);
  }

  if (!obj_pars.isParamValid("temperature"))
    return;

  // Nonisothermal: one enthalpy source, carrying the TOTAL mass rate at the injection temperature
  if (!obj_pars.isParamValid("fp") || !obj_pars.isParamValid("injection_temperature"))
    mooseError(nm,
               ": a nonisothermal well ('temperature' supplied) also needs 'fp' and "
               "'injection_temperature' so that the injected enthalpy can be computed");

  const PostprocessorName total_rate_pp = nm + "_q_total";
  addRatePostprocessor(total_rate_pp, rate, 1.0);

  const PostprocessorName t_in_pp = nm + "_T_in";
  addRatePostprocessor(t_in_pp, obj_pars.get<FunctionName>("injection_temperature"), 1.0);

  const std::string class_name = single_completion
                                     ? "PorousFlowPointEnthalpySourceFromPostprocessor"
                                     : "PorousFlowLineEnthalpySourceFromPostprocessor";
  InputParameters pars = _factory.getValidParams(class_name);
  pars.set<NonlinearVariableName>("variable") = obj_pars.get<VariableName>("temperature");
  if (single_completion)
    pars.set<Point>("point") = obj_pars.get<Point>("point");
  else
    setLineGeometry(pars);
  pars.set<PostprocessorName>("mass_flux") = total_rate_pp;
  pars.set<PostprocessorName>("T_in") = t_in_pp;
  pars.set<UserObjectName>("fp") = obj_pars.get<UserObjectName>("fp");
  pars.set<std::vector<VariableName>>("pressure") = {obj_pars.get<VariableName>("porepressure")};
  const std::string kernel_name = nm + "_enthalpy";
  _problem->addDiracKernel(class_name, kernel_name, pars);
}

bool
PorousFlowWellAction::finiteVolume() const
{
  bool fv = false;
  for (const auto i : make_range(_num_components))
  {
    const bool is_fv = _problem->getVariable(0, componentVariable(i)).isFV();
    if (i == 0)
      fv = is_fv;
    else if (is_fv != fv)
      mooseError(name(),
                 ": the fluid-component variables mix finite-element and finite-volume variables");
  }
  return fv;
}

void
PorousFlowWellAction::setLineGeometry(InputParameters & pars) const
{
  const InputParameters & obj_pars = getObjectParams();

  if (obj_pars.isParamValid("point_file"))
    pars.set<std::string>("point_file") = obj_pars.get<std::string>("point_file");
  if (obj_pars.isParamValid("line_base"))
    pars.set<std::vector<Real>>("line_base") = obj_pars.get<std::vector<Real>>("line_base");
  pars.set<Real>("line_length") = obj_pars.get<Real>("line_length");
  pars.set<RealVectorValue>("line_direction") = obj_pars.get<RealVectorValue>("line_direction");
}

void
PorousFlowWellAction::addConstrainedInjectionWell(const FunctionName & rate,
                                                  const std::vector<Real> & mass_fractions)
{
  const InputParameters & obj_pars = getObjectParams();
  const std::string nm = name();

  // The Peaceman well constant needs the radius of the hole, which point_file and line_base carry
  // as the first entry of each completion.  'point' has nowhere to put one
  if (obj_pars.isParamValid("point"))
    mooseError(nm,
               ": a well that respects 'bhp_max' is modelled by a Peaceman borehole, whose well "
               "constant depends on the radius of the hole.  Give the completions with "
               "'point_file' or 'line_base', whose first entry per completion is that radius, "
               "rather than with 'point'.");

  // [Wells] states rates positive for injection, the borehole states them positive for production
  const FunctionName injection_rate = nm + "_rate_injected";
  {
    InputParameters pars = _factory.getValidParams("ParsedFunction");
    pars.set<std::string>("expression") = "-q";
    pars.set<std::vector<std::string>>("symbol_names") = {"q"};
    pars.set<std::vector<std::string>>("symbol_values") = {rate};
    _problem->addFunction("ParsedFunction", injection_rate, pars);
  }

  // The same object for finite-element and finite-volume variables: it detects which
  const std::string class_name = "PorousFlowPeacemanBorehole";

  // Every borehole of this well shares one bottom-hole pressure, and each needs its own
  // PorousFlowSumQuantity, which accumulates and so cannot be shared
  auto boreholeParams = [&](const std::string & uo_name)
  {
    InputParameters pars = _factory.getValidParams(class_name);
    setLineGeometry(pars);
    pars.set<UserObjectName>("PorousFlowDictator") =
        obj_pars.get<UserObjectName>("PorousFlowDictator");
    pars.set<MooseEnum>("control") = "rate";
    pars.set<std::vector<VariableName>>("bottom_hole_pressure") = {bottomHolePressureVariable()};
    // The injected stream is carried at the density of the injection phase
    pars.set<unsigned int>("fluid_phase") = obj_pars.get<unsigned int>("injection_phase");
    // An injector weighted by the mobility of a phase absent from the reservoir has no injectivity
    // at all and can never begin to flow, so the total mobility is the only usable choice here
    pars.set<bool>("use_mobility") = true;
    pars.set<MooseEnum>("mobility_type") = "total";
    // The borehole wants the weight per unit volume of the fluid standing in the hole, which is
    // the wellbore density times gravity.  Both are the user's to give, since the density is that
    // of the injected stream rather than of anything in the formation, and the default of zero
    // leaves the wellbore with no gradient
    pars.set<RealVectorValue>("unit_weight") =
        obj_pars.get<Real>("wellbore_density") * obj_pars.get<RealVectorValue>("gravity");

    InputParameters uo_pars = _factory.getValidParams("PorousFlowSumQuantity");
    _problem->addUserObject("PorousFlowSumQuantity", uo_name, uo_pars);
    pars.set<UserObjectName>("SumQuantityUO") = uo_name;
    return pars;
  };

  // The borehole that owns the constraint takes no fluid component, so its flux is the well's
  // total, which is what 'rate' refers to.  It contributes no flux of its own: the per-component
  // boreholes below do that, and their fluxes sum to this one because the injected fractions do
  {
    InputParameters pars = boreholeParams(nm + "_total_uo");
    pars.set<NonlinearVariableName>("variable") = componentVariable(0);
    pars.set<bool>("compute_field_residuals") = false;
    pars.set<FunctionName>("rate") = injection_rate;
    pars.set<FunctionName>("bhp_max") = Moose::stringify(obj_pars.get<Real>("bhp_max"));
    _problem->addDiracKernel(class_name, nm + "_constraint", pars);
  }

  for (const auto i : make_range(_num_components))
  {
    if (mass_fractions[i] == 0.0)
      continue;

    InputParameters pars = boreholeParams(nm + "_uo_" + Moose::stringify(i));
    pars.set<NonlinearVariableName>("variable") = componentVariable(i);
    pars.set<unsigned int>("mass_fraction_component") = i;
    // Injected fluid carries the composition being injected, not the reservoir's
    pars.set<FunctionName>("injected_mass_fraction") = Moose::stringify(mass_fractions[i]);
    pars.set<bool>("compute_scalar_residuals") = false;
    _problem->addDiracKernel(class_name, nm + "_source_" + Moose::stringify(i), pars);
  }

  // No reporter is added for the pressure the well settled at: the scalar variable is itself
  // written to CSV, under the name bottomHolePressureVariable() gives it, so a postprocessor would
  // only duplicate it and would collide with the variable's own name
}

void
PorousFlowWellAction::addProductionWell()
{
  const InputParameters & obj_pars = getObjectParams();
  const std::string nm = name();

  if (!obj_pars.isParamValid("point_file") && !obj_pars.isParamValid("line_base"))
    mooseError(nm,
               ": a well's completions must be given by 'point_file' or 'line_base'.  [Wells] "
               "models a well as a point/line source in a 3D model; areal or radial injection "
               "across a boundary is a boundary condition, so use PorousFlowSink instead.");

  // A dt postprocessor, so that the cumulative mass reported by PorousFlowPlotQuantity can be
  // turned into a rate.  Hidden from the output: it is an implementation detail
  const PostprocessorName dt_pp = nm + "_dt";
  {
    const std::string class_name = "TimestepSize";
    InputParameters pars = _factory.getValidParams(class_name);
    pars.set<std::vector<OutputName>>("outputs") = {"none"};
    _problem->addPostprocessor(class_name, dt_pp, pars);
  }

  for (const auto i : make_range(_num_components))
  {
    // The line sink adds into its SumQuantity, so each kernel needs its own
    const UserObjectName sum_uo = nm + "_sum_" + Moose::stringify(i);
    {
      const std::string class_name = "PorousFlowSumQuantity";
      InputParameters pars = _factory.getValidParams(class_name);
      _problem->addUserObject(class_name, sum_uo, pars);
    }

    // Positive 'fluxes' are already flow out of the medium in PorousFlowPolyLineSink, which is
    // positive production: the natural-sense sign passes straight through
    {
      const std::string class_name = "PorousFlowPolyLineSink";
      InputParameters pars = _factory.getValidParams(class_name);
      pars.set<NonlinearVariableName>("variable") = componentVariable(i);
      pars.set<UserObjectName>("SumQuantityUO") = sum_uo;
      pars.set<UserObjectName>("PorousFlowDictator") =
          obj_pars.get<UserObjectName>("PorousFlowDictator");
      pars.set<unsigned int>("mass_fraction_component") = i;
      pars.set<std::vector<Real>>("p_or_t_vals") = obj_pars.get<std::vector<Real>>("p_or_t_vals");
      pars.set<std::vector<Real>>("fluxes") = obj_pars.get<std::vector<Real>>("fluxes");
      pars.set<MooseEnum>("function_of") = obj_pars.get<MooseEnum>("function_of");
      pars.set<unsigned int>("fluid_phase") = obj_pars.get<unsigned int>("fluid_phase");
      pars.set<bool>("use_mobility") = obj_pars.get<bool>("use_mobility");
      pars.set<bool>("use_relative_permeability") = obj_pars.get<bool>("use_relative_permeability");
      // A producer's rate comes from its deliverability curve rather than from a rate function, so
      // there is no rate to zero: its flux is multiplied by the gate AuxVariable that the gate task
      // built, whose name is likewise derived from name()
      if (isGated())
        pars.set<std::vector<VariableName>>("multiplying_var") = {VariableName(nm + "_gate_var")};
      setLineGeometry(pars);
      const std::string kernel_name = nm + "_sink_" + Moose::stringify(i);
      _problem->addDiracKernel(class_name, kernel_name, pars);
    }

    // PorousFlowPlotQuantity reports the mass produced during the timestep, in kg
    const PostprocessorName mass_pp = nm + "_mass_" + Moose::stringify(i);
    {
      const std::string class_name = "PorousFlowPlotQuantity";
      InputParameters pars = _factory.getValidParams(class_name);
      pars.set<UserObjectName>("uo") = sum_uo;
      _problem->addPostprocessor(class_name, mass_pp, pars);
    }
    {
      const std::string class_name = "CumulativeValuePostprocessor";
      InputParameters pars = _factory.getValidParams(class_name);
      pars.set<PostprocessorName>("postprocessor") = mass_pp;
      _problem->addPostprocessor(class_name, nm + "_cum_" + Moose::stringify(i), pars);
    }
    {
      const std::string class_name = "ParsedPostprocessor";
      InputParameters pars = _factory.getValidParams(class_name);
      pars.set<std::string>("expression") = "m / dt";
      pars.set<std::vector<PostprocessorName>>("pp_names") = {mass_pp, dt_pp};
      pars.set<std::vector<std::string>>("pp_symbols") = {"m", "dt"};
      _problem->addPostprocessor(class_name, nm + "_q_" + Moose::stringify(i), pars);
    }
  }
}

FunctionName
PorousFlowWellAction::gateFunction()
{
  const InputParameters & obj_pars = getObjectParams();
  const bool has_start = obj_pars.isParamValid("start_time");
  const bool has_end = obj_pars.isParamValid("end_time");
  const bool has_shut = obj_pars.isParamValid("shut_if");
  if (!has_start && !has_end && !has_shut)
    return "";

  std::vector<std::string> factors;
  std::vector<std::string> symbol_names;
  std::vector<std::string> symbol_values;

  if (has_start || has_end)
  {
    const Real never = 1.0E30;
    const Real start = has_start ? obj_pars.get<Real>("start_time") : -never;
    const Real end = has_end ? obj_pars.get<Real>("end_time") : never;
    if (end <= start)
      mooseError(name(), ": 'end_time' must be greater than 'start_time'");

    // The gate is sampled at the START of a step, at the step's END time t, and the resulting
    // constant factor applies across the whole of (t-dt, t].  The window is therefore half-open in
    // t: a step landing exactly on start_time lies entirely before the well opens, and a step
    // landing exactly on end_time is the last one during which the well flows.  With start_time and
    // end_time in the sync times, a fixed-rate well then delivers exactly rate * (end - start).
    factors.push_back("if(t > " + Moose::stringify(start) + " & t <= " + Moose::stringify(end) +
                      ", 1, 0)");
  }

  if (has_shut)
  {
    // ParsedFunction resolves a Postprocessor named in symbol_values, so the user's condition can
    // depend on the state of the simulation and not merely on t
    factors.push_back("if(s > 0.5, 0, 1)");
    symbol_names.push_back("s");
    symbol_values.push_back(obj_pars.get<FunctionName>("shut_if"));
  }

  std::string expression = factors[0];
  for (const auto i : make_range(std::size_t(1), factors.size()))
    expression += " * " + factors[i];

  const FunctionName gate_name = name() + "_gate";
  const std::string class_name = "ParsedFunction";
  InputParameters pars = _factory.getValidParams(class_name);
  pars.set<std::string>("expression") = expression;
  pars.set<std::vector<std::string>>("symbol_names") = symbol_names;
  pars.set<std::vector<std::string>>("symbol_values") = symbol_values;
  _problem->addFunction(class_name, gate_name, pars);

  return gate_name;
}

FunctionName
PorousFlowWellAction::gatedRate(const FunctionName & gate, const FunctionName & rate)
{
  if (gate.empty())
    return rate;

  const FunctionName gated_name = name() + "_rate_sched";
  const std::string class_name = "ParsedFunction";
  InputParameters pars = _factory.getValidParams(class_name);
  pars.set<std::string>("expression") = "g * r";
  pars.set<std::vector<std::string>>("symbol_names") = {"g", "r"};
  pars.set<std::vector<std::string>>("symbol_values") = {gate, rate};
  _problem->addFunction(class_name, gated_name, pars);

  return gated_name;
}

VariableName
PorousFlowWellAction::gateAuxVariable(const FunctionName & gate)
{
  if (gate.empty())
    return "";

  const VariableName var_name = name() + "_gate_var";
  {
    const std::string class_name = "MooseVariable";
    InputParameters pars = _factory.getValidParams(class_name);
    pars.set<MooseEnum>("family") = "LAGRANGE";
    pars.set<MooseEnum>("order") = "FIRST";
    _problem->addAuxVariable(class_name, var_name, pars);
  }
  {
    const std::string class_name = "FunctionAux";
    InputParameters pars = _factory.getValidParams(class_name);
    pars.set<AuxVariableName>("variable") = var_name;
    pars.set<FunctionName>("function") = gate;
    // The line sink reads this during the residual evaluation, so it must be current and frozen
    // across the nonlinear iterations of the step
    pars.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_BEGIN};
    _problem->addAuxKernel(class_name, name() + "_gate_aux", pars);
  }

  return var_name;
}
