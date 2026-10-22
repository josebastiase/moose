# Heat produced by a rate-controlled borehole, compared with the MRST geothermal module
# (mrst/bh_heat_mrst.m).
#
# A 20 x 20 x 30 m box of 5 m cubes holds water whose temperature increases with depth, from
# 293.15 K at the top to 353.15 K at the bottom.  A vertical borehole with three completions
# produces 1 kg/s of water.  PorousFlowPeacemanBorehole (control = rate) on pp solves for the
# bottom-hole pressure that delivers that rate, and a second PorousFlowPeacemanBorehole on the
# temperature, sharing that bottom-hole pressure, removes the enthalpy of the produced water.
# The heat rate depends on how the production is split between the completions, which are at
# different temperatures.
#
# This is the finite-volume version: heat_fe.i is the same model with finite-element variables.
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 6
  xmin = -10
  xmax = 10
  ymin = -10
  ymax = 10
  zmin = -15
  zmax = 15
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 -10'
[]

[Variables]
  [pp]
    type = MooseVariableFVReal
  []
  [temperature]
    type = MooseVariableFVReal
  []
  [bhp]
    family = SCALAR
    order = FIRST
    initial_condition = 1.0025E7
  []
[]

[ICs]
  [pp]
    type = FunctionIC
    variable = pp
    function = '1E7 - 1E4 * z'
  []
  [temperature]
    type = FunctionIC
    variable = temperature
    function = '293.15 + 2 * (15 - z)'
  []
[]

[FVKernels]
  [mass0]
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
  [flux0]
    type = FVPorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pp
  []
  [energy_dot]
    type = FVPorousFlowEnergyTimeDerivative
    variable = temperature
  []
  [heat_advection]
    type = FVPorousFlowHeatAdvection
    variable = temperature
  []
  [heat_conduction]
    type = FVPorousFlowHeatConduction
    variable = temperature
  []
[]

[DiracKernels]
  [bh_mass]
    type = PorousFlowPeacemanBorehole
    control = rate
    variable = pp
    SumQuantityUO = mass_uo
    point_file = borehole.bh
    unit_weight = '0 0 -1E4'
    use_mobility = true
    bottom_hole_pressure = bhp
    rate = 1
  []
  [bh_heat]
    type = PorousFlowPeacemanBorehole
    control = rate
    variable = temperature
    SumQuantityUO = heat_uo
    point_file = borehole.bh
    unit_weight = '0 0 -1E4'
    use_mobility = true
    use_enthalpy = true
    bottom_hole_pressure = bhp
    compute_scalar_residuals = false
  []
[]

[ScalarKernels]
  [bhp]
    type = NullScalarKernel
    variable = bhp
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp temperature'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [mass_uo]
    type = PorousFlowSumQuantity
  []
  [heat_uo]
    type = PorousFlowSumQuantity
  []
[]

[FluidProperties]
  [water]
    # density 1000 kg/m^3 at p = 0 and T = 293.15 K, as in MRST
    type = SimpleFluidProperties
    bulk_modulus = 2E9
    viscosity = 1E-3
    density0 = 1060.382826502196
    thermal_expansion = 2E-4
    cv = 4200
    cp = 4200
  []
[]

[Materials]
  [temperature]
    type = ADPorousFlowTemperature
    temperature = temperature
  []
  [ppss]
    type = ADPorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [massfrac]
    type = ADPorousFlowMassFraction
  []
  [water]
    type = ADPorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [porosity]
    type = ADPorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability]
    type = ADPorousFlowPermeabilityConst
    permeability = '1E-12 0 0 0 1E-12 0 0 0 1E-12'
  []
  [relperm]
    type = ADPorousFlowRelativePermeabilityConst
    phase = 0
  []
  [rock_heat]
    type = ADPorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1000
    density = 2700
  []
  [thermal_conductivity]
    # 2 * (1 - 0.1) + 0.6 * 0.1
    type = ADPorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '1.86 0 0 0 1.86 0 0 0 1.86'
  []
[]

[Postprocessors]
  [p_bot]
    type = ScalarVariable
    variable = bhp
  []
  [mass_out]
    type = PorousFlowPlotQuantity
    uo = mass_uo
    outputs = none
  []
  [heat_out]
    type = PorousFlowPlotQuantity
    uo = heat_uo
    outputs = none
  []
  [dt]
    type = TimestepSize
    outputs = none
  []
  [rate]
    type = ParsedPostprocessor
    expression = 'mass_out / dt'
    pp_names = 'mass_out dt'
  []
  [energy_rate]
    type = ParsedPostprocessor
    expression = 'heat_out / dt'
    pp_names = 'heat_out dt'
  []
  [T_top]
    type = PointValue
    variable = temperature
    point = '2.5 2.5 7.5'
  []
  [T_bot]
    type = PointValue
    variable = temperature
    point = '2.5 2.5 -2.5'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 100
  dt = 5
  # the energy residual is about 1E6 times the mass residual: 1E-3 is round-off for it
  nl_abs_tol = 1E-3
  nl_rel_tol = 1E-12
[]

[Outputs]
  csv = true
[]
