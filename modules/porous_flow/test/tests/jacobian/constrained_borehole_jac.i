# PorousFlowPeacemanBorehole with 2-phase, 3-components, with enthalpy, internal_energy, and thermal_conductivity
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = -1
  xmax = 1
  ymin = -1
  ymax = 1
  zmin = -1
  zmax = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [ppwater]
  []
  [ppgas]
  []
  [massfrac_ph0_sp0]
  []
  [massfrac_ph0_sp1]
  []
  [massfrac_ph1_sp0]
  []
  [massfrac_ph1_sp1]
  []
  [temp]
  []
  [bhp]
    family = SCALAR
    order = FIRST
    initial_condition = 0.5
  []
[]

[ScalarKernels]
  [bhp_null]
    # checkKernelCoverage cannot see a variable whose equation a DiracKernel assembles
    type = NullScalarKernel
    variable = bhp
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'temp ppwater ppgas massfrac_ph0_sp0 massfrac_ph0_sp1 massfrac_ph1_sp0 massfrac_ph1_sp1'
    number_fluid_phases = 2
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    m = 0.5
    alpha = 1
  []
  [dummy_outflow0]
    type = PorousFlowSumQuantity
  []
  [dummy_outflow1]
    type = PorousFlowSumQuantity
  []
  [dummy_outflow2]
    type = PorousFlowSumQuantity
  []
  [dummy_outflow3]
    type = PorousFlowSumQuantity
  []
  [dummy_outflow4]
    type = PorousFlowSumQuantity
  []
  [dummy_outflow5]
    type = PorousFlowSumQuantity
  []
  [dummy_outflow6]
    type = PorousFlowSumQuantity
  []
  [dummy_outflow7]
    type = PorousFlowSumQuantity
  []
[]

[ICs]
  [temp]
    type = RandomIC
    variable = temp
    min = 1
    max = 2
  []
  [ppwater]
    type = RandomIC
    variable = ppwater
    min = -1
    max = 0
  []
  [ppgas]
    type = RandomIC
    variable = ppgas
    min = 0
    max = 1
  []
  [massfrac_ph0_sp0]
    type = RandomIC
    variable = massfrac_ph0_sp0
    min = 0
    max = 1
  []
  [massfrac_ph0_sp1]
    type = RandomIC
    variable = massfrac_ph0_sp1
    min = 0
    max = 1
  []
  [massfrac_ph1_sp0]
    type = RandomIC
    variable = massfrac_ph1_sp0
    min = 0
    max = 1
  []
  [massfrac_ph1_sp1]
    type = RandomIC
    variable = massfrac_ph1_sp1
    min = 0
    max = 1
  []
[]

[Kernels]
  [dummy_temp]
    type = TimeDerivative
    variable = temp
  []
  [dummy_ppwater]
    type = TimeDerivative
    variable = ppwater
  []
  [dummy_ppgas]
    type = TimeDerivative
    variable = ppgas
  []
  [dummy_m00]
    type = TimeDerivative
    variable = massfrac_ph0_sp0
  []
  [dummy_m01]
    type = TimeDerivative
    variable = massfrac_ph0_sp1
  []
  [dummy_m10]
    type = TimeDerivative
    variable = massfrac_ph1_sp0
  []
  [dummy_m11]
    type = TimeDerivative
    variable = massfrac_ph1_sp1
  []
[]

[FluidProperties]
  [simple_fluid0]
    type = SimpleFluidProperties
    bulk_modulus = 1.5
    density0 = 1
    thermal_expansion = 0
    viscosity = 1
    cv = 1.1
  []
  [simple_fluid1]
    type = SimpleFluidProperties
    bulk_modulus = 0.5
    density0 = 0.5
    thermal_expansion = 0
    viscosity = 1.4
    cv = 1.8
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temp
  []
  [ppss]
    type = PorousFlow2PhasePP
    phase0_porepressure = ppwater
    phase1_porepressure = ppgas
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph0_sp1 massfrac_ph1_sp0 massfrac_ph1_sp1'
  []
  [simple_fluid0]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid0
    phase = 0
  []
  [simple_fluid1]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid1
    phase = 1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1 0 0 0 2 0 0 0 3'
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 3
    phase = 1
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '0.1 0.02 0.03 0.02 0.0 0.01 0.03 0.01 0.3'
  []
[]

[DiracKernels]
  [dirac0]
    type = PorousFlowPeacemanBorehole
    control = rate
    fluid_phase = 0
    variable = ppwater
    point_file = one_point.bh
    line_length = 1
    SumQuantityUO = dummy_outflow0
    unit_weight = '1 2 3'
    re_constant = 0.123
    use_mobility = true
    mass_fraction_component = 0
    bottom_hole_pressure = bhp
    rate = 1E-2
  []
[]

[Preconditioning]
  [check]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1
  end_time = 1
[]

[Outputs]
  file_base = line_sink01
[]
