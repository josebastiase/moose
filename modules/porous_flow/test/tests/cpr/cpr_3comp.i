# Single-phase, three-component pressure pulse in 1D solved with the CPR preconditioner (block
# size 3).  Same model as pressure_pulse/pressure_pulse_1d_3comp.i.
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0
  xmax = 100
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [pp]
    initial_condition = 2E6
  []
  [massfrac0]
    initial_condition = 0.1
  []
  [massfrac1]
    initial_condition = 0.3
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    variable = pp
    gravity = '0 0 0'
    fluid_component = 0
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = massfrac0
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    variable = massfrac0
    gravity = '0 0 0'
    fluid_component = 1
  []
  [mass2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 2
    variable = massfrac1
  []
  [flux2]
    type = PorousFlowAdvectiveFlux
    variable = massfrac1
    gravity = '0 0 0'
    fluid_component = 2
  []
[]

[UserObjects]
  [cpr_layout]
    type = CPRLayoutCheck
    pressure_variable = pp
  []
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp massfrac0 massfrac1'
    number_fluid_phases = 1
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    m = 0.5
    alpha = 1e-7
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2e9
    density0 = 1000
    thermal_expansion = 0
    viscosity = 1e-3
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [ppss]
    type = PorousFlow1PhaseP
    porepressure = pp
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac0 massfrac1'
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-15 0 0 0 1E-15 0 0 0 1E-15'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 0
    phase = 0
  []
[]

[BCs]
  [left]
    type = DirichletBC
    boundary = left
    value = 3E6
    variable = pp
  []
[]

[Preconditioning]
  [cpr]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -ksp_pc_side -pc_type -pc_cpr_block_size -pc_cpr_pressure_field'
    petsc_options_value = 'gmres right cpr 3 0'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1E3
  end_time = 1E4
  nl_abs_tol = 1E-7
  nl_rel_tol = 1E-10
  nl_max_its = 20
[]

[Postprocessors]
  [p000]
    type = PointValue
    variable = pp
    point = '0 0 0'
    execute_on = 'initial timestep_end'
  []
  [p010]
    type = PointValue
    variable = pp
    point = '10 0 0'
    execute_on = 'initial timestep_end'
  []
  [p020]
    type = PointValue
    variable = pp
    point = '20 0 0'
    execute_on = 'initial timestep_end'
  []
  [p030]
    type = PointValue
    variable = pp
    point = '30 0 0'
    execute_on = 'initial timestep_end'
  []
  [p040]
    type = PointValue
    variable = pp
    point = '40 0 0'
    execute_on = 'initial timestep_end'
  []
  [p050]
    type = PointValue
    variable = pp
    point = '50 0 0'
    execute_on = 'initial timestep_end'
  []
  [p060]
    type = PointValue
    variable = pp
    point = '60 0 0'
    execute_on = 'initial timestep_end'
  []
  [p070]
    type = PointValue
    variable = pp
    point = '70 0 0'
    execute_on = 'initial timestep_end'
  []
  [p080]
    type = PointValue
    variable = pp
    point = '80 0 0'
    execute_on = 'initial timestep_end'
  []
  [p090]
    type = PointValue
    variable = pp
    point = '90 0 0'
    execute_on = 'initial timestep_end'
  []
  [p100]
    type = PointValue
    variable = pp
    point = '100 0 0'
    execute_on = 'initial timestep_end'
  []
  [mf_0_010]
    type = PointValue
    variable = massfrac0
    point = '10 0 0'
    execute_on = 'timestep_end'
  []
  [mf_1_010]
    type = PointValue
    variable = massfrac1
    point = '10 0 0'
    execute_on = 'timestep_end'
  []
[]

[Outputs]
  file_base = cpr_3comp
  print_linear_residuals = true
  csv = true
[]
