# Finite-element counterpart of the single-phase rate-controlled case that the finite-volume
# Peaceman borehole was verified against MRST with (see compare.py for the reference).
#
# A 20 x 20 x 30 m box of 5 m cubes holds water at hydrostatic pressure.  A vertical borehole with
# three completions at element centres produces 1 kg/s until t = 50 s and then injects 0.5 kg/s.
# PorousFlowPeacemanBorehole (control = rate) solves for the bottom-hole pressure (the scalar bhp) that
# delivers that rate.
#
# The mesh is the finite-volume grid, so the pore volume is the same, but the porepressure lives
# at the element corners rather than the cell centres.  The borehole therefore sees the pressure
# interpolated to the completion, not a cell pressure, which is why the finite-element literature
# gives a different Peaceman re_constant (0.1594) from the finite-difference value (0.28).
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
  # Peaceman's finite-difference value.  compare.py also runs 0.1594, the value for rectangular
  # finite elements (Chen and Zhang, Int J Num Analysis and Modeling 3 (2008) 375-388, eqn 4.13)
  re_constant = 0.28
[]

[Variables]
  [pp]
  []
  [bhp]
    family = SCALAR
    order = FIRST
    # starting guess: hydrostatic pressure at the bottom completion
    initial_condition = 1.0025E7
  []
[]

[ICs]
  [pp]
    type = FunctionIC
    variable = pp
    function = '1E7 - 1E4 * z'
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
    fluid_component = 0
    variable = pp
    gravity = '0 0 -10'
  []
[]

[DiracKernels]
  [bh]
    type = PorousFlowPeacemanBorehole
    control = rate
    variable = pp
    SumQuantityUO = bh_outflow
    point_file = borehole.bh
    fluid_phase = 0
    unit_weight = '0 0 -1E4'
    use_mobility = true
    bottom_hole_pressure = bhp
    rate = rate
  []
[]

[ScalarKernels]
  [bhp]
    # The bhp equation is assembled by the borehole: this only tells MOOSE that bhp has one
    type = NullScalarKernel
    variable = bhp
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [bh_outflow]
    type = PorousFlowSumQuantity
  []
[]

[FluidProperties]
  [water]
    type = SimpleFluidProperties
    bulk_modulus = 2E9
    viscosity = 1E-3
    density0 = 1000
    thermal_expansion = 0
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-12 0 0 0 1E-12 0 0 0 1E-12'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  []
[]

[Functions]
  [rate]
    # kg/s, positive for production
    type = ParsedFunction
    expression = 'if(t <= 50, 1, -0.5)'
  []
[]

[Postprocessors]
  [p_bot]
    type = ScalarVariable
    variable = bhp
  []
  [mass]
    type = PorousFlowFluidMass
  []
  [mass_old]
    type = PorousFlowFluidMass
    execute_on = 'timestep_begin'
    outputs = none
  []
  [dt]
    type = TimestepSize
    outputs = none
  []
  [bh_rate]
    # mass rate leaving the porous medium, from the change in fluid mass
    type = ParsedPostprocessor
    expression = '(mass_old - mass) / dt'
    pp_names = 'mass mass_old dt'
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
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100'
  []
  nl_abs_tol = 1E-10
  nl_rel_tol = 1E-12
[]

[Outputs]
  csv = true
[]
