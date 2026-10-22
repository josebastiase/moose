# Finite-element counterpart of the two-phase rate-controlled gas-injection case that the
# finite-volume Peaceman borehole was verified against MRST with (see compare.py for the
# reference).
#
# Water and an immiscible CO2-like gas, with no capillary pressure.  The box of 1ph_rate.i starts
# full of water and gas is injected at 1 kg/s.  The bottom-hole pressure (the scalar bhp) is solved
# for.
#
# The gas phase is absent at the start, so a borehole weighted by the gas mobility could not
# inject at all.  mobility_type = total weights the injected gas by the total mobility of the
# fluid it displaces, which is the usual treatment of an injector.  The gas phase holds only the
# gas component, so this one borehole's flux is the whole well rate, and it both carries the flux
# and owns the rate constraint.  No water borehole is needed because an injector moves no water.
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
  [ppwater]
  []
  [sgas]
    initial_condition = 0
  []
  [bhp]
    family = SCALAR
    order = FIRST
    # starting guess: hydrostatic pressure at the bottom completion
    initial_condition = 1.0025E7
  []
[]

[AuxVariables]
  [massfrac_ph0_sp0]
    initial_condition = 1
  []
  [massfrac_ph1_sp0]
    initial_condition = 0
  []
[]

[ICs]
  [ppwater]
    type = FunctionIC
    variable = ppwater
    function = '1E7 - 1E4 * z'
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = ppwater
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = ppwater
    gravity = '0 0 -10'
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = sgas
    gravity = '0 0 -10'
  []
[]

[DiracKernels]
  [bh_gas]
    type = PorousFlowPeacemanBorehole
    control = rate
    variable = sgas
    SumQuantityUO = bh_gas_outflow
    point_file = borehole.bh
    fluid_phase = 1
    mass_fraction_component = 1
    unit_weight = '0 0 -1E4'
    use_mobility = true
    mobility_type = total
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
    porous_flow_vars = 'ppwater sgas'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
  [bh_gas_outflow]
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
  [gas]
    # CO2-like: density = 364 * exp(p / 2E7), about 600 kg/m^3 at 10 MPa
    type = SimpleFluidProperties
    bulk_modulus = 2E7
    viscosity = 5E-5
    density0 = 364
    thermal_expansion = 0
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = ppwater
    phase1_saturation = sgas
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [gas]
    type = PorousFlowSingleComponentFluid
    fp = gas
    phase = 1
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-12 0 0 0 1E-12 0 0 0 1E-12'
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
  []
[]

[Functions]
  [rate]
    # kg/s, positive for production
    type = ParsedFunction
    expression = '-1'
  []
[]

[Postprocessors]
  [p_bot]
    type = ScalarVariable
    variable = bhp
  []
  # Mass rates leaving the porous medium, from the change in the fluid mass of each component
  [water_mass]
    type = PorousFlowFluidMass
    fluid_component = 0
  []
  [water_mass_old]
    type = PorousFlowFluidMass
    fluid_component = 0
    execute_on = 'timestep_begin'
    outputs = none
  []
  [gas_mass]
    type = PorousFlowFluidMass
    fluid_component = 1
  []
  [gas_mass_old]
    type = PorousFlowFluidMass
    fluid_component = 1
    execute_on = 'timestep_begin'
    outputs = none
  []
  [dt]
    type = TimestepSize
    outputs = none
  []
  [water_rate]
    type = ParsedPostprocessor
    expression = '(water_mass_old - water_mass) / dt'
    pp_names = 'water_mass water_mass_old dt'
  []
  [gas_rate]
    type = ParsedPostprocessor
    expression = '(gas_mass_old - gas_mass) / dt'
    pp_names = 'gas_mass gas_mass_old dt'
  []
  [bh_rate]
    type = ParsedPostprocessor
    expression = 'water_rate + gas_rate'
    pp_names = 'water_rate gas_rate'
  []
  # Gas saturation at the top and bottom completions.  These are interpolated from the element
  # corners, so they are not the cell saturations the finite-volume reference reports
  [sgas_top]
    type = PointValue
    variable = sgas
    point = '2.5 2.5 7.5'
  []
  [sgas_bot]
    type = PointValue
    variable = sgas
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
