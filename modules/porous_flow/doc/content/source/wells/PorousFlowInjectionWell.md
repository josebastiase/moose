# PorousFlowInjectionWell

!syntax description /Wells/PorousFlowInjectionWell

A rate-controlled injection well, specified inside the [Wells](syntax/Wells/index.md) block.  The
well injects a fluid of fixed composition at a total mass rate given by the `rate` function.

The completions are described in either of two ways.  `point` gives the coordinates of a well
completed at a single point.  `point_file` or `line_base`, exactly as for
[PorousFlowProductionWell](PorousFlowProductionWell.md), instead describe a completion interval; with
`line_base` one completion is generated for each element the interval passes through, and those
completions follow the mesh if it is adapted or displaced.  The prescribed rate is divided between
the completions in proportion to the length of interval each represents, so the total injected is the
same however finely the interval is resolved.

The [PorousFlowWellAction](/Wells/index.md) expands this object into:

- one mass source per fluid component, each carrying $q_i = w_i q_\mathrm{total}$, where $w_i$ is the
  $i$th entry of `injected_mass_fractions`.  Components with a zero mass fraction are skipped.  These
  are [PorousFlowPointSourceFromPostprocessor](PorousFlowPointSourceFromPostprocessor.md) for a well
  completed at a single `point`, and
  [PorousFlowLineSourceFromPostprocessor](PorousFlowLineSourceFromPostprocessor.md) for one completed
  over an interval;
- when `temperature` is supplied, the matching
  [point](PorousFlowPointEnthalpySourceFromPostprocessor.md) or
  [line](PorousFlowLineEnthalpySourceFromPostprocessor.md) enthalpy source, carrying the *total* mass
  rate at the injection temperature;
- the rate postprocessors `<well>_q_<i>` and the cumulative-mass postprocessors `<well>_cum_<i>`.

Which source object is built follows from the geometry given, and is not chosen by the user.

Because the rate is prescribed rather than computed from a bottom-hole pressure, this well has no
injectivity limit and cannot fail to bootstrap: it is safe from a cold start, even when injecting
gas into a fully brine-saturated reservoir.  The cost is that the well has no feedback, so it will
pressurise the near-well cell without limit if the prescribed rate exceeds what the formation can
accept.  Supplying `bhp_max` gives that feedback: the well is then modelled by
[PorousFlowPeacemanBorehole](PorousFlowPeacemanBorehole.md) DiracKernels with `control = rate`, one
per fluid phase and component, which share the bottom-hole pressure, deliver `rate` until `bhp_max`
is reached and then hold it.

## Finite-volume variables

`[Wells]` works out whether the fluid-component variables are finite-element or finite-volume
variables, and builds the matching objects.  The point and line sources of a prescribed-rate well
need no material properties, so they are used as they are.  A well with `bhp_max` uses the same
[PorousFlowPeacemanBorehole](PorousFlowPeacemanBorehole.md) DiracKernels as with finite elements,
which detect the finite-volume variables themselves.  Their injecting completions inject
`injection_phase`, with its density and the total mobility of the cell, and the composition given by
`injected_mass_fractions`; completions where fluid flows out of the reservoir (crossflow) produce the
cell's fluid.

## Sign convention

`rate` is positive in the well's natural sense: a positive rate injects.  Every reported rate and
cumulative is likewise positive during injection.

## Shutting the well in

`start_time` and `end_time` shut the well in outside the operating window, and are added to the
simulation's sync times so that the solve lands exactly on them.  The well is shut by zeroing its
rate, so its rate and cumulative reporters correctly read zero while it is shut.

`shut_if` shuts the well in whenever a supplied `Function` evaluates to a nonzero (true) value.
Because a `ParsedFunction` resolves a `Postprocessor` named in its `symbol_values`, the condition
can depend on the state of the simulation rather than only on time:

```
[too_wet_fn]
  type = ParsedFunction
  expression = 'wc > 0.95'
  symbol_names = 'wc'
  symbol_values = 'my_water_cut_pp'
[]
```

The postprocessor is read from the end of the previous timestep, so the well shuts one timestep
after the condition is first met.  Unlike `start_time`/`end_time`, the crossing time is not known in
advance and so cannot be a sync time.  There is no hysteresis: a well whose trigger quantity hovers
on the threshold will chatter from step to step.

Rates are sampled at the *end* of each timestep, so a well opening at `start_time` first flows
during the step that ends at `start_time + dt`.  A `PiecewiseConstant` rate schedule should
therefore use `direction = left`.

!syntax parameters /Wells/PorousFlowInjectionWell

!syntax inputs /Wells/PorousFlowInjectionWell

!syntax children /Wells/PorousFlowInjectionWell
