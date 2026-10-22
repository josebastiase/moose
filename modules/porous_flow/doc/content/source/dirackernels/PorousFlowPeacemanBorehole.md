# PorousFlowPeacemanBorehole

A `PorousFlowPeacemanBorehole` is a special case of the general line sink in which a polyline (represented by a sequence of points) acts as a sink or source in the model.  Please see [sinks](sinks.md) for an extended discussion and examples.

!alert warning
The function given by [!param](/DiracKernels/PorousFlowPeacemanBorehole/bottom_p_or_t) is evaluated at the well bottom.  If a file is read in using [!param](/DiracKernels/PorousFlowPeacemanBorehole/point_file) to define the coordinates and weights of the PorousFlowPeacemanBorehole, the well bottom is assumed to be the last entry in this file and [!param](/DiracKernels/PorousFlowPeacemanBorehole/bottom_p_or_t) will be evaluated at the z-coordinate of the last entry in [!param](/DiracKernels/PorousFlowPeacemanBorehole/point_file).  It is an error if the first entry in the [!param](/DiracKernels/PorousFlowPeacemanBorehole/point_file) has a smaller z-coordinate than the last entry.

!alert note
The wellbore pressure along the borehole can be built in two ways.  By default, a single constant
[!param](/DiracKernels/PorousFlowPeacemanBorehole/unit_weight) (fluid density $\times$ gravity) is
used, exactly as described in [sinks](sinks.md).  If instead
[!param](/DiracKernels/PorousFlowPeacemanBorehole/unit_weight_fp) is supplied (together with
[!param](/DiracKernels/PorousFlowPeacemanBorehole/unit_weight_temperature) and
[!param](/DiracKernels/PorousFlowPeacemanBorehole/unit_weight_gravity)), the fluid unit weight is
instead computed at each borehole point from a temperature-dependent fluid density, and
integrated along the well &mdash; see [sinks](sinks.md) for the formula, its justification, and
the deliberate approximation made in its Jacobian (the residual always uses the fluid density
implied by the current nonlinear iterate's temperature, but the Jacobian does not differentiate
the wellbore pressure with respect to temperature, since that dependence couples this borehole
point's residual to temperature degrees of freedom in other elements along the well, which a
DiracKernel cannot assemble into).  `unit_weight` and `unit_weight_fp` are mutually exclusive,
and `unit_weight_fp` is not compatible with `function_of = temperature`.  These four new
parameters are deliberately not named `fp`/`gravity`/`temperature_variable` (even though that
would mirror convention elsewhere in PorousFlow) because those are common `[GlobalParams]` names
used by unrelated Darcy kernels or fluid-properties materials; an input file that sets `gravity`
or `fp` at the `[GlobalParams]` level would otherwise silently activate, or fail to validate,
this mode on every `PorousFlowPeacemanBorehole` in the input.

To report the flow rate at each individual borehole point (rather than only the well-wide total
available via [!param](/DiracKernels/PorousFlowPeacemanBorehole/SumQuantityUO)), set
[!param](/DiracKernels/PorousFlowPeacemanBorehole/PointFluxUO) to a
[`PorousFlowPointFluxQuantity`](PorousFlowPointFluxQuantity.md) UserObject, and read it out with a
[`PorousFlowPlotPointFluxQuantity`](PorousFlowPlotPointFluxQuantity.md) VectorPostprocessor.

## Bottom-hole pressure or rate control

[!param](/DiracKernels/PorousFlowPeacemanBorehole/control) chooses what is prescribed:

- `control = bhp` (the default) prescribes the bottom-hole pressure (or temperature), the function
  [!param](/DiracKernels/PorousFlowPeacemanBorehole/bottom_p_or_t), and
  [!param](/DiracKernels/PorousFlowPeacemanBorehole/character) makes the borehole a producer or an
  injector, as described in [sinks](sinks.md).
- `control = rate` prescribes the total mass rate,
  [!param](/DiracKernels/PorousFlowPeacemanBorehole/rate) (positive for production, negative for
  injection), and solves for the bottom-hole pressure, which is a first-order scalar variable given by
  [!param](/DiracKernels/PorousFlowPeacemanBorehole/bottom_hole_pressure).

With rate control the flux at each borehole point is the same Peaceman expression, with the borehole
pressure $P_{\mathrm{bot}} + \gamma \cdot (x_{i} - x^{\mathrm{bot}})$ built from the scalar variable
and [!param](/DiracKernels/PorousFlowPeacemanBorehole/unit_weight).  The bottom-hole pressure is
solved fully coupled with the PorousFlow variables, from the equation
\begin{equation}
\sum_{\mathrm{DiracKernels}} \sum_{i} Q_{i} - Q_{\mathrm{rate}}(t) = 0 \ .
\end{equation}
All the `PorousFlowPeacemanBorehole` DiracKernels that share the scalar variable are one borehole:
use one per fluid phase and component, each acting on that component's mass-balance variable.
`use_mobility = true` is needed so that the fluxes are mass fluxes.  A DiracKernel that should not
contribute to the rate (an energy borehole, for instance) sets `compute_scalar_residuals = false`.
The split of the rate between the borehole points is not prescribed: it follows the local
permeability and pressure difference, and the flow at each point may go in either direction
(crossflow), so `character` and `bottom_p_or_t` are not used.  The scalar variable needs a
[NullScalarKernel](NullScalarKernel.md), and its initial condition is the starting guess for the
bottom-hole pressure.  For injection of a phase that is not present in the reservoir, use
`mobility_type = total` and [!param](/DiracKernels/PorousFlowPeacemanBorehole/injected_mass_fraction).

With [!param](/DiracKernels/PorousFlowPeacemanBorehole/bhp_max), the borehole delivers the rate
only while the bottom-hole pressure stays below `bhp_max`; where it would need more, the bottom-hole
pressure is held at `bhp_max` and the borehole delivers whatever rate follows.  This is decided once
per timestep (at most [!param](/DiracKernels/PorousFlowPeacemanBorehole/max_constraint_swaps) times
per timestep), so the bottom-hole pressure may exceed `bhp_max` for the timestep in which the limit
is first reached.

!listing modules/porous_flow/test/tests/dirackernels/fv_borehole/fv_borehole_rate.i block=Variables DiracKernels ScalarKernels

## Finite-element and finite-volume variables

The same object works with finite-element and finite-volume PorousFlow models, and detects which
from its variable.  With finite-element variables it uses the nodal PorousFlow material properties.
With finite-volume variables it uses the AD material properties of the cell containing each borehole
point, and its residual and Jacobian contributions are computed by automatic differentiation.  The
default `re_constant = 0.28` is Peaceman's value for cell-centred grids, which is the appropriate one
for finite volumes.  `unit_weight_fp` is only available with finite-element variables.

## Heat

With `use_enthalpy = true` on the temperature variable, the borehole removes or adds the enthalpy
carried by the fluid.  With rate control it shares the bottom-hole pressure of the mass boreholes
and sets `compute_scalar_residuals = false`:

!listing modules/porous_flow/examples/borehole_mrst/heat_fv.i block=DiracKernels

## Comparison with MRST

The borehole has been compared with MRST 2023b [!citep](lie2019mrst) (`ad-blackoil` and, for heat,
`geothermal` modules).  The model is a $20 \times 20 \times 30$ m box of $5$ m cubic cells with
permeability $10^{-12}$ m$^{2}$, porosity $0.1$, gravity $10$ m s$^{-2}$ and a hydrostatic initial
porepressure.  A vertical borehole of radius $0.1$ m has three points at cell centres, $5$ m apart.
Both codes use the same Peaceman well constant (`re_constant = 0.28`) and a wellbore pressure
gradient of $10^{4}$ Pa m$^{-1}$, with $20$ time steps of $5$ s.  Water has density
$1000 \exp(P/2\times 10^{9})$ kg m$^{-3}$ and viscosity $10^{-3}$ Pa s.  The two-phase cases add an
immiscible CO$_{2}$-like gas (density $364 \exp(P/2\times 10^{7})$ kg m$^{-3}$, viscosity
$5\times 10^{-5}$ Pa s), Corey relative permeabilities with exponent $2$ and no capillary pressure.
In the heat case the temperature increases from $293.15$ K at the top to $353.15$ K at the bottom,
the water density also decreases with temperature (expansion coefficient $2\times 10^{-4}$ K$^{-1}$),
its internal energy is $4200 T$ J kg$^{-1}$, and the borehole produces $1$ kg s$^{-1}$.

!media media/porous_flow/borehole_vs_mrst.png style=width:100%;margin-left:10px caption=MOOSE with finite-element (FE) and finite-volume (FV) variables, and MRST.  The input files are in `modules/porous_flow/examples/borehole_mrst`.  id=borehole_vs_mrst.fig

With finite volumes, MOOSE and MRST solve the same discrete equations, and the results agree
closely ([borehole_mrst]).  With finite elements the porepressure lives at the element corners, not
at the cell centres where the borehole points are, so the borehole sees an interpolated pressure and
the Peaceman well constant derived for cell-centred grids does not apply exactly.  The finite-element
results are therefore close but not identical.

!table id=borehole_mrst caption=Largest difference from MRST over the 20 time steps.  With rate control the rate is prescribed, so the bottom-hole pressure is compared; with bottom-hole-pressure control the rate is compared.
| Case | Control | Compared | FE | FV |
| --- | --- | --- | --- | --- |
| Single phase | rate | bottom-hole pressure (Pa) | $1.2 \times 10^{4}$ | 0.37 |
| Single phase | bhp | rate (kg/s) | 0.39 | $4 \times 10^{-6}$ |
| Gas injection | rate | bottom-hole pressure (Pa) | $2.0 \times 10^{4}$ | 0.21 |
| Gas injection | bhp | gas rate (kg/s) | 0.26 | $2 \times 10^{-6}$ |
| Water and gas production | rate | bottom-hole pressure (Pa) | - | 11 |
| Water and gas production | bhp | gas rate (kg/s) | 0.23 | $2 \times 10^{-4}$ |
| Heat production | rate | bottom-hole pressure (Pa) | $1.2 \times 10^{4}$ | 24 |
| Heat production | rate | heat rate (W, of $1.35 \times 10^{6}$) | 256 | 287 |

!syntax parameters /DiracKernels/PorousFlowPeacemanBorehole

!syntax inputs /DiracKernels/PorousFlowPeacemanBorehole

!syntax children /DiracKernels/PorousFlowPeacemanBorehole
