% MRST (geothermal module) counterpart of the MOOSE heat cases (heat_fe.i and heat_fv.i): water
% with a vertical temperature gradient is produced at 1 kg/s through a vertical well with three
% completions, and the produced energy rate depends on how the flow is split between completions
% at different temperatures.  Run from this directory:
%   n = 1; bh_heat_mrst
% Writes ../mrst_reference/bh_heat_mrst_n<n>.csv with columns time, p_bot, rate (kg/s, positive =
% production), energy_rate (W, positive = production), T_top, T_bot (K, the cells of the top and
% bottom completions).
%
% Mapping to MOOSE (SimpleFluidProperties, PorousFlowMatrixInternalEnergy,
% PorousFlowThermalConductivityIdeal):
%  - density 1000*exp(p/2E9 - 2E-4*(T - 293.15)), ie density0 = 1000*exp(2E-4*293.15) and
%    thermal_expansion = 2E-4 in MOOSE.  rho = 1 and b = 1000*exp(p/2E9) make the surface volume
%    a mass, so the 'rate' control is a mass rate.
%  - internal energy 4200*T and enthalpy u + p/rho (cv = 4200 in MOOSE), viscosity 1E-3 Pa s
%  - rock density 2700, heat capacity 1000, thermal conductivity 2 (rock) and 0.6 (water), ie a
%    bulk conductivity 2*(1 - 0.1) + 0.6*0.1 = 1.86 in MOOSE
%  - the well exchanges no heat by conduction with the reservoir (WIth = 0): its heat flux is the
%    enthalpy carried by the produced water, as the MOOSE borehole with use_enthalpy = true
%  - the wellbore head is fixed at 1E4*dz, as MOOSE's unit_weight = (0, 0, -1E4)

if ~exist('n', 'var'), n = 1; end
mrstModule add ad-core ad-blackoil ad-props compositional geothermal

%% Grid: x, y in [-10, 10], MOOSE z in [-15, 15] (depth in [-15, 15])
G = cartGrid([4 4 6] * n, [20 20 30]);
G.nodes.coords = bsxfun(@minus, G.nodes.coords, [10 10 15]);
G = computeGeometry(G);
dx = 5 / n; dy = 5 / n; dz = 5 / n;

%% Rock, fluid
k = 1e-12;
rock = makeRock(G, k, 0.1);
rock = addThermalRockProps(rock, 'lambdaR', 2, 'rhoR', 2700, 'CpR', 1000);
gravity reset on
gravity([0 0 10]);
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1e-3, 'rho', 1, 'b', 1000, 'c', 1/2e9, 'pRef', 0);
fluid = addThermalFluidProps(fluid, 'Cp', 4200, 'lambdaF', 0.6, 'useEOS', false, ...
                             'cT', 2e-4, 'TRef', 293.15);
model = GeothermalModel(G, rock, fluid);
model.toleranceCNV = 1e-10;
model.toleranceMB = 1e-12;

%% Well: three completions at MOOSE (2.5, 2.5, z), z = 7.5, 2.5, -2.5, producing 1 kg/s
pts = [2.5 2.5 -7.5; 2.5 2.5 -2.5; 2.5 2.5 2.5];   % MRST coordinates (depth = -z)
ijk = floor(bsxfun(@rdivide, bsxfun(@plus, pts, [10 10 15]), [dx dy dz])) + 1;
cells = sub2ind(G.cartDims, ijk(:, 1), ijk(:, 2), ijk(:, 3));
assert(max(abs(G.cells.centroids(cells, :) - pts), [], 'all') < 1e-10);
rw = 0.1;
L = [2.5; 5; 2.5];
r0 = 0.28 * sqrt(dx^2 + dy^2) / 2;
WI = 2 * pi * k * L / log(r0 / rw);
W = addWell([], G, rock, cells, 'Type', 'rate', 'Val', -1, 'Radius', rw, 'WI', WI, ...
            'Comp_i', 1, 'refDepth', pts(end, 3), 'Name', 'bh');
W = addThermalWellProps(W, G, rock, fluid, 'T', 293.15, 'WIth', 0);

dt = 5 * ones(20, 1);
tEnd = cumsum(dt);
schedule = simpleSchedule(dt, 'W', W);

%% Initial state: hydrostatic pressure and T from 293.15 K at the top to 353.15 K at the bottom
depth = G.cells.centroids(:, 3);
state0 = initResSol(G, 1e7 + 1e4 * depth, 1);
state0.T = 293.15 + 2 * (15 + depth);

%% Wellbore head fixed at 1E4*dz
model = model.validateModel(struct('W', W));
model.FacilityModel.WellModels{1}.doUpdatePressureDrop = false;
state0.wellSol = initWellSolAD(W, model, state0);
state0.wellSol.cdp = 1e4 * W.dZ;

[wellSols, states] = simulateScheduleAD(state0, model, schedule);
dcdp = cellfun(@(ws) max(abs(ws.cdp - 1e4 * W.dZ)), wellSols);
assert(all(dcdp == 0), 'the wellbore head was not kept fixed');

%% Output
p_bot = cellfun(@(ws) ws.bhp, wellSols);
rate = -cellfun(@(ws) ws.qWs, wellSols);
energy_rate = -cellfun(@(ws) ws.effect, wellSols);
T_top = cellfun(@(s) s.T(cells(1)), states);
T_bot = cellfun(@(s) s.T(cells(end)), states);
T = table(tEnd, p_bot, rate, energy_rate, T_top, T_bot, ...
          'VariableNames', {'time', 'p_bot', 'rate', 'energy_rate', 'T_top', 'T_bot'});
writetable(T, sprintf('../mrst_reference/bh_heat_mrst_n%d.csv', n));
disp(T)
