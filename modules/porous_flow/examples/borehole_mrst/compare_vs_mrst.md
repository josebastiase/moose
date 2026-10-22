| implementation | case | runs? | comparison with MRST |
|---|---|---|---|
| FE re=0.28 | 1 phase, rate | runs: yes | max |MOOSE - MRST|: p_bot 1.22e+04 Pa |
| FE re=0.28 | 1 phase, bhp | runs: yes | max |MOOSE - MRST|: bh_rate 0.388 kg/s |
| FE re=0.28 | gas injection, rate | runs: yes | max |MOOSE - MRST|: p_bot 2.03e+04 Pa, sgas_top 0.00277, sgas_bot 0.00288 |
| FE re=0.28 | gas injection, bhp | runs: yes | max |MOOSE - MRST|: gas_rate 0.263 kg/s, sgas_top 0.00117, sgas_bot 0.00121 |
| FE re=0.28 | water+gas production, bhp | runs: yes | max |MOOSE - MRST|: water_rate 0.107 kg/s, gas_rate 0.234 kg/s, sgas_top 1.87e-05, sgas_bot 1.73e-05 |
| FE re=0.28 | heat production, rate | runs: yes | max |MOOSE - MRST|: p_bot 1.23e+04 Pa, energy_rate 256 W, T_top 0.000585 K, T_bot 0.000501 K |
| FE re=0.6 | 1 phase, rate | runs: yes | max |MOOSE - MRST|: p_bot 1.12e+03 Pa |
| FE re=0.6 | 1 phase, bhp | runs: yes | max |MOOSE - MRST|: bh_rate 0.0171 kg/s |
| FE re=0.6 | gas injection, rate | runs: yes | max |MOOSE - MRST|: p_bot 1.2e+03 Pa, sgas_top 0.00277, sgas_bot 0.00289 |
| FE re=0.6 | gas injection, bhp | runs: yes | max |MOOSE - MRST|: gas_rate 0.0136 kg/s, sgas_top 0.00118, sgas_bot 0.00123 |
| FE re=0.6 | water+gas production, bhp | runs: yes | max |MOOSE - MRST|: water_rate 0.0485 kg/s, gas_rate 0.106 kg/s, sgas_top 2.63e-05, sgas_bot 2.16e-05 |
| FE re=0.6 | heat production, rate | runs: yes | max |MOOSE - MRST|: p_bot 672 Pa, energy_rate 309 W, T_top 0.000586 K, T_bot 0.0005 K |
| FV | 1 phase, rate | runs: yes | max |MOOSE - MRST|: p_bot 0.367 Pa |
| FV | 1 phase, bhp | runs: yes | max |MOOSE - MRST|: bh_rate 4.04e-06 kg/s |
| FV | gas injection, rate | runs: yes | max |MOOSE - MRST|: p_bot 0.205 Pa, sgas_top 2.54e-08, sgas_bot 2.22e-08 |
| FV | gas injection, bhp | runs: yes | max |MOOSE - MRST|: gas_rate 1.52e-06 kg/s, sgas_top 2.85e-08, sgas_bot 2.01e-08 |
| FV | water+gas production, rate | runs: yes | max |MOOSE - MRST|: p_bot 11.1 Pa, water_rate 9.41e-08 kg/s, gas_rate 9.4e-08 kg/s, sgas_top 2.52e-08, sgas_bot 1.65e-08 |
| FV | water+gas production, bhp | runs: yes | max |MOOSE - MRST|: water_rate 8.46e-05 kg/s, gas_rate 0.000186 kg/s, sgas_top 4.55e-08, sgas_bot 3.78e-08 |
| FV | heat production, rate | runs: yes | max |MOOSE - MRST|: p_bot 23.5 Pa, energy_rate 287 W, T_top 1.55e-06 K, T_bot 1.48e-06 K |
| FV [Wells] | gas injection, injector at a prescribed rate | runs: yes | max |MOOSE - MRST|: gas_rate 9.99e-15 kg/s, sgas_top 0.000292, sgas_bot 0.000436 |
| FV [Wells] | gas injection, injector with bhp_max (not reached) | runs: yes | max |MOOSE - MRST|: inj_bhp 0.205 Pa, gas_rate 0 kg/s, sgas_top 2.54e-08, sgas_bot 2.22e-08 |
| FV [Wells] | production well | runs: **no** | /home/jose/projects/moose-wells/modules/porous_flow/examples/borehole_mrst/fv_checks/fv_wells_production.i:75.3: The following occurred in the Action 'prod' of type PorousFlowWellAction. |
