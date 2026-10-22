#!/usr/bin/env python3
"""
Compare the PorousFlow well/borehole objects with MRST 2023b on the MRST cases of idaholab/moose#34110
(3D box of 5 m cells, vertical well with three completions at cell centres, 20 steps of 5 s).

  FE   : PorousFlowPeacemanBorehole with finite-element variables (the inputs in this folder, by
         Chris Green), with re_constant 0.28 (Peaceman) and 0.6 (fitted)
  FV   : the same PorousFlowPeacemanBorehole with finite-volume variables
         (test/tests/dirackernels/fv_borehole)
  FV [Wells] : the [Wells] block in the finite-volume model (fv_checks/)

The MRST results are in mrst_reference/.  For rate control the bottom-hole pressure is compared, for
bottom-hole-pressure control the rates; the gas saturation in the cells of the top and bottom
completions is compared where both codes report it (for FE it is interpolated from the element
corners, so it is not exactly a cell value).  The heat case (heat_fe.i, heat_fv.i, MRST geothermal
module) compares the bottom-hole pressure, the produced heat rate and the completion temperatures.

usage: python3 compare_vs_mrst.py [path/to/porous_flow-opt]
Writes compare_vs_mrst.md.
"""
import csv
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PF = os.path.normpath(os.path.join(HERE, "..", ".."))
EXE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PF, "porous_flow-opt")
OUT = os.path.join(HERE, "compare_runs")
FVT = os.path.join(PF, "test", "tests", "dirackernels", "fv_borehole")
REF = os.path.join(HERE, "mrst_reference")


# label, case, input dir, input, cli, MRST file, compared columns (MOOSE name, MRST name)
P_BOT = [("p_bot", "p_bot")]
RATE1 = [("bh_rate", "rate")]
GAS = [("gas_rate", "gas_rate")]
WG = [("water_rate", "water_rate"), ("gas_rate", "gas_rate")]
SG = [("sgas_top", "sgas_top"), ("sgas_bot", "sgas_bot")]
HEAT = [("p_bot", "p_bot"), ("energy_rate", "energy_rate"), ("T_top", "T_top"), ("T_bot", "T_bot")]
CASES = []
for re_c in ("0.28", "0.6"):
    cli = [f"GlobalParams/re_constant={re_c}"]
    CASES += [
        (f"FE re={re_c}", "1 phase, rate", HERE, "1ph_rate.i", cli, "bh_rate_mrst_n1.csv", P_BOT),
        (f"FE re={re_c}", "1 phase, bhp", HERE, "1ph_bhp.i", cli, "bh_rate_mrst_bhp_n1.csv", RATE1),
        (f"FE re={re_c}", "gas injection, rate", HERE, "2ph_rate.i", cli, "bh_rate_mrst_2ph_inject_n1.csv", P_BOT + SG),
        (f"FE re={re_c}", "gas injection, bhp", HERE, "2ph_bhp.i", cli, "bh_rate_mrst_2ph_inject_bhp_n1.csv", GAS + SG),
        (f"FE re={re_c}", "water+gas production, bhp", HERE, "2ph_bhp_produce.i", cli, "bh_rate_mrst_2ph_produce_bhp_n1.csv", WG + SG),
        (f"FE re={re_c}", "heat production, rate", HERE, "heat_fe.i", cli, "bh_heat_mrst_n1.csv", HEAT),
    ]
CASES += [
    ("FV", "1 phase, rate", FVT, "fv_borehole_rate.i", [], "bh_rate_mrst_n1.csv", P_BOT),
    ("FV", "1 phase, bhp", FVT, "fv_borehole_bhp.i", [], "bh_rate_mrst_bhp_n1.csv", RATE1),
    ("FV", "gas injection, rate", FVT, "fv_borehole_2ph_rate.i", [], "bh_rate_mrst_2ph_inject_n1.csv", P_BOT + SG),
    ("FV", "gas injection, bhp", FVT, "fv_borehole_2ph_bhp.i", [], "bh_rate_mrst_2ph_inject_bhp_n1.csv", GAS + SG),
    ("FV", "water+gas production, rate", FVT, "fv_borehole_2ph_rate_produce.i", [], "bh_rate_mrst_2ph_produce_n1.csv", P_BOT + WG + SG),
    ("FV", "water+gas production, bhp", FVT, "fv_borehole_2ph_bhp_produce.i", [], "bh_rate_mrst_2ph_produce_bhp_n1.csv", WG + SG),
    ("FV", "heat production, rate", HERE, "heat_fv.i", [], "bh_heat_mrst_n1.csv", HEAT),
]
FVC = os.path.join(HERE, "fv_checks")
CASES += [
    ("FV [Wells]", "gas injection, injector at a prescribed rate", FVC, "fv_wells_inject.i", [], "bh_rate_mrst_2ph_inject_n1.csv", GAS + SG),
    ("FV [Wells]", "gas injection, injector with bhp_max (not reached)", FVC, "fv_wells_bhpmax.i", [], "bh_rate_mrst_2ph_inject_n1.csv", [("inj_bhp", "p_bot")] + GAS + SG),
    ("FV [Wells]", "production well", FVC, "fv_wells_production.i", [], None, []),
]


def read(path):
    with open(path) as f:
        return {round(float(r["time"]), 6): r for r in csv.DictReader(f) if float(r["time"]) > 0}


def first_error(text):
    m = re.search(r"\*\*\* ERROR \*\*\*\s*\n(?:.*\n){0,2}?(.*\S.*)\n(.*)", text)
    if not m:
        return "failed (see log)"
    msg = " ".join(x.strip() for x in m.groups() if x.strip())
    return re.sub(r"\x1b\[[0-9;]*m", "", msg)[:220]


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for impl, case, d, inp, cli, ref, cols in CASES:
        tag = re.sub(r"\W+", "_", f"{impl}_{case}")
        base = os.path.join(OUT, tag)
        p = subprocess.run([EXE, "-i", inp, f"Outputs/file_base={base}", "--color", "off"] + cli,
                           cwd=d, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
        text = p.stdout + p.stderr
        open(base + ".log", "w").write(text)
        failed = text.count("Solve Did NOT Converge!")
        if p.returncode != 0:
            rows.append((impl, case, "runs: **no**", first_error(text)))
            print(f"{impl:9s} {case:55s} ERROR", flush=True)
            continue
        if not ref:
            rows.append((impl, case, "runs: yes", "(no MRST counterpart)"))
            continue
        m, s = read(base + ".csv"), read(os.path.join(REF, ref))
        parts = []
        for mc, sc in cols:
            if mc not in next(iter(m.values())):
                continue
            worst = max(abs(float(m[t][mc]) - float(s[t][sc])) for t in s if t in m)
            unit = (" Pa" if mc in ("p_bot", "inj_bhp") else "" if mc.startswith("sgas") else
                    " W" if mc == "energy_rate" else " K" if mc.startswith("T_") else " kg/s")
            parts.append(f"{mc} {worst:.3g}{unit}")
        note = f"; {failed} failed solves" if failed else ""
        rows.append((impl, case, "runs: yes", "max |MOOSE - MRST|: " + ", ".join(parts) + note))
        print(f"{impl:9s} {case:55s} {', '.join(parts)}{note}", flush=True)

    lines = ["| implementation | case | runs? | comparison with MRST |", "|---|---|---|---|"]
    lines += [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in rows]
    table = "\n".join(lines)
    open(os.path.join(HERE, "compare_vs_mrst.md"), "w").write(table + "\n")
    print("\n" + table)


if __name__ == "__main__":
    main()
