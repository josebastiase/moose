#!/usr/bin/env python3
"""
Plot PorousFlowPeacemanBorehole with finite-element (FE) and finite-volume (FV) variables against
MRST 2023b for the cases of compare_vs_mrst.py (run that first: this reads compare_runs/).

usage: python3 plot_vs_mrst.py [output.png]
Default output: ../../doc/content/media/porous_flow/borehole_vs_mrst.png
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "compare_runs")
REF = os.path.join(HERE, "mrst_reference")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    HERE, "..", "..", "doc", "content", "media", "porous_flow", "borehole_vs_mrst.png")

# validated categorical palette (light surface)
FE, FV, MRST = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#2b2b29", "#6b6b66", "#e4e3df"

# title, FE run, FV run, MRST file, (MOOSE column, MRST column), y label, scale
PANELS = [
    ("1 phase, rate control", "FE_re_0_28_1_phase_rate", "FV_1_phase_rate", "bh_rate_mrst_n1",
     ("p_bot", "p_bot"), "bottom-hole pressure (MPa)", 1e-6),
    ("1 phase, pressure control", "FE_re_0_28_1_phase_bhp", "FV_1_phase_bhp", "bh_rate_mrst_bhp_n1",
     ("bh_rate", "rate"), "production rate (kg/s)", 1),
    ("gas injection, rate control", "FE_re_0_28_gas_injection_rate", "FV_gas_injection_rate",
     "bh_rate_mrst_2ph_inject_n1", ("p_bot", "p_bot"), "bottom-hole pressure (MPa)", 1e-6),
    ("water + gas production, pressure control", "FE_re_0_28_water_gas_production_bhp",
     "FV_water_gas_production_bhp", "bh_rate_mrst_2ph_produce_bhp_n1", ("gas_rate", "gas_rate"),
     "gas production rate (kg/s)", 1),
    ("heat production at 1 kg/s: pressure", "FE_re_0_28_heat_production_rate", "FV_heat_production_rate",
     "bh_heat_mrst_n1", ("p_bot", "p_bot"), "bottom-hole pressure (MPa)", 1e-6),
    ("heat production at 1 kg/s: heat rate", "FE_re_0_28_heat_production_rate", "FV_heat_production_rate",
     "bh_heat_mrst_n1", ("energy_rate", "energy_rate"), "heat production rate (MW)", 1e-6),
]


def series(path, col, scale):
    with open(path) as f:
        rows = [r for r in csv.DictReader(f) if float(r["time"]) > 0]
    return [float(r["time"]) for r in rows], [float(r[col]) * scale for r in rows]


def main():
    plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": GRID})
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), facecolor=SURFACE)
    for ax, (title, fe, fv, ref, (mc, sc), ylab, scale) in zip(axes.flat, PANELS):
        ax.set_facecolor(SURFACE)
        ax.grid(color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        t, y = series(os.path.join(RUNS, fe + ".csv"), mc, scale)
        ax.plot(t, y, color=FE, lw=2, label="MOOSE FE")
        t, y = series(os.path.join(RUNS, fv + ".csv"), mc, scale)
        ax.plot(t, y, color=FV, lw=2, label="MOOSE FV")
        t, y = series(os.path.join(REF, ref + ".csv"), sc, scale)
        ax.plot(t, y, ls="none", marker="o", ms=6, mfc=SURFACE, mec=MRST, mew=1.6, label="MRST")
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.set_ylabel(ylab)
        ax.set_xlabel("time (s)")
        ax.ticklabel_format(axis="y", useOffset=False)
        if mc == "energy_rate":
            # the three agree to 0.03%: a zoomed axis would exaggerate that
            ax.set_ylim(1.30, 1.40)
            ax.text(0.97, 0.06, "all within 0.03%", transform=ax.transAxes, ha="right", color=MUTED)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
    fig.savefig(OUT, dpi=150, facecolor=SURFACE)
    print("wrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
