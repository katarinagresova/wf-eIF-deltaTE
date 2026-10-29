"""The Ribo-seq QC panel of one experiment: three plots in a row, from ribokit's
<sample>.human.qc.tsv (length, region, phase, reads), a line or dot per library,
coloured by condition or by region / frame:

  1. footprint length distribution: % of a library's reads per read length
  2. P-site region (leader, CDS, trailer): % of a library's reads of that read
     length, per read length
  3. 5' end phase in the CDS (0, 1, 2): % of a library's reads of that read
     length, per read length

The read lengths shown are the window [lo, hi] (config read_lengths), the one
ribokit counts; reads outside it (and without a CDS) have no region.

Usage: plot_rpf_qc.py <out.pdf> <experiment> <lo> <hi> --lib <sample> <condition> <qc.tsv> ...
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

CONDITION_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
REGION_COLORS = {"leader": "#e08a2c", "CDS": "#2a6fb0", "trailer": "#8a5fa8"}
PHASE_COLORS = {0: "#deec3b", 1: "#44eaa8", 2: "#7363cc"}
MUTED, GRID = "#52514e", "#e6e6e3"


def read_libraries(libs, lo, hi):
    """All libraries' qc tables in one, within the read length window."""
    tables = []
    for sample, condition, path in libs:
        table = pd.read_csv(path, sep="\t")
        table["sample"], table["condition"] = sample, condition
        tables.append(table)
    qc = pd.concat(tables, ignore_index=True)
    return qc[qc["length"].between(lo, hi)]


def percent_per_length(qc, column):
    """% of a library's reads of one read length per value of `column` (NA: not counted)."""
    counts = qc.dropna(subset=[column]).groupby(["sample", "length", column])["reads"].sum().reset_index()
    counts["percent"] = 100 * counts["reads"] / counts.groupby(["sample", "length"])["reads"].transform("sum")
    return counts


def style(ax, lo, hi):
    ax.set_xlim(lo - 0.5, hi + 0.5)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("read length (nt)", color=MUTED, fontsize=10)
    ax.set_ylabel("% reads of that length", color=MUTED, fontsize=10)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def plot_lengths(ax, qc, colors, lo, hi):
    for condition, color in colors.items():
        reads = qc[qc["condition"] == condition]
        for j, (_, library) in enumerate(reads.groupby("sample")):
            per_length = library.groupby("length")["reads"].sum()
            ax.plot(per_length.index, 100 * per_length / per_length.sum(), color=color, linewidth=1, alpha=0.8,
                    label=condition if j == 0 else None)
    style(ax, lo, hi)
    ax.legend(frameon=False, fontsize=8)


def plot_regions(ax, qc, lo, hi):
    regions = percent_per_length(qc, "region")
    rng = np.random.default_rng(0)
    for region, color in REGION_COLORS.items():
        dots = regions[regions["region"] == region]
        ax.scatter(dots["length"] + rng.uniform(-0.15, 0.15, len(dots)), dots["percent"], s=12, color=color, alpha=0.8,
                   linewidths=0, label=region)
    ax.set_ylim(0, 100)
    style(ax, lo, hi)
    ax.legend(frameon=False, fontsize=8)


def plot_phases(ax, qc, lo, hi):
    phases = percent_per_length(qc, "phase")
    phases["phase"] = phases["phase"].astype(int)
    rng = np.random.default_rng(0)
    for phase, color in PHASE_COLORS.items():
        dots = phases[phases["phase"] == phase]
        # one column of dots per phase beside each read length, jittered a little
        ax.scatter(dots["length"] + (phase - 1) * 0.25 + rng.uniform(-0.06, 0.06, len(dots)), dots["percent"], s=12,
                   color=color, alpha=0.8, linewidths=0, label=f"phase {phase}")
    ax.set_ylim(0, 100)
    style(ax, lo, hi)
    ax.legend(frameon=False, fontsize=8)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out")
    parser.add_argument("experiment")
    parser.add_argument("lo", type=int)
    parser.add_argument("hi", type=int)
    parser.add_argument("--lib", nargs=3, action="append", required=True, metavar=("SAMPLE", "CONDITION", "QC_TSV"))
    args = parser.parse_args()

    qc = read_libraries(args.lib, args.lo, args.hi)
    conditions = sorted(qc["condition"].unique())
    colors = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)] for i, c in enumerate(conditions)}

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    plot_lengths(axes[0], qc, colors, args.lo, args.hi)
    plot_regions(axes[1], qc, args.lo, args.hi)
    plot_phases(axes[2], qc, args.lo, args.hi)
    fig.suptitle(args.experiment, fontsize=13)
    fig.savefig(args.out)


if __name__ == "__main__":
    main()
