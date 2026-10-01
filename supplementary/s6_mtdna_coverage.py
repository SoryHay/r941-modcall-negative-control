#!/usr/bin/env python3
"""S6 — mtDNA read depth per arm (a) and position of 6mA calls at threshold 0.98, read interiors (b).

USAGE: s6_mtdna_coverage.py <NT.bam> <PA.bam> <NT_chrM_calls.tsv.gz> <PA_chrM_calls.tsv.gz> <out_prefix>
"""
import subprocess
import sys
import gzip
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BAM = {"NT": sys.argv[1], "PA": sys.argv[2]}
CALLS = {"NT": sys.argv[3], "PA": sys.argv[4]}
OUT = sys.argv[5]
CHRM = "CP068254.1"
LEN = 16569
THRESHOLD = 0.98
INTERIOR = 250


INK, INK2, INK3 = "#15191A", "#3E4749", "#6C787A"
NT_C, PA_C = "#1F6FB2", "#C2571A"
RULE, GRID, DLOOP = "#B9C2BE", "#E4E8E5", "#F1F3F0"


def depth(bam):
    d = np.zeros(LEN + 1)
    p = subprocess.run(["samtools", "depth", "-a", "-r", CHRM, bam],
                       capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        d[int(f[1])] = int(f[2])
    return d[1:]


def calls(path):
    """Per-position modified and total adenine observations, interior positions only."""
    mod = defaultdict(int)
    cov = defaultdict(int)
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["canonical_base"] != "A":
                continue
            fp, L = int(r["forward_read_position"]), int(r["read_length"])
            if min(fp, L - 1 - fp) < INTERIOR:
                continue
            pos = int(r["ref_position"])
            if pos < 0:
                continue
            cov[pos] += 1
            p = 0.0 if r["inferred"].lower() == "true" else float(r["mod_qual"])
            if p > THRESHOLD:
                mod[pos] += 1
    return mod, cov


dep = {a: depth(b) for a, b in BAM.items()}
cal = {a: calls(p) for a, p in CALLS.items()}

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.4), sharex=True,
                               gridspec_kw={"height_ratios": [1, 1], "hspace": 0.12})
fig.patch.set_facecolor("white")
x = np.arange(1, LEN + 1)

for ax in (ax1, ax2):
    ax.axvspan(16024, LEN, color=DLOOP, zorder=0)
    ax.axvspan(1, 576, color=DLOOP, zorder=0)
    ax.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(RULE)
    ax.tick_params(colors=INK3, labelsize=9, length=3)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_xlim(0, LEN)

# ---- Panel A — depth --------------------------------------------------------
ax1.fill_between(x, dep["NT"], color=NT_C, alpha=0.20, lw=0, zorder=2)
ax1.plot(x, dep["NT"], color=NT_C, lw=1.1, zorder=4, label="Untreated (NT)")
ax1.fill_between(x, dep["PA"], color=PA_C, alpha=0.20, lw=0, zorder=3)
ax1.plot(x, dep["PA"], color=PA_C, lw=1.1, zorder=4, label="Palmitate (PA)")
ax1.set_ylabel("Read depth", fontsize=10, color=INK)
ax1.set_ylim(0, max(dep["NT"].max(), dep["PA"].max()) * 1.32)       # headroom for the labels
ax1.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper center", ncol=2)
# PvuII linearises the mitochondrial circle at a single site; reads begin there,
# so 1–2,073 is the terminal segment of the linear molecule and is under-covered.
ax1.axvline(2074, color=INK3, lw=0.9, ls=(0, (3, 3)), zorder=5)
ax1.annotate("PvuII site (2,074)\nmolecule linearised here",
             xy=(2074, ax1.get_ylim()[1] * 0.90), xytext=(2400, ax1.get_ylim()[1] * 0.90),
             fontsize=8, color=INK3, va="center",
             arrowprops=dict(arrowstyle="-", color=INK3, lw=0.8))

# ---- Panel B — 6mA calls ----------------------------------------------------
for arm, colour, marker, off in (("NT", NT_C, "o", 0), ("PA", PA_C, "s", 0)):
    mod, cov = cal[arm]
    pos = sorted(p for p, m in mod.items() if m > 0)
    frac = [100 * mod[p] / cov[p] for p in pos]
    size = [18 + 6 * mod[p] for p in pos]
    ax2.scatter(pos, frac, s=size, facecolor=colour, edgecolor="white", lw=0.8,
                marker=marker, alpha=0.85, zorder=4,
                label=f"{arm} — {sum(mod.values())} calls at {len(pos)} positions")
ax2.set_ylabel(f"6mA at threshold {THRESHOLD}\n(% of reads at that position)",
               fontsize=10, color=INK)
ax2.set_xlabel("Mitochondrial genome position (bp)", fontsize=10, color=INK)
ax2.set_ylim(0, max(ax2.get_ylim()[1], 1) * 1.15)
ax2.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper center", ncol=2)

for ax in (ax1, ax2):
    ax.text(16296, ax.get_ylim()[1] * 0.97, "D-loop", fontsize=8, color=INK3,
            ha="center", rotation=90, va="top")



fig.text(0.012, 0.005,
         f"Positions ≥{INTERIOR} bp from the nearest read end. Marker area scales with the "
         f"number of modified reads at that position. Shaded: the D-loop control region, "
         f"which carries the heavy-strand replication origin and both promoters.",
         fontsize=7.6, color=INK3)

for dpi in (300, 600):
    fig.savefig(f"{OUT}_{dpi}dpi.png", dpi=dpi, bbox_inches="tight", facecolor="white")
fig.savefig(f"{OUT}.pdf", bbox_inches="tight", facecolor="white")

for arm in ("NT", "PA"):
    mod, cov = cal[arm]
    tot = sum(mod.values())
    npos = len([p for p, m in mod.items() if m > 0])
    d = dep[arm]
    print(f"{arm}: depth median {np.median(d):.0f}, min {d.min():.0f}, max {d.max():.0f}; "
          f"{tot} 6mA calls at {npos} distinct positions")
    dl = [p for p in mod if (p >= 16024 or p <= 576) and mod[p] > 0]
    print(f"    of which in the D-loop: {sum(mod[p] for p in dl)} calls at {len(dl)} positions "
          f"({100*(576+ (LEN-16024))/LEN:.1f}% of the genome)")
print("wrote", OUT + "_{300,600}dpi.png and .pdf")
