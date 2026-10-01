"""Shared parameters, per-context table reader and figure style for the article's figures."""
import csv
import math
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
PARAMS = yaml.safe_load(open(ROOT / "config" / "params.yaml"))

INK, INK2, INK3 = "#15191A", "#3E4749", "#6C787A"
C1, C2 = "#1F6FB2", "#C2571A"
GRID, SURFACE, BAND, RULE = "#E4E8E5", "#FCFCFB", "#DCE3DE", "#B9C2BE"
WIDTH_IN = 135 / 25.4                      # journal text block, 135 mm

RC = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.labelsize": 7.2, "ytick.labelsize": 7.2,
    "axes.titlesize": 8.4, "axes.labelsize": 8,
}


def round_sig(x, n):
    return float(f"{x:.{n - 1}e}")


def requirement(mark, params=PARAMS):
    """Equation (1): the largest false-positive rate compatible with the target FDR."""
    r = params["requirement"]
    a, s, p = r["alpha"], r["sensitivity"], r["prevalence"][mark]
    return round_sig(a * p * s / ((1 - a) * (1 - p)), r["sig_digits"])


def read_fpr_table(path):
    """Rows of a per-5-mer FPR table (workflow/fpr_per_kmer_interior.awk output)."""
    return list(csv.DictReader(open(path), delimiter="\t"))


def rate_or_bound(row, t):
    """FPR at threshold t; a context with no false positive takes the 1/n bound."""
    f, n = float(row[f"FPR_{t}"]), int(row["n_obs"])
    return f if f > 0 else 1.0 / n


def q(f):
    return -10 * math.log10(f)


def tidy(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK3)
    ax.tick_params(colors=INK3, length=2.5, pad=2)
    ax.grid(True, axis=grid_axis, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def save(fig, out):
    for ext, dpi in ((".pdf", 600), ("_600dpi.png", 600), ("_300dpi.png", 300)):
        fig.savefig(f"{out}{ext}", dpi=dpi, bbox_inches="tight", facecolor=SURFACE)
