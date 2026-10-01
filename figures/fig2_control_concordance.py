#!/usr/bin/env python3
"""Fig. 2 — per-5-mer adenine FPR at 0.98 in two further amplified controls against the reference control."""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import (PARAMS, RC, INK, INK2, INK3, C1, C2, RULE, SURFACE, WIDTH_IN,
                     read_fpr_table, rate_or_bound, requirement, tidy, save)

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("ref", help="per-5-mer table, JAR cDNA run 1 (defines the allowlist)")
ap.add_argument("run2", help="per-5-mer table, JAR cDNA run 2")
ap.add_argument("jeg", help="per-5-mer table, JEG-3 cDNA")
ap.add_argument("out", help="output prefix")
args = ap.parse_args()

T = str(PARAMS["allowlist_threshold"])
FMAX = requirement("A")
plt.rcParams.update({**RC, "axes.titlesize": 8.2})


def arates(path):
    return {r["kmer"]: rate_or_bound(r, T) for r in read_fpr_table(path)
            if r["class"] == "A" and r[f"FPR_{T}"] != "NA"}


def spearman(a, b):
    ks = [k for k in a if k in b]
    ra = {k: i for i, k in enumerate(sorted(ks, key=lambda k: a[k]))}
    rb = {k: i for i, k in enumerate(sorted(ks, key=lambda k: b[k]))}
    n = len(ks)
    return 1 - 6 * sum((ra[k] - rb[k]) ** 2 for k in ks) / (n * (n * n - 1))


ref = arates(args.ref)
allow = {k for k, v in ref.items() if v <= FMAX}

fig, axes = plt.subplots(2, 1, figsize=(WIDTH_IN, 6.9), sharex=False, gridspec_kw={"hspace": 0.30})
fig.patch.set_facecolor(SURFACE)
panels = ((axes[0], args.run2, "a", "JAR cDNA run-2 — same line, different flow cell", "FPR, JAR cDNA run-2"),
          (axes[1], args.jeg, "b", "JEG-3 cDNA — different cell line", "FPR, JEG-3 cDNA"))

for ax, path, letter, subtitle, ylab in panels:
    other = arates(path)
    common = sorted(set(ref) & set(other))
    xs, ys = [ref[k] for k in common], [other[k] for k in common]
    lo, hi = min(min(xs), min(ys)) * 0.6, max(max(xs), max(ys)) * 1.6

    ax.set_facecolor(SURFACE)
    ax.plot([lo, hi], [lo, hi], color=RULE, lw=0.9, zorder=2)                 # identity, not a fit
    ax.axhline(FMAX, color=INK3, lw=0.7, ls=(0, (2, 2)), zorder=2)
    ax.axvline(FMAX, color=INK3, lw=0.7, ls=(0, (2, 2)), zorder=2)
    kept = [(x, y) for k, x, y in zip(common, xs, ys) if k in allow]
    drop = [(x, y) for k, x, y in zip(common, xs, ys) if k not in allow]
    ax.scatter([x for x, _ in kept], [y for _, y in kept], s=8, color=C1, lw=0, alpha=0.85, zorder=4,
               label=f"allowlist (n={len(kept)})")
    ax.scatter([x for x, _ in drop], [y for _, y in drop], s=15, marker="s", color=C2, lw=0, zorder=5,
               label=f"excluded (n={len(drop)})")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_title(letter, loc="left", color=INK, pad=13, fontweight="bold")
    ax.text(0.0, 1.015, subtitle, transform=ax.transAxes, ha="left", va="bottom", color=INK2, fontsize=7.4)
    ax.set_ylabel(ylab, color=INK)
    ax.set_xlabel("FPR, JAR cDNA run-1", color=INK)
    rho = spearman(ref, other)
    ax.text(0.035, 0.955, f"Spearman $\\rho$ = {rho:.3f}", transform=ax.transAxes,
            ha="left", va="top", color=INK2, fontsize=7.2)
    ax.legend(frameon=False, labelcolor=INK2, loc="lower right", handletextpad=0.3, fontsize=7,
              borderaxespad=0.4)
    tidy(ax, grid_axis="both")
    leak = sum(1 for k in common if k in allow and other[k] > FMAX)
    print(f"{letter}: rho={rho:.3f}  allowlisted contexts above requirement here: {leak}/{len(kept)}")

def tex_sci(x, digits):
    m, e = f"{x:.{digits}e}".split("e")
    return f"{m}\\times10^{{{int(e)}}}"


r = PARAMS["requirement"]
fig.text(0.5, 0.005,
         r"$\mathrm{FPR}_{\max}=\dfrac{\alpha\,p\,S}{(1-\alpha)(1-p)}=" + tex_sci(FMAX, 2) + "$"
         "\n"
         rf"at $\alpha={r['alpha']}$,  $p={tex_sci(r['prevalence']['A'], 1)}$,  $S={r['sensitivity']:g}$",
         ha="center", va="top", color=INK2, fontsize=8, linespacing=1.6)

save(fig, args.out)
print(f"{args.out}: allowlist n={len(allow)}, excluded n={len(ref) - len(allow)}")
