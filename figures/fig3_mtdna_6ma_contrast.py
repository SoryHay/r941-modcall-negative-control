#!/usr/bin/env python3
"""Fig. 3 — mitochondrial 6mA call rate per arm across the threshold sweep (a) and the PA/NT ratio (b)."""
import argparse
import csv
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import (PARAMS, RC, INK, INK2, INK3, C1, C2, BAND, SURFACE, WIDTH_IN,
                     read_fpr_table, rate_or_bound, requirement, tidy, save)

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("contrast", help="per-threshold 6mA counts, both arms (data/admsc_chrM_6mA_counts.tsv)")
ap.add_argument("background", help="control-arm chrM 6mA rate (data/control_chrM_6mA_background.tsv)")
ap.add_argument("ref", help="per-5-mer table of the reference control, for the qualifying thresholds")
ap.add_argument("out", help="output prefix")
args = ap.parse_args()

REQ = requirement("A")
plt.rcParams.update(RC)


def first_qualifying(rows, keep):
    """Lowest threshold at which the pooled FPR over the kept contexts meets the requirement."""
    for t in PARAMS["thresholds"]:
        fp = sum(int(r[f"fp_{t}"]) for r in rows if keep(r))
        n = sum(int(r["n_obs"]) for r in rows if keep(r))
        if fp / n <= REQ:
            return t


ref_rows = [r for r in read_fpr_table(args.ref) if r["class"] == "A"]
t_allow = str(PARAMS["allowlist_threshold"])
allow = {r["kmer"] for r in ref_rows if rate_or_bound(r, t_allow) <= REQ}
QUALIFY_ALLOWLIST = first_qualifying(ref_rows, lambda r: r["kmer"] in allow)
QUALIFY_ALL = first_qualifying(ref_rows, lambda r: True)

rows = sorted(csv.DictReader(open(args.contrast), delimiter="\t"), key=lambda r: float(r["threshold"]))
thr = [float(r["threshold"]) for r in rows]
x = list(range(len(thr)))
nt = [100 * int(r["nt_mod"]) / int(r["nt_n"]) for r in rows]
pa = [100 * int(r["pa_mod"]) / int(r["pa_n"]) for r in rows]
ratio = [p / n for p, n in zip(pa, nt)]
bgd = {float(r["threshold"]): float(r["fpr"]) for r in csv.DictReader(open(args.background), delimiter="\t")}
bg = [100 * bgd[t] for t in thr]

lo_r, hi_r = [], []                           # 95% CI on the rate ratio from the Poisson counts
for r, rr in zip(rows, ratio):
    se = (1.0 / int(r["nt_mod"]) + 1.0 / int(r["pa_mod"])) ** 0.5
    lo_r.append(rr * math.exp(-1.96 * se))
    hi_r.append(rr * math.exp(1.96 * se))

qstart = x[thr.index(QUALIFY_ALLOWLIST)] - 0.5

fig, axes = plt.subplots(2, 1, figsize=(WIDTH_IN, 6.9), sharex=True,
                         gridspec_kw={"hspace": 0.11, "height_ratios": [1.25, 1]})
fig.patch.set_facecolor(SURFACE)
for ax in axes:
    ax.set_facecolor(SURFACE)
    ax.axvspan(qstart, x[-1] + 0.5, color=BAND, zorder=1)       # thresholds where calls are admissible

ax = axes[0]
ax.plot(x, bg, color=INK3, lw=1.1, ls=(0, (4, 3)), zorder=3, label="caller background")
ax.plot(x, nt, color=C1, lw=1.4, marker="o", ms=4.5, zorder=5, label="untreated (NT)")
ax.plot(x, pa, color=C2, lw=1.4, marker="s", ms=4.2, zorder=5, label="palmitate (PA)")
ax.axhline(100 * REQ, color=INK3, lw=0.8, ls=(0, (1.5, 2)), zorder=2)
ax.set_yscale("log")
ax.set_ylabel("6mA calls, % of adenines", color=INK)
ax.text(x[0] - 0.06, nt[0] * 1.55, f"{nt[0]:.2f}%", color=C1, fontsize=6.9, ha="left")
ax.text(x[0] - 0.06, pa[0] * 0.55, f"{pa[0]:.2f}%", color=C2, fontsize=6.9, ha="left")
ax.text(x[0] + 0.08, 100 * REQ * 1.25, "required specificity for per-site 6mA",
        color=INK3, fontsize=6.6, ha="left", va="bottom")
for t in (QUALIFY_ALLOWLIST, QUALIFY_ALL):
    ax.axvline(x[thr.index(t)] - 0.5, color=INK3, lw=0.6, ls=(0, (1, 2.5)), zorder=2)
ax.legend(frameon=False, labelcolor=INK2, loc="upper right", fontsize=6.8,
          handletextpad=0.5, borderaxespad=0.3, labelspacing=0.35)
ax.set_title("a", loc="left", color=INK, pad=5, fontweight="bold")
tidy(ax)
for t, lab in ((QUALIFY_ALLOWLIST, "allowlist"), (QUALIFY_ALL, "all 256")):
    ax.text(x[thr.index(t)] - 0.58, ax.get_ylim()[0] * 1.35, lab, rotation=90,
            color=INK3, fontsize=6.5, ha="center", va="bottom")

ax = axes[1]
ax.axhline(1.0, color=INK3, lw=0.8, zorder=2)
ax.text(x[0], 1.02, "no effect", color=INK3, fontsize=6.6, ha="left", va="bottom")
ax.vlines(x, lo_r, hi_r, color=INK3, lw=0.9, zorder=4)
ax.plot(x, ratio, color=C2, lw=1.4, marker="D", ms=4.2, zorder=5)
for xi, rr in zip(x, ratio):
    ax.text(xi, min(lo_r) - 0.10, f"{rr:.2f}", color=INK2, fontsize=6.9, ha="center", va="top",
            fontweight="bold" if thr[xi] >= QUALIFY_ALLOWLIST else "normal")
ax.set_ylim(min(lo_r) - 0.22, max(1.28, max(hi_r) * 1.05))
ax.set_ylabel("palmitate / untreated", color=INK)
ax.set_xlabel("modified-base probability threshold", color=INK)
ax.set_xticks(x); ax.set_xticklabels([f"{t:g}" for t in thr])
ax.set_xlim(-0.5, x[-1] + 0.5)
ax.set_title("b", loc="left", color=INK, pad=5, fontweight="bold")
tidy(ax)

save(fig, args.out)
print(f"{args.out}: qualifying from {QUALIFY_ALLOWLIST} (allowlist, n={len(allow)}) and {QUALIFY_ALL} (all contexts)")
print("PA/NT ratios:", " ".join(f"{r:.3f}" for r in ratio))
