#!/usr/bin/env python3
"""Fig. 1 — per-5-mer specificity (Q) on the amplified control at threshold 0.98, contexts ranked."""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import (PARAMS, RC, INK, INK2, C1, C2, SURFACE, WIDTH_IN,
                     read_fpr_table, rate_or_bound, requirement, q, tidy, save)

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("table", help="per-5-mer FPR table of the reference control (JAR cDNA run 1)")
ap.add_argument("out", help="output prefix")
args = ap.parse_args()

T_MAIN, T_FAINT = "0.98", "0.7"
REQ = {m: requirement(m) for m in ("CpG", "CpH", "A")}
plt.rcParams.update(RC)


def load(path):
    """class -> [(Q at 0.98, Q at 0.7, no FP observed)], ranked by Q at 0.98."""
    d = {"CpG": [], "CpH": [], "A": []}
    for r in read_fpr_table(path):
        if r["class"] not in d or r[f"FPR_{T_MAIN}"] == "NA":
            continue
        d[r["class"]].append((q(rate_or_bound(r, T_MAIN)), q(rate_or_bound(r, T_FAINT)),
                              float(r[f"FPR_{T_MAIN}"]) == 0.0))
    for v in d.values():
        v.sort(key=lambda t: t[0])
    return d


def draw(ax, series):
    """series: (label, data, colour, requirement); requirement labels sit at the axis edge."""
    for label, data, col, req in series:
        xs = range(len(data))
        fil = [(x, y) for x, (y, _, s) in zip(xs, data) if not s]
        sat = [(x, y) for x, (y, _, s) in zip(xs, data) if s]
        ax.scatter([x for x, _ in fil], [y for _, y in fil], s=7, color=col, lw=0, zorder=4, label=label)
        ax.scatter([x for x, _ in sat], [y for _, y in sat], s=9, facecolors="none", edgecolors=col,
                   lw=0.6, zorder=4, label=f"{label}, no false positive observed")
        ax.scatter(list(xs), [y for _, y, _ in data], s=3.5, color=col, lw=0, alpha=0.18, zorder=3)
        ax.axhline(q(req), color=col, lw=0.8, ls=(0, (4, 3)), zorder=2)
    ax.set_ylabel("Q = −10·log₁₀(false-positive rate)\nat threshold 0.98", color=INK)
    ax.set_xlabel("5-mer contexts, ranked by specificity", color=INK)
    tidy(ax)
    for label, data, col, req in series:
        ax.text(0.995, q(req) + 0.6, f"required Q{q(req):.1f}", color=col, fontsize=6.9,
                ha="right", va="bottom", transform=ax.get_yaxis_transform(which="grid"))


d = load(args.table)
fig, axes = plt.subplots(2, 1, figsize=(WIDTH_IN, 6.9), gridspec_kw={"hspace": 0.34})
fig.patch.set_facecolor(SURFACE)
for ax in axes:
    ax.set_facecolor(SURFACE)

draw(axes[0], [("CpG", d["CpG"], C1, REQ["CpG"]), ("CpH", d["CpH"], C2, REQ["CpH"])])
axes[0].set_title("a", loc="left", color=INK, pad=5, fontweight="bold")
lo, hi = axes[0].get_ylim()
axes[0].set_ylim(lo, hi + 0.22 * (hi - lo))
axes[0].legend(frameon=False, labelcolor=INK2, loc="upper left", ncol=2, fontsize=6.6,
               handletextpad=0.3, columnspacing=1.0, borderaxespad=0.3, labelspacing=0.3)

draw(axes[1], [("A", d["A"], C1, REQ["A"])])
nfail = sum(1 for y, _, _ in d["A"] if y < q(REQ["A"]))
axes[1].set_title("b", loc="left", color=INK, pad=5, fontweight="bold")
axes[1].axvspan(-0.5, nfail - 0.5, color=C2, alpha=0.16, zorder=1)     # the rejected contexts
axes[1].legend(frameon=False, labelcolor=INK2, loc="lower right", fontsize=6.9,
               handletextpad=0.3, borderaxespad=0.4, labelspacing=0.3)

save(fig, args.out)
print(f"{args.out}: adenine contexts below requirement {nfail}/{len(d['A'])}")
