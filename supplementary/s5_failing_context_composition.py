#!/usr/bin/env python3
"""S5 — composition of the adenine contexts that fail equation (1) at 0.98 against those that pass (§3.2)."""
import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
from _common import PARAMS, read_fpr_table, rate_or_bound, requirement

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("table", help="per-5-mer FPR table of the reference control (JAR cDNA run 1)")
ap.add_argument("out", help="output prefix: <out>.positions.tsv, <out>.aggregate.tsv")
args = ap.parse_args()

T = str(PARAMS["allowlist_threshold"])
FMAX = requirement("A")


def two_prop_p(k1, n1, k2, n2):
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    return 1.0 if se == 0 else math.erfc(abs(k1 / n1 - k2 / n2) / se / math.sqrt(2))


def bh(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i], reverse=True)
    q, prev = [1.0] * len(ps), 1.0
    for rank, i in enumerate(order):
        prev = min(prev, ps[i] * len(ps) / (len(ps) - rank))
        q[i] = prev
    return q


def longest_run(s):
    best = run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


rows = [r for r in read_fpr_table(args.table) if r["class"] == "A" and r[f"FPR_{T}"] != "NA"]
fail = [r["kmer"] for r in rows if rate_or_bound(r, T) > FMAX]
ok = [r["kmer"] for r in rows if rate_or_bound(r, T) <= FMAX]

cells = [(pos, b, sum(k[pos - 1] == b for k in fail), sum(k[pos - 1] == b for k in ok))
         for pos in (1, 2, 4, 5) for b in "ACGT"]                     # position 3 is the called A
ps = [two_prop_p(kf, len(fail), kp, len(ok)) for _, _, kf, kp in cells]
with open(f"{args.out}.positions.tsv", "w") as fh:
    fh.write("position\tbase\tfailing\tn_failing\tpassing\tn_passing\tp\tq_BH\n")
    for (pos, b, kf, kp), p, qv in sorted(zip(cells, ps, bh(ps)), key=lambda x: x[1]):
        fh.write(f"{pos}\t{b}\t{kf}\t{len(fail)}\t{kp}\t{len(ok)}\t{p:.4f}\t{qv:.3f}\n")

features = {"longest_homopolymer_run": longest_run,
            "A_count": lambda k: k.count("A"),
            "A_count_flanks": lambda k: (k[:2] + k[3:]).count("A"),
            "GC_count": lambda k: k.count("G") + k.count("C"),
            "purine_count": lambda k: k.count("A") + k.count("G")}
with open(f"{args.out}.aggregate.tsv", "w") as fh:
    fh.write("feature\tmean_failing\tmean_passing\tabove_median_failing\tabove_median_passing\tp\n")
    for name, fn in features.items():
        vf, vp = [fn(k) for k in fail], [fn(k) for k in ok]
        med = sorted(vf + vp)[len(vf + vp) // 2]                     # proportion above the pooled median
        kf, kp = sum(v > med for v in vf), sum(v > med for v in vp)
        fh.write(f"{name}\t{sum(vf) / len(vf):.2f}\t{sum(vp) / len(vp):.2f}\t{kf / len(vf):.3f}\t"
                 f"{kp / len(vp):.3f}\t{two_prop_p(kf, len(vf), kp, len(vp)):.4f}\n")

print(f"{len(fail)} failing / {len(ok)} passing contexts; {len(cells)} position x base tests, BH across {len(cells)}")
