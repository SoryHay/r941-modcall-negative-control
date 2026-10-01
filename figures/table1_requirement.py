#!/usr/bin/env python3
"""Table 1 — pooled FPR per mark against equation (1), contexts qualifying, and the adenine allowlist gain (§3.1)."""
import argparse

from _common import PARAMS, read_fpr_table, rate_or_bound, requirement

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("table", help="per-5-mer FPR table of the reference control (JAR cDNA run 1)")
args = ap.parse_args()

rows = read_fpr_table(args.table)
T_SEL = str(PARAMS["allowlist_threshold"])


def pooled(rs, t):
    return sum(int(r[f"fp_{t}"]) for r in rs) / sum(int(r["n_obs"]) for r in rs)


print("mark\tprevalence\trequired_FPR\tn_obs\tFPR_0.7\tFPR_0.98\tqualifying\tallowlist_FPR_0.95")
for mark, label in (("CpG", "5mCpG"), ("CpH", "5mCpH"), ("A", "6mA")):
    rs = [r for r in rows if r["class"] == mark]
    req = requirement(mark)
    keep = [r for r in rs if rate_or_bound(r, T_SEL) <= req]
    allow = f"{pooled(keep, '0.95'):.2e}" if len(keep) < len(rs) else "not required"
    print(f"{label}\t{PARAMS['requirement']['prevalence'][mark]}\t{req:.2e}\t{sum(int(r['n_obs']) for r in rs)}\t"
          f"{pooled(rs, '0.7'):.2e}\t{pooled(rs, '0.98'):.2e}\t{len(keep)} / {len(rs)}\t{allow}")

a = [r for r in rows if r["class"] == "A"]
keep = [r for r in a if rate_or_bound(r, T_SEL) <= requirement("A")]
lost = 1 - sum(int(r["n_obs"]) for r in keep) / sum(int(r["n_obs"]) for r in a)
print(f"\n6mA at 0.7: {pooled(a, '0.7') / requirement('A'):.1f} x the admissible rate")
print(f"6mA at 0.95: all {pooled(a, '0.95'):.2e} -> allowlist {pooled(keep, '0.95'):.2e} "
      f"({pooled(a, '0.95') / pooled(keep, '0.95'):.1f}-fold), observations lost {100 * lost:.1f}%")
print("excluded contexts:", " ".join(sorted(r["kmer"] for r in a if r not in keep)))
