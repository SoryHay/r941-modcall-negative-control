#!/usr/bin/env python3
"""S2 — the 256 adenine contexts: allowlist status (from JAR cDNA run 1) and FPR at 0.98 in all three controls."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
from _common import PARAMS, read_fpr_table, rate_or_bound, requirement

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("ref", help="JAR cDNA run 1 (defines the allowlist)")
ap.add_argument("run2", help="JAR cDNA run 2")
ap.add_argument("jeg", help="JEG-3 cDNA")
ap.add_argument("out", help="output .tsv")
args = ap.parse_args()

T = str(PARAMS["allowlist_threshold"])
FMAX = requirement("A")
tables = {name: {r["kmer"]: r for r in read_fpr_table(p) if r["class"] == "A"}
          for name, p in (("JAR_run1", args.ref), ("JAR_run2", args.run2), ("JEG3", args.jeg))}
ref = tables["JAR_run1"]

with open(args.out, "w") as fh:
    fh.write("kmer\tstatus\t" + "\t".join(f"n_obs_{n}\tfp_{n}\tFPR_{n}\tabove_requirement_{n}" for n in tables) + "\n")
    for k in sorted(ref, key=lambda k: rate_or_bound(ref[k], T)):
        status = "allowlist" if rate_or_bound(ref[k], T) <= FMAX else "excluded"
        cells = []
        for t in tables.values():
            r = t[k]
            cells += [r["n_obs"], r[f"fp_{T}"], r[f"FPR_{T}"], str(rate_or_bound(r, T) > FMAX).lower()]
        fh.write(f"{k}\t{status}\t" + "\t".join(cells) + "\n")
excluded = {k for k in ref if rate_or_bound(ref[k], T) > FMAX}
for name in ("JAR_run2", "JEG3"):
    rates = {k: rate_or_bound(r, T) for k, r in tables[name].items()}
    worst = set(sorted(rates, key=rates.get, reverse=True)[:len(excluded)])
    print(f"{name}: {len(excluded & worst)} of {len(excluded)} excluded contexts are among its {len(excluded)} worst")
