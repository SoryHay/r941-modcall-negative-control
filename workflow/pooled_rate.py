#!/usr/bin/env python3
"""Pooled modification-call rate per threshold for one canonical base, read interiors only (modkit extract full input)."""
import argparse
import gzip
import sys
from pathlib import Path

import yaml

PARAMS = yaml.safe_load(open(Path(__file__).resolve().parents[1] / "config" / "params.yaml"))

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("calls", help="modkit 0.6.4 `extract full` table (.tsv or .tsv.gz)")
ap.add_argument("--base", default="A", choices=["A", "C"])
ap.add_argument("--context", choices=["CpG", "CpH"], help="for --base C: next base on the read is / is not G")
ap.add_argument("--chrom", help="restrict to one contig, e.g. CP068254.1 (chrM)")
args = ap.parse_args()

edge = PARAMS["read_end_exclusion_bp"]
thr = PARAMS["thresholds"]
opener = gzip.open if args.calls.endswith(".gz") else open
n, k = 0, [0] * len(thr)
with opener(args.calls, "rt") as fh:
    col = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[col["canonical_base"]] != args.base or (args.chrom and f[col["chrom"]] != args.chrom):
            continue
        if args.context and (f[col["query_kmer"]].upper()[3] == "G") != (args.context == "CpG"):
            continue
        pos, rlen = int(f[col["forward_read_position"]]), int(f[col["read_length"]])
        if min(pos, rlen - 1 - pos) < edge:
            continue
        p = 0.0 if f[col["inferred"]] == "true" else float(f[col["mod_qual"]])   # implicit canonical -> p = 0
        n += 1
        for i, t in enumerate(thr):
            if p > t:
                k[i] += 1

w = sys.stdout.write
w("threshold\tn_mod\tn_obs\tfpr\n")
for t, ki in zip(thr, k):
    w(f"{t}\t{ki}\t{n}\t{ki / n:.4e}\n")
