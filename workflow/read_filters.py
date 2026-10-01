#!/usr/bin/env python3
"""Read filters applied to both arms before extraction: primary, MAPQ, total soft-clip, shared mean-Q cutoff.

USAGE: read_filters.py <native.bam> <cdna.bam> <outdir>   (thresholds from config/params.yaml)
"""
import sys
import os
import math
from statistics import median

import pysam

native_bam, cdna_bam, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)

from pathlib import Path
import yaml
RF = yaml.safe_load(open(Path(__file__).resolve().parents[1] / "config" / "params.yaml"))["read_filters"]
MAPQ_MIN = RF["mapq_min"]
SOFTCLIP_MAX_BP = RF["softclip_max_bp"]


def mean_q(read):
    q = read.query_qualities
    if not q:
        return None
    return -10.0 * math.log10(sum(10.0 ** (-x / 10.0) for x in q) / len(q))


def softclip_bp(read):
    ct = read.cigartuples or []
    return sum(l for op, l in ct if op == 4)             # 4 == S


def is_primary(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary)


# ---- pass 1: per-arm mean-Q distribution over primary reads -----------------
stats = {}
for name, path in (("native", native_bam), ("cdna", cdna_bam)):
    qs = []
    total = 0
    primary = 0
    bam = pysam.AlignmentFile(path, "rb")
    for read in bam.fetch(until_eof=True):
        total += 1
        if not is_primary(read):
            continue
        primary += 1
        q = mean_q(read)
        if q is not None:
            qs.append(q)
    bam.close()
    stats[name] = {"total": total, "primary": primary, "medianQ": median(qs) if qs else 0.0}
    print(f"[pass1] {name}: total={total} primary={primary} medianQ={stats[name]['medianQ']:.3f}",
          flush=True)

Q_CUTOFF = max(stats["native"]["medianQ"], stats["cdna"]["medianQ"]) - RF["mean_q_below_higher_median"]
# Sensitivity only, not applied: the same rule anchored on the lower median.
Q_CUTOFF_LOW = min(stats["native"]["medianQ"], stats["cdna"]["medianQ"]) - 1.0
print(f"[pass1] shared mean-Q cutoff = {Q_CUTOFF:.3f} "
      f"(sensitivity-only alternative = {Q_CUTOFF_LOW:.3f})", flush=True)

# ---- pass 2: apply, count, write -------------------------------------------
rows = []
for name, path in (("native", native_bam), ("cdna", cdna_bam)):
    c = dict(total=0, primary=0, mapq10=0, mapq20=0, mapq30=0,
             after_mapq=0, after_softclip=0, after_q=0, sens_q_low=0)
    src = pysam.AlignmentFile(path, "rb")
    out_path = os.path.join(outdir, f"{name}.filt.bam")
    dst = pysam.AlignmentFile(out_path, "wb", template=src)
    for read in src.fetch(until_eof=True):
        c["total"] += 1
        if not is_primary(read):
            continue
        c["primary"] += 1
        if read.mapping_quality >= 10:
            c["mapq10"] += 1
        if read.mapping_quality >= 20:
            c["mapq20"] += 1
        if read.mapping_quality >= 30:
            c["mapq30"] += 1
        if read.mapping_quality < MAPQ_MIN:
            continue
        c["after_mapq"] += 1
        if softclip_bp(read) > SOFTCLIP_MAX_BP:
            continue
        c["after_softclip"] += 1
        q = mean_q(read)
        if q is not None and q >= Q_CUTOFF_LOW:
            c["sens_q_low"] += 1
        if q is None or q < Q_CUTOFF:
            continue
        c["after_q"] += 1
        dst.write(read)
    src.close()
    dst.close()
    pysam.index(out_path)
    c["medianQ"] = stats[name]["medianQ"]
    rows.append((name, c))
    print(f"[pass2] {name}: {c}", flush=True)

with open(os.path.join(outdir, "read_filter.tsv"), "w") as fh:
    fh.write(f"# shared mean-Q cutoff\t{Q_CUTOFF:.3f}\n")
    fh.write(f"# MAPQ_MIN\t{MAPQ_MIN}\n")
    fh.write(f"# SOFTCLIP_MAX_BP\t{SOFTCLIP_MAX_BP}\t(total per read)\n")
    fh.write(f"# meanQ cutoff sensitivity-only alternative\t{Q_CUTOFF_LOW:.3f}\tNOT APPLIED\n")
    fh.write("# end_reason filter\tnot applied (no adaptive sampling)\n")
    fh.write("arm\tstep\treads\tretained_pct_of_primary\n")
    for name, c in rows:
        p = c["primary"]
        fh.write(f"{name}\t0_total_records\t{c['total']}\t\n")
        fh.write(f"{name}\t1_primary\t{p}\t100.00\n")
        fh.write(f"{name}\t2_mapq_ge20\t{c['after_mapq']}\t{100*c['after_mapq']/p:.2f}\n")
        fh.write(f"{name}\t3_softclip_le300bp\t{c['after_softclip']}\t{100*c['after_softclip']/p:.2f}\n")
        fh.write(f"{name}\t4_meanQ_ge_cutoff\t{c['after_q']}\t{100*c['after_q']/p:.2f}\n")
        fh.write(f"{name}\tsens_meanQ_ge_lower\t{c['sens_q_low']}\t{100*c['sens_q_low']/p:.2f}\n")
        fh.write(f"{name}\tsens_mapq_ge10\t{c['mapq10']}\t{100*c['mapq10']/p:.2f}\n")
        fh.write(f"{name}\tsens_mapq_ge20\t{c['mapq20']}\t{100*c['mapq20']/p:.2f}\n")
        fh.write(f"{name}\tsens_mapq_ge30\t{c['mapq30']}\t{100*c['mapq30']/p:.2f}\n")
        fh.write(f"{name}\tmedian_meanQ\t{c['medianQ']:.3f}\t\n")
print("done")
