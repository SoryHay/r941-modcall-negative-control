#!/usr/bin/env python3
"""S7 — modkit and pysam must project every modification call to the same reference base across splice junctions.

USAGE: splice_projection_check.py <subset.bam> <modkit_extract_full.tsv> <out.tsv>   (exit 1 on any mismatch)
"""
import sys
import csv
from collections import defaultdict

import pysam

bam_path, modkit_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

# ---- 1. modkit's answer -----------------------------------------------------
modkit_calls = set()
modkit_rows = 0
inferred_rows = 0
with open(modkit_path) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for row in rdr:
        if row["inferred"].lower() == "true":
            inferred_rows += 1
            continue
        ref_pos = int(row["ref_position"])
        if ref_pos < 0:                      # soft-clipped / unaligned position
            continue
        modkit_rows += 1
        modkit_calls.add((row["read_id"], row["mod_code"], row["chrom"], ref_pos))

# ---- 2. pysam's independent answer ------------------------------------------
pysam_calls = set()
per_read = defaultdict(lambda: {"spliced": False, "n_calls": 0, "chrom": "", "len": 0})

bam = pysam.AlignmentFile(bam_path, "rb")
for read in bam.fetch(until_eof=True):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
    mods = read.modified_bases
    if not mods:
        continue
    q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}
    spliced = any(op == 3 for op, _ in (read.cigartuples or []))  # 3 == N
    info = per_read[read.query_name]
    info["spliced"] = spliced
    info["chrom"] = read.reference_name
    info["len"] = read.query_length
    for (_canon, _strand, mod_code), positions in mods.items():
        code = mod_code if isinstance(mod_code, str) else chr(mod_code)
        for qpos, _qual in positions:
            rpos = q2r.get(qpos)
            if rpos is None:                 # inside an insertion or soft clip
                continue
            pysam_calls.add((read.query_name, code, read.reference_name, rpos))
            info["n_calls"] += 1
bam.close()

# ---- 3. compare -------------------------------------------------------------
only_modkit = modkit_calls - pysam_calls
only_pysam = pysam_calls - modkit_calls
agree = modkit_calls & pysam_calls

spliced_reads = [r for r, i in per_read.items() if i["spliced"]]

with open(out_path, "w") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["reads_examined", len(per_read)])
    w.writerow(["reads_spanning_splice_junction", len(spliced_reads)])
    w.writerow(["chromosomes", len({i["chrom"] for i in per_read.values()})])
    w.writerow(["modkit_calls_compared", len(modkit_calls)])
    w.writerow(["modkit_rows_inferred_skipped", inferred_rows])
    w.writerow(["pysam_calls_compared", len(pysam_calls)])
    w.writerow(["agreeing_calls", len(agree)])
    w.writerow(["modkit_only", len(only_modkit)])
    w.writerow(["pysam_only", len(only_pysam)])
    verdict = "PASS" if (not only_modkit and not only_pysam) else "FAIL"
    w.writerow(["verdict", verdict])
    w.writerow([])
    w.writerow(["read_id", "chrom", "read_length", "spliced", "calls_projected"])
    for rid, i in sorted(per_read.items()):
        w.writerow([rid, i["chrom"], i["len"], i["spliced"], i["n_calls"]])
    if only_modkit or only_pysam:
        w.writerow([])
        w.writerow(["DISAGREEMENT", "read_id", "mod_code", "chrom", "ref_position"])
        for c in sorted(only_modkit)[:50]:
            w.writerow(["modkit_only", *c])
        for c in sorted(only_pysam)[:50]:
            w.writerow(["pysam_only", *c])

print(f"reads={len(per_read)} spliced={len(spliced_reads)} "
      f"chroms={len({i['chrom'] for i in per_read.values()})}")
print(f"modkit={len(modkit_calls)} pysam={len(pysam_calls)} agree={len(agree)} "
      f"modkit_only={len(only_modkit)} pysam_only={len(only_pysam)}")
print("VERDICT:", verdict)
sys.exit(0 if verdict == "PASS" else 1)
