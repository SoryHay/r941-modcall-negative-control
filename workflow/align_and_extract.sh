#!/usr/bin/env bash
# Guppy modbase BAMs (pass only) -> MM/ML tags -> T2T alignment -> per-read calls (modkit extract full).
# USAGE: align_and_extract.sh <pass_bam_dir> <T2T-CHM13v2.0.fa> <outdir> <genomic|splice> [threads]
# The read filters (read_filters.py) run between alignment and extraction; see workflow/README.md.
set -euo pipefail
IN=${1:?pass BAM dir}; FA=${2:?reference fasta}; OUT=${3:?outdir}; MODE=${4:?genomic|splice}; T=${5:-12}
MODKIT=${MODKIT:-modkit}
mkdir -p "$OUT"

case "$MODE" in
  genomic) MM2="-ax map-ont -y" ;;
  splice)  MM2="-ax splice -uf -k14 -y" ;;            # cDNA: 85% of reads span a junction
  *) echo "mode must be genomic or splice" >&2; exit 2 ;;
esac

samtools cat "$IN"/*.bam > "$OUT/merged.bam"
"$MODKIT" update-tags --mode implicit "$OUT/merged.bam" "$OUT/upd.bam"     # Guppy writes legacy Mm/Ml
rm -f "$OUT/merged.bam"
samtools fastq -T MM,ML "$OUT/upd.bam" \
  | minimap2 $MM2 -t "$T" "$FA" - 2>"$OUT/minimap2.log" \
  | samtools sort -@ 4 -o "$OUT/aln.bam" -
samtools index "$OUT/aln.bam"

if [ -n "${FILTERED_BAM:-}" ]; then IN_BAM=$FILTERED_BAM; else IN_BAM=$OUT/aln.bam; fi
"$MODKIT" extract full --reference "$FA" --mapped-only --kmer-size 5 --bgzf -t "$T" --io-threads 4 \
  "$IN_BAM" "$OUT/calls.tsv.gz"
