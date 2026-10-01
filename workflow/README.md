# Tier 2 — from raw signal to the tables in `data/`

The steps that produced the per-read calls. They need the raw signal (ENA, see the main README), a GPU for
basecalling, and T2T-CHM13v2.0 (chrM = `CP068254.1`). Parameters are read from `config/params.yaml`.

| step | command | version |
|---|---|---|
| 1. pod5 → fast5 (amplified controls only; Guppy 5 reads fast5) | `pod5 convert to_fast5 -r -t 4 --file-read-count 4000 -o fast5/ <pod5 dir>` | pod5 0.3.34 |
| 2. basecall + modified bases, pass reads only | `guppy_basecaller -i fast5/ -s out/ -c res_dna_r941_min_modbases-all-context_v001.cfg -d <rerio>/basecall_models -x cuda:0 --bam_out` (docker `genomicpariscentre/guppy-gpu:5.0.16`) | Guppy 5.0.16, Rerio all-context v001 |
| 3. tags, alignment | `align_and_extract.sh out/pass <fa> aln/ splice` (cDNA) or `genomic` | modkit 0.6.4, minimap2 2.29, samtools 1.19.2 |
| 4. read filters, both arms | `read_filters.py <native.bam> <cdna.bam> qc/` — primary, MAPQ ≥ 20, total soft-clip ≤ 300 bp, read mean-Q ≥ (higher median of the two arms − 1) = 9.145 | pysam |
| 5. per-read calls | `FILTERED_BAM=qc/cdna.filt.bam align_and_extract.sh …` → `modkit extract full --mapped-only --kmer-size 5` | modkit 0.6.4 |
| 6. per-5-mer tables (`data/per5mer_fpr/`) | `zcat calls.tsv.gz \| awk -v OUT=<prefix> -f fpr_per_kmer_interior.awk` | — |
| 7. read-end and read-length profile (S1) | `zcat calls.tsv.gz \| awk -v OUT=<prefix> -v Q1=663 -v Q2=748 -v Q3=967 -v Q4=1296 -f read_end_profile.awk` | — |
| 8. pooled rates (Fig. 3, S3, S4) | `pooled_rate.py calls.tsv.gz --base A [--chrom CP068254.1]`; `--base C --context CpG` for CpG | — |
| 9. splice projection check (S7) | `splice_projection_check.py subset.bam subset_calls.tsv out.tsv` | pysam |

**What was run on which reads.** JAR cDNA run 1 (FAS93941): all 36 pod5 files, 317,107 reads after the
filters of step 4. JAR cDNA run 2 (FAS90841) and JEG-3 cDNA (FAS87609): 8 of 72 and 8 of 63 pod5 files.
Native JAR genomic DNA (FAS94147): all reads basecalled, calls extracted from 10,005 reads
(`samtools view -s 42.01605`). Unrepaired JAR genomic DNA (FAS90851): every 8th fast5 file (239 of 1,909),
calls extracted from 9,888 reads. Dorado 0.3.3 + `NEMO_R9_6mA` (S3): JAR cDNA run 1, primary MAPQ ≥ 20
reads, same extraction. The file lists are given with the ENA accessions.

The reasoning behind each step (strand correction of reference k-mers, soft-masking, inferred rows as
observations with probability 0, the soft-clip rule) is in `docs/methods_notes.md`.
