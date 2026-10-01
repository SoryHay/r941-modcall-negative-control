# Methods notes

The processing decisions that change a number in the article, and the reason for each. Thresholds are in
`config/params.yaml`; the commands are in `workflow/README.md`.

## Reads

Only pass reads of Guppy 5.0.16 are used. Before extraction, both arms pass the same filters (`read_filters.py`):
primary alignments, MAPQ ≥ 20, total soft-clipping ≤ 300 bp per read, and read mean quality at or above the higher of
the two arms' median read quality minus 1 (9.145 here).

**Soft-clipping is limited in base pairs, not as a fraction of the read.** SQK-PCS111 leaves a primer of about 110 bp
at each end of every cDNA read (median total clip 212 bp, 95th percentile 255 bp, 98.4 % of reads ≤ 300 bp). With
a median control read of 779 bp, a rule such as "≤ 20 % of the read" would reject short reads for their primers and
remove most of the control arm while barely touching the genomic arm. A fixed allowance just above the control's 95th
percentile applies the same rule to both arms. The primers are soft-clipped by the aligner and calls inside a clip
carry no reference position, so they never enter a count.

**Read mean quality is computed from the mean error probability**, Q = −10 log₁₀(mean(10^(−q/10))), as ONT does, not
as the arithmetic mean of Phred scores, which overstates quality and would move the cutoff.

## Per-read calls and sequence context (`fpr_per_kmer_interior.awk`)

Calls come from `modkit extract full` (modkit 0.6.4) after `update-tags --mode implicit`.

**Every position is an observation, including those without an explicit probability.** In implicit mode, positions
not listed as modified are canonical; 96.7 % of rows are such inferred positions. They are observations of an
unmodified base that the caller did not call, so they enter the denominator with probability 0. Leaving them out
would shrink the denominator about thirty-fold and inflate every rate accordingly.

**A call is counted as modified when its probability exceeds the threshold** (p > t), at t = 0.5 … 0.99.

**Contexts are keyed on the read-derived 5-mer** (`query_kmer`), whose centre equals the called base in 100 % of rows.
The reference-derived 5-mer is also tabulated (`*.ref.tsv`) after two corrections: modkit writes it in reference (+)
orientation while the modification lies on the read strand, so it is reverse-complemented on minus-strand reads
(otherwise half of all observations fall in the complementary context); and the soft-masked reference is upper-cased
(otherwise `gcAAA` and `GCAAA` are counted as different contexts). Windows that run off a read end are padded with
N; those rows are dropped and their number reported.

**Contexts with fewer than 200 observations are reported as NA.** A context with no false positive gets Q from a
single-event pseudocount and is marked as a bound (`>x`).

## Read-end exclusion (`read_end_profile.awk`, supplement S1)

Distance to the nearest read end is d = min(position, read length − 1 − position) in read coordinates, where a
read-end effect on the signal would act. The adenine false-positive rate is not flat near read ends (S1), so
positions with d < 250 bp are removed from numerator and denominator in every arm.

## Splice junctions (`splice_projection_check.py`, supplement S7)

cDNA reads span introns, which minimap2 records as `N` operations; 85 % of primary control reads carry one. If the
projection of a modification onto the reference miscounted across `N`, every call after the first junction would
land on the wrong base while the output still looked valid. The projection was therefore re-derived independently
with pysam (`modified_bases` for read positions, `get_aligned_pairs` for the read-to-reference map) and compared with
modkit's `ref_position` on junction-spanning reads: 659 of 659 calls agree. Only calls with an explicit probability are
compared, since pysam does not list implied-canonical positions.
