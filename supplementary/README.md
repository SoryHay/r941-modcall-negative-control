# Supplementary results

Each item names the section of the article it supports and the command that produced it. All rates are
over read positions ≥ 250 bp from the nearest read end unless stated otherwise.

**S1 — false-positive rate by distance to the nearest read end, and by read length** (Methods,
read-end exclusion). `S1/S1_fpr_by_distance_to_read_end.tsv`, `S1/S1_fpr_by_read_length_quintile.tsv`:
JAR cDNA run 1, all read positions, per base class and threshold. The adenine rate peaks 50–100 bp from
a read end (3.83 × 10⁻³ at 0.7) and is flat beyond 250 bp; the outer-zone-to-interior ratio at 0.7 is
2.64. Produced by `workflow/read_end_profile.awk`.

**S2 — the 256 adenine contexts** (Fig. 1b, Fig. 2, Table 1, §3.3). `S2/S2_adenine_context_allowlist.tsv`:
allowlist status defined on JAR cDNA run 1 at threshold 0.98, and the observations, false positives and
rate of every context in all three controls. The full per-context tables for every base class and
threshold are in `../data/per5mer_fpr/`. `python supplementary/s2_allowlist.py …` also prints the
overlap of the 30 excluded contexts with the 30 worst of each other control (25 and 23 of 30).

**S3 — Dorado 0.3.3 + NEMO_R9_6mA** (§3.4). `S3/S3_nemo_6mA_rates.tsv`: adenine call rate per
threshold on JAR cDNA run 1 genome-wide and on chrM, and on chrM for the two AdMSC arms (aggregate
counts). Genome-wide rates run 5–17 % above the chrM ones. Produced by `workflow/pooled_rate.py`.

**S4 — omitting the repair step** (§3.5). `S4/S4_unrepaired_vs_repaired.tsv`: 6mA and CpG call rates
per threshold in unrepaired (FAS90851, 9,888 reads) and Fpg-repaired (FAS94147, 10,005 reads) JAR genomic
DNA. 6mA is 9.2- to 15.1-fold higher in unrepaired DNA; CpG falls to 0.80 at 0.7 and further at higher
thresholds. Produced by `workflow/pooled_rate.py` (`--base C --context CpG` for CpG).

**S5 — composition of the 30 failing adenine contexts** (§3.2). `S5/S5_failing_contexts.positions.tsv`:
base at each flanking position, failing against passing contexts, two-proportion test with
Benjamini–Hochberg correction across the 16 position × base comparisons (positions 1, 2, 4, 5; position 3
is the called adenine). `S5/S5_failing_contexts.aggregate.tsv`: homopolymer length, adenine, GC and
purine content. Produced by `s5_failing_context_composition.py`.

**S6 — mitochondrial coverage and 6mA call positions** (Methods, PvuII linearisation; §3.6).
`S6/S6_mtdna_coverage_6mA.{pdf,_300dpi.png,_600dpi.png}`: (a) read depth along chrM in the two AdMSC arms —
the segment 1–2,073 distal to the PvuII site is under-covered because the circle is linearised there;
(b) every 6mA call at threshold 0.98 in read interiors (Guppy 5.0.16 + Rerio, the calls of Fig. 3),
as the per cent of reads at that position. The calls are dispersed (NT 156 at 124 positions, PA 19 at
17); the D-loop carries 3 and 0. Produced by `s6_mtdna_coverage.py` from the AdMSC alignments, which are
released with a separate publication; only the figure is given here.

**S7 — splice-junction projection check** (Methods, alignment). `S7/S7_splice_projection_check.tsv`:
for 23 junction-spanning reads on 11 chromosomes, the reference position of each of 659 modification
calls derived independently by modkit and by pysam; all agree. Produced by
`workflow/splice_projection_check.py`.
