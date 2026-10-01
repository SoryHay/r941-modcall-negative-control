# r941-modcall-negative-control

Supplement to the article

> Gospodinova A., I. Tsvetkova, R. Tsvetankova, Y. Mariienko, K. Todorova, S. Hayrabedyan (2026)
> **A modification-free negative control defines the operating points at which 5mC and 6mA can be called
> from R9.4.1 nanopore data.** *Comptes rendus de l'Académie bulgare des Sciences* (Proceedings of the
> Bulgarian Academy of Sciences), in press. DOI and link: *to be added on publication*.

The published figures and Table 1 are in the article; they are not redistributed here. The scripts in
`figures/` regenerate them locally from the tables in `data/`.

This repository holds the code that reproduces the figures and analyses of the article, the tables behind
them, and supplementary results that did not fit the journal's limit of four figures and tables. It is
provided so that the logic of every number in the article can be checked.

## What the article does

PCR-amplified cDNA cannot carry 5mC or 6mA, so every modification call made on it is a false positive.
Three such controls (JAR cDNA on two flow cells, JEG-3 cDNA on a third; R9.4.1, SQK-PCS111) were basecalled
with Guppy 5.0.16 and the Rerio all-context model, and the false-positive rate was measured per 5-mer context
and compared with the rate each modification's prevalence allows:

FPR_max = α·p·S / [(1 − α)(1 − p)]   (equation 1; α = 0.05, S = 1)

## Contents

| | Item | Refers to |
|---|---|---|
| Fig. 1 | Per-5-mer false-positive rate, cytosine and adenine | §3.1, Table 1 |
| Fig. 2 | Per-context agreement between the three amplified controls | §3.3 |
| Fig. 3 | Mitochondrial 6mA contrast across the threshold sweep | §3.6 |
| **S1** | False-positive rate by distance from the nearest read end | §2, read-end exclusion |
| **S2** | Full per-5-mer false-positive tables, three controls, all thresholds | Fig. 1, Fig. 2, Table 1 |
| **S3** | Dorado 0.3.3 + NEMO: genome-wide and chrM rates by threshold | §3.4 |
| **S4** | Effect of omitting the repair step on 6mA and CpG calls, by threshold | §3.5 |
| **S5** | Sequence composition of the 30 failing adenine contexts (16 tests, Benjamini–Hochberg) | §3.2 |
| **S6** | mtDNA coverage and 6mA call positions (PvuII site, D-loop) | §2, §3.6 |
| **S7** | Splice-junction projection check | §2 |

## Layout

```
config/          paths and the article's parameters (no hard-coded paths in the code)
workflow/        raw signal -> per-read calls -> per-context tables
figures/         scripts that draw Fig. 1-3 from the tables in data/
supplementary/   S1-S7
data/            the small tables the figures and supplements are drawn from
docs/            notes on methods
```

## Reproducing the figures and Table 1

```bash
conda env create -f environment.yml && conda activate r941-modcall-negative-control
make            # -> figures/output/Fig1-3 (.pdf, 300/600 dpi .png) and Table1.tsv, locally
```

Parameters (α, S, prevalences, the 250 bp read-end exclusion, the allowlist threshold) are read from
`config/params.yaml`; none is repeated in the code. Checked 2026-09-30 against the published files:
Fig. 1 and Fig. 2 are pixel-identical; Fig. 3 differs in 266 of 2.6 million pixels by at most 1/255 in
anti-aliasing, because the rates are now computed from the counts rather than read back at six digits.
Every cell of Table 1 and the §3.1 figures (23.2-fold, 5.9-fold, 10.8 %, the 30 excluded contexts) are
reproduced exactly. Figures render in Times New Roman when it is installed, otherwise in DejaVu Serif.

## Data

`data/per5mer_fpr/` — per-5-mer false-positive tables for the three amplified controls, read interiors
(≥ 250 bp from a read end), contexts keyed on the read-derived k-mer. Columns: `kmer`, `class`
(CpG / CpH / A), `n_obs`, then for each threshold t: `fp_t` (calls with probability > t), `FPR_t`
(0 = no false positive observed; NA if `n_obs` < 200) and `Q_t` (−10 log₁₀ FPR; `>x` is a bound).
`data/control_chrM_6mA_background.tsv` — the reference control's adenine call rate on chrM, read
interiors (dashed curve of Fig. 3a). `data/admsc_chrM_6mA_counts.tsv` — 6mA calls and adenine
observations per threshold in the two AdMSC arms (aggregate counts only, as plotted in Fig. 3).

- **JAR and JEG-3 sequencing data** (the three amplified controls and the three genomic libraries):
  European Nucleotide Archive, project **[PRJEB127679](https://www.ebi.ac.uk/ena/browser/view/PRJEB127679)**
  (study ERP207156) — per library, the Guppy 5.0.16 + Rerio modified-base BAM that every number in the
  article was computed from, and the raw signal (fast5) of exactly the files analysed.
- **Adipose-derived mesenchymal stromal cell (AdMSC) data** (Fig. 3, S6): released with a separate
  publication; the link will be added here. Only aggregate values appear in this repository.

## Licence

Code: MIT (`LICENSE`). Tables in `data/` and `supplementary/`, and the S6 figure: CC BY 4.0 (`LICENSE-DATA`).
The article's own figures and Table 1 are not covered by these licences and are not included.

## Contact

Soren Hayrabedyan, Laboratory of Reproductive Omics Technologies, Institute of Biology and Immunology of
Reproduction "Acad. Kiril Bratanov", Bulgarian Academy of Sciences, Sofia. ORCID 0000-0002-2147-1982.
