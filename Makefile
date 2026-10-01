# Tier 1: the article's figures and Table 1 from the tables in data/ (seconds, no GPU).
PY  ?= python3
F   := data/per5mer_fpr
OUT := figures/output

all: figures table1 supplementary

figures: $(OUT)/Fig1.pdf $(OUT)/Fig2.pdf $(OUT)/Fig3.pdf

$(OUT)/Fig1.pdf: figures/fig1_context_specificity.py $(F)/JAR_cDNA_run1.tsv config/params.yaml
	@mkdir -p $(OUT); $(PY) $< $(F)/JAR_cDNA_run1.tsv $(OUT)/Fig1

$(OUT)/Fig2.pdf: figures/fig2_control_concordance.py $(F)/JAR_cDNA_run1.tsv $(F)/JAR_cDNA_run2.tsv $(F)/JEG3_cDNA.tsv config/params.yaml
	@mkdir -p $(OUT); $(PY) $< $(F)/JAR_cDNA_run1.tsv $(F)/JAR_cDNA_run2.tsv $(F)/JEG3_cDNA.tsv $(OUT)/Fig2

$(OUT)/Fig3.pdf: figures/fig3_mtdna_6ma_contrast.py data/admsc_chrM_6mA_counts.tsv data/control_chrM_6mA_background.tsv $(F)/JAR_cDNA_run1.tsv config/params.yaml
	@mkdir -p $(OUT); $(PY) $< data/admsc_chrM_6mA_counts.tsv data/control_chrM_6mA_background.tsv $(F)/JAR_cDNA_run1.tsv $(OUT)/Fig3

table1: figures/table1_requirement.py $(F)/JAR_cDNA_run1.tsv config/params.yaml
	@mkdir -p $(OUT); $(PY) $< $(F)/JAR_cDNA_run1.tsv | tee $(OUT)/Table1.tsv

supplementary: supplementary/s2_allowlist.py supplementary/s5_failing_context_composition.py
	$(PY) supplementary/s2_allowlist.py $(F)/JAR_cDNA_run1.tsv $(F)/JAR_cDNA_run2.tsv $(F)/JEG3_cDNA.tsv supplementary/S2/S2_adenine_context_allowlist.tsv
	$(PY) supplementary/s5_failing_context_composition.py $(F)/JAR_cDNA_run1.tsv supplementary/S5/S5_failing_contexts

clean:
	rm -rf $(OUT)

.PHONY: all figures table1 supplementary clean
