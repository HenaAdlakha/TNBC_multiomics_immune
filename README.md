# TNBC Multi-Omics Immune Analysis

## Overview
Unsupervised multi-omics factor analysis (MOFA2) identifying immune-cold and immune-hot 
subgroups in triple-negative breast cancer (TNBC), with independent validation in METABRIC 
and TIGIT as a candidate therapeutic target.

**Key findings:**
- Factor 4 from MOFA2 stratifies TNBC by immune infiltration (n=81 TCGA)
- Immune-cold patients show significantly worse overall survival (TCGA p=0.045; METABRIC p=0.00039)
- Cox HR=1.195 (95% CI 1.011–1.412, p=0.036) independent of age and NPI
- Immune-hot tumours show higher TIGIT, PD-1, PD-L1, CTLA4, LAG3 (all FDR<0.001 in METABRIC)
- CpG island enrichment confirms epigenetic basis (chi-square p=1.72×10⁻⁵⁶)

## Data Sources
- **TCGA BRCA:** Downloaded from [UCSC Xena](https://xenabrowser.net)
  - RNA-seq (n=1,218), methylation 450k (n=888), CNV GISTIC2 (n=1,080)
  - Clinical metadata and curated survival data
  - 68 pan-cancer immune signatures (Denise Wolf et al.)
- **METABRIC:** Downloaded from [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) (n=1,980)

Data files are not included in this repository due to size and data policy.
Download instructions are in each script header.

## Scripts

| Script | Description |
|---|---|
| `01_load_data.R` | Load all TCGA data files and initial sanity checks |
| `02_preprocessing.R` | TNBC filtering (ER−/PR−/HER2−), feature selection, QC plots |
| `03_mofa_setup.R` | Build MOFA2 object with 3 omics views |
| `04_mofa_train.R` | Train MOFA2 model (12 factors, seed=42) |
| `05_immune_characterization.R` | Correlate factors with 68 immune signatures; Wilcoxon tests; Paradigm pathway enrichment |
| `06_fixes.R` | Minor fixes applied during development |
| `06_survival_analysis.R` | Kaplan-Meier and Cox survival analysis in TCGA (n=55/81) |
| `07_full_cohort_preprocessing.R` | Preprocessing for full BRCA cohort (n=772) |
| `08_full_cohort_mofa.R` | MOFA2 on full BRCA cohort; subtype comparison |
| `09_TNBC_3omics_preprocessing.R` | Expanded TNBC model — RNA + methylation + CNV (n=81) |
| `10_TNBC_3omics_immune.R` | Immune characterisation of Factor 4 in n=81 model |
| `11_TNBC_3omics_survival.R` | Survival analysis in n=81 cohort |
| `12_METABRIC_validation.R` | Apply Factor 4 signature to METABRIC; immune group survival |
| `13_subtype_comparison.R` | Tumour size, grade and PAM50 subtype comparison |
| `14_drug_target.R` | Differential expression; CMAP gene lists; volcano plot |
| `15_checkpoint_expression.R` | Checkpoint gene expression TCGA + METABRIC validation |
| `16_methylation_epigenetic_features.R` | Annotate methylation weights; CpG island Fisher tests |
| `17_factor_stability.R` | Factor stability across 8 seeds and 4 factor numbers |
| `18_permutation_test.R` | Permutation test; full 12×68 correlation matrix with FDR |
| `19_purity_adjustment.R` | Tumour purity adjustment using PTPRC; partial correlations |
| `20_locked_metabric_validation.R` | Locked continuous Cox using signed TCGA weights |
| `21_subtype_comparison.R` | PAM50 subtype independence of Factor 4 in TCGA |
| `22_final_analysis.R` | Final Cox model; within-subgroup analyses; freeze analysis |
| `23_Tests.R` | PCA ANOVA tests and omics contribution analysis (AUC) |
| `24_TIGIT_target.R` | TIGIT survival analysis; combined immune + TIGIT stratification |
| `25_variance_omics.R` | Variance explained and single-omics predictive power per layer |
| `Untitled.R` | Scratch script used during development |

## Analysis Pipeline

Data download (Xena/cBioPortal)
↓
01–04: Load, filter to TNBC, MOFA2 training
↓
05–06: Immune characterisation and survival (TCGA, n=81)
↓
07–08: Full cohort context (n=772)
↓
09–11: Expanded TNBC model (n=81, 3 omics)
↓
12: METABRIC validation (n=320)
↓
13–16: Epigenetic features, checkpoints, clinical variables
↓
17–19: Robustness checks (stability, permutation, purity)
↓
20–22: Locked validation and final Cox model
↓
23–25: Statistical tests and omics contribution

## Requirements

**R version:** 4.5.3

**Key packages:**
```r
BiocManager::install(c("MOFA2", "limma",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"))

install.packages(c("tidyverse", "data.table", "ggplot2", "patchwork",
  "survival", "survminer", "ggforce", "pROC", "mediation"))
```

## Session Info
All analyses run on macOS (Apple Silicon, aarch64-apple-darwin20), R 4.5.3.

## Citation
Adlakha H, Arifin MZ. Multi-omics factor analysis identifies an immune-cold 
subgroup in triple-negative breast cancer with distinct epigenetic signatures 
and worse survival. 2026.

## Contact
Hena Adlakha — hena.adlakha1@ucdconnect.ie
