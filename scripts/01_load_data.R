# ── INSTALL PACKAGES ─────────────────────────────────────────
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("MOFA2", "limma"), update = FALSE)
#install.packages(c("tidyverse", "data.table", "pheatmap",
#                   "ggplot2", "RColorBrewer", "survival",
#                   "survminer", "ggpubr", "patchwork"))

# ── LIBRARIES ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(data.table)
library(MOFA2)

# ── PATHS ─────────────────────────────────────────────────────────────────────
base_dir <- "~/Documents/R Working Directory/BRCA_project"
raw_dir  <- file.path(base_dir, "data/raw")
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")
mod_dir  <- file.path(base_dir, "outputs/models")

# ── LOAD RAW FILES ────────────────────────────────────────────────────────────
cat("Loading RNA-seq...\n")
rna_raw <- fread(file.path(raw_dir, "HiSeqV2"), header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("sample") %>%
  as.matrix()

cat("Loading Methylation 450k...\n")
meth_raw <- fread(file.path(raw_dir, "HumanMethylation450"), header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("sample") %>%
  as.matrix()

cat("Loading CNV (GISTIC2 thresholded)...\n")
cnv_raw <- fread(file.path(raw_dir, "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"),
                 header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("Gene Symbol") %>%
  as.matrix()

cat("Loading RPPA...\n")
rppa_raw <- fread(file.path(raw_dir, "RPPA_RBN"), header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("Sample_description") %>%
  as.matrix()

cat("Loading clinical metadata...\n")
clinical <- fread(file.path(raw_dir, "TCGA.BRCA.sampleMap-BRCA_clinicalMatrix"),
                  header = TRUE)

cat("Loading survival data...\n")
survival <- fread(file.path(raw_dir, "survival-BRCA_survival.txt"),
                  header = TRUE)

cat("Loading MC3 mutations...\n")
mc3 <- fread(file.path(raw_dir, "BRCA_mc3_gene_level.txt"),
             header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("sample") %>%
  as.matrix()

cat("✓ All files loaded.\n")

# Sanity check
cat("\nDimensions:\n")
cat("RNA:      ", dim(rna_raw),  "\n")
cat("Methyl:   ", dim(meth_raw), "\n")
cat("CNV:      ", dim(cnv_raw),  "\n")
cat("RPPA:     ", dim(rppa_raw), "\n")
cat("Clinical: ", dim(clinical), "\n")
cat("Survival: ", dim(survival), "\n")
cat("MC3:      ", dim(mc3),      "\n")

# Save
save(rna_raw, meth_raw, cnv_raw, rppa_raw,
     clinical, survival, mc3, paradigm,
     file = file.path(proc_dir, "00_raw_loaded.RData"))
cat("Saved to processed/00_raw_loaded.RData\n")