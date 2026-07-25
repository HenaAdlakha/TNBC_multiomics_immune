# ── 07_full_cohort_preprocessing.R ───────────────────────────────────────────
# Full BRCA cohort rerun — all subtypes, no TNBC filter
# Goals:
#   1. Load all 1,200 BRCA samples across 4 omics
#   2. Align sample IDs
#   3. Feature selection
#   4. Add subtype labels as metadata
#   5. QC plots
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
raw_dir  <- file.path(base_dir, "data/raw")
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

# Helper for paths with spaces
f <- function(fn) path.expand(file.path(raw_dir, fn))

# ── 1. LOAD RAW DATA ──────────────────────────────────────────────────────────
cat("Loading RNA-seq...\n")
rna_raw <- fread(file = f("HiSeqV2"), header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("sample") %>%
  as.matrix()

cat("Loading methylation...\n")
meth_raw <- fread(file = f("HumanMethylation450"), header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("sample") %>%
  as.matrix()

cat("Loading CNV...\n")
cnv_raw <- fread(file = f("Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"),
                 header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("Gene Symbol") %>%
  as.matrix()

cat("Loading RPPA...\n")
rppa_raw <- fread(file = f("RPPA_RBN"), header = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("Sample_description") %>%
  as.matrix()

cat("Loading clinical...\n")
clinical <- fread(file = f("TCGA.BRCA.sampleMap-BRCA_clinicalMatrix"),
                  header = TRUE)

cat("Loading survival...\n")
survival <- fread(file = f("survival-BRCA_survival.txt"), header = TRUE)

# ── 2. TRIM SAMPLE IDs ────────────────────────────────────────────────────────
trim15 <- function(x) substr(x, 1, 15)

colnames(rna_raw)  <- trim15(colnames(rna_raw))
colnames(meth_raw) <- trim15(colnames(meth_raw))
colnames(cnv_raw)  <- trim15(colnames(cnv_raw))
colnames(rppa_raw) <- trim15(colnames(rppa_raw))

# ── 3. FIND COMMON SAMPLES ACROSS ALL 4 OMICS ─────────────────────────────────
common_all <- Reduce(intersect, list(
  colnames(rna_raw),
  colnames(meth_raw),
  colnames(cnv_raw),
  colnames(rppa_raw)
))
cat("Samples in all 4 omics:", length(common_all), "\n")

# Subset
rna_all  <- rna_raw[,  common_all]
meth_all <- meth_raw[, common_all]
cnv_all  <- cnv_raw[,  common_all]
rppa_all <- rppa_raw[, common_all]

# ── 4. CREATE SUBTYPE LABELS ──────────────────────────────────────────────────
cat("Creating subtype labels...\n")

subtype_df <- clinical %>%
  mutate(sampleID = trim15(sampleID)) %>%
  filter(sampleID %in% common_all) %>%
  mutate(
    subtype = case_when(
      breast_carcinoma_estrogen_receptor_status == "Negative" &
        breast_carcinoma_progesterone_receptor_status == "Negative" &
        lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"
      ~ "TNBC",
      lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive"
      ~ "HER2+",
      breast_carcinoma_estrogen_receptor_status == "Positive"
      ~ "ER+",
      TRUE ~ "Other"
    )
  ) %>%
  select(sampleID, subtype)

cat("Subtype distribution:\n")
print(table(subtype_df$subtype))

# ── 5. FEATURE FILTERING ──────────────────────────────────────────────────────
cat("\nFiltering features...\n")

# RNA: top 5000 by variance
rna_var <- apply(rna_all, 1, var, na.rm = TRUE)
rna_all <- rna_all[rowMeans(is.na(rna_all)) < 0.2, ]
rna_top <- rna_all[order(apply(rna_all, 1, var, na.rm=TRUE),
                         decreasing = TRUE)[1:5000], ]
cat("RNA features:", nrow(rna_top), "\n")

# Methylation: top 5000 by variance
meth_all <- meth_all[rowMeans(is.na(meth_all)) < 0.2, ]
meth_top <- meth_all[order(apply(meth_all, 1, var, na.rm=TRUE),
                           decreasing = TRUE)[1:5000], ]
cat("Methylation features:", nrow(meth_top), "\n")

# CNV: complete cases
cnv_clean <- cnv_all[complete.cases(cnv_all), ]
cat("CNV features:", nrow(cnv_clean), "\n")

# RPPA: complete cases
rppa_clean <- rppa_all[complete.cases(rppa_all), ]
cat("RPPA proteins:", nrow(rppa_clean), "\n")

# ── 6. QC PLOTS ───────────────────────────────────────────────────────────────

# Sample counts by subtype
p_subtype <- subtype_df %>%
  count(subtype) %>%
  ggplot(aes(x = reorder(subtype, -n), y = n, fill = subtype)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("TNBC"  = "#D85A30",
                               "HER2+" = "#534AB7",
                               "ER+"   = "#1D9E75",
                               "Other" = "#888780")) +
  geom_text(aes(label = n), vjust = -0.4, size = 4.5, fontface = "bold") +
  labs(title = "BRCA full cohort — sample counts by subtype",
       y = "Number of samples", x = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_subtype)
ggsave(file.path(fig_dir, "20_full_cohort_subtype_counts.png"),
       p_subtype, width = 7, height = 4, dpi = 150)

# PCA colored by subtype
plot_pca_subtype <- function(mat, title) {
  mat_t <- t(mat)
  mat_t[is.na(mat_t)] <- 0
  pca <- prcomp(mat_t, scale. = TRUE, center = TRUE)
  df  <- as.data.frame(pca$x[, 1:2])
  df$sample <- rownames(df)
  ve  <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  
  df <- df %>%
    left_join(subtype_df, by = c("sample" = "sampleID"))
  
  ggplot(df, aes(PC1, PC2, color = subtype)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("TNBC"  = "#D85A30",
                                  "HER2+" = "#534AB7",
                                  "ER+"   = "#1D9E75",
                                  "Other" = "#888780")) +
    labs(title = title,
         x = paste0("PC1 (", ve[1], "%)"),
         y = paste0("PC2 (", ve[2], "%)"),
         color = "Subtype") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

p_pca_rna  <- plot_pca_subtype(rna_top,  "RNA-seq — full BRCA cohort")
p_pca_meth <- plot_pca_subtype(meth_top, "Methylation — full BRCA cohort")

p_pca_full <- (p_pca_rna | p_pca_meth) +
  plot_annotation(
    title = "PCA by subtype — full BRCA cohort",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(p_pca_full)
ggsave(file.path(fig_dir, "21_full_cohort_pca.png"),
       p_pca_full, width = 12, height = 5, dpi = 150)

# ── 7. SAVE ───────────────────────────────────────────────────────────────────
save(rna_top, meth_top, cnv_clean, rppa_clean,
     common_all, subtype_df, clinical, survival,
     file = file.path(proc_dir, "07_full_cohort_clean.RData"))

cat("\n Full cohort preprocessing complete!\n")
cat("Samples:", length(common_all), "\n")
cat("Figures saved to outputs/figures/\n")
cat("Subtype distribution:\n")
print(table(subtype_df$subtype))