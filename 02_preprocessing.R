library(tidyverse)
library(data.table)
library(pheatmap)
library(ggplot2)
library(patchwork)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "00_raw_loaded.RData"))

# ── 1. IDENTIFY TNBC SAMPLES ──────────────────────────────────────────────────
tnbc_meta <- clinical %>%
  filter(
    breast_carcinoma_estrogen_receptor_status == "Negative",
    breast_carcinoma_progesterone_receptor_status == "Negative",
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"
  )

tnbc_ids <- tnbc_meta$sampleID
cat("TNBC samples found:", length(tnbc_ids), "\n")

# ── 2. HARMONISE SAMPLE IDs ───────────────────────────────────────────────────
trim15 <- function(x) substr(x, 1, 15)

colnames(rna_raw)  <- trim15(colnames(rna_raw))
colnames(meth_raw) <- trim15(colnames(meth_raw))
colnames(cnv_raw)  <- trim15(colnames(cnv_raw))
colnames(rppa_raw) <- trim15(colnames(rppa_raw))
colnames(mc3)      <- trim15(colnames(mc3))
tnbc_ids           <- trim15(tnbc_ids)

# Samples in all 4 main omics + TNBC
common_samples <- Reduce(intersect, list(
  colnames(rna_raw),
  colnames(meth_raw),
  colnames(cnv_raw),
  colnames(rppa_raw),
  tnbc_ids
))

cat("Samples in all 4 omics + TNBC:", length(common_samples), "\n")

# Subset
rna_tnbc  <- rna_raw[,  common_samples]
meth_tnbc <- meth_raw[, common_samples]
cnv_tnbc  <- cnv_raw[,  common_samples]
rppa_tnbc <- rppa_raw[, common_samples]

# MC3 separately — fewer samples
mc3_samples <- intersect(trim15(colnames(mc3)), common_samples)
mc3_tnbc    <- mc3[, mc3_samples]
cat("MC3 TNBC samples:", length(mc3_samples), "\n")

# ── 3. FEATURE FILTERING ──────────────────────────────────────────────────────

# RNA: remove >20% NA, top 5000 by variance
rna_tnbc <- rna_tnbc[rowMeans(is.na(rna_tnbc)) < 0.2, ]
rna_var  <- apply(rna_tnbc, 1, var, na.rm = TRUE)
rna_top  <- rna_tnbc[order(rna_var, decreasing = TRUE)[1:5000], ]
cat("RNA features:", nrow(rna_top), "\n")

# Methylation: remove >20% NA, top 5000 by variance
meth_tnbc <- meth_tnbc[rowMeans(is.na(meth_tnbc)) < 0.2, ]
meth_var  <- apply(meth_tnbc, 1, var, na.rm = TRUE)
meth_top  <- meth_tnbc[order(meth_var, decreasing = TRUE)[1:5000], ]
cat("Methylation features:", nrow(meth_top), "\n")

# CNV: complete cases only
cnv_clean  <- cnv_tnbc[complete.cases(cnv_tnbc), ]
cat("CNV features:", nrow(cnv_clean), "\n")

# RPPA: complete cases (~130 proteins, keep all)
rppa_clean <- rppa_tnbc[complete.cases(rppa_tnbc), ]
cat("RPPA proteins:", nrow(rppa_clean), "\n")

# ── 4. QC PLOTS ───────────────────────────────────────────────────────────────

# --- Sample counts ---
sample_counts <- data.frame(
  Layer   = c("RNA-seq", "Methylation", "CNV", "RPPA", "TNBC final"),
  Samples = c(ncol(rna_tnbc), ncol(meth_tnbc),
              ncol(cnv_tnbc), ncol(rppa_tnbc),
              length(common_samples))
)

p_counts <- ggplot(sample_counts, aes(x = Layer, y = Samples, fill = Layer)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("#1D9E75","#0F6E56",
                               "#534AB7","#D85A30","#085041")) +
  geom_text(aes(label = Samples), vjust = -0.4,
            size = 4.5, fontface = "bold") +
  labs(title = "TNBC sample counts across omics layers",
       y = "Number of samples", x = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_counts)
ggsave(file.path(fig_dir, "01_sample_counts.png"),
       p_counts, width = 7, height = 4, dpi = 150)

# --- PCA per omics ---
plot_pca <- function(mat, title, color = "#1D9E75") {
  mat_t <- t(mat)
  mat_t[is.na(mat_t)] <- 0
  pca  <- prcomp(mat_t, scale. = TRUE, center = TRUE)
  df   <- as.data.frame(pca$x[, 1:2])
  ve   <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  ggplot(df, aes(PC1, PC2)) +
    geom_point(alpha = 0.75, size = 2.5, color = color) +
    labs(title = title,
         x = paste0("PC1 (", ve[1], "%)"),
         y = paste0("PC2 (", ve[2], "%)")) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

p_rna  <- plot_pca(rna_top,    "RNA-seq",      "#1D9E75")
p_meth <- plot_pca(meth_top,   "Methylation",  "#534AB7")
p_cnv  <- plot_pca(cnv_clean,  "CNV",          "#D85A30")
p_rppa <- plot_pca(rppa_clean, "Protein/RPPA", "#0F6E56")

p_grid <- (p_rna | p_meth) / (p_cnv | p_rppa) +
  plot_annotation(
    title = "PCA per omics layer — TNBC samples",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(p_grid)
ggsave(file.path(fig_dir, "02_pca_per_omics.png"),
       p_grid, width = 10, height = 8, dpi = 150)

# --- Variance distribution ---
p_var <- data.frame(
  variance = c(rna_var[order(rna_var, decreasing=TRUE)[1:5000]],
               meth_var[order(meth_var, decreasing=TRUE)[1:5000]]),
  layer = rep(c("RNA-seq", "Methylation"), each = 5000)
) %>%
  ggplot(aes(x = variance, fill = layer)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("#1D9E75", "#534AB7")) +
  facet_wrap(~layer, scales = "free") +
  labs(title = "Variance distribution of selected features",
       x = "Variance", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

print(p_var)
ggsave(file.path(fig_dir, "03_variance_distribution.png"),
       p_var, width = 9, height = 4, dpi = 150)

# ── 5. SAVE ───────────────────────────────────────────────────────────────────
save(rna_top, meth_top, cnv_clean, rppa_clean,
     mc3_tnbc, mc3_samples,
     clinical, survival, tnbc_meta,
     common_samples,
     file = file.path(proc_dir, "01_TNBC_clean.RData"))

cat("\n Script 02 complete!\n")
cat("TNBC samples in final dataset:", length(common_samples), "\n")
cat("Figures saved to outputs/figures/\n")