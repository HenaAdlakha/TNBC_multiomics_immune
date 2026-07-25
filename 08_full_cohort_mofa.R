# ── 08_full_cohort_mofa.R ─────────────────────────────────────────────────────
# MOFA2 on full BRCA cohort (n=772, 3 omics: RNA + Methylation + CNV)
# Goals:
#   1. Build and train MOFA2 on full cohort
#   2. Check variance explained
#   3. Correlate factors with subtype labels
#   4. Identify immune factor in larger cohort
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(MOFA2)
library(ggplot2)
library(patchwork)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
raw_dir  <- file.path(base_dir, "data/raw")
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")
mod_dir  <- file.path(base_dir, "outputs/models")

load(file.path(proc_dir, "07_full_cohort_clean.RData"))

# ── 1. SUBSET TO 3-OMICS COMMON SAMPLES ──────────────────────────────────────
f <- function(fn) path.expand(file.path(raw_dir, fn))

# We need to reload rna/meth/cnv since 07 saved the RPPA-filtered versions
cat("Reloading raw data for 3-omics intersection...\n")
load(file.path(proc_dir, "00_raw_loaded.RData"))

trim15 <- function(x) substr(x, 1, 15)
colnames(rna_raw)  <- trim15(colnames(rna_raw))
colnames(meth_raw) <- trim15(colnames(meth_raw))
colnames(cnv_raw)  <- trim15(colnames(cnv_raw))

common_3omics <- Reduce(intersect, list(
  colnames(rna_raw),
  colnames(meth_raw),
  colnames(cnv_raw)
))
cat("Samples in 3-omics intersection:", length(common_3omics), "\n")

# Subset
rna_full  <- rna_raw[,  common_3omics]
meth_full <- meth_raw[, common_3omics]
cnv_full  <- cnv_raw[,  common_3omics]

# ── 2. FEATURE FILTERING ──────────────────────────────────────────────────────
cat("Filtering features...\n")

# RNA: remove >20% NA, top 5000 by variance
rna_full <- rna_full[rowMeans(is.na(rna_full)) < 0.2, ]
rna_var  <- apply(rna_full, 1, var, na.rm = TRUE)
rna_top  <- rna_full[order(rna_var, decreasing = TRUE)[1:5000], ]
cat("RNA features:", nrow(rna_top), "\n")

# Methylation: remove >20% NA, top 5000 by variance
cat("Computing methylation variance (this will take a few minutes)...\n")
meth_full <- meth_full[rowMeans(is.na(meth_full)) < 0.2, ]
meth_var  <- apply(meth_full, 1, var, na.rm = TRUE)
meth_top  <- meth_full[order(meth_var, decreasing = TRUE)[1:5000], ]
cat("Methylation features:", nrow(meth_top), "\n")

# CNV: complete cases
cnv_clean <- cnv_full[complete.cases(cnv_full), ]
cat("CNV features:", nrow(cnv_clean), "\n")

# ── 3. BUILD MOFA OBJECT ──────────────────────────────────────────────────────
cat("Building MOFA object...\n")

mofa_data_full <- list(
  RNA         = rna_top,
  Methylation = meth_top,
  CNV         = cnv_clean
)

mofa_full <- create_mofa(mofa_data_full)
plot_data_overview(mofa_full)

# ── 4. SET OPTIONS ────────────────────────────────────────────────────────────
data_opts  <- get_default_data_options(mofa_full)
model_opts <- get_default_model_options(mofa_full)
train_opts <- get_default_training_options(mofa_full)

# With 772 samples we can use more factors
model_opts$num_factors      <- 15
train_opts$seed             <- 42
train_opts$maxiter          <- 1000
train_opts$convergence_mode <- "fast"

mofa_full <- prepare_mofa(
  mofa_full,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

cat("MOFA object ready:\n")
cat("Views:   ", views_names(mofa_full), "\n")
cat("Samples: ", get_dimensions(mofa_full)[["N"]], "\n")
cat("Factors: ", model_opts$num_factors, "\n")

# ── 5. TRAIN ──────────────────────────────────────────────────────────────────
cat("\nStarting training — expected 20-40 minutes...\n")

mofa_trained_full <- run_mofa(
  mofa_full,
  use_basilisk = FALSE
)

# Save immediately
save(mofa_trained_full,
     file = file.path(proc_dir, "08_mofa_full_trained.RData"))
cat(" Model saved!\n")

# ── 6. FIRST LOOK ─────────────────────────────────────────────────────────────

# Variance explained
p_var <- plot_variance_explained(mofa_trained_full, max_r2 = 15)
print(p_var)
ggsave(file.path(fig_dir, "22_full_cohort_variance_explained.png"),
       p_var, width = 8, height = 5, dpi = 150)

p_var_total <- plot_variance_explained(mofa_trained_full, plot_total = TRUE)
print(p_var_total)
ggsave(file.path(fig_dir, "23_full_cohort_variance_total.png"),
       p_var_total, width = 6, height = 4, dpi = 150)

# Factor scores colored by subtype
factor_scores_full <- get_factors(mofa_trained_full, factors = "all")[[1]]

# Add subtype info
subtype_aligned <- subtype_df %>%
  filter(sampleID %in% rownames(factor_scores_full))

factor_df <- as.data.frame(factor_scores_full) %>%
  rownames_to_column("sample") %>%
  left_join(subtype_df, by = c("sample" = "sampleID"))

# Plot Factor 1 vs Factor 2 colored by subtype
p_f1f2 <- ggplot(factor_df, aes(Factor1, Factor2, color = subtype)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("TNBC"  = "#D85A30",
                                "HER2+" = "#534AB7",
                                "ER+"   = "#1D9E75",
                                "Other" = "#888780"),
                     na.value = "#888780") +
  labs(title    = "MOFA2 Factor 1 vs Factor 2 — full BRCA cohort",
       subtitle = "Colored by breast cancer subtype",
       color    = "Subtype") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_f1f2)
ggsave(file.path(fig_dir, "24_full_cohort_factor_scores.png"),
       p_f1f2, width = 8, height = 6, dpi = 150)

# ── 7. CORRELATE ALL FACTORS WITH SUBTYPE ─────────────────────────────────────
cat("\nCorrelating factors with TNBC subtype...\n")

# Binary TNBC label
tnbc_label <- ifelse(factor_df$subtype == "TNBC", 1, 0)
tnbc_label[is.na(tnbc_label)] <- 0

factor_cols <- grep("Factor", colnames(factor_df), value = TRUE)

tnbc_cors <- sapply(factor_cols, function(fc) {
  cor(factor_df[[fc]], tnbc_label,
      method = "pearson", use = "complete.obs")
})

cat("\n── Factor correlations with TNBC label ──\n")
tnbc_cor_df <- data.frame(
  factor = factor_cols,
  r      = round(tnbc_cors, 3)
) %>% arrange(desc(abs(r)))

print(tnbc_cor_df)

cat("\n Script 08 complete!\n")
cat("Check which factor correlates most strongly with TNBC\n")

#---
# What is Factor 1 driven by?
weights_f1 <- get_weights(mofa_trained_full,
                          factors = "Factor1",
                          as.data.frame = TRUE)

# Variance explained by Factor 1
var_exp_full <- calculate_variance_explained(mofa_trained_full)
cat("── Factor 1 variance explained per omics ──\n")
print(var_exp_full$r2_per_factor[[1]]["Factor1",])

# Top RNA genes
cat("\n── Top RNA genes in Factor 1 ──\n")
weights_f1 %>%
  filter(view == "RNA") %>%
  arrange(desc(abs(value))) %>%
  head(15) %>%
  print()

# Correlate Factor 1 with immune signatures in full cohort
cat("\n── Factor 1 vs immune signatures (full cohort) ──\n")
immune_t_full <- immune_t
common_full_imm <- intersect(rownames(factor_scores_full),
                             rownames(immune_t_full))
cat("Samples with immune scores:", length(common_full_imm), "\n")

key_sigs <- c("Tcell_21978456", "MHC1_21978456", "MHC2_21978456",
              "CD8A", "CD68", "IFNG_score_21050467",
              "PD1_PDL1_score", "Module4_TcellBcell_score")

f1_full <- factor_scores_full[common_full_imm, "Factor1"]

for(sig in key_sigs) {
  r <- round(cor(f1_full,
                 immune_t_full[common_full_imm, sig],
                 method = "pearson", use = "complete.obs"), 3)
  p <- round(cor.test(f1_full,
                      immune_t_full[common_full_imm, sig])$p.value, 4)
  cat(sig, ": r =", r, "  p =", p, "\n")
}

save(mofa_trained_full, factor_scores_full, 
     factor_df, tnbc_cor_df,
     file = file.path(proc_dir, "11_full_cohort_results.RData"))
cat(" Full cohort results saved!\n")