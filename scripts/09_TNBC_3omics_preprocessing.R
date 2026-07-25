# ── 09_TNBC_3omics_preprocessing.R ───────────────────────────────────────────
# TNBC rerun with RNA + CNV + Methylation (n=81)
# Drops RPPA to maximise sample size
# ─────────────────────────────────────────────────────────────────────────────
# Install mofapy2 in the current Python environment
#reticulate::py_install("mofapy2", pip = TRUE)

library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")
mod_dir  <- file.path(base_dir, "outputs/models")

# Load raw data
load(file.path(proc_dir, "00_raw_loaded.RData"))

trim15 <- function(x) substr(x, 1, 15)
colnames(rna_raw)  <- trim15(colnames(rna_raw))
colnames(meth_raw) <- trim15(colnames(meth_raw))
colnames(cnv_raw)  <- trim15(colnames(cnv_raw))

# ── 1. TNBC SAMPLES ───────────────────────────────────────────────────────────
tnbc_ids <- clinical %>%
  filter(
    breast_carcinoma_estrogen_receptor_status == "Negative",
    breast_carcinoma_progesterone_receptor_status == "Negative",
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"
  ) %>%
  pull(sampleID) %>%
  substr(1, 15)

cat("Total TNBC in clinical:", length(tnbc_ids), "\n")

# ── 2. 3-OMICS INTERSECTION ───────────────────────────────────────────────────
common_tnbc_3 <- Reduce(intersect, list(
  colnames(rna_raw),
  colnames(meth_raw),
  colnames(cnv_raw),
  tnbc_ids
))
cat("TNBC samples in 3 omics:", length(common_tnbc_3), "\n")

rna_tnbc3  <- rna_raw[,  common_tnbc_3]
meth_tnbc3 <- meth_raw[, common_tnbc_3]
cnv_tnbc3  <- cnv_raw[,  common_tnbc_3]

# ── 3. FEATURE FILTERING ──────────────────────────────────────────────────────
cat("Filtering features...\n")

rna_tnbc3  <- rna_tnbc3[rowMeans(is.na(rna_tnbc3)) < 0.2, ]
rna_var    <- apply(rna_tnbc3, 1, var, na.rm = TRUE)
rna_top3   <- rna_tnbc3[order(rna_var, decreasing = TRUE)[1:3000], ]
cat("RNA features:", nrow(rna_top3), "\n")

cat("Computing methylation variance...\n")
meth_tnbc3 <- meth_tnbc3[rowMeans(is.na(meth_tnbc3)) < 0.2, ]
meth_var   <- apply(meth_tnbc3, 1, var, na.rm = TRUE)
meth_top3  <- meth_tnbc3[order(meth_var, decreasing = TRUE)[1:3000], ]
cat("Methylation features:", nrow(meth_top3), "\n")

cnv_clean3 <- cnv_tnbc3[complete.cases(cnv_tnbc3), ]
cat("CNV features:", nrow(cnv_clean3), "\n")

# ── 4. MOFA2 SETUP ────────────────────────────────────────────────────────────
mofa_data_tnbc3 <- list(
  RNA         = rna_top3,
  Methylation = meth_top3,
  CNV         = cnv_clean3
)

mofa_tnbc3 <- create_mofa(mofa_data_tnbc3)
plot_data_overview(mofa_tnbc3)

data_opts3  <- get_default_data_options(mofa_tnbc3)
model_opts3 <- get_default_model_options(mofa_tnbc3)
train_opts3 <- get_default_training_options(mofa_tnbc3)

model_opts3$num_factors      <- 12
train_opts3$seed             <- 42
train_opts3$maxiter          <- 1000
train_opts3$convergence_mode <- "fast"

mofa_tnbc3 <- prepare_mofa(
  mofa_tnbc3,
  data_options     = data_opts3,
  model_options    = model_opts3,
  training_options = train_opts3
)

cat("Views:   ", views_names(mofa_tnbc3), "\n")
cat("Samples: ", get_dimensions(mofa_tnbc3)[["N"]], "\n")
cat("Factors: ", model_opts3$num_factors, "\n")

# ── 5. TRAIN ──────────────────────────────────────────────────────────────────
cat("\nTraining TNBC 3-omics model...\n")

mofa_trained_tnbc3 <- run_mofa(
  mofa_tnbc3,
  use_basilisk = FALSE
)

save(mofa_trained_tnbc3, common_tnbc_3,
     rna_top3, meth_top3, cnv_clean3,
     file = file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
cat(" TNBC 3-omics model saved!\n")

# ── 6. QUICK RESULTS ──────────────────────────────────────────────────────────
p_var <- plot_variance_explained(mofa_trained_tnbc3, max_r2 = 15)
print(p_var)
ggsave(file.path(fig_dir, "25_tnbc3_variance_explained.png"),
       p_var, width = 8, height = 5, dpi = 150)

p_var_total <- plot_variance_explained(mofa_trained_tnbc3, plot_total = TRUE)
print(p_var_total)
ggsave(file.path(fig_dir, "26_tnbc3_variance_total.png"),
       p_var_total, width = 6, height = 4, dpi = 150)

# Factor scores
factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
p_f12 <- as.data.frame(factor_scores3) %>%
  ggplot(aes(Factor1, Factor2)) +
  geom_point(alpha = 0.7, size = 2.5, color = "#D85A30") +
  labs(title    = "MOFA2 Factor 1 vs Factor 2 — TNBC (n=81)",
       subtitle = "3 omics: RNA + Methylation + CNV") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_f12)
ggsave(file.path(fig_dir, "27_tnbc3_factor_scores.png"),
       p_f12, width = 7, height = 6, dpi = 150)

cat("\n Script 09 complete!\n")