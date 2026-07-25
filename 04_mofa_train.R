library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
mod_dir  <- file.path(base_dir, "outputs/models")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "02_mofa_object_final.RData"))

# ── TRAIN ─────────────────────────────────────────────────────────────────────
cat("Starting MOFA2 training...\n")
cat("Expected time: 15-30 minutes\n")

mofa_trained <- run_mofa(
  mofa_obj2,
  outfile      = file.path(mod_dir, "mofa_tnbc_trained.hdf5"),
  use_basilisk = FALSE
)

cat("Training complete!\n")

# ── FIRST LOOK AT RESULTS ─────────────────────────────────────────────────────

# Variance explained per factor per omics — your key result plot
p_var1 <- plot_variance_explained(mofa_trained, max_r2 = 15)
print(p_var1)
ggsave(file.path(fig_dir, "05_variance_explained.png"),
       p_var1, width = 8, height = 5, dpi = 150)

# Total variance explained per omics layer
p_var2 <- plot_variance_explained(mofa_trained, plot_total = TRUE)
print(p_var2)
ggsave(file.path(fig_dir, "06_variance_explained_total.png"),
       p_var2, width = 6, height = 4, dpi = 150)

# Factor correlation plot — checks factors are not redundant
p_factors <- plot_factor_cor(mofa_trained)
print(p_factors)
ggsave(file.path(fig_dir, "07_factor_correlation.png"),
       p_factors, width = 6, height = 5, dpi = 150)

# ── SAVE ──────────────────────────────────────────────────────────────────────
save(mofa_trained,
     file = file.path(proc_dir, "03_mofa_trained.RData"))

cat("Model saved to outputs/models/mofa_tnbc_trained.hdf5\n")
cat("Figures saved to outputs/figures/\n")
cat("\nNext step: run 05_factor_characterization.R after training\n")