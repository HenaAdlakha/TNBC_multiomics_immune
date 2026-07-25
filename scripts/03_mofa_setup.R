library(tidyverse)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")
mod_dir  <- file.path(base_dir, "outputs/models")

load(file.path(proc_dir, "01_TNBC_clean.RData"))

# ── BUILD MOFA OBJECT ─────────────────────────────────────────────────────────
mofa_data <- list(
  RNA         = rna_top,
  Methylation = meth_top,
  CNV         = cnv_clean,
  Protein     = rppa_clean
)

mofa_obj <- create_mofa(mofa_data)

# Print summary — great to show in meeting
print(mofa_obj)

# This plot is your money slide for the meeting
p_overview <- plot_data_overview(mofa_obj)
print(p_overview)
ggsave(file.path(fig_dir, "04_mofa_data_overview.png"),
       p_overview, width = 8, height = 5, dpi = 150)

# ── SET OPTIONS ───────────────────────────────────────────────────────────────
data_opts  <- get_default_data_options(mofa_obj)
model_opts <- get_default_model_options(mofa_obj)
train_opts <- get_default_training_options(mofa_obj)

model_opts$num_factors          <- 15
train_opts$seed                 <- 42
train_opts$maxiter              <- 1000
train_opts$convergence_mode     <- "fast"

mofa_obj <- prepare_mofa(
  mofa_obj,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

cat("\nMOFA object ready!\n")
cat("Views:   ", views_names(mofa_obj), "\n")
cat("Samples: ", get_dimensions(mofa_obj)[["N"]], "\n")
cat("Factors: ", model_opts$num_factors, "\n")

save(mofa_obj,
     file = file.path(proc_dir, "02_mofa_object.RData"))
cat("Saved — ready to train\n")