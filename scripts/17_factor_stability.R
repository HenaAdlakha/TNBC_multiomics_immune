library(tidyverse)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "01_TNBC_clean.RData"))
load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))

# Reference Factor 4 scores
f4_ref <- get_factors(mofa_trained_tnbc3, 
                      factors = "all")[[1]][, "Factor4"]

# ── TEST ACROSS SEEDS ─────────────────────────────────────────────────────────
seeds <- c(1, 7, 21, 42, 99, 123, 456, 999)
n_factors <- 12

stability_results <- list()

for(seed in seeds) {
  cat("Training with seed:", seed, "\n")
  
  mofa_tmp <- create_mofa(list(
    RNA         = rna_top3,
    Methylation = meth_top3,
    CNV         = cnv_clean3
  ))
  
  data_opts  <- get_default_data_options(mofa_tmp)
  model_opts <- get_default_model_options(mofa_tmp)
  train_opts <- get_default_training_options(mofa_tmp)
  
  model_opts$num_factors      <- n_factors
  train_opts$seed             <- seed
  train_opts$maxiter          <- 1000
  train_opts$convergence_mode <- "fast"
  
  mofa_tmp <- prepare_mofa(mofa_tmp,
                           data_options     = data_opts,
                           model_options    = model_opts,
                           training_options = train_opts)
  
  mofa_trained_tmp <- run_mofa(mofa_tmp, use_basilisk = FALSE)
  
  # Get all factor scores
  scores_tmp <- get_factors(mofa_trained_tmp, factors = "all")[[1]]
  
  # Find factor most correlated with reference Factor 4
  common_s <- intersect(names(f4_ref), rownames(scores_tmp))
  cors <- apply(scores_tmp[common_s, ], 2, function(f) {
    cor(f4_ref[common_s], f, use = "complete.obs")
  })
  
  best_factor <- names(which.max(abs(cors)))
  best_r      <- max(abs(cors))
  
  stability_results[[as.character(seed)]] <- list(
    seed        = seed,
    best_factor = best_factor,
    max_r       = round(best_r, 3)
  )
  
  cat("  Best matching factor:", best_factor, 
      "| r =", round(best_r, 3), "\n")
}

# ── TEST ACROSS FACTOR NUMBERS ────────────────────────────────────────────────
factor_nums <- c(8, 10, 12, 15)
factor_stability <- list()

for(nf in factor_nums) {
  cat("Training with", nf, "factors...\n")
  
  mofa_tmp <- create_mofa(list(
    RNA         = rna_top3,
    Methylation = meth_top3,
    CNV         = cnv_clean3
  ))
  
  data_opts  <- get_default_data_options(mofa_tmp)
  model_opts <- get_default_model_options(mofa_tmp)
  train_opts <- get_default_training_options(mofa_tmp)
  
  model_opts$num_factors      <- nf
  train_opts$seed             <- 42
  train_opts$maxiter          <- 1000
  train_opts$convergence_mode <- "fast"
  
  mofa_tmp <- prepare_mofa(mofa_tmp,
                           data_options     = data_opts,
                           model_options    = model_opts,
                           training_options = train_opts)
  
  mofa_trained_tmp <- run_mofa(mofa_tmp, use_basilisk = FALSE)
  
  scores_tmp <- get_factors(mofa_trained_tmp, factors = "all")[[1]]
  common_s   <- intersect(names(f4_ref), rownames(scores_tmp))
  
  cors <- apply(scores_tmp[common_s, ], 2, function(f) {
    cor(f4_ref[common_s], f, use = "complete.obs")
  })
  
  best_factor <- names(which.max(abs(cors)))
  best_r      <- max(abs(cors))
  
  factor_stability[[as.character(nf)]] <- list(
    n_factors   = nf,
    best_factor = best_factor,
    max_r       = round(best_r, 3)
  )
  
  cat("  Best matching factor:", best_factor,
      "| r =", round(best_r, 3), "\n")
}

# ── SUMMARISE ─────────────────────────────────────────────────────────────────
seed_df <- do.call(rbind, lapply(stability_results, as.data.frame))
cat("\n── Stability across seeds ──\n")
print(seed_df)

factor_df <- do.call(rbind, lapply(factor_stability, as.data.frame))
cat("\n── Stability across factor numbers ──\n")
print(factor_df)

# Plot
p_seed <- ggplot(seed_df, aes(x = factor(seed), y = max_r)) +
  geom_col(fill = "#534AB7", width = 0.6) +
  geom_hline(yintercept = 0.8, linetype = "dashed",
             color = "#D85A30", linewidth = 0.8) +
  annotate("text", x = 1, y = 0.82,
           label = "r = 0.8 stability threshold",
           hjust = 0, size = 3.5, color = "#D85A30") +
  labs(title    = "Factor 4 stability across random seeds",
       subtitle = "r = max correlation with reference Factor 4 (seed=42)",
       x = "Random seed", y = "Max Pearson r") +
  ylim(0, 1) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_seed)
ggsave(file.path(fig_dir, "50_factor_stability_seeds.png"),
       p_seed, width = 8, height = 5, dpi = 150)

p_factors <- ggplot(factor_df, aes(x = factor(n_factors), y = max_r)) +
  geom_col(fill = "#1D9E75", width = 0.6) +
  geom_hline(yintercept = 0.8, linetype = "dashed",
             color = "#D85A30", linewidth = 0.8) +
  labs(title    = "Factor 4 stability across factor numbers",
       subtitle = "r = max correlation with reference Factor 4 (12 factors)",
       x = "Number of factors", y = "Max Pearson r") +
  ylim(0, 1) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_factors)
ggsave(file.path(fig_dir, "51_factor_stability_nfactors.png"),
       p_factors, width = 7, height = 5, dpi = 150)

save(seed_df, factor_df,
     file = file.path(proc_dir, "17_factor_stability.RData"))

cat("\n Stability analysis complete!\n")