library(tidyverse)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "01_TNBC_clean.RData"))

f <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/raw", fn))

# Load immune signatures
immune_sigs <- fread(
  file = f("TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena"),
  header = TRUE)

immune_mat <- immune_sigs %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  as.matrix()

immune_t <- t(immune_mat)
rownames(immune_t) <- substr(rownames(immune_t), 1, 15)

# Get factor scores
factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]

# Common samples
common_imm3 <- intersect(rownames(factor_scores3), rownames(immune_t))

# ── PREDEFINED IMMUNE COMPOSITE ───────────────────────────────────────────────
# Fix composite BEFORE looking at results
key_sigs <- c("Tcell_21978456", "MHC1_21978456",
              "CD8A", "IFNG_score_21050467")

composite <- rowMeans(
  immune_t[common_imm3, key_sigs],
  na.rm = TRUE
)

# ── OBSERVED: BEST FACTOR CORRELATION ─────────────────────────────────────────
obs_cors <- apply(factor_scores3[common_imm3, ], 2, function(f) {
  cor(f, composite, method = "pearson", use = "complete.obs")
})

obs_max <- max(abs(obs_cors))
obs_best_factor <- names(which.max(abs(obs_cors)))

cat("Observed correlations per factor:\n")
print(round(obs_cors, 3))
cat("\nBest factor:", obs_best_factor,
    "| r =", round(obs_cors[obs_best_factor], 3),
    "| max |r| =", round(obs_max, 3), "\n")

# ── PERMUTATION TEST ──────────────────────────────────────────────────────────
n_perm <- 1000
cat("\nRunning", n_perm, "permutations...\n")

set.seed(42)
perm_max <- numeric(n_perm)

for(i in 1:n_perm) {
  perm_composite <- sample(composite)
  names(perm_composite) <- names(composite)
  
  perm_cors <- apply(factor_scores3[common_imm3, ], 2, function(f) {
    cor(f, perm_composite, method = "pearson", use = "complete.obs")
  })
  
  perm_max[i] <- max(abs(perm_cors))
  
  if(i %% 100 == 0) cat("  Permutation", i, "\n")
}

# Empirical p-value
emp_p <- mean(perm_max >= obs_max)
cat("\n── Permutation test result ──\n")
cat("Observed max |r|:", round(obs_max, 3), "\n")
cat("Empirical p-value:", emp_p, "\n")
cat("(proportion of permutations with max |r| >=", round(obs_max, 3), ")\n")

# Plot null distribution
p_perm <- ggplot(data.frame(max_r = perm_max), aes(x = max_r)) +
  geom_histogram(bins = 50, fill = "#888780", alpha = 0.8) +
  geom_vline(xintercept = obs_max, color = "#D85A30",
             linewidth = 1, linetype = "solid") +
  annotate("text",
           x = obs_max + 0.01, y = 50,
           label = paste0("Observed r = ", round(obs_max, 3),
                          "\np = ", emp_p),
           hjust = 0, color = "#D85A30", size = 4) +
  labs(title    = "Permutation test: Factor 4 immune composite correlation",
       subtitle = paste0(n_perm, " permutations — maximum statistic null distribution"),
       x = "Max |r| across 12 factors (permuted)",
       y = "Count") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_perm)
ggsave(file.path(fig_dir, "52_permutation_test.png"),
       p_perm, width = 8, height = 5, dpi = 150)

# ── FULL 12x68 CORRELATION HEATMAP WITH FDR ───────────────────────────────────
cat("\nComputing full 12x68 correlation matrix...\n")

all_sigs <- rownames(immune_mat)
all_sigs <- intersect(all_sigs, colnames(immune_t))

cor_matrix <- matrix(NA, nrow = ncol(factor_scores3),
                     ncol = length(all_sigs))
rownames(cor_matrix) <- colnames(factor_scores3)
colnames(cor_matrix) <- all_sigs

for(sig in all_sigs) {
  for(fac in colnames(factor_scores3)) {
    cor_matrix[fac, sig] <- cor(
      factor_scores3[common_imm3, fac],
      immune_t[common_imm3, sig],
      use = "complete.obs"
    )
  }
}

# FDR correction across all cells
p_matrix <- matrix(NA, nrow = nrow(cor_matrix),
                   ncol = ncol(cor_matrix))
rownames(p_matrix) <- rownames(cor_matrix)
colnames(p_matrix) <- colnames(cor_matrix)

for(sig in all_sigs) {
  for(fac in colnames(factor_scores3)) {
    p_matrix[fac, sig] <- cor.test(
      factor_scores3[common_imm3, fac],
      immune_t[common_imm3, sig]
    )$p.value
  }
}

p_adj_matrix <- matrix(
  p.adjust(as.vector(p_matrix), method = "fdr"),
  nrow = nrow(p_matrix)
)
rownames(p_adj_matrix) <- rownames(p_matrix)
colnames(p_adj_matrix) <- colnames(p_matrix)

cat("\nFactor 4 significant immune signatures (FDR < 0.05):\n")
f4_sig <- which(p_adj_matrix["Factor4", ] < 0.05)
print(data.frame(
  signature = names(f4_sig),
  r         = round(cor_matrix["Factor4", names(f4_sig)], 3),
  p_adj     = round(p_adj_matrix["Factor4", names(f4_sig)], 4)
) %>% arrange(p_adj))

save(cor_matrix, p_matrix, p_adj_matrix,
     obs_max, perm_max, emp_p,
     file = file.path(proc_dir, "18_permutation_results.RData"))

cat("\n Permutation test complete!\n")