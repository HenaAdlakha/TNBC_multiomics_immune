# ── 10_TNBC_3omics_immune.R ──────────────────────────────────────────────────
# Immune characterization on expanded TNBC cohort (n=81, 3 omics)
# Goals:
#   1. Identify which factor captures immune signal in n=81 model
#   2. Correlate all factors with 68 immune signatures
#   3. Compare with original n=55 findings
#   4. Set up for survival analysis in script 11
# Input:  09_mofa_tnbc3_trained.RData, immune signatures
# Output: figures 28-30, immune factor identified
# ─────────────────────────────────────────────────────────────────────────────
# ── immune characterization on new n=81 model ──────────────────────────

library(tidyverse)
library(data.table)
library(MOFA2)
library(ggplot2)
library(survival)
library(survminer)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "01_TNBC_clean.RData"))

# Helper
f <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/raw", fn))

# ── 1. GET FACTOR SCORES ──────────────────────────────────────────────────────
factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
cat("Factor matrix dimensions:", dim(factor_scores3), "\n")

# ── 2. LOAD IMMUNE SIGNATURES ─────────────────────────────────────────────────
immune_sigs <- fread(
  file = f("TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena"),
  header = TRUE)

immune_mat <- immune_sigs %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  as.matrix()

immune_t <- t(immune_mat)
rownames(immune_t) <- substr(rownames(immune_t), 1, 15)

# Filter to n=81 TNBC samples
common_imm3 <- intersect(rownames(immune_t), common_tnbc_3)
cat("TNBC samples with immune scores:", length(common_imm3), "\n")
immune_tnbc3 <- immune_t[common_imm3, ]

# ── 3. CORRELATE ALL FACTORS WITH KEY IMMUNE SIGNATURES ──────────────────────
key_sigs <- c("Tcell_21978456", "MHC1_21978456", "MHC2_21978456",
              "CD8A", "CD68", "IFNG_score_21050467",
              "TGFB_score_21050467", "PD1_PDL1_score",
              "Module4_TcellBcell_score", "CD8_CD68_ratio")

cat("\n── Correlations of ALL factors with MHC1 ──\n")
for(i in 1:ncol(factor_scores3)) {
  common_f <- intersect(rownames(factor_scores3), common_imm3)
  r <- round(cor(factor_scores3[common_f, i],
                 immune_tnbc3[common_f, "MHC1_21978456"],
                 method = "pearson", use = "complete.obs"), 3)
  cat("Factor", i, ": r =", r, "\n")
}

cat("\n── Correlations of ALL factors with Tcell ──\n")
for(i in 1:ncol(factor_scores3)) {
  common_f <- intersect(rownames(factor_scores3), common_imm3)
  r <- round(cor(factor_scores3[common_f, i],
                 immune_tnbc3[common_f, "Tcell_21978456"],
                 method = "pearson", use = "complete.obs"), 3)
  cat("Factor", i, ": r =", r, "\n")
}

#-----
  # Check Factor 4 against all key signatures
  cat("── Factor 4 correlations with all immune signatures ──\n")
common_f <- intersect(rownames(factor_scores3), common_imm3)
f4 <- factor_scores3[common_f, "Factor4"]

cor_f4 <- data.frame(
  signature = key_sigs,
  r = sapply(key_sigs, function(s) {
    round(cor(f4, immune_tnbc3[common_f, s],
              method = "pearson", use = "complete.obs"), 3)
  }),
  p_value = sapply(key_sigs, function(s) {
    round(cor.test(f4, immune_tnbc3[common_f, s])$p.value, 4)
  })
) %>% arrange(r)

print(cor_f4)

# Also check Factor 10
cat("\n── Factor 10 correlations with all immune signatures ──\n")
f10 <- factor_scores3[common_f, "Factor10"]

cor_f10 <- data.frame(
  signature = key_sigs,
  r = sapply(key_sigs, function(s) {
    round(cor(f10, immune_tnbc3[common_f, s],
              method = "pearson", use = "complete.obs"), 3)
  }),
  p_value = sapply(key_sigs, function(s) {
    round(cor.test(f10, immune_tnbc3[common_f, s])$p.value, 4)
  })
) %>% arrange(r)

print(cor_f10)

# Also build composite immune score and check all factors
cat("\n── All factors vs composite immune score ──\n")
composite <- rowMeans(immune_tnbc3[common_f, 
                                   c("Tcell_21978456", "MHC1_21978456", 
                                     "CD8A", "IFNG_score_21050467")],
                      na.rm = TRUE)

for(i in 1:ncol(factor_scores3)) {
  r <- round(cor(factor_scores3[common_f, i], composite,
                 method = "pearson", use = "complete.obs"), 3)
  p <- round(cor.test(factor_scores3[common_f, i], composite)$p.value, 4)
  cat("Factor", i, ": r =", r, "  p =", p, "\n")
}