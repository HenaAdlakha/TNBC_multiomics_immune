library(tidyverse)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "01_TNBC_clean.RData"))

f <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/raw", fn))

# ── LOAD PURITY ESTIMATES ─────────────────────────────────────────────────────
# CD45/PTPRC as proxy for leukocyte fraction from RNA
# Also use clinical purity if available

factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
f4 <- factor_scores3[, "Factor4"]

# CD45 as leukocyte fraction proxy
ptprc_genes <- rownames(rna_top3)[grep("PTPRC|CD45",
                                       rownames(rna_top3),
                                       ignore.case = TRUE)]
cat("PTPRC/CD45 genes found:", ptprc_genes, "\n")

common_f <- intersect(names(f4), colnames(rna_top3))

if(length(ptprc_genes) > 0) {
  ptprc <- colMeans(rna_top3[ptprc_genes, common_f, drop=FALSE],
                    na.rm = TRUE)
} else {
  cat("PTPRC not in top variable genes — checking full RNA matrix\n")
  load(file.path(proc_dir, "00_raw_loaded.RData"))
  colnames(rna_raw) <- substr(colnames(rna_raw), 1, 15)
  ptprc_full <- rownames(rna_raw)[grep("PTPRC",
                                       rownames(rna_raw),
                                       ignore.case = TRUE)]
  cat("PTPRC in full matrix:", ptprc_full, "\n")
  ptprc <- rna_raw[ptprc_full[1], common_f]
}

# ── LOAD IMMUNE SIGNATURES ────────────────────────────────────────────────────
immune_sigs <- fread(
  file = f("TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena"),
  header = TRUE)

immune_mat <- immune_sigs %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  as.matrix()

immune_t <- t(immune_mat)
rownames(immune_t) <- substr(rownames(immune_t), 1, 15)

common_all <- Reduce(intersect, list(
  names(f4), names(ptprc),
  rownames(immune_t)
))
cat("Samples for purity analysis:", length(common_all), "\n")

tcell  <- immune_t[common_all, "Tcell_21978456"]
mhc1   <- immune_t[common_all, "MHC1_21978456"]
f4_sub <- f4[common_all]
pur    <- ptprc[common_all]

# ── CORRELATIONS ──────────────────────────────────────────────────────────────
cat("\n── Simple correlations ──\n")
cat("Factor 4 vs T cell:  r =",
    round(cor(f4_sub, tcell, use="complete.obs"), 3), "\n")
cat("Factor 4 vs MHC1:    r =",
    round(cor(f4_sub, mhc1,  use="complete.obs"), 3), "\n")
cat("Factor 4 vs PTPRC:   r =",
    round(cor(f4_sub, pur,   use="complete.obs"), 3), "\n")
cat("PTPRC vs T cell:     r =",
    round(cor(pur, tcell,    use="complete.obs"), 3), "\n")

# ── PARTIAL CORRELATIONS ──────────────────────────────────────────────────────
# Factor 4 vs T cell, controlling for PTPRC
purity_df <- data.frame(
  f4    = f4_sub,
  tcell = tcell,
  mhc1  = mhc1,
  pur   = pur
) %>% na.omit()

# Residualise F4 on purity
f4_resid   <- residuals(lm(f4  ~ pur, data = purity_df))
tcell_resid <- residuals(lm(tcell ~ pur, data = purity_df))
mhc1_resid  <- residuals(lm(mhc1  ~ pur, data = purity_df))

cat("\n── Partial correlations (controlling for PTPRC) ──\n")
cat("Factor 4 vs T cell (partial):  r =",
    round(cor(f4_resid, tcell_resid), 3), "\n")
cat("Factor 4 vs MHC1 (partial):    r =",
    round(cor(f4_resid, mhc1_resid),  3), "\n")

# ── REGRESSION MODELS ─────────────────────────────────────────────────────────
cat("\n── Regression: T cell ~ Factor4 + PTPRC ──\n")
m_tcell <- lm(tcell ~ f4 + pur, data = purity_df)
print(summary(m_tcell))

cat("\n── Regression: MHC1 ~ Factor4 + PTPRC ──\n")
m_mhc1 <- lm(mhc1 ~ f4 + pur, data = purity_df)
print(summary(m_mhc1))

# ── VISUALISE ─────────────────────────────────────────────────────────────────
p_purity <- ggplot(purity_df, aes(x = pur, y = f4)) +
  geom_point(alpha = 0.6, size = 2, color = "#534AB7") +
  geom_smooth(method = "lm", se = TRUE,
              color = "black", linewidth = 0.8) +
  labs(title    = "Factor 4 vs PTPRC (leukocyte fraction proxy)",
       subtitle = paste0("r = ", round(cor(purity_df$f4,
                                           purity_df$pur), 3)),
       x = "PTPRC expression (log2 RSEM)",
       y = "Factor 4 score") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_purity)
ggsave(file.path(fig_dir, "53_purity_factor4.png"),
       p_purity, width = 7, height = 5, dpi = 150)

save(purity_df, m_tcell, m_mhc1,
     file = file.path(proc_dir, "19_purity_results.RData"))

cat("\n Purity adjustment complete!\n")