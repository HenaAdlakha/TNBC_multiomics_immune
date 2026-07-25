library(tidyverse)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "16_methylation_features.RData"))

# ── BUILD CONTINGENCY TABLE ───────────────────────────────────────────────────
# Count CpG sites per island context per group
island_counts <- weights_annotated %>%
  filter(!is.na(relation_island)) %>%
  mutate(group = ifelse(value < 0, "Immune-cold", "Immune-hot")) %>%
  count(group, relation_island) %>%
  pivot_wider(names_from = relation_island, values_from = n, values_fill = 0)

cat("Contingency table:\n")
print(island_counts)

# Convert to matrix for tests
cont_mat <- island_counts %>%
  column_to_rownames("group") %>%
  as.matrix()

cat("\nContingency matrix:\n")
print(cont_mat)

# ── CHI-SQUARE TEST ───────────────────────────────────────────────────────────
chi_test <- chisq.test(cont_mat)
cat("\n── Chi-square test ──\n")
print(chi_test)

# ── FISHER TEST ───────────────────────────────────────────────────────────────
# Fisher exact works best on 2x2 — test each island context vs rest
cat("\n── Fisher exact: each context vs all others ──\n")

fisher_results <- data.frame(
  context   = colnames(cont_mat),
  cold_n    = cont_mat["Immune-cold", ],
  hot_n     = cont_mat["Immune-hot",  ],
  odds_ratio = NA,
  p_value    = NA
)

for(ctx in colnames(cont_mat)) {
  # 2x2: this context vs all others, cold vs hot
  mat_2x2 <- matrix(c(
    cont_mat["Immune-cold", ctx],
    cont_mat["Immune-hot",  ctx],
    sum(cont_mat["Immune-cold", colnames(cont_mat) != ctx]),
    sum(cont_mat["Immune-hot",  colnames(cont_mat) != ctx])
  ), nrow = 2,
  dimnames = list(
    c("Immune-cold", "Immune-hot"),
    c(ctx, "Other")
  ))
  
  ft <- fisher.test(mat_2x2)
  fisher_results$odds_ratio[fisher_results$context == ctx] <-
    round(ft$estimate, 3)
  fisher_results$p_value[fisher_results$context == ctx]   <-
    round(ft$p.value, 4)
}

fisher_results$p_adj <- round(
  p.adjust(fisher_results$p_value, method = "fdr"), 4
)
fisher_results$significant <- fisher_results$p_adj < 0.05
fisher_results$enriched_in <- ifelse(
  fisher_results$cold_n / sum(cont_mat["Immune-cold", ]) >
    fisher_results$hot_n  / sum(cont_mat["Immune-hot",  ]),
  "Immune-cold", "Immune-hot"
)

cat("\nFisher exact results (FDR corrected):\n")
print(fisher_results %>% arrange(p_adj))

# ── PLOT WITH SIGNIFICANCE ────────────────────────────────────────────────────
# Calculate proportions for fair comparison
island_prop <- weights_annotated %>%
  filter(!is.na(relation_island)) %>%
  mutate(group = ifelse(value < 0, "Immune-cold", "Immune-hot")) %>%
  count(group, relation_island) %>%
  group_by(group) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  left_join(
    fisher_results %>%
      dplyr::select(context, p_adj, significant),
    by = c("relation_island" = "context")
  )

# Add significance labels
island_prop <- island_prop %>%
  mutate(sig_label = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE          ~ ""
  ))

p_island_test <- ggplot(island_prop,
                        aes(x = relation_island, y = proportion,
                            fill = group)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Immune-cold" = "#D85A30",
                               "Immune-hot"  = "#1D9E75")) +
  geom_text(
    data = island_prop %>%
      filter(group == "Immune-hot") %>%
      distinct(relation_island, sig_label, p_adj),
    aes(x = relation_island,
        y = max(island_prop$proportion) + 0.02,
        label = sig_label),
    inherit.aes = FALSE,
    size = 5, color = "#444441"
  ) +
  labs(
    title    = "CpG island context — immune-cold vs immune-hot (Factor 4)",
    subtitle = paste0("Chi-square p = ",
                      formatC(chi_test$p.value,
                              format = "e", digits = 2),
                      " | * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001"),
    x        = "Relation to CpG island",
    y        = "Proportion of CpG sites",
    fill     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x   = element_text(angle = 45, hjust = 1)
  )

print(p_island_test)
ggsave(file.path(fig_dir, "59_cpg_island_fisher.png"),
       p_island_test, width = 9, height = 6, dpi = 150)

# Save
write.csv(fisher_results,
          file.path(base_dir, "outputs/cpg_island_fisher_results.csv"),
          row.names = FALSE)

cat("\n Done!\n")
cat("Chi-square p =",
    formatC(chi_test$p.value, format = "e", digits = 2), "\n")
cat("Significant contexts after FDR correction:\n")
print(fisher_results %>%
        filter(significant) %>%
        dplyr::select(context, enriched_in, odds_ratio, p_adj))


library(tidyverse)
library(ggplot2)
library(patchwork)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "01_TNBC_clean.RData"))
load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))

# Get immune groups from Factor 4
factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
f4 <- factor_scores3[, "Factor4"]

immune_groups <- data.frame(
  sample       = names(f4),
  immune_group = ifelse(f4 > median(f4), "Immune-cold", "Immune-hot")
)

# ── FUNCTION: PCA + ANOVA ─────────────────────────────────────────────────────
run_pca_anova <- function(mat, label, color, immune_groups) {
  
  # Align samples
  common_s <- intersect(colnames(mat), immune_groups$sample)
  mat_sub  <- mat[, common_s]
  groups   <- immune_groups %>% filter(sample %in% common_s)
  
  # PCA
  mat_t <- t(mat_sub)
  mat_t[is.na(mat_t)] <- 0
  pca  <- prcomp(mat_t, scale. = TRUE, center = TRUE)
  ve   <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  
  pca_df <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(groups, by = "sample")
  
  # ANOVA: does immune group predict PC1 and PC2?
  aov_pc1 <- summary(aov(PC1 ~ immune_group, data = pca_df))
  aov_pc2 <- summary(aov(PC2 ~ immune_group, data = pca_df))
  
  p1 <- round(aov_pc1[[1]]$`Pr(>F)`[1], 4)
  p2 <- round(aov_pc2[[1]]$`Pr(>F)`[1], 4)
  
  cat(label, "\n")
  cat("  PC1 ANOVA p =", p1, "\n")
  cat("  PC2 ANOVA p =", p2, "\n\n")
  
  # Plot colored by immune group with ANOVA p-values
  p <- ggplot(pca_df, aes(PC1, PC2, color = immune_group)) +
    geom_point(alpha = 0.8, size = 2.5) +
    scale_color_manual(values = c("Immune-cold" = "#D85A30",
                                  "Immune-hot"  = "#1D9E75")) +
    stat_ellipse(level = 0.95, linewidth = 0.8) +
    labs(
      title    = label,
      subtitle = paste0("PC1 ANOVA p=", p1,
                        " | PC2 ANOVA p=", p2),
      x        = paste0("PC1 (", ve[1], "%)"),
      y        = paste0("PC2 (", ve[2], "%)"),
      color    = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      legend.position = "bottom"
    )
  
  return(list(plot = p, p1 = p1, p2 = p2))
}

# ── RUN FOR EACH OMICS ────────────────────────────────────────────────────────
cat("── ANOVA results: immune group vs PCA axes ──\n\n")

res_rna  <- run_pca_anova(rna_top3,   "RNA-seq",     "#1D9E75", immune_groups)
res_meth <- run_pca_anova(meth_top3,  "Methylation", "#534AB7", immune_groups)
res_cnv  <- run_pca_anova(cnv_clean3, "CNV",         "#D85A30", immune_groups)

# ── COMBINED PLOT ─────────────────────────────────────────────────────────────
p_combined <- (res_rna$plot | res_meth$plot) /
  (res_cnv$plot | plot_spacer()) +
  plot_annotation(
    title    = "PCA per omics layer — colored by immune group",
    subtitle = "ANOVA p-values test whether immune-cold vs immune-hot differ on each PC",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  )

print(p_combined)
ggsave(file.path(fig_dir, "64_PCA_ANOVA_immune_groups.png"),
       p_combined, width = 12, height = 10, dpi = 150)

# ── SUMMARY TABLE ─────────────────────────────────────────────────────────────
anova_summary <- data.frame(
  Omics    = c("RNA-seq", "Methylation", "CNV"),
  PC1_p    = c(res_rna$p1, res_meth$p1, res_cnv$p1),
  PC2_p    = c(res_rna$p2, res_meth$p2, res_cnv$p2)
) %>%
  mutate(
    PC1_sig = ifelse(PC1_p < 0.05, "*", "ns"),
    PC2_sig = ifelse(PC2_p < 0.05, "*", "ns")
  )

cat("\n── Summary ──\n")
print(anova_summary)

write.csv(anova_summary,
          file.path(base_dir, "outputs/PCA_ANOVA_results.csv"),
          row.names = FALSE)

cat("\n PCA ANOVA complete!\n")

# Save ANOVA results
anova_summary_final <- data.frame(
  Omics   = c("RNA-seq", "Methylation", "CNV"),
  PC1_p   = c(0.2537, 0.0387, 0.0582),
  PC2_p   = c(0.6511, 0.0385, 0.3058),
  PC1_sig = c("ns", "*", "trend"),
  PC2_sig = c("ns", "*", "ns")
)

write.csv(anova_summary_final,
          file.path(base_dir, "outputs/PCA_ANOVA_final.csv"),
          row.names = FALSE)

cat(" PCA ANOVA results saved\n")