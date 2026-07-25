# ── 05_immune_characterization.R ─────────────────────────────────────────────
# Step 1 of post-MOFA analysis
# Goals:
#   1. Correlate Factor 6 with 68 validated immune signatures
#   2. Wilcoxon tests: immune-cold vs immune-hot groups
#   3. Extract Factor 6 feature weights (top genes + CpGs)
#   4. Correlate Factor 6 with 1387 Paradigm pathway scores
#   5. Generate supplementary gene table
# Input:  03_mofa_trained.RData, 01_TNBC_clean.RData, raw immune/pathway files
# Output: figures 11-13, supplementary_table_factor6_genes.csv,
#         04_step1_immune_analysis.RData
# ─────────────────────────────────────────────────────────────────────────────
# ── STEP 1: CIBERSORT + IMMUNE SIGNATURES + FACTOR 6 CHARACTERIZATION ────────

# ── 1. TRANSPOSE IMMUNE SIGS AND FILTER TO BRCA TNBC ─────────────────────────
cat("Processing immune signatures...\n")

# Rows are signatures, columns are samples — need to transpose
immune_mat <- immune_sigs %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  as.matrix()

# Transpose so samples are rows
immune_t <- t(immune_mat)
cat("Immune matrix dimensions (samples x signatures):", dim(immune_t), "\n")

# Trim sample IDs to 15 chars to match our TNBC samples
rownames(immune_t) <- substr(rownames(immune_t), 1, 15)

# Filter to our TNBC samples
common_imm <- intersect(rownames(immune_t), common_samples)
cat("TNBC samples with immune scores:", length(common_imm), "\n")

immune_tnbc <- immune_t[common_imm, ]

# ── 2. GET FACTOR SCORES ──────────────────────────────────────────────────────
factor_scores <- get_factors(mofa_trained, factors = "all")[[1]]
f6 <- factor_scores[common_imm, "Factor6"]

# ── 3. KEY IMMUNE SIGNATURES TO TEST ─────────────────────────────────────────
key_sigs <- c(
  "Tcell_21978456",
  "Bcell_21978456", 
  "MHC1_21978456",
  "MHC2_21978456",
  "IFNG_score_21050467",
  "TGFB_score_21050467",
  "PD1_PDL1_score",
  "CD8A",
  "CD68",
  "CD8_CD68_ratio",
  "Module3_IFN_score",
  "Module4_TcellBcell_score"
)

# Check all exist
key_sigs <- intersect(key_sigs, colnames(immune_tnbc))
cat("Signatures found:", length(key_sigs), "\n")

# ── 4. CORRELATE FACTOR 6 WITH EACH IMMUNE SIGNATURE ─────────────────────────
cat("\n── Factor 6 correlations with immune signatures ──\n")
cor_results <- data.frame(
  signature = key_sigs,
  r         = sapply(key_sigs, function(s) {
    cor(f6, immune_tnbc[, s], method = "pearson", use = "complete.obs")
  }),
  p_value   = sapply(key_sigs, function(s) {
    cor.test(f6, immune_tnbc[, s], method = "pearson")$p.value
  })
) %>%
  mutate(
    r         = round(r, 3),
    p_value   = round(p_value, 4),
    significant = p_value < 0.05
  ) %>%
  arrange(r)

print(cor_results)

# ── 5. WILCOXON TESTS: IMMUNE-COLD VS IMMUNE-HOT ─────────────────────────────
cat("\n── Wilcoxon tests: immune-cold vs immune-hot ──\n")

# Define groups by Factor 6 median split
f6_median <- median(f6)
groups <- ifelse(f6 > f6_median, "Immune-cold", "Immune-hot")
cat("Immune-cold n =", sum(groups == "Immune-cold"), "\n")
cat("Immune-hot  n =", sum(groups == "Immune-hot"),  "\n")

wilcox_results <- data.frame(
  signature = key_sigs,
  p_value   = sapply(key_sigs, function(s) {
    wilcox.test(
      immune_tnbc[groups == "Immune-cold", s],
      immune_tnbc[groups == "Immune-hot",  s]
    )$p.value
  }),
  median_cold = sapply(key_sigs, function(s) {
    round(median(immune_tnbc[groups == "Immune-cold", s], na.rm = TRUE), 3)
  }),
  median_hot = sapply(key_sigs, function(s) {
    round(median(immune_tnbc[groups == "Immune-hot", s], na.rm = TRUE), 3)
  })
) %>%
  mutate(
    p_value     = round(p_value, 4),
    significant = p_value < 0.05,
    direction   = ifelse(median_cold < median_hot, "↓ in cold", "↑ in cold")
  ) %>%
  arrange(p_value)

print(wilcox_results)

# ── 6. VISUALIZATION ──────────────────────────────────────────────────────────

# Correlation barplot
p_cor <- cor_results %>%
  mutate(signature = factor(signature, levels = signature)) %>%
  ggplot(aes(x = signature, y = r, 
             fill = ifelse(significant, "significant", "not significant"))) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("significant"     = "#D85A30",
                               "not significant" = "#B4B2A9")) +
  coord_flip() +
  labs(title    = "Factor 6 correlation with immune signatures — TNBC",
       subtitle = "Negative r = immune-cold phenotype",
       x = NULL, y = "Pearson r", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "bottom")

print(p_cor)
ggsave(file.path(fig_dir, "11_factor6_immune_correlations.png"),
       p_cor, width = 9, height = 6, dpi = 150)

# Boxplots for significant signatures
sig_sigs <- wilcox_results %>% filter(significant) %>% pull(signature)

if(length(sig_sigs) > 0) {
  plot_list <- lapply(sig_sigs, function(s) {
    df <- data.frame(
      score = immune_tnbc[, s],
      group = groups
    )
    p_val <- wilcox_results %>% filter(signature == s) %>% pull(p_value)
    
    ggplot(df, aes(x = group, y = score, fill = group)) +
      geom_boxplot(alpha = 0.8, outlier.size = 1.5) +
      scale_fill_manual(values = c("Immune-cold" = "#D85A30",
                                   "Immune-hot"  = "#1D9E75")) +
      labs(title = s,
           subtitle = paste0("Wilcoxon p = ", p_val),
           x = NULL, y = "Score") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 10))
  })
  
  p_boxes <- wrap_plots(plot_list, ncol = 3)
  print(p_boxes)
  ggsave(file.path(fig_dir, "12_immune_boxplots.png"),
         p_boxes, width = 12, height = 8, dpi = 150)
} else {
  cat("No signatures significant at p<0.05 — showing top 3 by p-value\n")
}

# ── 7. FACTOR 6 WEIGHT EXTRACTION ─────────────────────────────────────────────
cat("\n── Extracting Factor 6 weights ──\n")

weights <- get_weights(mofa_trained, factors = "Factor6", as.data.frame = TRUE)

# Top RNA weights
rna_weights <- weights %>%
  filter(view == "RNA") %>%
  arrange(desc(abs(value))) %>%
  head(30)

cat("Top 10 RNA genes driving Factor 6:\n")
print(rna_weights[1:10, c("feature", "value")])

# Top methylation weights
meth_weights <- weights %>%
  filter(view == "Methylation") %>%
  arrange(desc(abs(value))) %>%
  head(30)

cat("\nTop 10 CpG sites driving Factor 6:\n")
print(meth_weights[1:10, c("feature", "value")])

# ── 8. SUPPLEMENTARY TABLE ────────────────────────────────────────────────────
supp_table <- weights %>%
  filter(view == "RNA") %>%
  arrange(desc(abs(value))) %>%
  head(50) %>%
  select(Gene = feature, 
         Factor6_weight = value,
         View = view) %>%
  mutate(
    Factor6_weight = round(Factor6_weight, 4),
    Direction = ifelse(Factor6_weight > 0, "Positive", "Negative")
  )

write.csv(supp_table,
          file = file.path(base_dir, "outputs/supplementary_table_factor6_genes.csv"),
          row.names = FALSE)

cat("\n Supplementary table saved!\n")

# ── 9. PARADIGM PATHWAY CORRELATION ───────────────────────────────────────────
cat("\nProcessing Paradigm pathway scores...\n")

# Transpose paradigm (rows = pathways, cols = samples)
paradigm_mat <- paradigm_z %>%
  as.data.frame() %>%
  column_to_rownames("sample") %>%
  as.matrix()

paradigm_t <- t(paradigm_mat)
rownames(paradigm_t) <- substr(rownames(paradigm_t), 1, 15)

# Filter to TNBC
common_par <- intersect(rownames(paradigm_t), common_samples)
cat("TNBC samples with pathway scores:", length(common_par), "\n")

paradigm_tnbc <- paradigm_t[common_par, ]
f6_par <- factor_scores[common_par, "Factor6"]

# Correlate all 1387 pathways with Factor 6
cat("Correlating Factor 6 with 1387 pathways...\n")
pathway_cors <- sapply(colnames(paradigm_tnbc), function(p) {
  cor(f6_par, paradigm_tnbc[, p], 
      method = "pearson", use = "complete.obs")
})

# Top positively and negatively correlated pathways
pathway_df <- data.frame(
  pathway = names(pathway_cors),
  r       = round(pathway_cors, 3)
) %>% arrange(r)

cat("\nTop 10 pathways NEGATIVELY correlated with Factor 6 (immune-cold):\n")
print(head(pathway_df, 10))

cat("\nTop 10 pathways POSITIVELY correlated with Factor 6:\n")
print(tail(pathway_df, 10))

# Plot top pathways
p_pathways <- pathway_df %>%
  slice(c(1:10, (nrow(.)-9):nrow(.))) %>%
  mutate(pathway = factor(pathway, levels = pathway),
         direction = ifelse(r < 0, "Negative", "Positive")) %>%
  ggplot(aes(x = pathway, y = r, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Negative" = "#1D9E75",
                               "Positive" = "#D85A30")) +
  coord_flip() +
  labs(title    = "Top pathways correlated with Factor 6 — TNBC",
       subtitle = "Negative = depleted in immune-cold tumors",
       x = NULL, y = "Pearson r with Factor 6", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

print(p_pathways)
ggsave(file.path(fig_dir, "13_factor6_pathway_correlations.png"),
       p_pathways, width = 10, height = 8, dpi = 150)

# Save everything
save(immune_tnbc, immune_subtype, paradigm_tnbc,
     cor_results, wilcox_results, pathway_df,
     supp_table, groups,
     file = file.path(proc_dir, "04_step1_immune_analysis.RData"))

cat("\n Step 1 complete!\n")
cat("Figures saved to outputs/figures/\n")
cat("Supplementary table saved to outputs/\n")