library(tidyverse)
library(ggplot2)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "01_TNBC_clean.RData"))

# ── 1. VARIANCE EXPLAINED BY FACTOR 4 PER OMICS ───────────────────────────────
var_exp <- calculate_variance_explained(mofa_trained_tnbc3)

# Factor 4 specifically
f4_var <- var_exp$r2_per_factor[[1]]["Factor4", ]
cat("── Variance explained by Factor 4 per omics ──\n")
print(round(f4_var, 4))

# All factors variance for context
cat("\n── Total variance explained per omics (all factors) ──\n")
print(round(colSums(var_exp$r2_per_factor[[1]]), 3))

# ── 2. FACTOR 4 WEIGHT MAGNITUDE PER OMICS ────────────────────────────────────
weights_f4 <- get_weights(mofa_trained_tnbc3,
                          factors = "Factor4",
                          as.data.frame = TRUE)

weight_summary <- weights_f4 %>%
  group_by(view) %>%
  summarise(
    n_features      = n(),
    mean_abs_weight = round(mean(abs(value)), 5),
    max_abs_weight  = round(max(abs(value)), 5),
    top5_mean       = round(mean(sort(abs(value), decreasing=TRUE)[1:5]), 5)
  ) %>%
  arrange(desc(mean_abs_weight))

cat("\n── Factor 4 weight magnitude per omics ──\n")
print(weight_summary)

# ── 3. PARTIAL R² APPROACH ────────────────────────────────────────────────────
# Get Factor 4 scores
factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
f4 <- factor_scores3[, "Factor4"]

# Get immune group
immune_group <- ifelse(f4 > median(f4), 1, 0)  # 1=cold, 0=hot

# For each omics: how much does a PCA summary of that omics 
# predict immune group independently?
# Fix predict conflict
predict_immune <- function(mat, label) {
  common_s <- intersect(names(f4), colnames(mat))
  mat_sub  <- mat[, common_s]
  mat_sub[is.na(mat_sub)] <- 0
  
  pca  <- prcomp(t(mat_sub), scale. = TRUE, center = TRUE)
  pcs  <- as.data.frame(pca$x[, 1:min(5, ncol(pca$x))])
  grp  <- immune_group[common_s]
  
  df   <- cbind(pcs, group = grp)
  mod  <- glm(group ~ ., data = df, family = binomial)
  
  null_dev  <- mod$null.deviance
  res_dev   <- mod$deviance
  pseudo_r2 <- round(1 - (res_dev / null_dev), 4)
  
  # Use stats::predict explicitly
  pred <- stats::predict(mod, type = "response")
  auc  <- round(as.numeric(
    pROC::auc(pROC::roc(grp, pred, quiet = TRUE))
  ), 3)
  
  cat(label, "— Pseudo R² =", pseudo_r2, "| AUC =", auc, "\n")
  return(data.frame(omics = label, pseudo_r2 = pseudo_r2, auc = auc))
}

cat("\n── Predictive power of each omics for immune group ──\n")
res_rna  <- predict_immune(rna_top3,   "RNA-seq")
res_meth <- predict_immune(meth_top3,  "Methylation")
res_cnv  <- predict_immune(cnv_clean3, "CNV")

pred_df <- bind_rows(res_rna, res_meth, res_cnv)
print(pred_df)
# ── 4. VISUALISE ──────────────────────────────────────────────────────────────
# A: Variance explained by Factor 4 per omics
var_df <- data.frame(
  omics    = names(f4_var),
  variance = as.numeric(f4_var)
)

p_var <- ggplot(var_df, aes(x = reorder(omics, variance),
                            y = variance, fill = omics)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("RNA"         = "#1D9E75",
                               "Methylation" = "#534AB7",
                               "CNV"         = "#D85A30")) +
  geom_text(aes(label = paste0(round(variance, 2), "%")),
            hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  labs(
    title    = "Variance explained by Factor 4 per omics layer",
    subtitle = "Contribution of each molecular layer to the immune axis",
    x = NULL, y = "Variance explained (%)"
  ) +
  xlim(NA, max(var_df$variance) * 1.2) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_var)
ggsave(file.path(fig_dir, "65_factor4_variance_per_omics.png"),
       p_var, width = 8, height = 4, dpi = 150)

# B: AUC — predictive power
p_auc <- ggplot(pred_df, aes(x = reorder(omics, auc),
                             y = auc, fill = omics)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("RNA-seq"     = "#1D9E75",
                               "Methylation" = "#534AB7",
                               "CNV"         = "#D85A30")) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "#888780", linewidth = 0.8) +
  geom_text(aes(label = paste0("AUC=", auc)),
            hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  ylim(0, 1) +
  labs(
    title    = "Predictive power of each omics for immune group",
    subtitle = "AUC from logistic regression using top 5 PCs per omics",
    x = NULL, y = "AUC (0.5 = random)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_auc)
ggsave(file.path(fig_dir, "66_omics_predictive_power.png"),
       p_auc, width = 8, height = 4, dpi = 150)

# C: Weight distribution per omics
p_weights <- ggplot(weights_f4, aes(x = abs(value), fill = view)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  scale_fill_manual(values = c("RNA"         = "#1D9E75",
                               "Methylation" = "#534AB7",
                               "CNV"         = "#D85A30")) +
  facet_wrap(~view, scales = "free", ncol = 1) +
  labs(
    title    = "Factor 4 weight distribution per omics layer",
    subtitle = "Larger weights = stronger contribution to immune axis",
    x = "|Weight|", y = "Count", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "none")

print(p_weights)
ggsave(file.path(fig_dir, "67_factor4_weight_distributions.png"),
       p_weights, width = 7, height = 8, dpi = 150)

cat("\n Omics contribution analysis complete!\n")


# Plots
# Plot variance explained
p_var <- data.frame(
  omics    = c("RNA-seq", "Methylation", "CNV"),
  variance = c(2.13, 1.75, 3.48)
) %>%
  ggplot(aes(x = reorder(omics, variance), y = variance, fill = omics)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("RNA-seq"     = "#1D9E75",
                               "Methylation" = "#534AB7",
                               "CNV"         = "#D85A30")) +
  geom_text(aes(label = paste0(round(variance, 2), "%")),
            hjust = -0.1, size = 4.5, fontface = "bold") +
  coord_flip() +
  expand_limits(x = c(0, 5)) +
  labs(
    title    = "Variance explained by Factor 4 per omics layer",
    subtitle = "CNV contributes most to the immune axis",
    x = NULL, y = "Variance explained (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_var)
ggsave(file.path(fig_dir, "65_factor4_variance_per_omics.png"),
       p_var, width = 8, height = 4, dpi = 150)

# Plot AUC
p_auc <- pred_df %>%
  ggplot(aes(x = reorder(omics, auc), y = auc, fill = omics)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("RNA-seq"     = "#1D9E75",
                               "Methylation" = "#534AB7",
                               "CNV"         = "#D85A30")) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "#888780", linewidth = 0.8) +
  geom_hline(yintercept = 0.7, linetype = "dotted",
             color = "#444441", linewidth = 0.6) +
  geom_text(aes(label = paste0("AUC = ", auc)),
            hjust = -0.1, size = 4.5, fontface = "bold") +
  coord_flip() +
  ylim(0, 1.1) +
  annotate("text", x = 0.6, y = 0.51,
           label = "random", hjust = 0,
           size = 3, color = "#888780") +
  annotate("text", x = 0.6, y = 0.71,
           label = "AUC = 0.7 threshold",
           hjust = 0, size = 3, color = "#444441") +
  labs(
    title    = "Single-omics predictive power for immune group",
    subtitle = "AUC from logistic regression using top 5 PCs per omics layer",
    x = NULL, y = "AUC (area under ROC curve)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_auc)
ggsave(file.path(fig_dir, "66_omics_predictive_power.png"),
       p_auc, width = 8, height = 4, dpi = 150)

# Combined figure — variance + AUC side by side
library(patchwork)

p_var <- data.frame(
  omics    = c("RNA-seq", "Methylation", "CNV"),
  variance = c(2.13, 1.75, 3.48)
) %>%
  ggplot(aes(x = reorder(omics, variance),
             y = variance, fill = omics)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("RNA-seq"     = "#1D9E75",
                               "Methylation" = "#534AB7",
                               "CNV"         = "#D85A30")) +
  geom_text(aes(label = paste0(round(variance, 2), "%")),
            hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  expand_limits(y = 5) +
  labs(
    title = "Variance explained by Factor 4",
    x = NULL, y = "Variance explained (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 12))

p_auc <- pred_df %>%
  ggplot(aes(x = reorder(omics, auc), y = auc, fill = omics)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("RNA-seq"     = "#1D9E75",
                               "Methylation" = "#534AB7",
                               "CNV"         = "#D85A30")) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "#888780", linewidth = 0.7) +
  geom_hline(yintercept = 0.7, linetype = "dotted",
             color = "#444441", linewidth = 0.6) +
  geom_text(aes(label = paste0("AUC=", auc)),
            hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  ylim(0, 1.1) +
  annotate("text", x = 0.55, y = 0.52,
           label = "random", size = 3, color = "#888780") +
  annotate("text", x = 0.55, y = 0.72,
           label = "AUC=0.7", size = 3, color = "#444441") +
  labs(
    title = "Single-omics predictive power",
    x = NULL, y = "AUC"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 12))

p_combined <- p_var | p_auc +
  plot_annotation(
    title    = "Omics contribution to the immune-cold axis (Factor 4)",
    subtitle = "CNV dominant — all three layers independently informative",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  )

print(p_combined)
ggsave(file.path(fig_dir, "65_omics_contribution.png"),
       p_combined, width = 12, height = 5, dpi = 150)

cat(" Done!\n")
cat("CNV: AUC=0.815 — strongest contributor\n")
cat("RNA: AUC=0.780 — strong second\n")
cat("Methylation: AUC=0.703 — meaningful, weakest\n")
cat("No single omics achieves AUC>0.85 — validates multi-omics approach\n")
