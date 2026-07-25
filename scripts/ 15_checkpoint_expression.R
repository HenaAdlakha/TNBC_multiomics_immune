# ── 15_checkpoint_expression.R ────────────────────────────────────────────────
# Goals:
#   1. Checkpoint gene expression in TCGA TNBC (n=81)
#      by immune-cold vs immune-hot groups
#   2. Validate in METABRIC TNBC (n=320)
#   3. Assess immunotherapy sensitivity features of immune-hot subgroup
#
# Input:  01_TNBC_clean.RData, 09_mofa_tnbc3_trained.RData,
#         12_metabric_validation.RData, METABRIC expression matrix
# Output: figures 46-47, results tables
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(ggplot2)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

fm <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/metabric", fn))

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
load(file.path(proc_dir, "01_TNBC_clean.RData"))
load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "12_metabric_validation.RData"))

# ── 1. TCGA: CHECKPOINT EXPRESSION BY IMMUNE GROUP ────────────────────────────
cat("── TCGA checkpoint analysis ──\n")

factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
f4 <- factor_scores3[, "Factor4"]

checkpoint_genes <- c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "SIGLEC15")

# Confirm genes present
present_tcga <- intersect(checkpoint_genes, rownames(rna_top))
cat("Checkpoint genes in TCGA RNA matrix:", present_tcga, "\n")

# Build dataframe
common_ck <- intersect(names(f4), colnames(rna_top))

ck_df <- as.data.frame(t(rna_top[present_tcga, common_ck])) %>%
  rownames_to_column("sample") %>%
  mutate(
    factor4      = f4[sample],
    immune_group = ifelse(factor4 > median(f4[common_ck]),
                          "Immune-cold", "Immune-hot")
  )

cat("Group sizes — TCGA:\n")
print(table(ck_df$immune_group))

# Wilcoxon tests
results_tcga_ck <- data.frame(
  gene      = present_tcga,
  mean_cold = sapply(present_tcga, function(g)
    round(mean(ck_df[[g]][ck_df$immune_group == "Immune-cold"],
               na.rm = TRUE), 3)),
  mean_hot  = sapply(present_tcga, function(g)
    round(mean(ck_df[[g]][ck_df$immune_group == "Immune-hot"],
               na.rm = TRUE), 3)),
  p_value   = sapply(present_tcga, function(g)
    round(wilcox.test(ck_df[[g]] ~ ck_df$immune_group)$p.value, 4))
) %>%
  mutate(
    p_adj       = round(p.adjust(p_value, method = "fdr"), 4),
    direction   = ifelse(mean_cold < mean_hot,
                         "Higher in immune-hot",
                         "Higher in immune-cold"),
    significant = p_adj < 0.05,
    cohort      = "TCGA"
  ) %>%
  arrange(p_value)

cat("\n── TCGA checkpoint results ──\n")
print(results_tcga_ck)

# Plot TCGA
ck_long_tcga <- ck_df %>%
  dplyr::select(sample, immune_group, all_of(present_tcga)) %>%
  pivot_longer(cols = all_of(present_tcga),
               names_to = "gene",
               values_to = "expression")

p_ck_tcga <- ggplot(ck_long_tcga,
                    aes(x = immune_group, y = expression,
                        fill = immune_group)) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("Immune-cold" = "#D85A30",
                               "Immune-hot"  = "#1D9E75")) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(
    title    = "Checkpoint gene expression by immune group — TCGA TNBC",
    subtitle = paste0("n=", nrow(ck_df), " — discovery cohort"),
    x = NULL, y = "Expression (log2 RSEM)", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "bottom",
        strip.text      = element_text(face = "bold"))

print(p_ck_tcga)
ggsave(file.path(fig_dir, "46_checkpoint_TCGA.png"),
       p_ck_tcga, width = 10, height = 7, dpi = 150)

# ── 2. METABRIC: CHECKPOINT VALIDATION ────────────────────────────────────────
cat("\n── METABRIC checkpoint validation ──\n")

# Load METABRIC expression
expr_raw <- fread(
  file = fm("data_mrna_illumina_microarray.txt"),
  header = TRUE)

expr_mat <- expr_raw %>%
  as.data.frame() %>%
  filter(!duplicated(Hugo_Symbol)) %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

# Get immune group samples
cold_meta <- surv_meta %>%
  filter(immune_group == "Immune-cold") %>%
  pull(sample_id)
hot_meta  <- surv_meta %>%
  filter(immune_group == "Immune-hot") %>%
  pull(sample_id)

cold_meta <- intersect(cold_meta, colnames(expr_mat))
hot_meta  <- intersect(hot_meta,  colnames(expr_mat))

cat("METABRIC group sizes:\n")
cat("Immune-cold:", length(cold_meta), "\n")
cat("Immune-hot: ", length(hot_meta),  "\n")

# Confirm genes present
present_meta <- intersect(checkpoint_genes, rownames(expr_mat))
cat("Checkpoint genes in METABRIC:", present_meta, "\n")

# Wilcoxon tests
results_meta_ck <- data.frame(
  gene      = present_meta,
  mean_cold = sapply(present_meta, function(g)
    round(mean(expr_mat[g, cold_meta], na.rm = TRUE), 3)),
  mean_hot  = sapply(present_meta, function(g)
    round(mean(expr_mat[g, hot_meta],  na.rm = TRUE), 3)),
  p_value   = sapply(present_meta, function(g)
    round(wilcox.test(expr_mat[g, cold_meta],
                      expr_mat[g, hot_meta])$p.value, 4))
) %>%
  mutate(
    p_adj       = round(p.adjust(p_value, method = "fdr"), 4),
    direction   = ifelse(mean_cold < mean_hot,
                         "Higher in immune-hot",
                         "Higher in immune-cold"),
    significant = p_adj < 0.05,
    cohort      = "METABRIC"
  ) %>%
  arrange(p_value)

cat("\n── METABRIC checkpoint results ──\n")
print(results_meta_ck)

# Plot METABRIC
ck_long_meta <- as.data.frame(
  t(expr_mat[present_meta, c(cold_meta, hot_meta)])
) %>%
  rownames_to_column("sample_id") %>%
  left_join(surv_meta %>%
              dplyr::select(sample_id, immune_group),
            by = "sample_id") %>%
  pivot_longer(cols = all_of(present_meta),
               names_to  = "gene",
               values_to = "expression")

p_ck_meta <- ggplot(ck_long_meta,
                    aes(x = immune_group, y = expression,
                        fill = immune_group)) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.size = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.8) +
  scale_fill_manual(values = c("Immune-cold" = "#D85A30",
                               "Immune-hot"  = "#1D9E75")) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(
    title    = "Checkpoint gene expression — METABRIC TNBC validation",
    subtitle = paste0("n=", length(cold_meta) + length(hot_meta),
                      " — independent cohort"),
    x = NULL, y = "Expression", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "bottom",
        strip.text      =