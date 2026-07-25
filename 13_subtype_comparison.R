# ── 13_subtype_comparison.R ───────────────────────────────────────────────────
# Goals:
#   1. Formal subtype comparison — immune-cold score across TNBC, HER2+, ER+
#   2. Tumour size and grade comparison — immune-cold vs immune-hot METABRIC
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(ggplot2)
library(patchwork)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "08_mofa_full_trained.RData"))
load(file.path(proc_dir, "07_full_cohort_clean.RData"))
load(file.path(proc_dir, "12_metabric_validation.RData"))

f <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/metabric", fn))

# ── 1. SUBTYPE COMPARISON IN FULL TCGA COHORT ─────────────────────────────────
cat("── Subtype comparison ──\n")

# Get Factor 1 scores from full cohort model
factor_scores_full <- get_factors(mofa_trained_full, factors = "all")[[1]]
f1_full <- factor_scores_full[, "Factor1"]

# Build dataframe with subtype labels
subtype_score_df <- data.frame(
  sample  = names(f1_full),
  factor1 = f1_full
) %>%
  left_join(subtype_df, by = c("sample" = "sampleID")) %>%
  filter(subtype %in% c("TNBC", "HER2+", "ER+")) %>%
  mutate(subtype = factor(subtype, levels = c("TNBC", "HER2+", "ER+")))

cat("Sample counts per subtype:\n")
print(table(subtype_score_df$subtype))

# Kruskal-Wallis test — are scores different across subtypes?
kw_test <- kruskal.test(factor1 ~ subtype, data = subtype_score_df)
cat("\nKruskal-Wallis test:\n")
print(kw_test)

# Pairwise Wilcoxon tests
pairwise <- pairwise.wilcox.test(
  subtype_score_df$factor1,
  subtype_score_df$subtype,
  p.adjust.method = "bonferroni"
)
cat("\nPairwise Wilcoxon (Bonferroni corrected):\n")
print(pairwise)

# Plot — boxplot per subtype
p_subtype <- ggplot(subtype_score_df,
                    aes(x = subtype, y = factor1, fill = subtype)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.5, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("TNBC"  = "#D85A30",
                               "HER2+" = "#534AB7",
                               "ER+"   = "#1D9E75")) +
  labs(
    title    = "MOFA2 Factor 1 score by breast cancer subtype",
    subtitle = paste0("Kruskal-Wallis p = ",
                      formatC(kw_test$p.value, format = "e", digits = 2)),
    x        = "Subtype",
    y        = "Factor 1 score (immune-cold axis)",
    fill     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "none"
  ) +
  annotate("text", x = 2, y = max(subtype_score_df$factor1) * 0.95,
           label = "* TNBC significantly lower than ER+ and HER2+",
           size = 3.5, color = "#5F5E5A")

print(p_subtype)
ggsave(file.path(fig_dir, "38_subtype_comparison.png"),
       p_subtype, width = 7, height = 5, dpi = 150)

# ── 2. TUMOUR SIZE AND GRADE — METABRIC ───────────────────────────────────────
cat("\n── Tumour size and grade comparison ──\n")

# Load clinical sample data
clin_sam <- fread(file = f("data_clinical_sample.txt"), header = TRUE)
colnames(clin_sam)[1] <- "patient_id"

# Merge with immune groups
clinical_immune <- surv_meta %>%
  left_join(
    clin_sam %>%
      select(sample_id     = `Sample Identifier`,
             tumour_size   = `Tumor Size`,
             grade         = `Neoplasm Histologic Grade`),
    by = "sample_id"
  ) %>%
  mutate(
    tumour_size = as.numeric(tumour_size),
    grade       = as.numeric(grade)
  ) %>%
  filter(!is.na(tumour_size), !is.na(grade))

cat("Samples with tumour size + grade:", nrow(clinical_immune), "\n")

# Summary statistics
cat("\n── Tumour size by immune group ──\n")
clinical_immune %>%
  group_by(immune_group) %>%
  summarise(
    n      = n(),
    mean   = round(mean(tumour_size, na.rm=TRUE), 1),
    median = round(median(tumour_size, na.rm=TRUE), 1),
    sd     = round(sd(tumour_size, na.rm=TRUE), 1)
  ) %>% print()

# Wilcoxon test for tumour size
wt_size <- wilcox.test(
  tumour_size ~ immune_group,
  data = clinical_immune
)
cat("\nWilcoxon test — tumour size: p =",
    round(wt_size$p.value, 4), "\n")

# Summary statistics for grade
cat("\n── Grade by immune group ──\n")
clinical_immune %>%
  group_by(immune_group) %>%
  summarise(
    n      = n(),
    mean   = round(mean(grade, na.rm=TRUE), 2),
    median = median(grade, na.rm=TRUE)
  ) %>% print()

# Wilcoxon test for grade
wt_grade <- wilcox.test(
  grade ~ immune_group,
  data = clinical_immune
)
cat("\nWilcoxon test — grade: p =",
    round(wt_grade$p.value, 4), "\n")

# Plot tumour size
p_size <- ggplot(clinical_immune,
                 aes(x = immune_group, y = tumour_size, fill = immune_group)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.5, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("Immune-cold" = "#D85A30",
                               "Immune-hot"  = "#1D9E75")) +
  labs(
    title    = "Tumour size by immune group — METABRIC TNBC",
    subtitle = paste0("Wilcoxon p = ", round(wt_size$p.value, 4)),
    x        = NULL,
    y        = "Tumour size (mm)",
    fill     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "none")

# Plot grade
p_grade <- ggplot(clinical_immune,
                  aes(x = immune_group, y = grade, fill = immune_group)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.5, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("Immune-cold" = "#D85A30",
                               "Immune-hot"  = "#1D9E75")) +
  labs(
    title    = "Tumour grade by immune group — METABRIC TNBC",
    subtitle = paste0("Wilcoxon p = ", round(wt_grade$p.value, 4)),
    x        = NULL,
    y        = "Histologic grade",
    fill     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "none")

p_clinical <- p_size | p_grade
print(p_clinical)
ggsave(file.path(fig_dir, "39_clinical_variables.png"),
       p_clinical, width = 10, height = 5, dpi = 150)

# ── 3. SAVE ───────────────────────────────────────────────────────────────────
save(subtype_score_df, kw_test, pairwise,
     clinical_immune, wt_size, wt_grade,
     file = file.path(proc_dir, "13_subtype_clinical.RData"))

cat("\n Script 13 complete!\n")
cat("Figures 38-39 saved to outputs/figures/\n")