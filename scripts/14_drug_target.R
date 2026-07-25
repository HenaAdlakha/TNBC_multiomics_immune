# ── 14_drug_target.R ──────────────────────────────────────────────────────────
# Goals:
#   1. Chemotherapy response in immune-cold vs immune-hot METABRIC TNBC
#   2. Identify top upregulated genes in immune-cold as drug targets
#   3. Prepare gene signature for CMAP query
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(ggplot2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

fm <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/metabric", fn))

load(file.path(proc_dir, "12_metabric_validation.RData"))

# ── 1. LOAD TREATMENT DATA ────────────────────────────────────────────────────
clin_pat <- fread(file = fm("data_clinical_patient.txt"), header = TRUE)
colnames(clin_pat)[1] <- "patient_id"

cat("Treatment columns available:\n")
print(colnames(clin_pat)[grep("chemo|therapy|radio|treat",
                              colnames(clin_pat), ignore.case = TRUE)])

# ── 2. MERGE TREATMENT WITH IMMUNE GROUPS ─────────────────────────────────────
chemo_df <- surv_meta %>%
  left_join(
    clin_pat %>%
      select(patient_id,
             chemo   = Chemotherapy,
             hormone = `Hormone Therapy`,
             radio   = `Radio Therapy`),
    by = "patient_id"
  )

cat("\nChemotherapy distribution:\n")
print(table(chemo_df$chemo, useNA = "always"))

cat("\nChemotherapy by immune group:\n")
print(table(chemo_df$immune_group, chemo_df$chemo))

# ── 3. SURVIVAL BY CHEMO IN IMMUNE-COLD ───────────────────────────────────────
# Does chemotherapy benefit immune-cold patients?
chemo_cold <- chemo_df %>%
  filter(immune_group == "Immune-cold",
         !is.na(chemo))

cat("\nImmune-cold patients:\n")
print(table(chemo_cold$chemo))

km_chemo_cold <- survfit(Surv(OS_time, OS_event) ~ chemo,
                         data = chemo_cold)

p_chemo_cold <- ggsurvplot(
  km_chemo_cold,
  data              = chemo_cold,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#534AB7", "#D85A30"),
  title             = "OS in immune-cold TNBC — by chemotherapy",
  subtitle          = "METABRIC — does chemo benefit immune-cold patients?",
  xlab              = "Time (months)",
  ylab              = "Survival probability",
  legend.labs       = c("No chemo", "Chemo"),
  ggtheme           = theme_minimal(base_size = 13)
)

png(file.path(fig_dir, "40_chemo_immune_cold.png"),
    width = 9, height = 7, units = "in", res = 150)
print(p_chemo_cold)
dev.off()
print(p_chemo_cold)

# ── 4. SURVIVAL BY CHEMO IN IMMUNE-HOT ────────────────────────────────────────
# Compare — does chemo benefit immune-hot patients more?
chemo_hot <- chemo_df %>%
  filter(immune_group == "Immune-hot",
         !is.na(chemo))

km_chemo_hot <- survfit(Surv(OS_time, OS_event) ~ chemo,
                        data = chemo_hot)

p_chemo_hot <- ggsurvplot(
  km_chemo_hot,
  data              = chemo_hot,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#534AB7", "#D85A30"),
  title             = "OS in immune-hot TNBC — by chemotherapy",
  subtitle          = "METABRIC — does chemo benefit immune-hot patients?",
  xlab              = "Time (months)",
  ylab              = "Survival probability",
  legend.labs       = c("No chemo", "Chemo"),
  ggtheme           = theme_minimal(base_size = 13)
)

png(file.path(fig_dir, "41_chemo_immune_hot.png"),
    width = 9, height = 7, units = "in", res = 150)
print(p_chemo_hot)
dev.off()
print(p_chemo_hot)

# ── 5. INTERACTION TEST ───────────────────────────────────────────────────────
# Does immune group MODIFY the effect of chemotherapy?
# This is the key scientific question
cat("\n── Interaction test: immune group x chemotherapy ──\n")

cox_interaction <- coxph(
  Surv(OS_time, OS_event) ~ immune_group * chemo,
  data = chemo_df %>% filter(!is.na(chemo))
)
print(summary(cox_interaction))

# ── 6. IDENTIFY TOP UPREGULATED GENES IN IMMUNE-COLD ─────────────────────────
# Load METABRIC expression and immune groups
clin_sam <- fread(file = fm("data_clinical_sample.txt"), header = TRUE)
colnames(clin_sam)[1] <- "patient_id"

expr_raw <- fread(
  file = path.expand("~/Documents/R Working Directory/BRCA_project/data/metabric/data_mrna_illumina_microarray.txt"),
  header = TRUE)

expr_mat <- expr_raw %>%
  as.data.frame() %>%
  filter(!duplicated(Hugo_Symbol)) %>%
  select(-Entrez_Gene_Id) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

# Get cold vs hot sample IDs
cold_samples <- surv_meta %>%
  filter(immune_group == "Immune-cold") %>%
  pull(sample_id)

hot_samples <- surv_meta %>%
  filter(immune_group == "Immune-hot") %>%
  pull(sample_id)

cold_present <- intersect(cold_samples, colnames(expr_mat))
hot_present  <- intersect(hot_samples,  colnames(expr_mat))

cat("\nCold samples with expression:", length(cold_present), "\n")
cat("Hot samples with expression:",  length(hot_present),  "\n")

# Differential expression — Wilcoxon per gene
cat("Running differential expression (this may take 2-3 minutes)...\n")

de_results <- data.frame(
  gene    = rownames(expr_mat),
  mean_cold = rowMeans(expr_mat[, cold_present], na.rm = TRUE),
  mean_hot  = rowMeans(expr_mat[, hot_present],  na.rm = TRUE)
) %>%
  mutate(
    log2FC  = mean_cold - mean_hot,
    p_value = sapply(rownames(expr_mat), function(g) {
      tryCatch(
        wilcox.test(expr_mat[g, cold_present],
                    expr_mat[g, hot_present])$p.value,
        error = function(e) NA
      )
    }),
    p_adj   = p.adjust(p_value, method = "fdr")
  ) %>%
  arrange(p_adj)

cat("\nTop 20 upregulated genes in immune-cold TNBC:\n")
de_results %>%
  filter(log2FC > 0) %>%
  head(20) %>%
  select(gene, mean_cold, mean_hot, log2FC, p_value, p_adj) %>%
  print()

cat("\nTop 20 upregulated genes in immune-hot TNBC:\n")
de_results %>%
  filter(log2FC < 0) %>%
  head(20) %>%
  select(gene, mean_cold, mean_hot, log2FC, p_value, p_adj) %>%
  print()

# ── 7. PREPARE CMAP SIGNATURE ─────────────────────────────────────────────────
# Top 150 up and down genes for CMAP query
cmap_up <- de_results %>%
  filter(log2FC > 0, !is.na(p_adj)) %>%
  arrange(p_adj) %>%
  head(150) %>%
  pull(gene)

cmap_down <- de_results %>%
  filter(log2FC < 0, !is.na(p_adj)) %>%
  arrange(p_adj) %>%
  head(150) %>%
  pull(gene)

cat("\n── CMAP query signature ──\n")
cat("Upregulated in immune-cold (top 10):\n")
print(head(cmap_up, 10))
cat("\nDownregulated in immune-cold (top 10):\n")
print(head(cmap_down, 10))

# Save for CMAP upload
write.csv(
  data.frame(gene = cmap_up, direction = "up"),
  file.path(base_dir, "outputs/CMAP_upregulated_immune_cold.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = cmap_down, direction = "down"),
  file.path(base_dir, "outputs/CMAP_downregulated_immune_cold.csv"),
  row.names = FALSE
)

# Full DE results table
write.csv(
  de_results %>% filter(p_adj < 0.05) %>% head(200),
  file.path(base_dir, "outputs/DE_immune_cold_vs_hot_METABRIC.csv"),
  row.names = FALSE
)

# ── 8. VOLCANO PLOT ───────────────────────────────────────────────────────────
p_volcano <- de_results %>%
  filter(!is.na(p_adj)) %>%
  mutate(
    sig      = p_adj < 0.05 & abs(log2FC) > 0.5,
    direction = case_when(
      log2FC > 0.5  & p_adj < 0.05 ~ "Up in immune-cold",
      log2FC < -0.5 & p_adj < 0.05 ~ "Up in immune-hot",
      TRUE ~ "Not significant"
    )
  ) %>%
  ggplot(aes(x = log2FC, y = -log10(p_adj), color = direction)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c(
    "Up in immune-cold" = "#D85A30",
    "Up in immune-hot"  = "#1D9E75",
    "Not significant"   = "#B4B2A9"
  )) +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed", color = "#888780", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "#888780", linewidth = 0.4) +
  labs(
    title    = "Differential expression: immune-cold vs immune-hot TNBC",
    subtitle = "METABRIC cohort — Wilcoxon FDR corrected",
    x        = "Mean expression difference (cold - hot)",
    y        = "-log10(FDR)",
    color    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_volcano)
ggsave(file.path(fig_dir, "42_volcano_plot.png"),
       p_volcano, width = 9, height = 6, dpi = 150)

# ── 9. SAVE ───────────────────────────────────────────────────────────────────
save(chemo_df, de_results, cmap_up, cmap_down,
     file = file.path(proc_dir, "14_drug_target.RData"))

cat("\n Script 14 complete!\n")
cat("Figures 40-42 saved\n")
cat("CMAP files saved to outputs/\n")
cat("\nNext step: upload CMAP files to https://clue.io\n")