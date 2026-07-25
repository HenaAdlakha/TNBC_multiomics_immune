# ── TIGIT as potential drug target ────────────────────────────────────────────
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
load(file.path(proc_dir, "20_locked_validation.RData"))

# Load METABRIC expression
expr_raw <- fread(file = fm("data_mrna_illumina_microarray.txt"),
                  header = TRUE)

expr_mat <- expr_raw %>%
  as.data.frame() %>%
  filter(!duplicated(Hugo_Symbol)) %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

clin_pat <- fread(file = fm("data_clinical_patient.txt"), header = TRUE)
colnames(clin_pat)[1] <- "patient_id"

clin_sam <- fread(file = fm("data_clinical_sample.txt"), header = TRUE)
colnames(clin_sam)[1] <- "patient_id"

# ── 1. TIGIT EXPRESSION IN TNBC ───────────────────────────────────────────────
cat("TIGIT in METABRIC expression matrix:\n")
cat("Present:", "TIGIT" %in% rownames(expr_mat), "\n")

tnbc_samples <- surv_meta$sample_id
tigit_expr   <- expr_mat["TIGIT", intersect(tnbc_samples,
                                            colnames(expr_mat))]

cat("TIGIT expression summary:\n")
print(summary(tigit_expr))

# ── 2. TIGIT vs IMMUNE GROUP ───────────────────────────────────────────────────
tigit_df <- data.frame(
  sample_id  = names(tigit_expr),
  tigit      = tigit_expr
) %>%
  left_join(surv_meta %>%
              dplyr::select(sample_id, immune_group,
                            OS_time, OS_event,
                            RFS_time, RFS_event),
            by = "sample_id") %>%
  left_join(surv_locked %>%
              dplyr::select(sample_id, score_std),
            by = "sample_id") %>%
  left_join(
    clin_pat %>%
      dplyr::select(patient_id,
                    age = `Age at Diagnosis`,
                    npi = `Nottingham prognostic index`),
    by = c("sample_id" = "patient_id")
  ) %>%
  mutate(
    age        = as.numeric(age),
    npi        = as.numeric(npi),
    tigit_std  = as.numeric(scale(tigit)),
    tigit_group = ifelse(tigit > median(tigit, na.rm=TRUE),
                         "TIGIT-high", "TIGIT-low")
  )

cat("\nGroup sizes:\n")
print(table(tigit_df$tigit_group))
cat("Events:", sum(tigit_df$OS_event, na.rm=TRUE), "\n")

# Wilcoxon: TIGIT by immune group
wt_tigit <- wilcox.test(tigit ~ immune_group, data = tigit_df)
cat("\nTIGIT expression: immune-cold vs immune-hot\n")
cat("Wilcoxon p =", round(wt_tigit$p.value, 4), "\n")

tigit_df %>%
  group_by(immune_group) %>%
  summarise(
    n           = n(),
    median_tigit = round(median(tigit, na.rm=TRUE), 3),
    mean_tigit   = round(mean(tigit, na.rm=TRUE), 3)
  ) %>% print()

# ── 3. KM: OVERALL SURVIVAL BY TIGIT EXPRESSION ───────────────────────────────
km_tigit_os <- survfit(Surv(OS_time, OS_event) ~ tigit_group,
                       data = tigit_df)

p_km_tigit_os <- ggsurvplot(
  km_tigit_os,
  data              = tigit_df,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#534AB7", "#D85A30"),
  title    = "OS by TIGIT expression — METABRIC TNBC",
  subtitle = "TIGIT-high = potential anti-TIGIT responders",
  xlab     = "Time (months)",
  ylab     = "Survival probability",
  legend.labs = c("TIGIT-high", "TIGIT-low"),
  ggtheme  = theme_minimal(base_size = 13)
)

png(file.path(fig_dir, "60_KM_TIGIT_OS.png"),
    width = 9, height = 7, units = "in", res = 150)
print(p_km_tigit_os)
dev.off()
print(p_km_tigit_os)

# ── 4. KM: RFS BY TIGIT ───────────────────────────────────────────────────────
km_tigit_rfs <- survfit(Surv(RFS_time, RFS_event) ~ tigit_group,
                        data = tigit_df)

p_km_tigit_rfs <- ggsurvplot(
  km_tigit_rfs,
  data              = tigit_df,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#534AB7", "#D85A30"),
  title    = "RFS by TIGIT expression — METABRIC TNBC",
  xlab     = "Time (months)",
  ylab     = "Survival probability",
  legend.labs = c("TIGIT-high", "TIGIT-low"),
  ggtheme  = theme_minimal(base_size = 13)
)

png(file.path(fig_dir, "61_KM_TIGIT_RFS.png"),
    width = 9, height = 7, units = "in", res = 150)
print(p_km_tigit_rfs)
dev.off()
print(p_km_tigit_rfs)

# ── 5. COX: TIGIT AS CONTINUOUS PREDICTOR ─────────────────────────────────────
cox_tigit <- coxph(
  Surv(OS_time, OS_event) ~ tigit_std + age + npi,
  data = tigit_df
)
cat("\n── Cox: TIGIT continuous score ──\n")
print(summary(cox_tigit))

tigit_coef <- summary(cox_tigit)$coefficients["tigit_std", ]
cat("\nTIGIT HR per SD =", round(exp(tigit_coef["coef"]), 3), "\n")
cat("95% CI =",
    round(exp(tigit_coef["coef"] - 1.96*tigit_coef["se(coef)"]), 3),
    "–",
    round(exp(tigit_coef["coef"] + 1.96*tigit_coef["se(coef)"]), 3), "\n")
cat("p =", round(tigit_coef["Pr(>|z|)"], 4), "\n")

# ── 6. TIGIT CORRELATION WITH IMMUNE-COLD SCORE ────────────────────────────────
common_ts <- intersect(tigit_df$sample_id,
                       surv_locked$sample_id)

tigit_score_df <- tigit_df %>%
  filter(sample_id %in% common_ts) %>%
  left_join(surv_locked %>%
              dplyr::select(sample_id, locked_score),
            by = "sample_id")

r_tigit_score <- cor(tigit_score_df$tigit,
                     tigit_score_df$locked_score,
                     use = "complete.obs")
cat("\nTIGIT vs immune-cold score: r =",
    round(r_tigit_score, 3), "\n")
cat("Negative = TIGIT higher in immune-hot\n")

p_tigit_score <- ggplot(tigit_score_df,
                        aes(x = locked_score, y = tigit)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#534AB7") +
  geom_smooth(method = "lm", se = TRUE,
              color = "black", linewidth = 0.8) +
  labs(
    title    = "TIGIT expression vs immune-cold score — METABRIC TNBC",
    subtitle = paste0("r = ", round(r_tigit_score, 3),
                      " | n = ", nrow(tigit_score_df)),
    x        = "Immune-cold score",
    y        = "TIGIT expression"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_tigit_score)
ggsave(file.path(fig_dir, "62_TIGIT_vs_immune_score.png"),
       p_tigit_score, width = 7, height = 5, dpi = 150)

# ── 7. TIGIT + IMMUNE-COLD SCORE COMBINED ─────────────────────────────────────
# The key clinical question: does TIGIT high + immune-hot = best prognosis?
# And TIGIT low + immune-cold = worst prognosis?

tigit_combined_df <- tigit_df %>%
  filter(!is.na(immune_group)) %>%
  mutate(
    combined_group = case_when(
      immune_group == "Immune-hot" & tigit_group == "TIGIT-high"
      ~ "Hot + TIGIT-high\n(potential responders)",
      immune_group == "Immune-hot" & tigit_group == "TIGIT-low"
      ~ "Hot + TIGIT-low",
      immune_group == "Immune-cold" & tigit_group == "TIGIT-high"
      ~ "Cold + TIGIT-high",
      immune_group == "Immune-cold" & tigit_group == "TIGIT-low"
      ~ "Cold + TIGIT-low\n(worst prognosis?)"
    )
  )

cat("\nCombined group sizes:\n")
print(table(tigit_combined_df$combined_group))

km_combined <- survfit(
  Surv(OS_time, OS_event) ~ combined_group,
  data = tigit_combined_df
)

p_km_combined <- ggsurvplot(
  km_combined,
  data              = tigit_combined_df,
  pval              = TRUE,
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.3,
  palette           = c("#1D9E75", "#85B7EB",
                        "#F0997B", "#D85A30"),
  title    = "OS by immune group and TIGIT — METABRIC TNBC",
  subtitle = "Combined stratification: immune status + TIGIT expression",
  xlab     = "Time (months)",
  ylab     = "Survival probability",
  ggtheme  = theme_minimal(base_size = 12)
)

png(file.path(fig_dir, "63_KM_TIGIT_combined.png"),
    width = 10, height = 8, units = "in", res = 150)
print(p_km_combined)
dev.off()
print(p_km_combined)

# ── 8. SAVE ───────────────────────────────────────────────────────────────────
save(tigit_df, tigit_combined_df, cox_tigit,
     file = file.path(proc_dir, "24_TIGIT_analysis.RData"))

cat("\n TIGIT analysis complete!\n")
cat("Figures 60-63 saved\n")