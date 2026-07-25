# ── 11_TNBC_3omics_survival.R ─────────────────────────────────────────────────
# Survival analysis on expanded TNBC cohort (n=81)
# Using Factor 4 as immune-cold signature
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(MOFA2)
library(survival)
library(survminer)
library(ggplot2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

f <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/raw", fn))

load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "01_TNBC_clean.RData"))

# ── 1. PREPARE SURVIVAL DATA ──────────────────────────────────────────────────
factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
f4 <- factor_scores3[, "Factor4"]

survival$sample <- substr(survival$sample, 1, 15)
common_surv3 <- intersect(names(f4), survival$sample)
cat("TNBC samples with survival data:", length(common_surv3), "\n")

surv_df3 <- survival %>%
  filter(sample %in% common_surv3) %>%
  mutate(
    Factor4      = f4[sample],
    immune_group = ifelse(Factor4 > median(f4[common_surv3]),
                          "Immune-cold", "Immune-hot")
  )

cat("Group sizes:\n")
print(table(surv_df3$immune_group))
cat("Total events (deaths):", sum(surv_df3$OS), "\n")

# ── 2. KAPLAN-MEIER: OVERALL SURVIVAL ─────────────────────────────────────────
km_fit3 <- survfit(Surv(OS.time, OS) ~ immune_group, data = surv_df3)

p_km3 <- ggsurvplot(
  km_fit3,
  data              = surv_df3,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#D85A30", "#1D9E75"),
  title             = "Overall survival: immune-cold vs immune-hot TNBC (n=81)",
  subtitle          = "Factor 4 median split — 3 omics model",
  xlab              = "Time (days)",
  ylab              = "Survival probability",
  legend.labs       = c("Immune-cold", "Immune-hot"),
  ggtheme           = theme_minimal(base_size = 13)
)
print(p_km3)
ggsave(file.path(fig_dir, "28_KM_OS_n81.png"),
       print(p_km3), width = 9, height = 7, dpi = 150)

# ── 3. KAPLAN-MEIER: PFI ──────────────────────────────────────────────────────
km_pfi3 <- survfit(Surv(PFI.time, PFI) ~ immune_group, data = surv_df3)

p_pfi3 <- ggsurvplot(
  km_pfi3,
  data              = surv_df3,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#D85A30", "#1D9E75"),
  title             = "Progression-free interval: immune-cold vs immune-hot TNBC (n=81)",
  xlab              = "Time (days)",
  ylab              = "Survival probability",
  legend.labs       = c("Immune-cold", "Immune-hot"),
  ggtheme           = theme_minimal(base_size = 13)
)
print(p_pfi3)
ggsave(file.path(fig_dir, "29_KM_PFI_n81.png"),
       print(p_pfi3), width = 9, height = 7, dpi = 150)

# ── 4. MHC1 SURVIVAL IN n=81 ──────────────────────────────────────────────────
# Reload immune scores for n=81 samples
immune_sigs <- fread(
  file = f("TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena"),
  header = TRUE)

immune_mat <- immune_sigs %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  as.matrix()

immune_t <- t(immune_mat)
rownames(immune_t) <- substr(rownames(immune_t), 1, 15)

common_imm3 <- intersect(rownames(immune_t), common_tnbc_3)
immune_tnbc3 <- immune_t[common_imm3, ]

# MHC1 survival
surv_imm3 <- surv_df3 %>%
  filter(sample %in% rownames(immune_tnbc3))

mhc1_scores <- immune_tnbc3[surv_imm3$sample, "MHC1_21978456"]
surv_imm3$mhc1_group <- ifelse(mhc1_scores > median(mhc1_scores, na.rm=TRUE),
                               "MHC1-high", "MHC1-low")

km_mhc1_3 <- survfit(Surv(OS.time, OS) ~ mhc1_group, data = surv_imm3)

p_mhc1_3 <- ggsurvplot(
  km_mhc1_3,
  data              = surv_imm3,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#534AB7", "#D85A30"),
  title             = "OS by MHC1 expression — TNBC (n=81)",
  xlab              = "Time (days)",
  ylab              = "Survival probability",
  ggtheme           = theme_minimal(base_size = 13)
)
print(p_mhc1_3)
ggsave(file.path(fig_dir, "30_KM_MHC1_n81.png"),
       print(p_mhc1_3), width = 9, height = 7, dpi = 150)

# ── 5. COX MODEL ──────────────────────────────────────────────────────────────
surv_df3_cox <- surv_df3 %>%
  left_join(
    clinical %>%
      mutate(sampleID = substr(sampleID, 1, 15)) %>%
      select(sampleID,
             age   = Age_at_Initial_Pathologic_Diagnosis_nature2012,
             stage = AJCC_Stage_nature2012),
    by = c("sample" = "sampleID")
  ) %>%
  mutate(
    age          = as.numeric(age),
    stage_binary = ifelse(grepl("III|IV", as.character(stage)),
                          "Late", "Early"),
    stage_binary = as.factor(stage_binary),
    group        = ifelse(immune_group == "Immune-cold", 1, 0)
  )

cox3 <- coxph(
  Surv(OS.time, OS) ~ group + age + stage_binary,
  data = surv_df3_cox
)

cat("\n── Cox model (n=81) ──\n")
print(summary(cox3))

# ── 6. SAVE ───────────────────────────────────────────────────────────────────
save(surv_df3, surv_df3_cox, cox3,
     file = file.path(proc_dir, "10_survival_n81.RData"))

cat("\n Survival analysis n=81 complete!\n")
cat("Check figures 28-30 in outputs/figures/\n")