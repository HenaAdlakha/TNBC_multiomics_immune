# ── 06_survival_analysis.R ───────────────────────────────────────────────────
# Goals:
#   1. Kaplan-Meier survival curves: immune-cold vs immune-hot
#   2. Cox proportional hazards model
#   3. Survival per individual immune signature
#   4. Check group balance
# Input:  03_mofa_trained.RData, 01_TNBC_clean.RData, 
#         04_step1_immune_analysis.RData
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
raw_dir  <- file.path(base_dir, "data/raw")
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
load(file.path(proc_dir, "03_mofa_trained.RData"))
load(file.path(proc_dir, "01_TNBC_clean.RData"))
load(file.path(proc_dir, "04_step1_immune_analysis.RData"))

# ── 1. PREPARE SURVIVAL DATA ──────────────────────────────────────────────────
cat("Preparing survival data...\n")

# Check column names
cat("Survival columns:\n")
print(colnames(survival))

# Trim sample IDs
survival$sample <- substr(survival$sample, 1, 15)

# Get Factor 6 scores
factor_scores <- get_factors(mofa_trained, factors = "all")[[1]]
f6 <- factor_scores[, "Factor6"]

# Align with survival data
common_surv <- intersect(names(f6), survival$sample)
cat("TNBC samples with survival data:", length(common_surv), "\n")

# Build survival dataframe
surv_df <- survival %>%
  filter(sample %in% common_surv) %>%
  mutate(
    Factor6     = f6[sample],
    immune_group = ifelse(Factor6 > median(f6[common_surv]),
                          "Immune-cold", "Immune-hot")
  )

cat("Group sizes:\n")
print(table(surv_df$immune_group))

# Check survival column names — adjust if different
cat("\nSurvival data preview:\n")
print(head(surv_df[, c("sample", "OS", "OS.time", "immune_group")], 3))

# ── 2. KAPLAN-MEIER: OVERALL SURVIVAL ─────────────────────────────────────────
cat("Running Kaplan-Meier analysis...\n")

# Create survival object
surv_obj <- Surv(time  = surv_df$OS.time,
                 event = surv_df$OS)

# Fit KM curves
km_fit <- survfit(surv_obj ~ immune_group, data = surv_df)

# Plot
p_km_os <- ggsurvplot(
  km_fit,
  data          = surv_df,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.25,
  palette       = c("#D85A30", "#1D9E75"),
  title         = "Overall survival: immune-cold vs immune-hot TNBC",
  xlab          = "Time (days)",
  ylab          = "Survival probability",
  legend.labs   = c("Immune-cold", "Immune-hot"),
  legend.title  = "Group",
  ggtheme       = theme_minimal(base_size = 13)
)

print(p_km_os)

ggsave(file.path(fig_dir, "14_KM_overall_survival.png"),
       print(p_km_os),
       width = 9, height = 7, dpi = 150)

# ── 3. KAPLAN-MEIER: PROGRESSION FREE INTERVAL ────────────────────────────────
km_fit_pfi <- survfit(
  Surv(PFI.time, PFI) ~ immune_group,
  data = surv_df
)

p_km_pfi <- ggsurvplot(
  km_fit_pfi,
  data          = surv_df,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.25,
  palette       = c("#D85A30", "#1D9E75"),
  title         = "Progression-free interval: immune-cold vs immune-hot TNBC",
  xlab          = "Time (days)",
  ylab          = "Survival probability",
  legend.labs   = c("Immune-cold", "Immune-hot"),
  legend.title  = "Group",
  ggtheme       = theme_minimal(base_size = 13)
)

print(p_km_pfi)
ggsave(file.path(fig_dir, "15_KM_PFI.png"),
       print(p_km_pfi),
       width = 9, height = 7, dpi = 150)

# ── 4. COX PROPORTIONAL HAZARDS MODEL ─────────────────────────────────────────
cat("\nRunning Cox model...\n")

# Add clinical covariates
surv_df_cox <- surv_df %>%
  left_join(
    clinical %>%
      mutate(sampleID = substr(sampleID, 1, 15)) %>%
      select(sampleID,
             age    = Age_at_Initial_Pathologic_Diagnosis_nature2012,
             stage  = AJCC_Stage_nature2012),
    by = c("sample" = "sampleID")
  ) %>%
  mutate(
    age   = as.numeric(age),
    stage = as.factor(stage),
    group = ifelse(immune_group == "Immune-cold", 1, 0)
  )

# Simple Cox — Factor 6 group only
cox_simple <- coxph(
  Surv(OS.time, OS) ~ group,
  data = surv_df_cox
)
cat("\n── Simple Cox model ──\n")
print(summary(cox_simple))

# Adjusted Cox — controlling for age and stage
cox_adjusted <- coxph(
  Surv(OS.time, OS) ~ group + age + stage,
  data = surv_df_cox
)
cat("\n── Adjusted Cox model (age + stage) ──\n")
print(summary(cox_adjusted))

# Forest plot
p_forest <- ggforest(cox_adjusted,
                     data  = surv_df_cox,
                     main  = "Cox model: OS in TNBC\n(adjusted for age and stage)")
print(p_forest)
ggsave(file.path(fig_dir, "16_cox_forest_plot.png"),
       p_forest, width = 9, height = 6, dpi = 150)

# ── 5. SURVIVAL PER INDIVIDUAL IMMUNE SIGNATURE ───────────────────────────────
cat("\nSurvival analysis per immune signature...\n")

# Use the 3 significant signatures from Step 1
sig_sigs <- c("MHC1_21978456", "TGFB_score_21050467", "MHC2_21978456")

# Align immune scores with survival samples
immune_surv <- immune_tnbc[intersect(rownames(immune_tnbc), common_surv), ]
surv_sub    <- surv_df %>% filter(sample %in% rownames(immune_surv))

km_plots <- lapply(sig_sigs, function(sig) {
  scores <- immune_surv[surv_sub$sample, sig]
  surv_sub$sig_group <- ifelse(scores > median(scores, na.rm = TRUE),
                               "High", "Low")
  fit <- survfit(Surv(OS.time, OS) ~ sig_group, data = surv_sub)
  
  ggsurvplot(
    fit,
    data        = surv_sub,
    pval        = TRUE,
    conf.int    = TRUE,
    palette     = c("#534AB7", "#D85A30"),
    title       = paste0(sig, "\nOS by high vs low"),
    xlab        = "Time (days)",
    ylab        = "Survival probability",
    legend.labs = c("High", "Low"),
    ggtheme     = theme_minimal(base_size = 11)
  )
})

# Save each individually
for(i in seq_along(sig_sigs)) {
  ggsave(
    file.path(fig_dir, paste0("17_KM_", sig_sigs[i], ".png")),
    print(km_plots[[i]]),
    width = 7, height = 6, dpi = 150
  )
}

# ── 6. SAVE ───────────────────────────────────────────────────────────────────
save(surv_df, surv_df_cox, cox_simple, cox_adjusted,
     file = file.path(proc_dir, "05_survival_results.RData"))

cat("\n Survival analysis complete!\n")
cat("Figures 14-17 saved to outputs/figures/\n")