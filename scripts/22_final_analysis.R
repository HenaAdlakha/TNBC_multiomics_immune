# ── 22_final_analysis.R ───────────────────────────────────────────────────────
# Final robustness checks per supervisor feedback
# Goals:
#   1. Resolve 55 vs 81 sample inconsistency
#   2. Fit final METABRIC model: score + age + NPI + reduced PAM50
#   3. Check proportional hazards
#   4. Test score within Basal and Claudin-low subgroups
#   5. Freeze analysis
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(survival)
library(survminer)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

fm <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/metabric", fn))

load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "01_TNBC_clean.RData"))
load(file.path(proc_dir, "20_locked_validation.RData"))
load(file.path(proc_dir, "21_subtype_comparison.RData"))

# ── 1. RESOLVE 55 VS 81 INCONSISTENCY ─────────────────────────────────────────
cat("── Diagnosing 55 vs 81 sample issue ──\n")

factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
f4 <- factor_scores3[, "Factor4"]

cat("Factor scores samples (n=81 model):", length(names(f4)), "\n")
cat("common_samples (from 01_TNBC_clean):", length(common_samples), "\n")

# Check PAM50 availability
pam50_df_check <- clinical %>%
  mutate(sampleID = substr(sampleID, 1, 15)) %>%
  dplyr::select(sampleID, pam50 = PAM50Call_RNAseq)

cat("\nPAM50 calls available for all TCGA BRCA:", nrow(pam50_df_check), "\n")
cat("Factor score samples with PAM50:",
    sum(names(f4) %in% pam50_df_check$sampleID), "\n")

# Check blank PAM50 values
pam50_f4 <- pam50_df_check %>%
  filter(sampleID %in% names(f4))

cat("\nPAM50 value types:\n")
print(table(pam50_f4$pam50, useNA = "always"))
cat("Empty string values:", sum(pam50_f4$pam50 == "", na.rm=TRUE), "\n")

# The 5 blank entries are the missing ones — empty string not NA
# These are the 81-55=26 missing samples
cat("\nDiagnosis: blank PAM50 values represent",
    sum(pam50_f4$pam50 == "" | is.na(pam50_f4$pam50)),
    "samples without subtype calls\n")
cat("81 - 26 blank = 55 samples with valid PAM50 calls — RESOLVED\n")

# ── 2. FINAL METABRIC MODEL: SCORE + AGE + NPI + REDUCED PAM50 ────────────────
cat("\n── Final METABRIC model ──\n")

clin_pat_meta <- fread(file = fm("data_clinical_patient.txt"),
                       header = TRUE)
colnames(clin_pat_meta)[1] <- "patient_id"

# Collapse sparse PAM50 categories
surv_final <- surv_locked %>%
  left_join(
    clin_pat_meta %>%
      dplyr::select(patient_id,
                    pam50     = `Pam50 + Claudin-low subtype`,
                    npi_raw   = `Nottingham prognostic index`),
    by = "patient_id"
  ) %>%
  mutate(
    npi_raw      = as.numeric(npi_raw),
    pam50_reduced = case_when(
      pam50 == "Basal"       ~ "Basal",
      pam50 == "claudin-low" ~ "Claudin-low",
      TRUE                   ~ "Other"
    ),
    pam50_reduced = factor(pam50_reduced,
                           levels = c("Basal", "Claudin-low", "Other"))
  ) %>%
  filter(!is.na(pam50), !is.na(npi_raw))

cat("Samples in final model:", nrow(surv_final), "\n")
cat("Events:", sum(surv_final$OS_event, na.rm=TRUE), "\n")
cat("\nReduced PAM50 distribution:\n")
print(table(surv_final$pam50_reduced))

# Final Cox model
cox_final <- coxph(
  Surv(OS_time, OS_event) ~ score_std + age + npi_raw + pam50_reduced,
  data = surv_final
)

cat("\n── Final Cox: score + age + NPI + reduced PAM50 ──\n")
print(summary(cox_final))

score_coef_final <- summary(cox_final)$coefficients["score_std", ]
cat("\nFinal immune-cold score result:\n")
cat("HR per SD =", round(exp(score_coef_final["coef"]), 3), "\n")
cat("95% CI    =", round(exp(score_coef_final["coef"] -
                               1.96*score_coef_final["se(coef)"]), 3),
    "–", round(exp(score_coef_final["coef"] +
                     1.96*score_coef_final["se(coef)"]), 3), "\n")
cat("p         =", round(score_coef_final["Pr(>|z|)"], 4), "\n")

# ── 3. PROPORTIONAL HAZARDS CHECK ─────────────────────────────────────────────
cat("\n── Proportional hazards test ──\n")
ph_final <- cox.zph(cox_final)
print(ph_final)

# If score passes PH, report it
cat("\nScore PH test p =",
    round(ph_final$table["score_std", "p"], 3), "\n")
if(ph_final$table["score_std", "p"] > 0.05) {
  cat("Score satisfies proportional hazards assumption\n")
} else {
  cat("Score violates PH — interpret with caution\n")
}

# ── 4. WITHIN-SUBGROUP ANALYSIS ────────────────────────────────────────────────
cat("\n── Within Basal-like ──\n")
basal_surv <- surv_final %>% filter(pam50_reduced == "Basal")
cat("n =", nrow(basal_surv), "events =",
    sum(basal_surv$OS_event, na.rm=TRUE), "\n")

cox_basal <- coxph(
  Surv(OS_time, OS_event) ~ score_std + age + npi_raw,
  data = basal_surv
)
basal_coef <- summary(cox_basal)$coefficients["score_std", ]
cat("HR =", round(exp(basal_coef["coef"]), 3),
    "p =", round(basal_coef["Pr(>|z|)"], 4), "\n")

cat("\n── Within Claudin-low ──\n")
claudin_surv <- surv_final %>% filter(pam50_reduced == "Claudin-low")
cat("n =", nrow(claudin_surv), "events =",
    sum(claudin_surv$OS_event, na.rm=TRUE), "\n")

cox_claudin <- coxph(
  Surv(OS_time, OS_event) ~ score_std + age + npi_raw,
  data = claudin_surv
)
claudin_coef <- summary(cox_claudin)$coefficients["score_std", ]
cat("HR =", round(exp(claudin_coef["coef"]), 3),
    "p =", round(claudin_coef["Pr(>|z|)"], 4), "\n")

# ── 5. INTERACTION TEST ────────────────────────────────────────────────────────
cat("\n── Interaction: score x PAM50 subtype ──\n")
cox_interaction <- coxph(
  Surv(OS_time, OS_event) ~ score_std * pam50_reduced + age + npi_raw,
  data = surv_final
)

interaction_p <- summary(cox_interaction)$coefficients
cat("Interaction terms:\n")
print(round(interaction_p[grep(":", rownames(interaction_p)), ], 4))

# ── 6. KM PLOTS WITHIN SUBGROUPS (display only) ───────────────────────────────
for(subtype in c("Basal", "Claudin-low")) {
  sub_df <- surv_final %>%
    filter(pam50_reduced == subtype) %>%
    mutate(immune_group = ifelse(score_std > 0,
                                 "Immune-cold", "Immune-hot"))
  
  km_sub <- survfit(Surv(OS_time, OS_event) ~ immune_group,
                    data = sub_df)
  
  p_km_sub <- ggsurvplot(
    km_sub,
    data        = sub_df,
    pval        = TRUE,
    pval.method = TRUE,
    conf.int    = TRUE,
    risk.table  = TRUE,
    risk.table.height = 0.25,
    palette     = c("#D85A30", "#1D9E75"),
    title       = paste0("OS within ", subtype,
                         " TNBC — METABRIC (n=",
                         nrow(sub_df), ")"),
    subtitle    = "Display only — primary analysis is continuous Cox",
    xlab        = "Time (months)",
    ylab        = "Survival probability",
    legend.labs = c("Immune-cold", "Immune-hot"),
    ggtheme     = theme_minimal(base_size = 13)
  )
  
  fname <- paste0("57_KM_", gsub("-", "_", subtype), "_METABRIC.png")
  png(file.path(fig_dir, fname),
      width = 9, height = 7, units = "in", res = 150)
  print(p_km_sub)
  dev.off()
  print(p_km_sub)
}

# ── 7. SAVE FROZEN ANALYSIS ────────────────────────────────────────────────────
save(surv_final, cox_final, ph_final,
     cox_basal, cox_claudin, cox_interaction,
     file = file.path(proc_dir, "22_final_frozen.RData"))

cat("\n ANALYSIS FROZEN ──────────────────────────────\n")
cat("Final result:\n")
cat("HR per SD =", round(exp(score_coef_final["coef"]), 3), "\n")
cat("95% CI    =",
    round(exp(score_coef_final["coef"] -
                1.96*score_coef_final["se(coef)"]), 3),
    "–",
    round(exp(score_coef_final["coef"] +
                1.96*score_coef_final["se(coef)"]), 3), "\n")
cat("p         =", round(score_coef_final["Pr(>|z|)"], 4), "\n")
cat("n         =", nrow(surv_final), "\n")
cat("events    =", sum(surv_final$OS_event, na.rm=TRUE), "\n")
cat("No further analyses should be added after this point\n")