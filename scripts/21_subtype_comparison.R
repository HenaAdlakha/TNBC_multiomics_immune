# ── 21_subtype_comparison.R ───────────────────────────────────────────────────
# Goals:
#   1. Test Factor 4 association with PAM50 subtypes within TNBC
#   2. Test whether Factor 4 predicts survival WITHIN PAM50 subtype
#   3. This is the key novelty question: does Factor 4 add information
#      beyond established subtype classification?
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(survival)
library(survminer)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "01_TNBC_clean.RData"))

# ── 1. GET PAM50 CALLS FOR TNBC SAMPLES ───────────────────────────────────────
cat("PAM50 distribution in full cohort:\n")
print(table(clinical$PAM50Call_RNAseq, useNA = "always"))

pam50_df <- clinical %>%
  mutate(sampleID = substr(sampleID, 1, 15)) %>%
  dplyr::select(sampleID, pam50 = PAM50Call_RNAseq) %>%
  filter(sampleID %in% common_samples)

cat("\nPAM50 distribution in TNBC samples (n=81):\n")
print(table(pam50_df$pam50, useNA = "always"))

# ── 2. FACTOR 4 BY PAM50 SUBTYPE ──────────────────────────────────────────────
factor_scores3 <- get_factors(mofa_trained_tnbc3, factors = "all")[[1]]
f4 <- factor_scores3[, "Factor4"]

pam50_factor_df <- data.frame(
  sample = names(f4),
  f4     = f4
) %>%
  left_join(pam50_df, by = c("sample" = "sampleID")) %>%
  filter(!is.na(pam50))

cat("\nSamples with PAM50 calls:", nrow(pam50_factor_df), "\n")
cat("PAM50 distribution:\n")
print(table(pam50_factor_df$pam50))

# Kruskal-Wallis: does Factor 4 vary by PAM50?
kw_pam50 <- kruskal.test(f4 ~ pam50, data = pam50_factor_df)
cat("\nKruskal-Wallis: Factor 4 across PAM50 subtypes:\n")
print(kw_pam50)

# Plot
p_pam50 <- ggplot(pam50_factor_df,
                  aes(x = pam50, y = f4, fill = pam50)) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  labs(
    title    = "Factor 4 score by PAM50 subtype — TNBC",
    subtitle = paste0("Kruskal-Wallis p = ",
                      round(kw_pam50$p.value, 3)),
    x = "PAM50 subtype", y = "Factor 4 score", fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "none")

print(p_pam50)
ggsave(file.path(fig_dir, "55_factor4_by_pam50.png"),
       p_pam50, width = 7, height = 5, dpi = 150)

# ── 3. DOES FACTOR 4 PREDICT SURVIVAL WITHIN PAM50 SUBTYPE? ──────────────────
# Merge with survival
survival$sample <- substr(survival$sample, 1, 15)

surv_pam50 <- data.frame(
  sample = names(f4),
  f4     = f4
) %>%
  left_join(pam50_df,  by = c("sample" = "sampleID")) %>%
  left_join(survival %>%
              dplyr::select(sample, OS, OS.time, PFI, PFI.time),
            by = "sample") %>%
  filter(!is.na(pam50), !is.na(OS.time)) %>%
  mutate(
    immune_group = ifelse(f4 > median(f4), "Immune-cold", "Immune-hot"),
    f4_std       = as.numeric(scale(f4))
  )

cat("\nSamples with PAM50 + survival:", nrow(surv_pam50), "\n")
cat("Events:", sum(surv_pam50$OS, na.rm=TRUE), "\n")

# Cox: Factor 4 + PAM50 subtype
cox_pam50 <- coxph(
  Surv(OS.time, OS) ~ f4_std + pam50,
  data = surv_pam50
)
cat("\n── Cox: Factor 4 + PAM50 subtype ──\n")
print(summary(cox_pam50))

# Key question: is Factor 4 significant after adjusting for PAM50?
f4_coef <- summary(cox_pam50)$coefficients["f4_std", ]
cat("\nFactor 4 after adjusting for PAM50:\n")
cat("HR =", round(exp(f4_coef["coef"]), 3), "\n")
cat("p  =", round(f4_coef["Pr(>|z|)"], 4), "\n")

# ── 4. SURVIVAL WITHIN BASAL-LIKE ONLY ────────────────────────────────────────
# Most TNBC is Basal-like — does Factor 4 stratify within Basal-like?
basal_df <- surv_pam50 %>%
  filter(pam50 == "Basal")

cat("\nBasal-like TNBC samples:", nrow(basal_df), "\n")
cat("Events in Basal-like:", sum(basal_df$OS, na.rm=TRUE), "\n")

if(nrow(basal_df) >= 20 & sum(basal_df$OS, na.rm=TRUE) >= 5) {
  km_basal <- survfit(Surv(OS.time, OS) ~ immune_group,
                      data = basal_df)
  
  p_km_basal <- ggsurvplot(
    km_basal,
    data        = basal_df,
    pval        = TRUE,
    pval.method = TRUE,
    conf.int    = TRUE,
    risk.table  = TRUE,
    risk.table.height = 0.25,
    palette     = c("#D85A30", "#1D9E75"),
    title       = "OS within Basal-like TNBC — immune-cold vs immune-hot",
    subtitle    = "Factor 4 stratification within PAM50 Basal-like subtype",
    xlab        = "Time (days)",
    ylab        = "Survival probability",
    legend.labs = c("Immune-cold", "Immune-hot"),
    ggtheme     = theme_minimal(base_size = 13)
  )
  
  png(file.path(fig_dir, "56_KM_basal_like.png"),
      width = 9, height = 7, units = "in", res = 150)
  print(p_km_basal)
  dev.off()
  print(p_km_basal)
}

# ── 5. METABRIC PAM50 COMPARISON ──────────────────────────────────────────────
# METABRIC has PAM50 in clinical data
fm <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/metabric", fn))

clin_pat_meta <- fread(file = fm("data_clinical_patient.txt"),
                       header = TRUE)
colnames(clin_pat_meta)[1] <- "patient_id"

cat("\nMETABRIC PAM50 columns:\n")
print(colnames(clin_pat_meta)[grep("PAM50|subtype|Pam",
                                   colnames(clin_pat_meta),
                                   ignore.case = TRUE)])

load(file.path(proc_dir, "20_locked_validation.RData"))

# Merge locked score with PAM50
surv_locked_pam50 <- surv_locked %>%
  left_join(
    clin_pat_meta %>%
      dplyr::select(patient_id,
                    pam50 = `Pam50 + Claudin-low subtype`),
    by = "patient_id"
  )

cat("\nMETABRIC PAM50 distribution in TNBC:\n")
print(table(surv_locked_pam50$pam50, useNA = "always"))

# Cox within METABRIC: locked score + PAM50
cox_meta_pam50 <- coxph(
  Surv(OS_time, OS_event) ~ score_std + pam50 + age,
  data = surv_locked_pam50 %>% filter(!is.na(pam50))
)

cat("\n── METABRIC Cox: locked score + PAM50 + age ──\n")
print(summary(cox_meta_pam50))

score_coef <- summary(cox_meta_pam50)$coefficients["score_std", ]
cat("\nLocked score after adjusting for PAM50 + age:\n")
cat("HR =", round(exp(score_coef["coef"]), 3), "\n")
cat("p  =", round(score_coef["Pr(>|z|)"], 4), "\n")

# ── 6. SAVE ───────────────────────────────────────────────────────────────────
save(pam50_factor_df, surv_pam50, cox_pam50,
     surv_locked_pam50, cox_meta_pam50,
     file = file.path(proc_dir, "21_subtype_comparison.RData"))

cat("\n Script 21 complete!\n")