# Fix 1: Continuous Factor 6 score in Cox
cox_continuous <- coxph(
  Surv(PFI.time, PFI) ~ Factor6,
  data = surv_df
)
cat("── Continuous Factor 6 Cox (PFI) ──\n")
print(summary(cox_continuous))

# Fix 2: Simplified stage binary Cox
surv_df_cox2 <- surv_df_cox %>%
  mutate(
    stage_binary = ifelse(grepl("III|IV", as.character(stage)),
                          "Late", "Early"),
    stage_binary = as.factor(stage_binary)
  )

cox_simple2 <- coxph(
  Surv(PFI.time, PFI) ~ group + age + stage_binary,
  data = surv_df_cox2
)
cat("\n── Simplified Cox (PFI, age, binary stage) ──\n")
print(summary(cox_simple2))

# Fix 3: MHC1 high vs low survival — your strongest immune signature
immune_surv <- immune_tnbc[intersect(rownames(immune_tnbc), common_surv), ]
surv_sub <- surv_df %>% filter(sample %in% rownames(immune_surv))
mhc1_scores <- immune_surv[surv_sub$sample, "MHC1_21978456"]

surv_sub$mhc1_group <- ifelse(mhc1_scores > median(mhc1_scores, na.rm=TRUE),
                              "MHC1-high", "MHC1-low")

km_mhc1 <- survfit(Surv(OS.time, OS) ~ mhc1_group, data = surv_sub)

p_mhc1 <- ggsurvplot(
  km_mhc1,
  data        = surv_sub,
  pval        = TRUE,
  pval.method = TRUE,
  conf.int    = TRUE,
  risk.table  = TRUE,
  risk.table.height = 0.25,
  palette     = c("#534AB7", "#D85A30"),
  title       = "OS by MHC1 expression — TNBC",
  xlab        = "Time (days)",
  ylab        = "Survival probability",
  ggtheme     = theme_minimal(base_size = 13)
)
print(p_mhc1)
ggsave(file.path(fig_dir, "18_KM_MHC1.png"),
       print(p_mhc1), width = 9, height = 7, dpi = 150)

mhc2_scores <- immune_surv[surv_sub$sample, "MHC2_21978456"]
surv_sub$mhc2_group <- ifelse(mhc2_scores > median(mhc2_scores, na.rm=TRUE),
                              "MHC2-high", "MHC2-low")

km_mhc2 <- survfit(Surv(OS.time, OS) ~ mhc2_group, data = surv_sub)

p_mhc2 <- ggsurvplot(
  km_mhc2,
  data        = surv_sub,
  pval        = TRUE,
  pval.method = TRUE,
  conf.int    = TRUE,
  risk.table  = TRUE,
  risk.table.height = 0.25,
  palette     = c("#534AB7", "#D85A30"),
  title       = "OS by MHC2 expression — TNBC",
  xlab        = "Time (days)",
  ylab        = "Survival probability",
  ggtheme     = theme_minimal(base_size = 13)
)
print(p_mhc2)
ggsave(file.path(fig_dir, "19_KM_MHC2.png"),
       print(p_mhc2), width = 9, height = 7, dpi = 150)