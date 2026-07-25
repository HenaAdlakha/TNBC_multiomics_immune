# ── 12_METABRIC_validation.R ──────────────────────────────────────────────────
# Independent validation of immune-cold TNBC signature in METABRIC
# Goals:
#   1. Filter METABRIC to TNBC samples
#   2. Apply Factor 4 immune-cold gene signature
#   3. Kaplan-Meier survival analysis
#   4. Compare with TCGA findings
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(ggplot2)
library(MOFA2)

base_dir  <- "~/Documents/R Working Directory/BRCA_project"
meta_dir  <- file.path(base_dir, "data/metabric")
proc_dir  <- file.path(base_dir, "data/processed")
fig_dir   <- file.path(base_dir, "outputs/figures")

fm <- function(fn) path.expand(file.path(meta_dir, fn))

# ── 1. PEEK AT FILES ──────────────────────────────────────────────────────────
cat("── Clinical patient columns ──\n")
clin_pat <- fread(file = fm("data_clinical_patient.txt"), 
                  header = TRUE, nrows = 3)
print(colnames(clin_pat))

cat("\n── Clinical sample columns ──\n")
clin_sam <- fread(file = fm("data_clinical_sample.txt"),
                  header = TRUE, nrows = 3)
print(colnames(clin_sam))

cat("\n── Expression file preview ──\n")
expr_peek <- fread(file = fm("data_mrna_illumina_microarray.txt"),
                   header = TRUE, nrows = 3)
cat("Dimensions:", dim(expr_peek), "\n")
print(colnames(expr_peek)[1:5])

# ── 2. LOAD ALL METABRIC DATA ─────────────────────────────────────────────────
cat("Loading METABRIC data...\n")

# Load full clinical files
clin_pat <- fread(file = fm("data_clinical_patient.txt"), header = TRUE)
clin_sam <- fread(file = fm("data_clinical_sample.txt"), header = TRUE)

# Fix column names
colnames(clin_pat)[1] <- "patient_id"
colnames(clin_sam)[1] <- "patient_id"

cat("Clinical patients:", nrow(clin_pat), "\n")
cat("Clinical samples:", nrow(clin_sam), "\n")

# Load expression
cat("Loading expression data...\n")
expr_raw <- fread(file = fm("data_mrna_illumina_microarray.txt"),
                  header = TRUE)

cat("Expression dimensions:", dim(expr_raw), "\n")

# ── 3. IDENTIFY TNBC SAMPLES ──────────────────────────────────────────────────
cat("\nER Status values:\n")
print(table(clin_sam$`ER Status`))
cat("\nPR Status values:\n")
print(table(clin_sam$`PR Status`))
cat("\nHER2 Status values:\n")
print(table(clin_sam$`HER2 Status`))

# Filter TNBC
tnbc_meta <- clin_sam %>%
  filter(
    `ER Status`   == "Negative",
    `PR Status`   == "Negative",
    `HER2 Status` == "Negative"
  )

cat("\nTNBC samples in METABRIC:", nrow(tnbc_meta), "\n")

# ── 4. PREPARE EXPRESSION MATRIX ──────────────────────────────────────────────
# Set gene symbols as rownames
expr_mat <- expr_raw %>%
  as.data.frame() %>%
  filter(!duplicated(Hugo_Symbol)) %>%
  select(-Entrez_Gene_Id) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

cat("Expression matrix:", dim(expr_mat), "\n")

# Filter to TNBC samples
tnbc_samples_meta <- intersect(tnbc_meta$`Sample Identifier`,
                               colnames(expr_mat))
cat("TNBC samples with expression:", length(tnbc_samples_meta), "\n")

expr_tnbc <- expr_mat[, tnbc_samples_meta]

# ── 5. LOAD FACTOR 4 GENE SIGNATURE FROM TCGA ─────────────────────────────────
load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))

weights_f4 <- get_weights(mofa_trained_tnbc3,
                          factors = "Factor4",
                          as.data.frame = TRUE)

# Get top RNA genes — positive and negative weights separately
top_f4_positive <- weights_f4 %>%
  filter(view == "RNA") %>%
  arrange(desc(value)) %>%
  head(50) %>%
  pull(feature) %>%
  gsub("_RNA", "", .)

top_f4_negative <- weights_f4 %>%
  filter(view == "RNA") %>%
  arrange(value) %>%
  head(50) %>%
  pull(feature) %>%
  gsub("_RNA", "", .)

cat("\nTop positive Factor 4 genes (immune-hot markers):\n")
print(head(top_f4_positive, 10))

cat("\nTop negative Factor 4 genes (immune-cold markers):\n")
print(head(top_f4_negative, 10))

# ── 6. CALCULATE IMMUNE-COLD SCORE IN METABRIC ────────────────────────────────
# Immune-cold score = mean of negative weight genes - mean of positive weight genes
# Higher score = more immune-cold

# Check which genes are present in METABRIC
pos_present <- intersect(top_f4_positive, rownames(expr_tnbc))
neg_present  <- intersect(top_f4_negative, rownames(expr_tnbc))

cat("\nPositive genes found in METABRIC:", length(pos_present), "/50\n")
cat("Negative genes found in METABRIC:", length(neg_present), "/50\n")

# Calculate scores
pos_score <- colMeans(expr_tnbc[pos_present, ], na.rm = TRUE)
neg_score  <- colMeans(expr_tnbc[neg_present, ], na.rm = TRUE)

# Immune-cold score: high = immune-cold
# pos_score = positive weight genes = expressed in immune-cold tumors
# neg_score = negative weight genes = expressed in immune-hot tumors
# Therefore: immune-cold score = pos_score - neg_score
immune_cold_score <- pos_score - neg_score

cat("\nImmune-cold score distribution:\n")
print(summary(immune_cold_score))

# ── 7. ALSO CALCULATE IFN-GAMMA SCORE ─────────────────────────────────────────
# IFN-gamma pathway genes — directly test the primary mechanism
ifng_genes <- c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10",
                "CXCL11", "IDO1", "HLA-A", "HLA-B", "HLA-C",
                "B2M", "TAP1", "TAP2", "PSMB9", "PSMB8")

ifng_present <- intersect(ifng_genes, rownames(expr_tnbc))
cat("\nIFN-gamma genes found:", length(ifng_present), "/", length(ifng_genes), "\n")

ifng_score <- colMeans(expr_tnbc[ifng_present, ], na.rm = TRUE)

# ── 8. MERGE WITH SURVIVAL DATA ───────────────────────────────────────────────
# Match sample to patient IDs
sample_to_patient <- clin_sam %>%
  select(patient_id, sample_id = `Sample Identifier`)

surv_meta <- data.frame(
  sample_id         = tnbc_samples_meta,
  immune_cold_score = immune_cold_score[tnbc_samples_meta],
  ifng_score        = ifng_score[tnbc_samples_meta]
) %>%
  left_join(sample_to_patient, by = "sample_id") %>%
  left_join(clin_pat %>%
              select(patient_id,
                     OS_months = `Overall Survival (Months)`,
                     OS_status = `Overall Survival Status`,
                     RFS_months = `Relapse Free Status (Months)`,
                     RFS_status = `Relapse Free Status`),
            by = "patient_id") %>%
  mutate(
    OS_time   = as.numeric(OS_months),
    OS_event  = ifelse(OS_status == "1:DECEASED", 1, 0),
    RFS_time  = as.numeric(RFS_months),
    RFS_event = ifelse(RFS_status == "1:Recurred", 1, 0),
    immune_group = ifelse(immune_cold_score > median(immune_cold_score,
                                                     na.rm = TRUE),
                          "Immune-cold", "Immune-hot"),
    ifng_group   = ifelse(ifng_score < median(ifng_score, na.rm = TRUE),
                          "IFNg-low", "IFNg-high")
  )

cat("\nSurvival data:\n")
cat("Samples:", nrow(surv_meta), "\n")
cat("Events (deaths):", sum(surv_meta$OS_event, na.rm = TRUE), "\n")
cat("Group sizes:\n")
print(table(surv_meta$immune_group))

# ── 9. KAPLAN-MEIER: IMMUNE-COLD SCORE ────────────────────────────────────────
km_meta <- survfit(Surv(OS_time, OS_event) ~ immune_group,
                   data = surv_meta)

p_km_meta <- ggsurvplot(
  km_meta,
  data              = surv_meta,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#D85A30", "#1D9E75"),
  title             = "OS: immune-cold vs immune-hot — METABRIC TNBC",
  subtitle          = "Independent validation cohort",
  xlab              = "Time (months)",
  ylab              = "Survival probability",
  legend.labs       = c("Immune-cold", "Immune-hot"),
  ggtheme           = theme_minimal(base_size = 13)
)

# Save correctly
png(file.path(fig_dir, "31_METABRIC_KM_immune_cold.png"),
    width = 9, height = 7, units = "in", res = 150)
print(p_km_meta)
dev.off()
print(p_km_meta)

# ── 10. KAPLAN-MEIER: IFN-GAMMA SCORE ─────────────────────────────────────────
km_ifng <- survfit(Surv(OS_time, OS_event) ~ ifng_group,
                   data = surv_meta)

p_km_ifng <- ggsurvplot(
  km_ifng,
  data              = surv_meta,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#534AB7", "#D85A30"),
  title             = "OS by IFN-gamma score — METABRIC TNBC",
  subtitle          = "IFNg-low = immune-cold phenotype",
  xlab              = "Time (months)",
  ylab              = "Survival probability",
  legend.labs       = c("IFNg-high", "IFNg-low"),
  ggtheme           = theme_minimal(base_size = 13)
)

png(file.path(fig_dir, "32_METABRIC_KM_IFNg.png"),
    width = 9, height = 7, units = "in", res = 150)
print(p_km_ifng)
dev.off()
print(p_km_ifng)

# ── 11. RELAPSE FREE SURVIVAL ──────────────────────────────────────────────────
km_rfs <- survfit(Surv(RFS_time, RFS_event) ~ immune_group,
                  data = surv_meta)

p_km_rfs <- ggsurvplot(
  km_rfs,
  data              = surv_meta,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#D85A30", "#1D9E75"),
  title             = "Relapse-free survival: immune-cold vs immune-hot — METABRIC TNBC",
  xlab              = "Time (months)",
  ylab              = "Survival probability",
  legend.labs       = c("Immune-cold", "Immune-hot"),
  ggtheme           = theme_minimal(base_size = 13)
)

png(file.path(fig_dir, "33_METABRIC_KM_RFS.png"),
    width = 9, height = 7, units = "in", res = 150)
print(p_km_rfs)
dev.off()
print(p_km_rfs)

# ── 12. COX MODEL ─────────────────────────────────────────────────────────────
# Corrected Cox — using NPI instead of Tumor Stage
surv_meta_cox <- surv_meta %>%
  left_join(
    clin_pat %>%
      select(patient_id,
             age = `Age at Diagnosis`,
             npi = `Nottingham prognostic index`),
    by = "patient_id"
  ) %>%
  mutate(
    age   = as.numeric(age),
    npi   = as.numeric(npi),
    group = ifelse(immune_group == "Immune-cold", 1, 0)
  )

cox_meta <- coxph(
  Surv(OS_time, OS_event) ~ group + age + npi,
  data = surv_meta_cox
)

cat("\n── Cox model — METABRIC (corrected) ──\n")
print(summary(cox_meta))

save(surv_meta, surv_meta_cox, cox_meta,
     immune_cold_score, ifng_score,
     pos_present, neg_present,
     file = file.path(proc_dir, "12_metabric_validation.RData"))

cat("\n METABRIC validation complete!\n")

cat("\n── Cox model — METABRIC ──\n")
print(summary(cox_meta))

# ── 13. SAVE ──────────────────────────────────────────────────────────────────
save(surv_meta, surv_meta_cox, cox_meta,
     immune_cold_score, ifng_score,
     pos_present, neg_present,
     file = file.path(proc_dir, "12_metabric_validation.RData"))

cat("\n METABRIC validation complete!\n")
cat("Figures 31-33 saved to outputs/figures/\n")

# Fix Cox model - use Nottingham prognostic index instead of stage
surv_meta_cox <- surv_meta %>%
  left_join(
    clin_pat %>%
      select(patient_id,
             age = `Age at Diagnosis`,
             npi = `Nottingham prognostic index`),
    by = "patient_id"
  ) %>%
  mutate(
    age   = as.numeric(age),
    npi   = as.numeric(npi),
    group = ifelse(immune_group == "Immune-cold", 1, 0)
  )

cox_meta <- coxph(
  Surv(OS_time, OS_event) ~ group + age + npi,
  data = surv_meta_cox
)

cat("\n── Cox model — METABRIC ──\n")
print(summary(cox_meta))

# Save
save(surv_meta, surv_meta_cox, cox_meta,
     immune_cold_score, ifng_score,
     pos_present, neg_present,
     file = file.path(proc_dir, "12_metabric_validation.RData"))

cat("\n METABRIC validation complete!\n")

# Hazard Forest plot
library(survminer)

p_forest <- ggforest(
  cox_meta,
  data     = surv_meta_cox,
  main     = "Cox model: OS in METABRIC TNBC (n=320)",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.8
)

png(file.path(fig_dir, "37_METABRIC_cox_forest.png"),
    width = 10, height = 5, units = "in", res = 150)
print(p_forest)
dev.off()
print(p_forest)