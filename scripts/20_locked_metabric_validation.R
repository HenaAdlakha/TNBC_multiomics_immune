library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

fm <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/metabric", fn))

load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))
load(file.path(proc_dir, "12_metabric_validation.RData"))

# ── LOCKED SCORE: SIGNED TCGA WEIGHTS ─────────────────────────────────────────
# Use ALL Factor 4 RNA weights — not selected top 150
# Score = sum(weight_g * standardised_expression_g)
cat("Extracting all Factor 4 RNA weights from TCGA...\n")

weights_f4_rna <- get_weights(mofa_trained_tnbc3,
                              factors = "Factor4",
                              as.data.frame = TRUE) %>%
  filter(view == "RNA") %>%
  mutate(gene = gsub("_RNA", "", feature)) %>%
  dplyr::select(gene, weight = value)

cat("Total RNA weights:", nrow(weights_f4_rna), "\n")

# Load METABRIC expression
expr_raw <- fread(file = fm("data_mrna_illumina_microarray.txt"),
                  header = TRUE)

expr_mat <- expr_raw %>%
  as.data.frame() %>%
  filter(!duplicated(Hugo_Symbol)) %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

# Filter to TNBC samples
tnbc_samples_meta <- surv_meta$sample_id
expr_tnbc <- expr_mat[, intersect(tnbc_samples_meta, colnames(expr_mat))]

# Genes available in both TCGA weights and METABRIC
common_genes <- intersect(weights_f4_rna$gene, rownames(expr_tnbc))
cat("Genes available for scoring:", length(common_genes), "/",
    nrow(weights_f4_rna), "\n")

weights_sub <- weights_f4_rna %>% filter(gene %in% common_genes)

# Standardise METABRIC expression per gene
expr_std <- t(scale(t(expr_tnbc[common_genes, ])))

# Calculate signed weighted score
score_locked <- as.vector(
  t(weights_sub$weight) %*% expr_std[common_genes, ]
)
names(score_locked) <- colnames(expr_tnbc)

cat("Score distribution:\n")
print(summary(score_locked))

# ── VERIFY DIRECTION ──────────────────────────────────────────────────────────
f <- function(fn) path.expand(
  file.path("~/Documents/R Working Directory/BRCA_project/data/raw", fn))

immune_sigs <- fread(
  file = f("TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena"),
  header = TRUE)

immune_mat <- immune_sigs %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  as.matrix()

immune_t <- t(immune_mat)
rownames(immune_t) <- substr(rownames(immune_t), 1, 15)

common_ver <- intersect(names(score_locked), rownames(immune_t))
tcell_meta <- immune_t[common_ver, "Tcell_21978456"]

r_tcell <- cor(score_locked[common_ver], tcell_meta,
               use = "complete.obs")
cat("\nScore vs T cell in METABRIC: r =", round(r_tcell, 3), "\n")
cat("Negative = immune-cold score correctly oriented\n")

# Flip if needed
if(r_tcell > 0) {
  cat("Score is positively correlated — flipping\n")
  score_locked <- -score_locked
}

# ── BUILD SURVIVAL DATAFRAME ───────────────────────────────────────────────────
clin_pat <- fread(file = fm("data_clinical_patient.txt"), header = TRUE)
colnames(clin_pat)[1] <- "patient_id"

clin_sam <- fread(file = fm("data_clinical_sample.txt"), header = TRUE)
colnames(clin_sam)[1] <- "patient_id"

sample_to_patient <- clin_sam %>%
  dplyr::select(patient_id, sample_id = `Sample Identifier`)

surv_locked <- data.frame(
  sample_id    = names(score_locked),
  locked_score = score_locked
) %>%
  left_join(sample_to_patient, by = "sample_id") %>%
  left_join(clin_pat %>%
              dplyr::select(
                patient_id,
                OS_months  = `Overall Survival (Months)`,
                OS_status  = `Overall Survival Status`,
                RFS_months = `Relapse Free Status (Months)`,
                RFS_status = `Relapse Free Status`,
                age        = `Age at Diagnosis`,
                npi        = `Nottingham prognostic index`
              ),
            by = "patient_id") %>%
  mutate(
    OS_time   = as.numeric(OS_months),
    OS_event  = ifelse(OS_status == "1:DECEASED", 1, 0),
    RFS_time  = as.numeric(RFS_months),
    RFS_event = ifelse(RFS_status == "1:Recurred", 1, 0),
    age       = as.numeric(age),
    npi       = as.numeric(npi),
    score_std = as.numeric(scale(locked_score))
  ) %>%
  filter(!is.na(OS_time), !is.na(OS_event))

cat("\nSamples:", nrow(surv_locked), "\n")
cat("Events:", sum(surv_locked$OS_event, na.rm=TRUE), "\n")

# ── PRIMARY COX: CONTINUOUS SCORE ─────────────────────────────────────────────
cat("\n── Primary Cox model: continuous locked score ──\n")
cox_continuous <- coxph(
  Surv(OS_time, OS_event) ~ score_std + age + npi,
  data = surv_locked
)
print(summary(cox_continuous))

# Test proportional hazards assumption
ph_test <- cox.zph(cox_continuous)
cat("\nProportional hazards test:\n")
print(ph_test)

# ── SECONDARY: KM CURVE (MEDIAN SPLIT FOR DISPLAY ONLY) ──────────────────────
surv_locked <- surv_locked %>%
  mutate(immune_group = ifelse(
    locked_score > median(locked_score, na.rm=TRUE),
    "Immune-cold", "Immune-hot"
  ))

cat("\nGroup sizes:\n")
print(table(surv_locked$immune_group))

km_locked <- survfit(Surv(OS_time, OS_event) ~ immune_group,
                     data = surv_locked)

p_km_locked <- ggsurvplot(
  km_locked,
  data              = surv_locked,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = c("#D85A30", "#1D9E75"),
  title    = "OS: locked continuous score — METABRIC TNBC",
  subtitle = "Signed weighted TCGA Factor 4 projection (display only — primary analysis is continuous Cox)",
  xlab     = "Time (months)",
  ylab     = "Survival probability",
  legend.labs = c("Immune-cold", "Immune-hot"),
  ggtheme  = theme_minimal(base_size = 13)
)

png(file.path(fig_dir, "54_locked_METABRIC_KM.png"),
    width = 9, height = 7, units = "in", res = 150)
print(p_km_locked)
dev.off()
print(p_km_locked)

# RFS
cox_rfs <- coxph(
  Surv(RFS_time, RFS_event) ~ score_std + age + npi,
  data = surv_locked %>% filter(!is.na(RFS_time), !is.na(RFS_event))
)
cat("\n── Cox RFS: continuous locked score ──\n")
print(summary(cox_rfs))

save(surv_locked, cox_continuous, cox_rfs, score_locked,
     file = file.path(proc_dir, "20_locked_validation.RData"))

cat("\n Locked METABRIC validation complete!\n")
cat("Primary result: HR per SD =",
    round(exp(coef(cox_continuous)["score_std"]), 3),
    "p =", round(summary(cox_continuous)$coefficients["score_std", 5], 4),
    "\n")