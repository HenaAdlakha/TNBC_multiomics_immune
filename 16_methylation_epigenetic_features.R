# ── 16_methylation_epigenetic_features.R ──────────────────────────────────────
# Goals:
#   1. Annotate top Factor 4 methylation weights to genes
#   2. Identify immune gene promoter methylation
#   3. Compare methylation at immune loci between groups
#   4. Address epigenetic features part of research question
#
# Input:  01_TNBC_clean.RData, 09_mofa_tnbc3_trained.RData
# Output: figures 48-49, epigenetic feature table
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(ggplot2)
library(MOFA2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

# ── INSTALL ANNOTATION PACKAGE ────────────────────────────────────────────────
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                      quietly = TRUE)) {
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
}

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
load(file.path(proc_dir, "01_TNBC_clean.RData"))
load(file.path(proc_dir, "09_mofa_tnbc3_trained.RData"))

# ── 1. GET METHYLATION WEIGHTS FROM FACTOR 4 ──────────────────────────────────
cat("Extracting Factor 4 methylation weights...\n")

weights_f4_meth <- get_weights(mofa_trained_tnbc3,
                               factors = "Factor4",
                               as.data.frame = TRUE) %>%
  filter(view == "Methylation") %>%
  mutate(probe_id = feature)

cat("Total methylation features in Factor 4:", nrow(weights_f4_meth), "\n")

# ── 2. ANNOTATE CpG PROBES ────────────────────────────────────────────────────
cat("Loading 450k annotation...\n")

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno_df <- as.data.frame(anno) %>%
  rownames_to_column("probe_id") %>%
  dplyr::select(
    probe_id,
    chr,
    pos,
    gene             = UCSC_RefGene_Name,
    gene_group       = UCSC_RefGene_Group,
    relation_island  = Relation_to_Island,
    regulatory       = Regulatory_Feature_Group
  ) %>%
  mutate(
    gene_clean = sapply(strsplit(gene, ";"), function(x) {
      paste(unique(trimws(x)), collapse = ";")
    })
  )

cat("Annotation loaded:", nrow(anno_df), "probes\n")

# Join weights with annotation
weights_annotated <- weights_f4_meth %>%
  left_join(anno_df, by = "probe_id") %>%
  arrange(value)

cat("Probes annotated:", sum(!is.na(weights_annotated$gene)), "\n")

# ── 3. TOP NEGATIVE WEIGHTS — immune-cold associated methylation ───────────────
cat("\n── Top CpGs driving immune-cold (negative weights) ──\n")
top_cold_meth <- weights_annotated %>%
  filter(!is.na(gene), gene != "") %>%
  head(30)
print(top_cold_meth %>%
        dplyr::select(probe_id, gene_clean, value,
                      relation_island, gene_group, chr))

# ── 4. TOP POSITIVE WEIGHTS — immune-hot associated methylation ───────────────
cat("\n── Top CpGs driving immune-hot (positive weights) ──\n")
top_hot_meth <- weights_annotated %>%
  filter(!is.na(gene), gene != "") %>%
  arrange(desc(value)) %>%
  head(30)
print(top_hot_meth %>%
        dplyr::select(probe_id, gene_clean, value,
                      relation_island, gene_group, chr))

# ── 5. CHECK IMMUNE GENES SPECIFICALLY ────────────────────────────────────────
immune_genes_list <- c(
  "CD8A", "CD3D", "CD3E", "GZMB", "PRF1",
  "IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10",
  "HLA-A", "HLA-B", "HLA-C", "B2M",
  "PDCD1", "CD274", "LAG3", "TIGIT", "CTLA4",
  "TAP1", "TAP2", "PSMB8", "PSMB9"
)

cat("\n── Immune genes in Factor 4 methylation weights ──\n")
immune_meth <- weights_annotated %>%
  filter(sapply(gene_clean, function(g) {
    any(sapply(immune_genes_list, function(ig) grepl(ig, g,
                                                     fixed = TRUE)))
  })) %>%
  arrange(value)

if(nrow(immune_meth) > 0) {
  cat("Found", nrow(immune_meth), "immune gene CpG sites\n")
  print(immune_meth %>%
          dplyr::select(probe_id, gene_clean, value,
                        relation_island, gene_group) %>%
          head(20))
} else {
  cat("No immune genes found in top methylation weights\n")
}

# ── 6. PROMOTER CpGs SPECIFICALLY ─────────────────────────────────────────────
cat("\n── Promoter CpGs in Factor 4 weights ──\n")
promoter_meth <- weights_annotated %>%
  filter(grepl("TSS|Promoter|5'UTR", gene_group, ignore.case = TRUE),
         !is.na(gene), gene != "") %>%
  arrange(value)

cat("Promoter CpGs:", nrow(promoter_meth), "\n")
cat("\nTop 20 promoter CpGs (most immune-cold associated):\n")
print(promoter_meth %>%
        dplyr::select(probe_id, gene_clean, value,
                      relation_island, gene_group) %>%
        head(20))

# ── 7. METHYLATION DIFFERENCES AT IMMUNE LOCI ─────────────────────────────────
# Get Factor 4 groups
f4 <- get_factors(mofa_trained_tnbc3,
                  factors = "all")[[1]][, "Factor4"]

common_meth <- intersect(names(f4), colnames(meth_top))

immune_group_meth <- ifelse(
  f4[common_meth] > median(f4[common_meth]),
  "Immune-cold", "Immune-hot"
)

# If immune genes found in weights, test their methylation
if(nrow(immune_meth) > 0) {
  immune_probes <- intersect(immune_meth$probe_id, rownames(meth_top))
  
  if(length(immune_probes) > 0) {
    cat("\n── Methylation at immune loci: cold vs hot ──\n")
    
    meth_immune_df <- as.data.frame(
      t(meth_top[immune_probes, common_meth, drop = FALSE])
    ) %>%
      rownames_to_column("sample") %>%
      mutate(immune_group = immune_group_meth[sample]) %>%
      pivot_longer(cols = all_of(immune_probes),
                   names_to  = "probe",
                   values_to = "beta")
    
    # Add gene names
    meth_immune_df <- meth_immune_df %>%
      left_join(immune_meth %>%
                  dplyr::select(probe_id, gene_clean),
                by = c("probe" = "probe_id"))
    
    # Wilcoxon per probe
    meth_tests <- meth_immune_df %>%
      group_by(probe, gene_clean) %>%
      summarise(
        mean_cold = round(mean(beta[immune_group == "Immune-cold"],
                               na.rm = TRUE), 3),
        mean_hot  = round(mean(beta[immune_group == "Immune-hot"],
                               na.rm = TRUE), 3),
        p_value   = tryCatch(
          round(wilcox.test(
            beta[immune_group == "Immune-cold"],
            beta[immune_group == "Immune-hot"]
          )$p.value, 4),
          error = function(e) NA
        ),
        .groups = "drop"
      ) %>%
      mutate(
        direction   = ifelse(mean_cold > mean_hot,
                             "Hypermethylated in cold",
                             "Hypomethylated in cold"),
        significant = p_value < 0.05
      ) %>%
      arrange(p_value)
    
    print(meth_tests)
    
    # Plot
    p_meth_immune <- ggplot(
      meth_immune_df,
      aes(x = immune_group, y = beta, fill = immune_group)
    ) +
      geom_boxplot(alpha = 0.8, width = 0.5, outlier.size = 1) +
      scale_fill_manual(values = c("Immune-cold" = "#D85A30",
                                   "Immune-hot"  = "#1D9E75")) +
      facet_wrap(~gene_clean, scales = "free_y") +
      labs(
        title    = "Methylation at immune gene loci — TNBC",
        subtitle = "Beta values: higher = more methylated",
        x = NULL, y = "DNA methylation (beta value)",
        fill = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title      = element_text(face = "bold"),
            legend.position = "bottom",
            strip.text      = element_text(face = "bold", size = 9))
    
    print(p_meth_immune)
    ggsave(file.path(fig_dir, "48_immune_loci_methylation.png"),
           p_meth_immune, width = 10, height = 7, dpi = 150)
  }
}

# ── 8. RELATION TO CpG ISLAND ─────────────────────────────────────────────────
cat("\n── CpG island relation for top weights ──\n")
island_summary <- weights_annotated %>%
  filter(!is.na(relation_island)) %>%
  mutate(cold_hot = ifelse(value < 0, "Immune-cold", "Immune-hot")) %>%
  group_by(cold_hot, relation_island) %>%
  summarise(n = n(), mean_weight = round(mean(abs(value)), 4),
            .groups = "drop") %>%
  arrange(cold_hot, desc(n))

print(island_summary)

p_island <- ggplot(island_summary,
                   aes(x = relation_island, y = n, fill = cold_hot)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Immune-cold" = "#D85A30",
                               "Immune-hot"  = "#1D9E75")) +
  labs(
    title = "CpG island context of Factor 4 methylation weights",
    x = "Relation to CpG island",
    y = "Number of CpG sites",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 45, hjust = 1))

print(p_island)
ggsave(file.path(fig_dir, "49_cpg_island_context.png"),
       p_island, width = 8, height = 5, dpi = 150)

# ── 9. SAVE ───────────────────────────────────────────────────────────────────
write.csv(
  weights_annotated %>%
    filter(!is.na(gene), gene != "") %>%
    dplyr::select(probe_id, gene_clean, value,
                  relation_island, gene_group, chr) %>%
    head(100),
  file.path(base_dir, "outputs/factor4_methylation_annotated.csv"),
  row.names = FALSE
)

save(weights_annotated, promoter_meth, immune_meth,
     file = file.path(proc_dir, "16_methylation_features.RData"))

cat("\n Script 16 complete!\n")
cat("Figures 48-49 saved\n")
cat("Annotated methylation weights saved to outputs/\n")