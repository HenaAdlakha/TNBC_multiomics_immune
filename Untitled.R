library(tidyverse)
library(ggplot2)

base_dir <- "~/Documents/R Working Directory/BRCA_project"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "outputs/figures")

load(file.path(proc_dir, "22_final_frozen.RData"))

# ── 1. FOREST PLOT: ALL COX RESULTS ───────────────────────────────────────────

# Extract results from each model
extract_score_result <- function(cox_model, label, n, events) {
  coef <- summary(cox_model)$coefficients["score_std", ]
  data.frame(
    model   = label,
    n       = n,
    events  = events,
    hr      = round(exp(coef["coef"]), 3),
    lower   = round(exp(coef["coef"] - 1.96*coef["se(coef)"]), 3),
    upper   = round(exp(coef["coef"] + 1.96*coef["se(coef)"]), 3),
    p       = round(coef["Pr(>|z|)"], 4)
  )
}

# Build results dataframe
forest_df <- bind_rows(
  extract_score_result(
    cox_final,
    "Full model\n(score + age + NPI + PAM50)",
    nrow(surv_final),
    sum(surv_final$OS_event, na.rm=TRUE)
  ),
  extract_score_result(
    cox_basal,
    "Within Basal-like",
    nrow(surv_final %>% filter(pam50_reduced == "Basal")),
    sum(surv_final$OS_event[surv_final$pam50_reduced == "Basal"],
        na.rm=TRUE)
  ),
  extract_score_result(
    cox_claudin,
    "Within Claudin-low",
    nrow(surv_final %>% filter(pam50_reduced == "Claudin-low")),
    sum(surv_final$OS_event[surv_final$pam50_reduced == "Claudin-low"],
        na.rm=TRUE)
  )
) %>%
  mutate(
    label     = paste0(model, "\nn=", n, ", events=", events),
    sig       = p < 0.05,
    p_label   = ifelse(p < 0.001, "p<0.001",
                       ifelse(p < 0.01,  paste0("p=", sprintf("%.3f", p)),
                              paste0("p=", sprintf("%.3f", p)))),
    ci_label  = paste0("HR=", hr, " (", lower, "–", upper, ")")
  )

print(forest_df %>% dplyr::select(model, n, events, hr, lower, upper, p))

# Forest plot
p_forest <- ggplot(forest_df,
                   aes(x = hr, y = reorder(label, hr),
                       color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "#888780", linewidth = 0.8) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                width = 0.2, linewidth = 1,
                orientation = "y") +
  geom_point(size = 4) +
  scale_color_manual(values = c("TRUE"  = "#D85A30",
                                "FALSE" = "#888780"),
                     labels = c("TRUE"  = "p < 0.05",
                                "FALSE" = "p >= 0.05")) +
  geom_text(aes(label = paste0(ci_label, "\n", p_label)),
            x         = max(forest_df$upper) + 0.05,
            hjust     = 0,
            size      = 3.5,
            color     = "#444441") +
  scale_x_continuous(limits = c(0.5,
                                max(forest_df$upper) + 0.7)) +
  labs(
    title    = "Immune-cold score: Cox regression — METABRIC TNBC",
    subtitle = "Hazard ratio per 1 SD increase in immune-cold score",
    x        = "Hazard ratio (95% CI)",
    y        = NULL,
    color    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "bottom",
    axis.text.y      = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

print(p_forest)
ggsave(file.path(fig_dir, "58_forest_plot_final.png"),
       p_forest, width = 11, height = 6, dpi = 150)