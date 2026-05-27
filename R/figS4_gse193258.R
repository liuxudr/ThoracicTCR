# Supp Fig S4 — GSE193258 cohort-exclusion diagnostic (5 data panels)
#
# Real inputs:
#   data/metrics/gse193258_diagnostic_per_assembly.csv
#   data/metrics/gse193258_diagnostic_summary.csv
#
# Layout (text-only "verdict" panel removed — exclusion rationale now lives in
# the figure caption / plot_annotation(caption), not in the panel grid):
#   "ABC
#    DDE"
#   Panel D widened (per-sample V-only and V+J rates are the key evidence).
#
# Output: paper/figures/figS4_gse193258.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(tidyr)
})

asm <- read_csv(file.path(METRICS_DIR,
                          "gse193258_diagnostic_per_assembly.csv"),
                show_col_types = FALSE)
ssm <- read_csv(file.path(METRICS_DIR,
                          "gse193258_diagnostic_summary.csv"),
                show_col_types = FALSE)

asm <- asm |>
  mutate(completeness = factor(completeness,
                               levels = c("V+J", "V-only", "J-only", "None")),
         # Abbreviate SRR to last-5 form so the x-axis label fits.
         label2 = paste0(short_id(sample_id), "\n",
                         cohort, " — ",
                         ifelse(label == "Problematic (excluded)",
                                "Problematic", "Good")))

# ----- A. Stacked-bar completeness (representative samples) -------------------
asm_count <- asm |>
  count(sample_id, label2, completeness, name = "n") |>
  group_by(sample_id) |>
  mutate(frac = n / sum(n)) |>
  ungroup()

p_A <- asm_count |>
  ggplot(aes(x = label2, y = frac, fill = completeness)) +
  geom_col(position = "stack", color = "black",
           linewidth = 0.2, width = 0.7) +
  scale_fill_manual(values = c(`V+J` = npg_pal[3],
                               `V-only` = npg_pal[1],
                               `J-only` = npg_pal[2],
                               `None` = "grey70"),
                    name = "Assembly") +
  scale_y_continuous(labels = label_percent(), expand = c(0, 0),
                     limits = c(0, 1.001)) +
  labs(x = NULL, y = "Fraction of assemblies",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(size = 5.5, lineheight = 1,
                                   angle = 45, hjust = 1),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        legend.position = "right",
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0))

# ----- B. CDR3 nt length overlay (truncated vs full) -------------------------
p_B <- asm |>
  filter(!is.na(cdr3_nt_len), cdr3_nt_len > 0, cdr3_nt_len < 120) |>
  ggplot(aes(x = cdr3_nt_len, color = label, fill = label)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  scale_color_manual(values = c("Good cohort" = npg_pal[3],
                                "Problematic (excluded)" = npg_pal[1]),
                     name = NULL) +
  scale_fill_manual(values = c("Good cohort" = npg_pal[3],
                               "Problematic (excluded)" = npg_pal[1]),
                    guide = "none") +
  labs(x = "CDR3 nucleotide length",
       y = "Density",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.72, 0.78),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- C. Per-sample V+J completeness rate ----------------------------------
p_C <- ssm |>
  mutate(sample_short = short_id(sample_id),
         sample_short = reorder(sample_short, pct_V_plus_J)) |>
  ggplot(aes(x = pct_V_plus_J / 100, y = sample_short, fill = label)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_vline(xintercept = 0.10, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", pct_V_plus_J)),
            hjust = -0.1, size = 2.4) +
  scale_fill_manual(values = c("Good cohort" = npg_pal[3],
                               "Problematic (excluded)" = npg_pal[1]),
                    name = NULL) +
  scale_x_continuous(labels = label_percent(),
                     limits = c(0, 1.18),
                     expand = c(0, 0)) +
  labs(x = "V+J fraction",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "top",
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 5.5),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0))

# ----- D. Per-sample partial-assembly (V-only) rate -------------------------
p_D <- ssm |>
  mutate(pct_V_only = n_V_only / n_assemblies * 100,
         sample_short = short_id(sample_id),
         sample_short = reorder(sample_short, pct_V_only)) |>
  ggplot(aes(x = pct_V_only / 100, y = sample_short, fill = label)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_V_only)),
            hjust = -0.1, size = 2.4) +
  scale_fill_manual(values = c("Good cohort" = npg_pal[3],
                               "Problematic (excluded)" = npg_pal[1]),
                    guide = "none") +
  scale_x_continuous(labels = label_percent(),
                     limits = c(0, 1.18),
                     expand = c(0, 0)) +
  labs(x = "V-only fraction",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0))

# ----- E. Number of assemblies per sample (workload comparison) -------------
p_E <- ssm |>
  mutate(sample_short = short_id(sample_id),
         sample_short = reorder(sample_short, n_assemblies)) |>
  ggplot(aes(x = n_assemblies, y = sample_short, fill = label)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = n_assemblies), hjust = -0.15, size = 2.4) +
  scale_fill_manual(values = c("Good cohort" = npg_pal[3],
                               "Problematic (excluded)" = npg_pal[1]),
                    guide = "none") +
  scale_x_continuous(labels = label_comma(),
                     expand = expansion(mult = c(0, 0.22))) +
  labs(x = "TRUST4 assemblies",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0))

# ----- Verdict numbers (now exposed as plot caption, not as a text-only panel)
n_bad <- ssm |> filter(label == "Problematic (excluded)") |> nrow()
n_gse <- ssm |> filter(cohort == "GSE193258") |> nrow()
med_pct_bad <- median(ssm$pct_V_plus_J[ssm$cohort == "GSE193258"], na.rm = TRUE)
med_pct_good <- median(ssm$pct_V_plus_J[ssm$cohort != "GSE193258"], na.rm = TRUE)

verdict_caption <- sprintf(
  paste("Exclusion verdict: %d / %d GSE193258 samples have <10%% V+J assemblies.",
        "Median V+J = %.1f%% (GSE193258) vs %.1f%% (other cohorts).",
        "V-only assemblies dominate and CDR3 nt lengths fall in a narrow,",
        "short band, consistent with incomplete CDR3 reads.\n",
        "Cohort-level QC fails; cohort removed before atlas integration.",
        sep = " "),
  n_bad, n_gse, med_pct_bad, med_pct_good)

# 5-panel design: drop the text-only "verdict" panel, widen panel D (key
# evidence) to span 2 cells. Each panel now shows data, no free prose.
design <- "
ABC
DDE
"
figS4 <- p_A + p_B + p_C + p_D + p_E +
  plot_layout(design = design) +
  plot_annotation(
    caption = verdict_caption,
    theme = theme(plot.caption = element_text(size = 6, hjust = 0,
                                              lineheight = 1.1,
                                              face = "italic",
                                              color = "grey20"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(figS4, "FigS4",
         width_mm = FIG_W_DOUBLE, height_mm = 170)
