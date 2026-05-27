# Supp Fig S6 — Rarefaction & sampling depth analysis (5 panels)
#
# Real inputs:
#   data/metrics/rarefaction.csv     (sample × fraction × mean ± SD unique clones)
#   data/metrics/atlas_diversity.csv (chao1, n_clones)
#
# Layout:
#   "AAB
#    CDE"
#
# Output: paper/figures/figS6_rarefaction.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(tidyr)
})

rar <- read_csv(file.path(METRICS_DIR, "rarefaction.csv"),
                show_col_types = FALSE)
div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)

# ----- A. Per-sample rarefaction with ±1 SD ribbon ---------------------------
p_A <- rar |>
  ggplot(aes(x = fraction, y = unique_clones,
             group = sample_id, color = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = pmax(unique_clones - unique_clones_sd, 0),
                  ymax = unique_clones + unique_clones_sd),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 1.1, alpha = 0.85) +
  scale_color_npg(name = NULL) +
  scale_fill_npg(guide = "none") +
  scale_x_continuous(labels = label_percent()) +
  scale_y_log10(labels = label_comma()) +
  labs(x = "Subsampling fraction",
       y = "Unique TRB clones (log10)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.18, 0.78),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- B. Predicted plateau: Chao1 vs observed n_clones ----------------------
plateau <- div |>
  mutate(sample_short = short_id(sample_id),
         sample_short = reorder(sample_short, chao1 / n_clones))

p_B <- plateau |>
  ggplot(aes(x = chao1 / n_clones, y = sample_short, fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40",
             linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", chao1 / n_clones)),
            hjust = -0.1, size = 2.3) +
  scale_fill_npg(guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Chao1 / observed",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6))

# ----- C. Per-cohort mean rarefaction (smoothed) -----------------------------
rar_cohort <- rar |>
  group_by(cohort, fraction) |>
  summarise(mean_clones = mean(unique_clones),
            sd_clones   = sd(unique_clones),
            .groups = "drop")

p_C <- rar_cohort |>
  ggplot(aes(x = fraction, y = mean_clones,
             color = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = pmax(mean_clones - sd_clones, 0),
                  ymax = mean_clones + sd_clones),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_color_npg(name = NULL) +
  scale_fill_npg(guide = "none") +
  scale_x_continuous(labels = label_percent()) +
  scale_y_log10(labels = label_comma()) +
  labs(x = "Subsampling fraction",
       y = "Mean unique clones (log10)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.22, 0.78),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- D. Unique clones at 50% subsample per sample -------------------------
rar50 <- rar |>
  filter(abs(fraction - 0.5) < 0.01) |>
  arrange(cohort, unique_clones) |>
  mutate(sample_id = factor(sample_id, levels = sample_id))

p_D <- rar50 |>
  mutate(sample_short = factor(short_id(as.character(sample_id)),
                               levels = short_id(levels(sample_id)))) |>
  ggplot(aes(x = sample_short, y = unique_clones, fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_errorbar(aes(ymin = pmax(unique_clones - unique_clones_sd, 0),
                    ymax = unique_clones + unique_clones_sd),
                width = 0.25, linewidth = 0.3) +
  scale_fill_npg(guide = "none") +
  scale_y_log10(labels = label_comma()) +
  labs(x = NULL,
       y = "Clones at 50% subsample (log10)",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

# ----- E. Saturation index = 1 - (clones@50% / clones@100%) ----------------
sat <- rar |>
  filter(fraction %in% c(0.5, 1.0)) |>
  select(sample_id, cohort, fraction, unique_clones) |>
  pivot_wider(names_from = fraction, values_from = unique_clones,
              names_prefix = "f_") |>
  mutate(sat_index = f_0.5 / f_1) |>
  arrange(cohort, sat_index) |>
  mutate(sample_id = factor(sample_id, levels = sample_id))

p_E <- sat |>
  mutate(sample_short = factor(short_id(as.character(sample_id)),
                               levels = short_id(levels(sample_id)))) |>
  ggplot(aes(x = sat_index, y = sample_short, fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", sat_index)),
            hjust = -0.1, size = 2.3) +
  scale_fill_npg(guide = "none") +
  scale_x_continuous(limits = c(0, 1.18), expand = c(0, 0)) +
  labs(x = "Clones@50% / Clones@100%",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold", hjust = 0))

design <- "
AAB
CDE
"
figS6 <- p_A + p_B + p_C + p_D + p_E +
  plot_layout(design = design, heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(figS6, "FigS6",
         width_mm = FIG_W_DOUBLE, height_mm = 170)
