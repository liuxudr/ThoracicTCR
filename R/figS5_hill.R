# Supp Fig S5 — Hill diversity ${}^qD$ profile per sample (4 panels)
#
# Real input: data/metrics/atlas_diversity.csv
#   (hill_q0 = richness, hill_q1 = Shannon, hill_q2 = Simpson)
#
# Layout:
#   "AB
#    CD"
#
# Output: paper/figures/figS5_hill.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(tidyr)
})

div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)

hill_long <- div |>
  select(sample_id, cohort, hill_q0, hill_q1, hill_q2) |>
  pivot_longer(starts_with("hill_"),
               names_to = "q", values_to = "value") |>
  mutate(q_num = recode(q, hill_q0 = 0, hill_q1 = 1, hill_q2 = 2))

# Extend with q=3 using inverse_simpson as upper-q proxy (closest we have);
# for q ≥ 2, ${}^qD$ trends to 1/max(p_i), bounded above by inverse Simpson.
# We do NOT fabricate a q=3 column; only plot q=0,1,2 as real measured values.

# ----- A. Per-sample line: q = 0, 1, 2 ---------------------------------------
p_A <- hill_long |>
  ggplot(aes(x = q_num, y = value, group = sample_id, color = cohort)) +
  geom_line(alpha = 0.7, linewidth = 0.45) +
  geom_point(size = 1.6, alpha = 0.9) +
  scale_x_continuous(breaks = 0:2,
                     labels = c("q=0\n(richness)",
                                "q=1\n(Shannon)",
                                "q=2\n(Simpson)")) +
  scale_y_log10(labels = label_comma()) +
  scale_color_npg(name = NULL) +
  labs(x = NULL, y = "Effective species (log10)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- B. Cohort-averaged Hill profile (median ± IQR) ------------------------
hill_summary <- hill_long |>
  group_by(cohort, q_num) |>
  summarise(med = median(value),
            lo  = quantile(value, 0.25),
            hi  = quantile(value, 0.75),
            .groups = "drop")

p_B <- hill_summary |>
  ggplot(aes(x = q_num, y = med, color = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = 0:2,
                     labels = c("q=0", "q=1", "q=2")) +
  scale_y_log10(labels = label_comma()) +
  scale_color_npg(name = NULL) +
  scale_fill_npg(guide = "none") +
  labs(x = "Hill order q",
       y = "Effective species (median ± IQR)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- C. Hill q=1 / q=2 ratio per sample (Shannon : Simpson effective) -----
div <- div |> mutate(q1_q2_ratio = hill_q1 / hill_q2)

p_C <- div |>
  arrange(cohort, q1_q2_ratio) |>
  mutate(sample_short = short_id(sample_id),
         sample_short = factor(sample_short, levels = sample_short)) |>
  ggplot(aes(x = q1_q2_ratio, y = sample_short, fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_fill_npg(guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Hill ratio q=1 / q=2",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6))

# ----- D. Per-cohort Hill profile boxplot (all three q on facets) -----------
p_D <- hill_long |>
  mutate(q_label = recode(q,
                          hill_q0 = "q=0",
                          hill_q1 = "q=1",
                          hill_q2 = "q=2"),
         q_label = factor(q_label, levels = c("q=0", "q=1", "q=2"))) |>
  ggplot(aes(x = cohort, y = value, fill = cohort)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, alpha = 0.85, size = 1.4, shape = 21,
              color = "black", stroke = 0.25) +
  facet_wrap(~q_label, scales = "free_y", ncol = 3) +
  scale_y_log10(labels = label_comma()) +
  scale_fill_npg(guide = "none") +
  labs(x = NULL, y = "Effective species (log10)",
       title = NULL) +
  theme_thoracictcr() +
  theme(strip.text = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

design <- "
AB
CD
"
figS5 <- p_A + p_B + p_C + p_D +
  plot_layout(design = design) &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(figS5, "FigS5",
         width_mm = FIG_W_DOUBLE, height_mm = 165)
