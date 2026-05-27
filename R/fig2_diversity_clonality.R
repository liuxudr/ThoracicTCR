# Figure 2 — Diversity & clonal architecture (5 panels, 2 over 3)
#
# Real input: data/metrics/atlas_diversity.csv
#
# Layout:
#   "ABB
#    CDE"
#
# Output: paper/figures/fig2_diversity_clonality.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ggridges)
  library(scales)
})

div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)

# ----- A. Shannon ridge per cohort -------------------------------------------
p_A <- div |>
  ggplot(aes(x = shannon, y = cohort, fill = cohort)) +
  geom_density_ridges(alpha = 0.6, scale = 0.95,
                      jittered_points = TRUE,
                      point_size = 1.2,
                      point_alpha = 0.85,
                      position = position_points_jitter(width = 0, height = 0.05),
                      linewidth = 0.3) +
  scale_fill_npg(guide = "none") +
  labs(x = "Shannon entropy", y = NULL,
       title = NULL) +
  theme_thoracictcr()

# ----- B. Inverse Simpson per cohort -----------------------------------------
p_B <- div |>
  ggplot(aes(x = cohort, y = inverse_simpson, fill = cohort)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.85, size = 1.6, shape = 21,
              color = "black", stroke = 0.25) +
  scale_y_log10(labels = label_comma()) +
  scale_fill_npg(guide = "none") +
  labs(x = NULL, y = "Inverse Simpson (1/D)",
       title = NULL) +
  theme_thoracictcr()

# ----- C. Top-10 fraction vs D50 (oligoclonality) ----------------------------
p_C <- div |>
  ggplot(aes(x = d50, y = top10_frac, color = cohort)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_npg(name = NULL) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(1, 3, 5, 10, 20)) +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1)) +
  labs(x = "D50 (clones to 50% reads, log1p)",
       y = "Top-10 clone fraction",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- D. Pielou evenness density --------------------------------------------
p_D <- div |>
  ggplot(aes(x = pielou_evenness, fill = cohort)) +
  geom_density(alpha = 0.45, linewidth = 0.3, color = "black") +
  scale_fill_npg(name = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Pielou evenness J", y = "Density",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.22, 0.6),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- E. Chao1 vs n_clones (richness saturation) ----------------------------
p_E <- div |>
  ggplot(aes(x = n_clones, y = chao1, color = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.3) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_npg(guide = "none") +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  labs(x = "Observed clones (S)",
       y = "Chao1 richness",
       title = NULL) +
  theme_thoracictcr()

design <- "
ABB
CDE
"
fig2 <- p_A + p_B + p_C + p_D + p_E +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(fig2, "Fig2",
         width_mm = FIG_W_DOUBLE, height_mm = 135)
