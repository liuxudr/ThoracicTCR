# Figure 1 — Pan-thoracic TCR atlas landscape (7 panels, narrative layout)
#
# Real inputs:
#   data/metrics/atlas_diversity.csv  (12 productive samples, 17 metric cols)
#
# Layout (patchwork design):
#   "AAB
#    CDE
#    FGG"
#
# Output: paper/figures/fig1_landscape.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ggrepel)
  library(scales)
})

div_csv <- file.path(METRICS_DIR, "atlas_diversity.csv")
stopifnot(file.exists(div_csv))
div <- read_csv(div_csv, show_col_types = FALSE)

# ----- A. Cohort composition --------------------------------------------------
p_A <- div |>
  count(cohort, name = "n_samples") |>
  ggplot(aes(x = reorder(cohort, -n_samples), y = n_samples, fill = cohort)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.25) +
  geom_text(aes(label = n_samples), vjust = -0.3, size = 2.6, fontface = "bold") +
  scale_fill_npg(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Samples", title = NULL) +
  theme_thoracictcr()

# ----- B. Repertoire depth (n_clones, log10) ---------------------------------
p_B <- div |>
  ggplot(aes(x = cohort, y = n_clones, fill = cohort)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.85, size = 1.6, shape = 21,
              color = "black", stroke = 0.25) +
  scale_y_log10(labels = label_comma()) +
  scale_fill_npg(guide = "none") +
  labs(x = NULL, y = "Productive TRB clones",
       title = NULL) +
  theme_thoracictcr()

# ----- C. Shannon diversity (violin) ------------------------------------------
p_C <- div |>
  ggplot(aes(x = cohort, y = shannon, fill = cohort)) +
  geom_violin(alpha = 0.45, scale = "width", width = 0.85,
              color = "black", linewidth = 0.25) +
  geom_jitter(width = 0.12, alpha = 0.85, size = 1.4, shape = 21,
              color = "black", stroke = 0.25) +
  scale_fill_npg(guide = "none") +
  labs(x = NULL, y = "Shannon entropy",
       title = NULL) +
  theme_thoracictcr()

# ----- D. Clonality vs depth scatter -----------------------------------------
p_D <- div |>
  ggplot(aes(x = n_clones, y = clonality, color = cohort)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_npg(name = NULL) +
  scale_x_log10(labels = label_comma()) +
  labs(x = "Productive TRB clones",
       y = "Clonality (1 - H/log S)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.78, 0.82),
        legend.background = element_rect(fill = alpha("white", 0.6),
                                         color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- E. Hill q profile (q=0,1,2) per cohort --------------------------------
hill_long <- div |>
  select(sample_id, cohort, hill_q0, hill_q1, hill_q2) |>
  pivot_longer(starts_with("hill_"), names_to = "q", values_to = "value") |>
  mutate(q_num = recode(q, hill_q0 = 0, hill_q1 = 1, hill_q2 = 2))

p_E <- hill_long |>
  ggplot(aes(x = q_num, y = value, group = sample_id, color = cohort)) +
  geom_line(alpha = 0.45, linewidth = 0.4) +
  stat_summary(aes(group = cohort), fun = median, geom = "line",
               linewidth = 1.1) +
  stat_summary(aes(group = cohort), fun = median, geom = "point",
               size = 2.2, shape = 18) +
  scale_y_log10(labels = label_comma()) +
  scale_color_npg(guide = "none") +
  scale_x_continuous(breaks = 0:2,
                     labels = c("q=0", "q=1", "q=2")) +
  labs(x = "Hill order",
       y = "Effective species",
       title = NULL) +
  theme_thoracictcr()

# ----- F. D50 distribution per cohort ----------------------------------------
p_F <- div |>
  ggplot(aes(x = cohort, y = d50, fill = cohort)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.85, size = 1.6, shape = 21,
              color = "black", stroke = 0.25) +
  scale_fill_npg(guide = "none") +
  labs(x = NULL, y = "D50 clones",
       title = NULL) +
  theme_thoracictcr()

# ----- G. Gini coefficient violin (oligoclonality summary) -------------------
p_G <- div |>
  ggplot(aes(x = cohort, y = gini, fill = cohort)) +
  geom_violin(alpha = 0.45, scale = "width", width = 0.85,
              color = "black", linewidth = 0.25) +
  geom_jitter(width = 0.1, alpha = 0.85, size = 1.4, shape = 21,
              color = "black", stroke = 0.25) +
  scale_fill_npg(guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Gini coefficient",
       title = NULL) +
  theme_thoracictcr()

# ----- Compose with design layout --------------------------------------------
design <- "
AAABB
CCDDE
FFGGG
"
fig1 <- p_A + p_B + p_C + p_D + p_E + p_F + p_G +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(fig1, "Fig1", width_mm = FIG_W_DOUBLE, height_mm = 175)
