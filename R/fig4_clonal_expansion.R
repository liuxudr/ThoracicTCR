# Figure 4 — Clonal expansion landscape (6 panels)
#
# Real inputs:
#   data/metrics/clone_size_distribution.csv  (full rank list per sample)
#   data/metrics/clone_size_summary.csv       (per-sample summary)
#   data/metrics/atlas_diversity.csv          (d50, n_clones)
#
# Layout:
#   "AAB
#    CDE
#    FFF"
#
# Output: paper/figures/fig4_clonal_expansion.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(tidyr)
})

dist <- read_csv(file.path(METRICS_DIR, "clone_size_distribution.csv"),
                 show_col_types = FALSE)
csum <- read_csv(file.path(METRICS_DIR, "clone_size_summary.csv"),
                 show_col_types = FALSE)
div  <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                 show_col_types = FALSE)

# Sample order: by cohort then n_clones; short_id used on axes to avoid overlap
sample_levels <- csum |> arrange(cohort, desc(n_clones)) |> pull(sample_id)
sample_levels_short <- short_id(sample_levels)

# ----- A. Rank-frequency log-log per sample ----------------------------------
p_A <- dist |>
  mutate(sample_id = factor(sample_id, levels = sample_levels)) |>
  ggplot(aes(x = rank, y = freq, color = cohort, group = sample_id)) +
  geom_line(alpha = 0.85, linewidth = 0.45) +
  scale_x_log10() +
  scale_y_log10(labels = label_percent()) +
  scale_color_npg(name = NULL) +
  labs(x = "Clone rank (log10)",
       y = "Clone frequency (log10)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.18, 0.22),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- B. Top-1 / Top-10 / Top-100 grouped bar per sample --------------------
top_long <- csum |>
  select(sample_id, cohort, top1_frac, top10_frac, top100_frac) |>
  pivot_longer(starts_with("top"), names_to = "metric", values_to = "frac") |>
  mutate(
    metric = recode(metric,
                    top1_frac = "Top-1",
                    top10_frac = "Top-10",
                    top100_frac = "Top-100"),
    metric = factor(metric, levels = c("Top-1", "Top-10", "Top-100")),
    sample_short = factor(short_id(sample_id), levels = sample_levels_short)
  )

p_B <- top_long |>
  ggplot(aes(x = sample_short, y = frac, fill = metric)) +
  geom_col(position = position_dodge(0.8), width = 0.75,
           color = "black", linewidth = 0.15) +
  scale_fill_manual(values = c(`Top-1` = npg_pal[1],
                               `Top-10` = npg_pal[3],
                               `Top-100` = npg_pal[5]),
                    name = NULL) +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1.02),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Cumulative read fraction",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5.5),
        legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ----- C. CDR3-aa length density per cohort ---------------------------------
p_C <- dist |>
  filter(!is.na(cdr3_len), cdr3_len >= 5, cdr3_len <= 28) |>
  ggplot(aes(x = cdr3_len, color = cohort, fill = cohort, weight = count)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  scale_color_npg(name = NULL) +
  scale_fill_npg(guide = "none") +
  labs(x = "CDR3 aa length",
       y = "Density (read-weighted)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.78, 0.78),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- D. Hyperexpanded count per sample (bubble) ----------------------------
p_D <- csum |>
  mutate(sample_short = factor(short_id(sample_id),
                               levels = sample_levels_short)) |>
  ggplot(aes(x = sample_short, y = n_hyperexpanded,
             color = cohort, size = n_clones)) +
  geom_point(alpha = 0.85) +
  scale_color_npg(name = NULL) +
  scale_size_continuous(range = c(2, 6), name = "n_clones",
                        breaks = c(10, 100, 500)) +
  labs(x = NULL, y = "Hyperexpanded clones (>1%)",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5.5),
        legend.position = "right",
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))

# ----- E. D50 boxplot per cohort ---------------------------------------------
p_E <- div |>
  ggplot(aes(x = cohort, y = d50, fill = cohort)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.85, size = 1.6, shape = 21,
              color = "black", stroke = 0.25) +
  scale_fill_npg(guide = "none") +
  labs(x = NULL, y = "D50", title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

# ----- F. Stacked bar: clones above expansion thresholds ---------------------
# For each sample, count clones with freq > 0.01, > 0.05, > 0.10 (mutually
# exclusive bins for a stacked bar)
thr_long <- dist |>
  mutate(bin = cut(freq,
                   breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
                   labels = c("< 1%", "1–5%", "5–10%", "> 10%"),
                   right = FALSE)) |>
  count(sample_id, cohort, bin, name = "n") |>
  group_by(sample_id) |>
  mutate(frac = n / sum(n)) |>
  ungroup() |>
  mutate(sample_short = factor(short_id(sample_id),
                               levels = sample_levels_short),
         bin = factor(bin, levels = c("< 1%", "1–5%", "5–10%", "> 10%")))

p_F <- thr_long |>
  ggplot(aes(x = sample_short, y = frac, fill = bin)) +
  geom_col(position = "stack", color = "black", linewidth = 0.15, width = 0.75) +
  scale_fill_manual(values = c(`< 1%` = "grey80",
                               `1–5%` = npg_pal[4],
                               `5–10%` = npg_pal[2],
                               `> 10%` = npg_pal[1]),
                    name = "Clone freq") +
  scale_y_continuous(labels = label_percent(), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of clones",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

design <- "
AAB
CDE
FFF
"
fig4 <- p_A + p_B + p_C + p_D + p_E + p_F +
  plot_layout(design = design, heights = c(1, 1, 0.85)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(fig4, "Fig4",
         width_mm = FIG_W_DOUBLE, height_mm = 200)
