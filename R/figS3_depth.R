# Supp Fig S3 — Sequencing depth vs repertoire depth (5 panels)
#
# Real input: data/metrics/atlas_diversity.csv (n_clones, total_reads, cohort)
#
# Layout:
#   "AAB
#    CDE"
#
# Output: paper/figures/figS3_depth.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(tidyr)
})

div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)

cohorts <- sort(unique(div$cohort))

# Spearman + Pearson across all 12 samples
spear_all <- cor(div$total_reads, div$n_clones, method = "spearman")
pear_all  <- cor(log10(div$total_reads + 1),
                 log10(div$n_clones  + 1), method = "pearson")

# Per-cohort summary
per_cohort_cor <- div |>
  group_by(cohort) |>
  summarise(
    n = dplyr::n(),
    spearman = if (n() >= 3) round(cor(total_reads, n_clones, method = "spearman"), 3) else NA_real_,
    pearson_log = if (n() >= 3) round(cor(log10(total_reads + 1),
                                          log10(n_clones + 1),
                                          method = "pearson"), 3) else NA_real_,
    slope_log   = if (n() >= 3) round(coef(lm(log10(n_clones + 1) ~
                                              log10(total_reads + 1)))[2], 3) else NA_real_,
    .groups = "drop"
  )

# ----- A. Big scatter: n_clones vs total_reads colored by cohort + global fit
p_A <- div |>
  ggplot(aes(x = total_reads, y = n_clones, color = cohort)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "grey40", fill = "grey80",
              linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = max(div$total_reads), y = min(div$n_clones),
           hjust = 1, vjust = 0,
           label = sprintf("Spearman ρ = %.3f\nPearson r (log10) = %.3f\nn = %d",
                           spear_all, pear_all, nrow(div)),
           size = 2.6, fontface = "italic") +
  scale_color_npg(name = NULL) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  labs(x = "TRUST4 mapped reads (log10)",
       y = "Productive TRB clones (log10)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.18, 0.78),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- B. Per-cohort log-log linear fit --------------------------------------
p_B <- div |>
  ggplot(aes(x = total_reads, y = n_clones, color = cohort, fill = cohort)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, alpha = 0.18) +
  scale_color_npg(guide = "none") +
  scale_fill_npg(guide = "none") +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  facet_wrap(~cohort, scales = "free", ncol = 1) +
  labs(x = "TRUST4 mapped reads",
       y = "Productive TRB clones",
       title = NULL) +
  theme_thoracictcr() +
  theme(strip.text = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 5.5))

# ----- C. Residuals per sample (log-log fit) ---------------------------------
fit <- lm(log10(n_clones + 1) ~ log10(total_reads + 1), data = div)
div$resid <- residuals(fit)

p_C <- div |>
  arrange(cohort, resid) |>
  mutate(sample_short = short_id(sample_id),
         sample_short = factor(sample_short, levels = sample_short)) |>
  ggplot(aes(x = sample_short, y = resid, fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.65) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  scale_fill_npg(guide = "none") +
  labs(x = NULL,
       y = "Residual (log10 n_clones)",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

# ----- D. Spearman + Pearson per cohort (forest-style points) -----------------
cor_long <- per_cohort_cor |>
  pivot_longer(c(spearman, pearson_log),
               names_to = "stat", values_to = "value") |>
  mutate(stat = recode(stat,
                       spearman = "Spearman ρ",
                       pearson_log = "Pearson r (log10)"))

p_D <- cor_long |>
  ggplot(aes(x = value, y = cohort, color = stat)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_point(size = 2.6, alpha = 0.9,
             position = position_dodge(0.4)) +
  scale_color_manual(values = c("Spearman ρ" = npg_pal[1],
                                "Pearson r (log10)" = npg_pal[3]),
                     name = NULL) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Correlation coefficient", y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ----- E. Regression summary table (slope of log10–log10 fit) ----------------
p_E <- per_cohort_cor |>
  mutate(label = sprintf("n=%d\nslope=%.2f", n, slope_log)) |>
  ggplot(aes(x = cohort, y = slope_log, fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.6) +
  geom_text(aes(label = label),
            vjust = -0.2, size = 2.4, lineheight = 0.9,
            fontface = "bold") +
  scale_fill_npg(guide = "none") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL,
       y = "log10–log10 slope",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

design <- "
AAB
CDE
"
figS3 <- p_A + p_B + p_C + p_D + p_E +
  plot_layout(design = design, heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(figS3, "FigS3",
         width_mm = FIG_W_DOUBLE, height_mm = 175)
