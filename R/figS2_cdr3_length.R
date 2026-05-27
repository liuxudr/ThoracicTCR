# Supp Fig S2 — CDR3-aa length distribution by cohort (4 panels)
#
# Real input: data/metrics/clone_size_distribution.csv
#
# Layout:
#   "ABC
#    DDD"
#
# Output: paper/figures/figS2_cdr3_length.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
})

dist <- read_csv(file.path(METRICS_DIR, "clone_size_distribution.csv"),
                 show_col_types = FALSE)

dist <- dist |> filter(!is.na(cdr3_len), cdr3_len >= 5, cdr3_len <= 28)
cohorts <- sort(unique(dist$cohort))

palette_map <- setNames(pal_npg("nrc")(length(cohorts)), cohorts)

make_cohort_panel <- function(coh, letter) {
  d <- dist |> filter(cohort == coh)
  ggplot(d, aes(x = cdr3_len, weight = count)) +
    geom_histogram(binwidth = 1, fill = palette_map[[coh]],
                   color = "black", linewidth = 0.2) +
    scale_x_continuous(breaks = seq(6, 28, 2)) +
    scale_y_continuous(labels = label_comma()) +
    labs(x = "CDR3 aa length",
         y = "Read count",
         title = sprintf("%s. %s (n=%d clones)",
                         letter, coh, length(unique(d$cdr3_aa)))) +
    theme_thoracictcr()
}

letters3 <- LETTERS[seq_along(cohorts)]
panels <- mapply(make_cohort_panel, cohorts, letters3, SIMPLIFY = FALSE)

# Panel D — cross-cohort overlay (density, read-weighted) -------------------
p_D <- ggplot(dist,
              aes(x = cdr3_len, color = cohort, fill = cohort,
                  weight = count)) +
  geom_density(alpha = 0.25, linewidth = 0.55) +
  scale_color_npg(name = NULL) +
  scale_fill_npg(guide = "none") +
  scale_x_continuous(breaks = seq(6, 28, 2)) +
  labs(x = "CDR3 aa length",
       y = "Density (read-weighted)",
       title = sprintf("%s. Cross-cohort overlay",
                       LETTERS[length(cohorts) + 1])) +
  theme_thoracictcr() +
  theme(legend.position = c(0.82, 0.78),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# Build with patchwork — design depending on cohort count -------------------
design <- if (length(cohorts) == 3) "ABC\nDDD" else "AB\nCD"
plots <- c(panels, list(p_D))
names(plots) <- c(letters3, LETTERS[length(cohorts) + 1])

figS2 <- Reduce(`+`, plots) +
  plot_layout(design = design, heights = c(1, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(figS2, "FigS2",
         width_mm = FIG_W_DOUBLE, height_mm = 135)
