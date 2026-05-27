# Supp Fig S8 — Cohort attrition flow (CONSORT-style, single panel)
#
# Real numbers from configs/atlas_manifest.csv, configs/geo_thoracic_manifest.csv,
# data/metrics/atlas_diversity.csv, data/metrics/gse193258_diagnostic_summary.csv
#
# Output: paper/figures/figS8_attrition.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages(library(tibble))

man  <- read_csv(file.path(REPO_ROOT, "configs", "atlas_manifest.csv"),
                 show_col_types = FALSE)
geo  <- read_csv(file.path(REPO_ROOT, "configs", "geo_thoracic_manifest.csv"),
                 show_col_types = FALSE)
div  <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                 show_col_types = FALSE)
diag <- read_csv(file.path(METRICS_DIR, "gse193258_diagnostic_summary.csv"),
                 show_col_types = FALSE)

# Real counts
n_ena  <- nrow(geo)
n_atlas <- nrow(man)
n_atlas_minus_excluded_cohort <- man |>
  filter(queried_accession != "GSE193258") |> nrow()
n_prod  <- nrow(div)
n_gse193_atlas <- man |> filter(queried_accession == "GSE193258") |> nrow()
n_excluded_query <- n_ena - n_atlas

# Layout
boxes <- tribble(
  ~id, ~x, ~y, ~w,   ~h, ~stage,        ~label,
  1,   2,  10, 3.4, 1.15, "Query",
    sprintf("ENA query: %d RNA-Seq runs\n4 cohorts: GSE126044, GSE135222,\nGSE145370, GSE193258", n_ena),
  2,   2,  8,  3.4, 1.15, "Selected",
    sprintf("Atlas selection: %d runs\n(cohort cap & productivity rank)", n_atlas),
  3,   2,  6,  3.4, 1.15, "Cohort QC",
    sprintf("Cohort-level assembly QC\nGSE193258 excluded\n(%d / %d <10%% V+J)",
            diag |> filter(cohort == "GSE193258", pct_V_plus_J < 10) |> nrow(),
            diag |> filter(cohort == "GSE193258") |> nrow()),
  4,   2,  4,  3.4, 1.15, "Download QC", "Download + gzip QC\nincomplete / corrupt files excluded",
  5,   2,  2,  3.4, 1.15, "Productive",
    sprintf("Final productive atlas: %d samples\n3 retained cohorts", n_prod)
)

stage_colors <- c(Query = npg_pal[1], Selected = npg_pal[2],
                  `Cohort QC` = npg_pal[3], `Download QC` = npg_pal[4],
                  Productive = npg_pal[5])

arrows <- tribble(
  ~x, ~y,   ~xend, ~yend,
  2,  9.42, 2,     8.58,
  2,  7.42, 2,     6.58,
  2,  5.42, 2,     4.58,
  2,  3.42, 2,     2.58
)

excl <- tribble(
  ~x, ~y, ~label,
  6,  9,  sprintf("Excluded:\n%d runs not retained\n(cohort cap & productivity)", n_excluded_query),
  6,  7,  sprintf("Excluded:\n%d GSE193258 atlas runs\n(<6%% V+J assemblies)", n_gse193_atlas),
  6,  5,  sprintf("Excluded:\n%d incomplete download\nor corrupt gzip",
                  n_atlas_minus_excluded_cohort - n_prod),
  6,  3,  "Retained:\nproductive TRB output\nin all remaining samples"
)

ex_arrows <- tribble(
  ~x,  ~y, ~xend, ~yend,
  3.7, 9, 4.5,    9,
  3.7, 7, 4.5,    7,
  3.7, 5, 4.5,    5,
  3.7, 3, 4.5,    3
)

p <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = x - w / 2, xmax = x + w / 2,
                ymin = y - h / 2, ymax = y + h / 2,
                fill = stage),
            color = "black", linewidth = 0.4, alpha = 0.85) +
  geom_text(data = boxes, aes(x = x, y = y, label = label),
            size = 2.5, lineheight = 1) +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
               linewidth = 0.55, color = "black") +
  geom_segment(data = ex_arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
               linewidth = 0.4, color = "grey45") +
  geom_text(data = excl, aes(x = x, y = y, label = label),
            size = 2.3, lineheight = 1, hjust = 0,
            color = "grey20") +
  scale_fill_manual(values = stage_colors, name = "Stage") +
  scale_x_continuous(limits = c(0, 9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1, 11), expand = c(0, 0)) +
  labs(title = NULL,
       x = NULL, y = NULL) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

save_fig(p, "FigS8",
         width_mm = FIG_W_DOUBLE, height_mm = 165)
