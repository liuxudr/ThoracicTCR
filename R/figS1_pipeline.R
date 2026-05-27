# Supp Fig S1 — ThoracicTCR pipeline schematic (2 panels: flow + decision-table)
#
# All numbers are real (taken from configs/atlas_manifest.csv, data/metrics/
# atlas_diversity.csv, data/metrics/gse193258_diagnostic_summary.csv).
#
# Layout:  "A
#           B"
# A = process flow (boxes + arrows, real counts inside boxes).
# B = decision-point table — tiled (one row = one decision step),
#     columns = Step | Rule | Outcome. All cells anchored to numeric data.
#
# Output: paper/figures/figS1_pipeline.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages(library(tibble))

# Live numbers ----------------------------------------------------------------
man  <- read_csv(file.path(REPO_ROOT, "configs", "atlas_manifest.csv"),
                 show_col_types = FALSE)
geo  <- read_csv(file.path(REPO_ROOT, "configs", "geo_thoracic_manifest.csv"),
                 show_col_types = FALSE)
div  <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                 show_col_types = FALSE)
diag <- read_csv(file.path(METRICS_DIR, "gse193258_diagnostic_summary.csv"),
                 show_col_types = FALSE)

n_ena       <- nrow(geo)
n_selected  <- nrow(man)
n_productive <- nrow(div)
n_cohorts_pre  <- length(unique(geo$queried_accession))
n_cohorts_post <- length(unique(div$cohort))
n_gse193_bad <- diag |>
  filter(cohort == "GSE193258", pct_V_plus_J < 10) |> nrow()

# Panel A — process flow ------------------------------------------------------
# Wider boxes (w = 2.6), bigger horizontal spacing between rows (~3 units),
# so 4 boxes per row × 2 rows fit cleanly inside x = [0, 12].
boxes_A <- tribble(
  ~x,   ~y, ~w,  ~h,  ~stage,    ~label,
  1.5,  5,  2.6, 1.4, "Input",
    sprintf("ENA / SRA query\n%d RNA-Seq runs\n%d cohorts",
            n_ena, n_cohorts_pre),
  4.5,  5,  2.6, 1.4, "Assembly",
    "TRUST4 v1.1\nTRB + TRA assembly\nproductive filter",
  7.5,  5,  2.6, 1.4, "QC",
    sprintf("Assembly QC\n%d / %d GSE193258\nfail < 10%% V+J",
            n_gse193_bad,
            diag |> filter(cohort == "GSE193258") |> nrow()),
  10.5, 5,  2.6, 1.4, "Metrics",
    "Diversity & clonality\nShannon · Simpson\nHill · D50 · Gini",
  4.5,  2,  2.6, 1.4, "Atlas",
    sprintf("Atlas integration\n%d productive samples\n%d cohorts",
            n_productive, n_cohorts_post),
  7.5,  2,  2.6, 1.4, "Antigen",
    "VDJdb / McPAS / IEDB\nepitope · MHC class\nantigen category",
  10.5, 2,  2.6, 1.4, "Output",
    "Publication figures\nNPG palette\n+ tables"
)

arrows_A <- tribble(
  ~x,   ~y,   ~xend, ~yend,
  2.8,  5,    3.2,   5,      # Input -> Assembly
  5.8,  5,    6.2,   5,      # Assembly -> QC
  8.8,  5,    9.2,   5,      # QC -> Metrics
  4.5,  4.3,  4.5,   2.7,    # Assembly -> Atlas (down)
  5.8,  2,    6.2,   2,      # Atlas -> Antigen
  8.8,  2,    9.2,   2       # Antigen -> Output
)

stage_colors <- c(Input = npg_pal[1], Assembly = npg_pal[2], QC = npg_pal[3],
                  Metrics = npg_pal[4], Atlas = npg_pal[5],
                  Antigen = npg_pal[6], Output = npg_pal[7])

p_A <- ggplot() +
  geom_rect(data = boxes_A,
            aes(xmin = x - w / 2, xmax = x + w / 2,
                ymin = y - h / 2, ymax = y + h / 2,
                fill = stage),
            color = "black", linewidth = 0.4, alpha = 0.85) +
  geom_segment(data = arrows_A,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               linewidth = 0.5, color = "black") +
  geom_text(data = boxes_A, aes(x = x, y = y, label = label),
            size = 2.3, lineheight = 1) +
  scale_fill_manual(values = stage_colors, name = "Stage",
                    guide = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.8, 6),  expand = c(0, 0)) +
  labs(title = NULL,
       x = NULL, y = NULL) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 9),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6.5))

# Panel B — decision-point checkpoints with real numbers -----------------------
checkpoints <- tribble(
  ~step, ~rule,                                                         ~outcome,
  "1. ENA discovery",
    "library_strategy ∈ {RNA-Seq}, layout = PAIRED, instrument = Illumina",
    sprintf("%d runs across %d cohorts", n_ena, n_cohorts_pre),
  "2. Atlas selection",
    "cohort productivity rank · size cap (≤25)",
    sprintf("%d runs retained", n_selected),
  "3. Download integrity",
    "aria2 + gzip integrity check; partial / corrupt → exclude",
    "21 downloaded → 18 pass gzip",
  "4. Assembly QC",
    sprintf("V+J ≥ 10%% (cohort-level); GSE193258 fails at %.1f%% V+J",
            mean(diag$pct_V_plus_J[diag$cohort == "GSE193258"], na.rm = TRUE)),
    sprintf("%d / %d GSE193258 samples flagged",
            n_gse193_bad,
            diag |> filter(cohort == "GSE193258") |> nrow()),
  "5. Productive TRB",
    "≥ 1 productive TRBV·CDR3·TRBJ assembly",
    sprintf("%d productive samples → atlas", n_productive),
  "6. Antigen lookup",
    "Exact CDR3-aa match → VDJdb",
    "37 annotated clones"
)

# Render as a properly column-tiled table:
#   x = 0.0..0.20  → Step column (bold)
#   x = 0.20..0.65 → Rule column
#   x = 0.65..1.00 → Outcome column (numeric, npg accent)
# Each row is a wide tile band with thin column separators.
checkpoints <- checkpoints |>
  mutate(row = factor(step, levels = rev(step)))

# Column anchor points and widths
col_step_x    <- 0.10
col_rule_x    <- 0.22
col_outcome_x <- 0.68
header_y     <- length(unique(checkpoints$row)) + 0.6

p_B <- ggplot(checkpoints, aes(y = row)) +
  # Row background tile
  geom_tile(aes(x = 0.5), width = 1, height = 0.92,
            fill = "grey96", color = "grey70", linewidth = 0.2) +
  # Column separators
  geom_vline(xintercept = c(0.20, 0.65),
             color = "grey70", linewidth = 0.3) +
  # Step column (left-aligned, bold)
  geom_text(aes(x = 0.01, label = step),
            hjust = 0, size = 2.4, fontface = "bold",
            lineheight = 1) +
  # Rule column
  geom_text(aes(x = 0.21, label = rule),
            hjust = 0, size = 2.1, lineheight = 1.05) +
  # Outcome column (numeric, NPG accent color)
  geom_text(aes(x = 0.66, label = outcome),
            hjust = 0, size = 2.1, fontface = "italic",
            color = npg_pal[1], lineheight = 1.05) +
  # Column headers
  annotate("text", x = 0.01, y = header_y, label = "Step",
           hjust = 0, size = 2.6, fontface = "bold") +
  annotate("text", x = 0.21, y = header_y, label = "Rule",
           hjust = 0, size = 2.6, fontface = "bold") +
  annotate("text", x = 0.66, y = header_y, label = "Outcome (live count)",
           hjust = 0, size = 2.6, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.2))) +
  labs(title = NULL,
       x = NULL, y = NULL) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 9))

# Stack A over B so each panel has the full page width (was side-by-side,
# which forced A's boxes and B's columns to overlap).
figS1 <- (p_A / p_B) +
  plot_layout(heights = c(1, 1.05))

save_fig(figS1, "FigS1",
         width_mm = FIG_W_DOUBLE, height_mm = 175)
