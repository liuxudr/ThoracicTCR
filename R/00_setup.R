# R setup for ThoracicTCR SCI figures
# All publication figures (Fig 1–6 + supplementary) should source() this file first.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Per-figure modules load these as needed (kept out of mandatory list to allow
# bootstrapping on systems where they're missing):
#   library(ComplexHeatmap)
#   library(circlize)
#   library(survminer)
#   library(forestplot)
#   library(gtsummary)
#   library(immunarch)
#   library(ggseqlogo)
#   library(ggridges)
#   library(ggrepel)

# Paths
REPO_ROOT <- normalizePath(
  Sys.getenv("THORACICTCR_REPO", file.path(dirname(sys.frame(1)$ofile), "..")),
  mustWork = FALSE
)
FIG_DIR <- file.path(REPO_ROOT, "paper", "figures")
TBL_DIR <- file.path(REPO_ROOT, "paper", "tables")
DATA_DIR <- file.path(REPO_ROOT, "data")
METRICS_DIR <- file.path(DATA_DIR, "metrics")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

# NPG (Nature Publishing Group) palette via ggsci
npg_pal <- pal_npg("nrc")(10)
thoracic_cancer_colors <- c(
  LUAD = npg_pal[1],
  LUSC = npg_pal[2],
  ESCA = npg_pal[3],
  MESO = npg_pal[4],
  THYM = npg_pal[5]
)

# Lancet palette (used for survival panels)
lancet_pal <- pal_lancet("lanonc")(9)

# JCO palette (used for forest plots)
jco_pal <- pal_jco("default")(10)

# Project theme
theme_thoracictcr <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text       = element_text(color = "black"),
      axis.title      = element_text(color = "black", face = "bold"),
      axis.line       = element_line(color = "black", linewidth = 0.3),
      axis.ticks      = element_line(color = "black", linewidth = 0.3),
      legend.position = "right",
      legend.title    = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5),
      plot.tag         = element_text(face = "bold", size = base_size + 2)
    )
}
theme_set(theme_thoracictcr())

# Standard SCI figure widths (mm) for single/double-column journals
FIG_W_SINGLE <- 89
FIG_W_DOUBLE <- 183

# Helper to abbreviate SRR / ERR / DRR sample IDs to the last 5 digits.
# Keeps a one-letter prefix so accession class is still visible (S/E/D).
# Non-SRR strings (already short, e.g. cohort names) are returned unchanged.
short_id <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^[SED]RR[0-9]+$", x),
         paste0(substr(x, 1, 1), substr(x, nchar(x) - 4, nchar(x))),
         x)
}

# Helper to save 600-dpi PDF + PNG side-by-side
save_fig <- function(plot, name, width_mm = FIG_W_DOUBLE, height_mm = 120) {
  pdf_path <- file.path(FIG_DIR, paste0(name, ".pdf"))
  png_path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width_mm, height = height_mm, units = "mm",
         device = cairo_pdf)
  ggsave(png_path, plot, width = width_mm, height = height_mm, units = "mm",
         dpi = 600)
  message("Saved: ", pdf_path)
  message("Saved: ", png_path)
}

message("ThoracicTCR R setup loaded. Palette: ", paste(names(thoracic_cancer_colors), collapse = ", "))
