# Figure 5 — Antigen recognition landscape (7 panels)
#
# Real inputs:
#   data/atlas/vdjdb_matches.csv          (37 antigen-annotated clones)
#   data/metrics/atlas_diversity.csv      (n_clones per sample, cohort)
#   data/metrics/clone_size_distribution.csv  (counts for annotated CDR3s)
#
# Layout (asymmetric — give panel G extra vertical space):
#   "AB
#    CD
#    EF
#    GG"
#
# Output: paper/figures/fig5_antigen.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(stringr)
  library(tidyr)
  library(ggrepel)
})

vdj <- read_csv(file.path(REPO_ROOT, "data", "atlas", "vdjdb_matches.csv"),
                show_col_types = FALSE)
div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)
dist <- read_csv(file.path(METRICS_DIR, "clone_size_distribution.csv"),
                 show_col_types = FALSE)

# Join cohort onto vdj
vdj <- vdj |> left_join(div |> select(sample_id, cohort), by = "sample_id")

first_field <- function(x) {
  vapply(strsplit(as.character(x), ";"), function(v) {
    v <- v[!is.na(v) & nchar(v) > 0]
    if (length(v) == 0) NA_character_ else v[1]
  }, character(1))
}

vdj <- vdj |>
  mutate(species_first = first_field(species),
         mhc_first     = first_field(mhc_class),
         antigen_first = first_field(antigens),
         epitope_first = first_field(epitopes))

classify <- function(sp, ag) {
  sp <- toupper(sp); ag <- toupper(ag)
  if (is.na(sp)) sp <- ""
  if (is.na(ag)) ag <- ""
  viral_kw <- c("CMV", "EBV", "INFLUENZA", "SARS", "HIV", "MCPYV",
                "MCMV", "HCV", "HBV", "YFV", "HSV", "HPV", "DENGUE",
                "BZLF", "VACCINIA")
  tumor_kw <- c("HOMOSAPIENS", "MELANOMA", "NY-ESO", "MART", "MLANA",
                "MAGE", "WT1", "TUMOR", "CANCER")
  auto_kw  <- c("INS", "DIABETES", "MYELIN", "AUTO")
  if (any(sapply(viral_kw, function(k) grepl(k, sp) | grepl(k, ag)))) return("Viral")
  if (any(sapply(tumor_kw, function(k) grepl(k, sp) | grepl(k, ag)))) return("Tumor/self")
  if (any(sapply(auto_kw,  function(k) grepl(k, sp) | grepl(k, ag)))) return("Autoimmune")
  "Other"
}

vdj$category <- mapply(classify, vdj$species_first, vdj$antigen_first)

sample_meta <- div |> distinct(sample_id, cohort)

# ----- A. Total VDJdb hits per sample ----------------------------------------
hits_per_sample <- vdj |>
  count(sample_id, cohort, name = "n_hits")
hits_with_zero <- sample_meta |>
  left_join(hits_per_sample, by = c("sample_id", "cohort")) |>
  mutate(n_hits = tidyr::replace_na(n_hits, 0L)) |>
  arrange(cohort, desc(n_hits))

p_A <- hits_with_zero |>
  mutate(sample_short = short_id(sample_id),
         sample_short = factor(sample_short, levels = sample_short)) |>
  ggplot(aes(x = sample_short, y = n_hits, fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = n_hits), vjust = -0.3, size = 2.4) +
  scale_fill_npg(name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "VDJdb hits",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "top",
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- B. Antigen species pie / stacked bar (single bar) ----------------------
species_count <- vdj |>
  count(species_first, name = "n") |>
  arrange(desc(n)) |>
  mutate(label = paste0(species_first, " (", n, ")"))

p_B <- species_count |>
  ggplot(aes(x = "", y = n, fill = species_first)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  scale_fill_npg(name = "Species") +
  labs(x = NULL, y = NULL,
       title = NULL) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ----- C. MHC class breakdown per cohort -------------------------------------
mhc_cohort <- vdj |>
  count(cohort, mhc_first, name = "n") |>
  group_by(cohort) |>
  mutate(frac = n / sum(n)) |>
  ungroup()

p_C <- mhc_cohort |>
  ggplot(aes(x = cohort, y = frac, fill = mhc_first)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2, width = 0.65) +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            size = 2.6, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(MHCI = npg_pal[1], MHCII = npg_pal[3],
                               `MHCI,MHCII` = npg_pal[5]),
                    name = "MHC") +
  scale_y_continuous(labels = label_percent(), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of hits",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

# ----- D. Top-15 antigen-gene horizontal bar ---------------------------------
ag_top <- vdj |>
  filter(!is.na(antigen_first), nchar(antigen_first) > 0) |>
  count(antigen_first, name = "n") |>
  arrange(desc(n)) |>
  slice_head(n = 15)

p_D <- ag_top |>
  ggplot(aes(x = n, y = reorder(antigen_first, n))) +
  geom_col(fill = npg_pal[2], color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Annotated clones",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6))

# ----- E. Category breakdown per cohort (viral / tumor / autoimmune / other) -
cat_long <- vdj |>
  count(cohort, category, name = "n") |>
  group_by(cohort) |>
  mutate(frac = n / sum(n)) |>
  ungroup() |>
  mutate(category = factor(category,
                           levels = c("Viral", "Tumor/self",
                                      "Autoimmune", "Other")))

p_E <- cat_long |>
  ggplot(aes(x = cohort, y = frac, fill = category)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2, width = 0.65) +
  scale_fill_manual(values = c(`Viral` = npg_pal[1],
                               `Tumor/self` = npg_pal[3],
                               `Autoimmune` = npg_pal[4],
                               `Other` = "grey70"),
                    name = NULL) +
  scale_y_continuous(labels = label_percent(), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of hits",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

# ----- F. VDJdb hits vs n_clones scatter -------------------------------------
hit_vs_depth <- hits_with_zero |>
  left_join(div |> select(sample_id, n_clones), by = "sample_id")

p_F <- hit_vs_depth |>
  ggplot(aes(x = n_clones, y = n_hits, color = cohort)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "grey40", linetype = "dashed", linewidth = 0.4) +
  scale_color_npg(name = NULL) +
  scale_x_log10(labels = label_comma()) +
  labs(x = "Productive TRB clones (log10)",
       y = "VDJdb hits",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.20, 0.82),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- G. Annotated CDR3 list — labeled scatter (clone size vs category) -----
# Pull clone count for each annotated CDR3
dist_join <- dist |>
  select(sample_id, cdr3_aa, count, freq, cohort)
vdj_annotated <- vdj |>
  left_join(dist_join, by = c("sample_id", "cdr3_aa", "cohort")) |>
  filter(!is.na(count)) |>
  arrange(desc(count))

# Label top 12 by count
lbl <- vdj_annotated |>
  slice_head(n = 12) |>
  mutate(plot_label = sprintf("%s\n[%s | %s]",
                              cdr3_aa, antigen_first, species_first))

p_G <- vdj_annotated |>
  ggplot(aes(x = count, y = freq, color = category, shape = cohort)) +
  geom_point(size = 2.4, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = lbl,
    aes(label = plot_label),
    size = 2.2, lineheight = 0.95,
    box.padding = 0.4, point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.25,
    max.overlaps = 20,
    seed = 42, show.legend = FALSE
  ) +
  scale_color_manual(values = c(`Viral` = npg_pal[1],
                                `Tumor/self` = npg_pal[3],
                                `Autoimmune` = npg_pal[4],
                                `Other` = "grey50"),
                     name = "Category") +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_percent()) +
  labs(x = "Clone count (log10)",
       y = "Clone frequency in repertoire (log10)",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7, face = "bold"))

design <- "
AB
CD
EF
GG
"
fig5 <- p_A + p_B + p_C + p_D + p_E + p_F + p_G +
  plot_layout(design = design, heights = c(1, 1, 1, 1.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(fig5, "Fig5",
         width_mm = FIG_W_DOUBLE, height_mm = 235)
