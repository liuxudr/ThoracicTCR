# Supp Fig S7 — VDJdb annotation breakdown (6 panels)
#
# Real input: data/atlas/vdjdb_matches.csv  +  data/metrics/atlas_diversity.csv
#
# Layout:
#   "ABC
#    DEF"
#
# Output: paper/figures/figS7_vdjdb_detail.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(tidyr)
})

vdj <- read_csv(file.path(REPO_ROOT, "data", "atlas", "vdjdb_matches.csv"),
                show_col_types = FALSE)
div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)

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

sample_order <- div |>
  arrange(cohort, desc(n_clones)) |>
  pull(sample_id)
vdj_samples <- intersect(sample_order, unique(vdj$sample_id))
vdj_samples_short <- short_id(vdj_samples)

# ----- A. Per-sample species stacked bar ------------------------------------
sp_long <- vdj |>
  count(sample_id, cohort, species_first, name = "n") |>
  mutate(sample_short = factor(short_id(sample_id),
                               levels = vdj_samples_short))

p_A <- sp_long |>
  ggplot(aes(x = sample_short, y = n, fill = species_first)) +
  geom_col(position = "stack", color = "black",
           linewidth = 0.2, width = 0.7) +
  scale_fill_npg(name = "Species") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "VDJdb hits",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 6),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- B. Per-sample MHC class stacked --------------------------------------
mhc_long <- vdj |>
  count(sample_id, cohort, mhc_first, name = "n") |>
  mutate(sample_short = factor(short_id(sample_id),
                               levels = vdj_samples_short))

p_B <- mhc_long |>
  ggplot(aes(x = sample_short, y = n, fill = mhc_first)) +
  geom_col(position = "stack", color = "black",
           linewidth = 0.2, width = 0.7) +
  scale_fill_manual(values = c(MHCI = npg_pal[1], MHCII = npg_pal[3],
                               `MHCI,MHCII` = npg_pal[5]),
                    name = "MHC") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "VDJdb hits",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 6),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- C. MHC-I top epitopes -------------------------------------------------
mhci_top <- vdj |>
  filter(mhc_first == "MHCI") |>
  count(epitope_first, name = "n") |>
  arrange(desc(n)) |>
  slice_head(n = 10)

p_C <- mhci_top |>
  ggplot(aes(x = n, y = reorder(epitope_first, n))) +
  geom_col(fill = npg_pal[1], color = "black",
           linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Hits", y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6, family = "mono"))

# ----- D. MHC-II top epitopes -----------------------------------------------
mhcii_top <- vdj |>
  filter(mhc_first == "MHCII") |>
  count(epitope_first, name = "n") |>
  arrange(desc(n)) |>
  slice_head(n = 10)

# Handle case when there are no MHC-II hits
if (nrow(mhcii_top) == 0) {
  mhcii_top <- tibble(epitope_first = "(no MHC-II hits)", n = 0)
}

p_D <- mhcii_top |>
  ggplot(aes(x = n, y = reorder(epitope_first, n))) +
  geom_col(fill = npg_pal[3], color = "black",
           linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Hits", y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6, family = "mono"))

# ----- E. Hit-count distribution (n_hits per CDR3) --------------------------
# Per CDR3 number of epitope IDs hit — uses n_hits column directly
p_E <- vdj |>
  ggplot(aes(x = n_hits)) +
  geom_histogram(binwidth = 1, fill = npg_pal[2],
                 color = "black", linewidth = 0.2) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "VDJdb hits per CDR3",
       y = "CDR3 clones",
       title = NULL) +
  theme_thoracictcr()

# ----- F. Per-cohort antigen-category breakdown -----------------------------
cat_long <- vdj |>
  count(cohort, category, name = "n") |>
  group_by(cohort) |>
  mutate(frac = n / sum(n)) |>
  ungroup() |>
  mutate(category = factor(category,
                           levels = c("Viral", "Tumor/self",
                                      "Autoimmune", "Other")))

p_F <- cat_long |>
  ggplot(aes(x = cohort, y = frac, fill = category)) +
  geom_col(position = "stack", color = "black",
           linewidth = 0.2, width = 0.65) +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            size = 2.4, color = "white", fontface = "bold") +
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
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
        plot.title = element_text(size = 9, face = "bold", hjust = 0))

design <- "
ABC
DEF
"
figS7 <- p_A + p_B + p_C + p_D + p_E + p_F +
  plot_layout(design = design) &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(figS7, "FigS7",
         width_mm = FIG_W_DOUBLE, height_mm = 175)
