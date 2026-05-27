# Supp Table S2 — Full diversity metrics dump (12 productive samples)
#
# Caption: Full per-sample diversity / clonality metric dump for the 12
#   productive atlas samples (Shannon, Simpson, inverse Simpson, Pielou
#   evenness, Gini, Chao1, Hill q=0/1/2, top-N fractions, D50, clonality).
#
# Input:  data/metrics/atlas_diversity.csv
# Output: paper/tables/TableS2.tsv

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))

div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)

ts2 <- div |>
  select(sample_id, cohort, library_layout, instrument_model,
         n_clones, total_reads, shannon, simpson, inverse_simpson,
         pielou_evenness, gini, chao1,
         hill_q0, hill_q1, hill_q2,
         top1_frac, top10_frac, top100_frac,
         d50, clonality) |>
  arrange(cohort, desc(n_clones)) |>
  mutate(across(c(shannon, simpson, inverse_simpson, pielou_evenness, gini,
                  chao1, hill_q0, hill_q1, hill_q2,
                  top1_frac, top10_frac, top100_frac, clonality),
                ~ round(.x, 4)))

write_tsv(ts2, file.path(TBL_DIR, "TableS2.tsv"))
txt <- knitr::kable(ts2, format = "simple")
writeLines(c("Supplementary Table S2. Full diversity metrics dump.", "", txt),
           file.path(TBL_DIR, "TableS2.txt"))
message("Wrote: ", file.path(TBL_DIR, "TableS2.tsv"))
