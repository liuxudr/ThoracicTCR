# Supp Table S6 — GSE193258 exclusion evidence
#
# Caption: Per-sample diagnostic numbers used to justify excluding GSE193258
#   from the productive TRB atlas.
#
# Input:  data/metrics/gse193258_diagnostic_summary.csv
# Output: paper/tables/TableS6.tsv

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))

ds <- read_csv(file.path(METRICS_DIR, "gse193258_diagnostic_summary.csv"),
               show_col_types = FALSE)

n_problem <- ds |> filter(label == "Problematic (excluded)") |> nrow()
n_total_problem <- ds |> filter(grepl("GSE193258", cohort)) |> nrow()

ts6 <- ds |>
  mutate(interpretation = case_when(
    cohort == "GSE193258" & pct_V_plus_J < 10 ~
      sprintf("Supports exclusion: %.1f%% V+J assemblies (< 10%% threshold); %d V-only vs %d V+J",
              pct_V_plus_J, n_V_only, n_V_plus_J),
    cohort != "GSE193258" & pct_V_plus_J > 50 ~
      sprintf("Reference: well-behaved cohort, %.1f%% V+J assemblies", pct_V_plus_J),
    TRUE ~ "—"
  )) |>
  arrange(cohort, sample_id)

write_tsv(ts6, file.path(TBL_DIR, "TableS6.tsv"))
txt <- knitr::kable(ts6, format = "simple")
hdr <- sprintf("Supplementary Table S6. GSE193258 exclusion evidence (%d/%d GSE193258 samples have <10%% V+J assemblies).",
               n_problem, n_total_problem)
writeLines(c(hdr, "", txt),
           file.path(TBL_DIR, "TableS6.txt"))
message("Wrote: ", file.path(TBL_DIR, "TableS6.tsv"))
