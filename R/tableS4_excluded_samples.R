# Supp Table S4 — Excluded samples and reasons
#
# Caption: All atlas samples excluded from the final productive set, with
#   evidence notes. GSE193258 cohort fails assembly QC (<6% V+J assemblies);
#   additional exclusions due to incomplete download or corrupt gzip integrity
#   or no productive TRB output.
#
# Inputs: configs/atlas_manifest.csv, data/raw/fastq/, data/repertoire/trust4/,
#         data/metrics/atlas_diversity.csv
# Output: paper/tables/TableS4.tsv

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))

# Reuse the S1 logic but keep only excluded
man <- read_csv(file.path(REPO_ROOT, "configs", "atlas_manifest.csv"),
                show_col_types = FALSE)
div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)
diag_sum <- read_csv(file.path(METRICS_DIR,
                               "gse193258_diagnostic_summary.csv"),
                     show_col_types = FALSE)

FQ_DIR <- file.path(REPO_ROOT, "data", "raw", "fastq")
T4_DIR <- file.path(REPO_ROOT, "data", "repertoire", "trust4", "geo_pilot")

check_fastq <- function(sid, read) {
  base <- file.path(FQ_DIR, paste0(sid, "_R", read, ".fastq.gz"))
  partial <- paste0(base, ".aria2")
  if (file.exists(base) && file.info(base)$size > 1000) return("OK")
  if (file.exists(partial)) return("PARTIAL")
  if (file.exists(base)) return("CORRUPT")
  "MISSING"
}
check_trust4 <- function(sid) {
  rpt <- file.path(T4_DIR, sid, paste0(sid, "_report.tsv"))
  if (file.exists(rpt) && file.info(rpt)$size > 0) "yes" else "no"
}
n_clones_lookup <- setNames(div$n_clones, div$sample_id)

gse193_evidence <- function(sid) {
  row <- diag_sum |> filter(sample_id == sid)
  if (nrow(row) == 0) return("V-only assemblies dominate (cohort-level QC)")
  sprintf("V+J=%.1f%% (n_assemblies=%d); median CDR3 nt=%.0f",
          row$pct_V_plus_J[1], row$n_assemblies[1],
          row$median_cdr3_nt_len[1])
}

excluded <- man |>
  rowwise() |>
  mutate(
    sid = run_accession,
    cohort = queried_accession,
    fastq_R1_status = check_fastq(sid, 1),
    fastq_R2_status = check_fastq(sid, 2),
    trust4_run      = check_trust4(sid),
    n_clones        = unname(n_clones_lookup[sid]),
    in_productive   = !is.na(n_clones) && cohort != "GSE193258",
    exclusion_reason = case_when(
      cohort == "GSE193258"
        ~ "GSE193258 cohort excluded (assembly QC failure)",
      fastq_R1_status %in% c("PARTIAL", "MISSING") |
        fastq_R2_status %in% c("PARTIAL", "MISSING")
        ~ "fastq incomplete or missing",
      fastq_R1_status == "CORRUPT" | fastq_R2_status == "CORRUPT"
        ~ "fastq corrupt (gzip integrity)",
      trust4_run == "no"
        ~ "TRUST4 not run",
      is.na(n_clones)
        ~ "no productive TRB",
      TRUE ~ ""
    ),
    evidence_note = case_when(
      cohort == "GSE193258" ~ gse193_evidence(sid),
      fastq_R1_status %in% c("PARTIAL", "MISSING") |
        fastq_R2_status %in% c("PARTIAL", "MISSING")
        ~ paste0("R1=", fastq_R1_status, "; R2=", fastq_R2_status),
      fastq_R1_status == "CORRUPT" | fastq_R2_status == "CORRUPT"
        ~ "gzip -t failed on one or both reads",
      trust4_run == "no" ~ "no <sample>_report.tsv produced",
      is.na(n_clones) ~ "TRUST4 report empty after productive TRB filter",
      TRUE ~ ""
    )
  ) |>
  ungroup() |>
  filter(exclusion_reason != "" | !in_productive) |>
  filter(exclusion_reason != "") |>
  transmute(sample_id = sid, cohort, exclusion_reason, evidence_note) |>
  arrange(cohort, sample_id)

write_tsv(excluded, file.path(TBL_DIR, "TableS4.tsv"))
txt <- knitr::kable(excluded, format = "simple")
writeLines(c("Supplementary Table S4. Excluded atlas samples with reasons.",
             "", txt),
           file.path(TBL_DIR, "TableS4.txt"))
message("Wrote: ", file.path(TBL_DIR, "TableS4.tsv"))
