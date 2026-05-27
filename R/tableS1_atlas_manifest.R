# Supp Table S1 — Atlas manifest with QC status (25 samples)
#
# Caption: All 25 selected atlas samples with download / QC / TRUST4 status
#   and exclusion reasons.
#
# Inputs:
#   configs/atlas_manifest.csv
#   data/raw/fastq/                     (presence of fastq files)
#   data/repertoire/trust4/geo_pilot/   (presence of TRUST4 output)
#   data/metrics/atlas_diversity.csv    (productive n_clones)
#
# Output: paper/tables/TableS1.tsv

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))

man <- read_csv(file.path(REPO_ROOT, "configs", "atlas_manifest.csv"),
                show_col_types = FALSE)
div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
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

ts1 <- man |>
  rowwise() |>
  mutate(
    fastq_R1_status = check_fastq(run_accession, 1),
    fastq_R2_status = check_fastq(run_accession, 2),
    trust4_run      = check_trust4(run_accession),
    n_clones        = unname(n_clones_lookup[run_accession]),
    n_clones        = if (is.na(n_clones)) NA_real_ else n_clones,
    exclusion_reason = case_when(
      queried_accession == "GSE193258" ~ "GSE193258 cohort excluded (V-only assemblies)",
      fastq_R1_status %in% c("PARTIAL", "MISSING") |
        fastq_R2_status %in% c("PARTIAL", "MISSING")
                        ~ "fastq incomplete or missing",
      fastq_R1_status == "CORRUPT" | fastq_R2_status == "CORRUPT"
                        ~ "fastq corrupt (failed gzip integrity)",
      trust4_run == "no" ~ "TRUST4 not run",
      is.na(n_clones)   ~ "no productive TRB output",
      TRUE              ~ ""
    )
  ) |>
  ungroup() |>
  transmute(sample_id = run_accession,
            cohort    = queried_accession,
            library_layout, instrument = instrument_model,
            fastq_R1_status, fastq_R2_status,
            trust4_run, n_clones, exclusion_reason) |>
  arrange(cohort, sample_id)

write_tsv(ts1, file.path(TBL_DIR, "TableS1.tsv"))

txt <- knitr::kable(ts1, format = "simple")
writeLines(c("Supplementary Table S1. Atlas manifest with QC status.", "", txt),
           file.path(TBL_DIR, "TableS1.txt"))

message("Wrote: ", file.path(TBL_DIR, "TableS1.tsv"))
