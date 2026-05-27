# Supp Table S3 — Full VDJdb match table with renamed columns
#
# Caption: All VDJdb matches (CDR3-aa exact lookup) across the productive
#   atlas, with column names made human-readable.
#
# Input:  data/atlas/vdjdb_matches.csv
# Output: paper/tables/TableS3.tsv

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))

vdj <- read_csv(file.path(REPO_ROOT, "data", "atlas", "vdjdb_matches.csv"),
                show_col_types = FALSE)

ts3 <- vdj |>
  transmute(sample_id,
            CDR3_aa     = cdr3_aa,
            n_hits      = n_hits,
            epitopes    = epitopes,
            antigens    = antigens,
            species     = species,
            MHC_class   = mhc_class) |>
  arrange(sample_id, desc(n_hits))

write_tsv(ts3, file.path(TBL_DIR, "TableS3.tsv"))
txt <- knitr::kable(ts3, format = "simple")
writeLines(c("Supplementary Table S3. Full VDJdb matches.", "", txt),
           file.path(TBL_DIR, "TableS3.txt"))
message("Wrote: ", file.path(TBL_DIR, "TableS3.tsv"))
