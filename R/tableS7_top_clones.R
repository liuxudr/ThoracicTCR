# Main Table 3 — Top 20 expanded TRB clones across atlas + VDJdb annotation
#
# Caption: Top 20 expanded productive TRB clones across the atlas
#   (ranked by count) joined with their VDJdb antigen / epitope annotation.
#
# Inputs:
#   data/repertoire/trust4/geo_pilot/<sample>/<sample>_report.tsv
#   data/metrics/atlas_diversity.csv  (defines the 12 productive samples)
#   data/atlas/vdjdb_matches.csv
#
# Outputs:
#   paper/tables/TableS7.tsv
#   paper/tables/TableS7.txt

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))

div <- read_csv(file.path(METRICS_DIR, "atlas_diversity.csv"),
                show_col_types = FALSE)
vdj <- read_csv(file.path(REPO_ROOT, "data", "atlas", "vdjdb_matches.csv"),
                show_col_types = FALSE)

TRUST4_DIR <- file.path(REPO_ROOT, "data", "repertoire", "trust4", "geo_pilot")

is_trb <- function(v, j) {
  vc <- !is.na(v) & grepl("^TRBV", v)
  jc <- !is.na(j) & grepl("^TRBJ", j)
  vc | jc
}

read_one <- function(sid, cohort) {
  rpt <- file.path(TRUST4_DIR, sid, paste0(sid, "_report.tsv"))
  if (!file.exists(rpt)) return(NULL)
  # report.tsv header begins with '#count' — read with col names manually
  d <- tryCatch(
    read_tsv(rpt, comment = "", show_col_types = FALSE,
             col_names = c("count", "frequency", "cdr3_nt", "cdr3_aa",
                           "V", "D", "J", "C", "cid", "cid_full_length"),
             skip = 1),
    error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- d |>
    filter(!is.na(cdr3_aa), cdr3_aa != "out_of_frame", cdr3_aa != "",
           !grepl("\\*", cdr3_aa),
           !grepl("_", cdr3_aa),
           is_trb(V, J))
  if (nrow(d) == 0) return(NULL)
  d$sample_id <- sid; d$cohort <- cohort
  d
}

all_clones <- purrr::pmap_dfr(list(div$sample_id, div$cohort), read_one)

if (is.null(all_clones) || nrow(all_clones) == 0)
  stop("No productive TRB clones loaded from TRUST4 reports.")

top20 <- all_clones |>
  arrange(desc(count)) |>
  slice_head(n = 20) |>
  mutate(v_call = sub(",.*", "", V), j_call = sub(",.*", "", J),
         v_call = sub("\\*.*", "", v_call), j_call = sub("\\*.*", "", j_call))

vdj_join <- vdj |>
  group_by(cdr3_aa, sample_id) |>
  summarise(vdjdb_epitope = paste(unique(unlist(strsplit(paste(epitopes, collapse = ";"), ";"))),
                                  collapse = ";"),
            vdjdb_antigen = paste(unique(unlist(strsplit(paste(antigens, collapse = ";"), ";"))),
                                  collapse = ";"),
            vdjdb_species = paste(unique(unlist(strsplit(paste(species, collapse = ";"), ";"))),
                                  collapse = ";"),
            vdjdb_mhc_class = paste(unique(unlist(strsplit(paste(mhc_class, collapse = ";"), ";"))),
                                    collapse = ";"),
            .groups = "drop")

t3 <- top20 |>
  left_join(vdj_join, by = c("cdr3_aa", "sample_id")) |>
  select(sample_id, cohort, junction_aa = cdr3_aa, v_call, j_call,
         count, frequency, vdjdb_epitope, vdjdb_antigen, vdjdb_species,
         vdjdb_mhc_class) |>
  mutate(frequency = signif(frequency, 4),
         vdjdb_epitope = tidyr::replace_na(vdjdb_epitope, "-"),
         vdjdb_antigen = tidyr::replace_na(vdjdb_antigen, "-"),
         vdjdb_species = tidyr::replace_na(vdjdb_species, "-"),
         vdjdb_mhc_class = tidyr::replace_na(vdjdb_mhc_class, "-"))

write_tsv(t3, file.path(TBL_DIR, "TableS7.tsv"))

txt <- knitr::kable(t3, format = "simple")
writeLines(c("Table 3. Top 20 expanded productive TRB clones with VDJdb annotation.",
             "", txt),
           file.path(TBL_DIR, "TableS7.txt"))

message("Wrote: ", file.path(TBL_DIR, "TableS7.tsv"))
message("Wrote: ", file.path(TBL_DIR, "TableS7.txt"))
