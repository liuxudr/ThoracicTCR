# Figure 6 — Cross-cohort sharing & motif analysis (5 panels)
#
# Real inputs:
#   data/metrics/shared_clones.csv          (pairwise jaccard, kmer overlaps)
#   data/metrics/top_shared_motifs.csv      (top k-mers, k = 3,4,5)
#   data/metrics/cumulative_shared_by_cohort.csv
#   data/metrics/clone_size_distribution.csv  (per-cohort unique motifs)
#
# Layout:
#   "AAB
#    CDE"
#
# Output: paper/figures/fig6_shared_motifs.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(stringr)
  library(tidyr)
})

shared <- read_csv(file.path(METRICS_DIR, "shared_clones.csv"),
                   show_col_types = FALSE)
motifs <- read_csv(file.path(METRICS_DIR, "top_shared_motifs.csv"),
                   show_col_types = FALSE)
cumshare <- read_csv(file.path(METRICS_DIR,
                               "cumulative_shared_by_cohort.csv"),
                     show_col_types = FALSE)
dist <- read_csv(file.path(METRICS_DIR, "clone_size_distribution.csv"),
                 show_col_types = FALSE)

cohorts <- sort(unique(c(shared$cohort_a, shared$cohort_b)))

# ----- A. Cohort × cohort Jaccard heatmap (symmetric, diagonal = 1) ---------
jacc_long <- shared |>
  select(cohort_a, cohort_b, jaccard) |>
  bind_rows(shared |>
              select(cohort_a = cohort_b, cohort_b = cohort_a, jaccard)) |>
  bind_rows(tibble(cohort_a = cohorts, cohort_b = cohorts, jaccard = 1))

p_A <- jacc_long |>
  mutate(cohort_a = factor(cohort_a, levels = cohorts),
         cohort_b = factor(cohort_b, levels = cohorts)) |>
  ggplot(aes(x = cohort_a, y = cohort_b, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", jaccard)),
            size = 2.6, fontface = "bold",
            color = "black") +
  scale_fill_gradient(low = "white", high = npg_pal[1],
                      limits = c(0, 1), name = "Jaccard") +
  labs(x = NULL, y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ----- B. Top-10 shared k-mer motifs across ≥ 2 cohorts ----------------------
motif_top <- motifs |>
  filter(n_cohorts >= 2) |>
  group_by(k) |>
  arrange(desc(n_samples), .by_group = TRUE) |>
  slice_head(n = 10) |>
  ungroup()

p_B <- motif_top |>
  ggplot(aes(x = n_samples, y = reorder(motif_kmer, n_samples),
             fill = factor(k))) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  facet_wrap(~ paste0("k = ", k), scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(`3` = npg_pal[1], `4` = npg_pal[3], `5` = npg_pal[5]),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Samples carrying motif",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6, family = "mono"),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0))

# ----- C. k = 3,4,5-mer overlap percent grouped bar by cohort pair ----------
kmer_long <- shared |>
  select(cohort_a, cohort_b, kmer3_overlap_pct, kmer4_overlap_pct,
         kmer5_overlap_pct) |>
  pivot_longer(starts_with("kmer"),
               names_to = "k", values_to = "overlap_pct") |>
  mutate(k = recode(k, kmer3_overlap_pct = "k=3",
                    kmer4_overlap_pct = "k=4",
                    kmer5_overlap_pct = "k=5"),
         pair = paste0(sub("^GSE", "", cohort_a), " · ",
                       sub("^GSE", "", cohort_b)))

p_C <- kmer_long |>
  ggplot(aes(x = pair, y = overlap_pct, fill = k)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           color = "black", linewidth = 0.15) +
  scale_fill_manual(values = c(`k=3` = npg_pal[1], `k=4` = npg_pal[3],
                               `k=5` = npg_pal[5]),
                    name = NULL) +
  labs(x = NULL, y = "k-mer overlap (%)",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8),
        legend.position = "top",
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ----- D. Cumulative shared clones vs # cohorts -------------------------------
p_D <- cumshare |>
  ggplot(aes(x = factor(n_cohorts_min), y = n_clones_shared)) +
  geom_col(fill = npg_pal[2], color = "black", linewidth = 0.2, width = 0.6) +
  geom_text(aes(label = n_clones_shared),
            vjust = -0.3, size = 2.6, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Clones shared by ≥ N cohorts",
       y = "Clone count",
       title = NULL) +
  theme_thoracictcr()

# ----- E. Per-cohort unique-motif count -------------------------------------
# Compute the number of distinct k=4 k-mers per cohort (within-cohort uniqueness)
get_kmers <- function(seq, k) {
  seq <- as.character(seq)
  if (is.na(seq) || nchar(seq) < k) return(character(0))
  vapply(seq_len(nchar(seq) - k + 1),
         function(i) substr(seq, i, i + k - 1),
         character(1))
}
k <- 4
cdr3_per_cohort <- dist |>
  distinct(cohort, cdr3_aa)

kmer_per_cohort <- cdr3_per_cohort |>
  rowwise() |>
  mutate(kmers = list(get_kmers(cdr3_aa, k))) |>
  ungroup() |>
  tidyr::unnest(cols = kmers)

cohort_kmer_sets <- split(kmer_per_cohort$kmers, kmer_per_cohort$cohort)
cohort_kmer_sets <- lapply(cohort_kmer_sets, unique)

unique_per_cohort <- tibble(
  cohort = names(cohort_kmer_sets),
  n_unique_kmers = sapply(names(cohort_kmer_sets), function(c) {
    others <- unique(unlist(cohort_kmer_sets[names(cohort_kmer_sets) != c]))
    length(setdiff(cohort_kmer_sets[[c]], others))
  }),
  n_total_kmers = sapply(cohort_kmer_sets, length)
)

p_E <- unique_per_cohort |>
  pivot_longer(c(n_total_kmers, n_unique_kmers),
               names_to = "set", values_to = "n") |>
  mutate(set = recode(set,
                      n_total_kmers = "Total k-mers",
                      n_unique_kmers = "Cohort-unique")) |>
  ggplot(aes(x = cohort, y = n, fill = set)) +
  geom_col(position = position_dodge(0.8), width = 0.7,
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = n),
            position = position_dodge(0.8),
            vjust = -0.3, size = 2.3) +
  scale_fill_manual(values = c(`Total k-mers` = "grey75",
                               `Cohort-unique` = npg_pal[4]),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = paste0("k=", k, " k-mer count"),
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = "top",
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

design <- "
AAB
CDE
"
fig6 <- p_A + p_B + p_C + p_D + p_E +
  plot_layout(design = design, heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(fig6, "Fig6",
         width_mm = FIG_W_DOUBLE, height_mm = 150)
