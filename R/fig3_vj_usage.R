# Figure 3 — V/J gene usage (8 panels)
#
# Real inputs:
#   data/metrics/vj_usage_long.csv
#   data/metrics/vj_pair_counts.csv
#
# Layout (4 rows x 2 cols, but G+H share a wide block):
#   "AB
#    CD
#    EF
#    GH"
#
# Output: paper/figures/fig3_vj_usage.{pdf,png}

source(file.path(Sys.getenv("THORACICTCR_REPO", "."), "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(scales)
  library(stringr)
  library(tidyr)
})

usage <- read_csv(file.path(METRICS_DIR, "vj_usage_long.csv"),
                  show_col_types = FALSE)
pairs <- read_csv(file.path(METRICS_DIR, "vj_pair_counts.csv"),
                  show_col_types = FALSE)

# Compute per-(sample, gene) frequency from counts (some samples may have zero
# for a gene → join with full sample×gene grid)
usage_v <- usage |> filter(gene_type == "V")
usage_j <- usage |> filter(gene_type == "J")

sample_meta <- usage |> distinct(sample_id, cohort)

# ---------------- A. Top-30 TRBV heatmap (sample x gene) ---------------------
top_v <- usage_v |>
  group_by(gene_name) |>
  summarise(total = sum(count), .groups = "drop") |>
  arrange(desc(total)) |> slice_head(n = 30) |> pull(gene_name)

mat_v <- usage_v |>
  filter(gene_name %in% top_v) |>
  tidyr::complete(sample_id = unique(usage_v$sample_id),
                  gene_name = top_v,
                  fill = list(frequency = 0)) |>
  left_join(sample_meta, by = "sample_id") |>
  mutate(gene_name = factor(gene_name, levels = top_v),
         sample_id = factor(sample_id, levels = sample_meta |>
                              arrange(cohort, sample_id) |> pull(sample_id)))

p_A <- mat_v |>
  mutate(sample_id = factor(short_id(as.character(sample_id)),
                            levels = short_id(levels(sample_id)))) |>
  ggplot(aes(x = gene_name, y = sample_id, fill = frequency)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient(low = "white", high = npg_pal[1],
                      labels = label_percent(accuracy = 1),
                      name = "Freq") +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ---------------- B. Top-12 TRBJ heatmap -------------------------------------
top_j <- usage_j |>
  group_by(gene_name) |>
  summarise(total = sum(count), .groups = "drop") |>
  arrange(desc(total)) |> slice_head(n = 12) |> pull(gene_name)

mat_j <- usage_j |>
  filter(gene_name %in% top_j) |>
  tidyr::complete(sample_id = unique(usage_j$sample_id),
                  gene_name = top_j,
                  fill = list(frequency = 0)) |>
  left_join(sample_meta, by = "sample_id") |>
  mutate(gene_name = factor(gene_name, levels = top_j),
         sample_id = factor(sample_id, levels = sample_meta |>
                              arrange(cohort, sample_id) |> pull(sample_id)))

p_B <- mat_j |>
  mutate(sample_id = factor(short_id(as.character(sample_id)),
                            levels = short_id(levels(sample_id)))) |>
  ggplot(aes(x = gene_name, y = sample_id, fill = frequency)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient(low = "white", high = npg_pal[2],
                      labels = label_percent(accuracy = 1),
                      name = "Freq") +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ---------------- C. TRBV rank curve per cohort -------------------------------
v_rank <- usage_v |>
  group_by(cohort, gene_name) |>
  summarise(freq = sum(count) / sum(sample_meta$cohort == first(cohort)),
            .groups = "drop") |>
  group_by(cohort) |>
  arrange(desc(freq), .by_group = TRUE) |>
  mutate(rank = row_number()) |>
  ungroup()

p_C <- v_rank |>
  filter(rank <= 40) |>
  ggplot(aes(x = rank, y = freq, color = cohort, group = cohort)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1, alpha = 0.85) +
  scale_color_npg(name = NULL) +
  scale_y_continuous(labels = label_number()) +
  labs(x = "TRBV rank", y = "Mean count / sample",
       title = NULL) +
  theme_thoracictcr() +
  theme(legend.position = c(0.75, 0.82),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 6))

# ---------------- D. TRBV Shannon-of-V-usage per sample ----------------------
v_shannon <- usage_v |>
  group_by(sample_id, cohort) |>
  summarise(
    total = sum(count),
    H_v = {
      f <- count[count > 0] / sum(count)
      -sum(f * log(f))
    },
    .groups = "drop"
  )

p_D <- v_shannon |>
  ggplot(aes(x = cohort, y = H_v, fill = cohort)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.85, size = 1.6, shape = 21,
              color = "black", stroke = 0.25) +
  scale_fill_npg(guide = "none") +
  labs(x = NULL, y = "Shannon over TRBV",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

# ---------------- E. Top-15 V-J pair frequency -------------------------------
pair_top <- pairs |>
  group_by(v_gene, j_gene) |>
  summarise(total = sum(count), .groups = "drop") |>
  arrange(desc(total)) |>
  slice_head(n = 15) |>
  mutate(pair = paste(v_gene, j_gene, sep = " · "),
         pair = factor(pair, levels = rev(pair)))

p_E <- pair_top |>
  ggplot(aes(x = total, y = pair)) +
  geom_col(fill = npg_pal[4], color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = total), hjust = -0.15, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = label_comma()) +
  labs(x = "Total count", y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 6))

# ---------------- F. V→J heatmap (top V × top J) -----------------------------
top_v_pairs <- pairs |>
  group_by(v_gene) |>
  summarise(s = sum(count), .groups = "drop") |>
  arrange(desc(s)) |> slice_head(n = 15) |> pull(v_gene)
top_j_pairs <- pairs |>
  group_by(j_gene) |>
  summarise(s = sum(count), .groups = "drop") |>
  arrange(desc(s)) |> slice_head(n = 8) |> pull(j_gene)

pair_mat <- pairs |>
  filter(v_gene %in% top_v_pairs, j_gene %in% top_j_pairs) |>
  group_by(v_gene, j_gene) |>
  summarise(count = sum(count), .groups = "drop") |>
  tidyr::complete(v_gene = top_v_pairs, j_gene = top_j_pairs,
                  fill = list(count = 0)) |>
  mutate(v_gene = factor(v_gene, levels = top_v_pairs),
         j_gene = factor(j_gene, levels = top_j_pairs))

p_F <- pair_mat |>
  ggplot(aes(x = j_gene, y = v_gene, fill = log10(count + 1))) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "white", high = npg_pal[5],
                      name = "log10\n(count+1)") +
  labs(x = "TRBJ", y = "TRBV",
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 5.5),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# ---------------- G. Convergent V genes per cohort ---------------------------
# Convergent: V genes present in ≥ 50% of samples in a cohort
samples_per_cohort <- sample_meta |>
  count(cohort, name = "n_samples")

v_presence <- usage_v |>
  filter(count > 0) |>
  distinct(sample_id, cohort, gene_name) |>
  count(cohort, gene_name, name = "n_present") |>
  left_join(samples_per_cohort, by = "cohort") |>
  mutate(presence_pct = n_present / n_samples)

conv_top <- v_presence |>
  filter(presence_pct >= 0.5) |>
  group_by(cohort) |>
  arrange(desc(presence_pct), .by_group = TRUE) |>
  slice_head(n = 10) |>
  ungroup()

p_G <- conv_top |>
  ggplot(aes(x = presence_pct, y = reorder(gene_name, presence_pct),
             fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.75) +
  facet_wrap(~cohort, scales = "free_y", ncol = 3,
             labeller = as_labeller(function(x) sub("^GSE", "", x))) +
  scale_fill_npg(guide = "none") +
  scale_x_continuous(labels = label_percent(accuracy = 1),
                     limits = c(0, 1.05),
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = c(0, 0.5, 1)) +
  labs(x = "Samples carrying the V gene (cohort label = GSE suffix)",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 5.8),
        axis.text.x = element_text(size = 5.5),
        strip.text = element_text(size = 7, face = "bold"),
        panel.spacing.x = unit(2, "mm"))

# ---------------- H. Cohort-specific TRBV (top 5 per cohort, log2-enrich) -----
# Compute per-cohort mean frequency, then log2 ratio cohort / mean(other cohorts)
v_cohort_freq <- usage_v |>
  group_by(cohort, gene_name) |>
  summarise(mean_freq = mean(frequency), .groups = "drop") |>
  tidyr::complete(cohort = unique(sample_meta$cohort),
                  gene_name = unique(usage_v$gene_name),
                  fill = list(mean_freq = 0))

v_specificity <- v_cohort_freq |>
  group_by(gene_name) |>
  mutate(
    other_mean = (sum(mean_freq) - mean_freq) / (n() - 1),
    log2fc = log2((mean_freq + 1e-4) / (other_mean + 1e-4))
  ) |>
  ungroup()

spec_top <- v_specificity |>
  filter(mean_freq >= 0.005) |>
  group_by(cohort) |>
  slice_max(log2fc, n = 5, with_ties = FALSE) |>
  ungroup()

p_H <- spec_top |>
  ggplot(aes(x = log2fc, y = reorder(gene_name, log2fc),
             fill = cohort)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.75) +
  facet_wrap(~cohort, scales = "free_y", ncol = 3,
             labeller = as_labeller(function(x) sub("^GSE", "", x))) +
  scale_fill_npg(guide = "none") +
  labs(x = "log2(cohort / mean other) — cohort = GSE suffix",
       y = NULL,
       title = NULL) +
  theme_thoracictcr() +
  theme(axis.text.y = element_text(size = 5.8),
        axis.text.x = element_text(size = 5.5),
        strip.text = element_text(size = 7, face = "bold"),
        panel.spacing.x = unit(2, "mm"))

design <- "
AB
CD
EF
GH
"
fig3 <- p_A + p_B + p_C + p_D + p_E + p_F + p_G + p_H +
  plot_layout(design = design, heights = c(1.1, 1, 1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_fig(fig3, "Fig3",
         width_mm = FIG_W_DOUBLE, height_mm = 240)
