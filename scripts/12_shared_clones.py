#!/usr/bin/env python
"""Cross-cohort CDR3 sharing analysis.

Outputs:
  - data/metrics/shared_clones.csv (cohort_a, cohort_b, n_shared_exact, jaccard,
                                    kmer3_overlap_pct, kmer4_overlap_pct, kmer5_overlap_pct)
  - data/metrics/top_shared_motifs.csv (motif_kmer, k, n_cohorts, n_samples, count)
  - data/metrics/top_shared_clones.csv (cdr3_aa, n_samples, cohorts, total_count)
  - data/metrics/cohort_pairs_kmer_matrix.csv (long format for k=3,4,5)
"""
from __future__ import annotations

import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from thoracictcr.repertoire.repertoire_io import load_trust4_clones  # noqa: E402

TRUST4_DIR = REPO / "data" / "repertoire" / "trust4" / "geo_pilot"
DIVERSITY_CSV = REPO / "data" / "metrics" / "atlas_diversity.csv"
OUT_PAIRS = REPO / "data" / "metrics" / "shared_clones.csv"
OUT_MOTIFS = REPO / "data" / "metrics" / "top_shared_motifs.csv"
OUT_TOP_CLONES = REPO / "data" / "metrics" / "top_shared_clones.csv"


def kmers(seq: str, k: int) -> set[str]:
    if not isinstance(seq, str) or len(seq) < k:
        return set()
    # strip canonical C/F termini for biological motifs (keep stems)
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def main() -> int:
    div = pd.read_csv(DIVERSITY_CSV)
    sample_cohort = dict(zip(div.sample_id, div.cohort))

    cohort_cdr3: dict[str, set[str]] = defaultdict(set)
    cohort_kmers: dict[str, dict[int, set[str]]] = defaultdict(lambda: {3: set(), 4: set(), 5: set()})
    sample_cdr3: dict[str, set[str]] = {}
    cdr3_to_samples: dict[str, set[str]] = defaultdict(set)
    cdr3_to_count: dict[str, int] = defaultdict(int)

    for sample_id, cohort in sample_cohort.items():
        sdir = TRUST4_DIR / sample_id
        clones = load_trust4_clones(sdir, chain="TRB", productive_only=True)
        if clones is None or len(clones) == 0:
            continue
        seqs = set(clones["junction_aa"].dropna().astype(str).tolist())
        sample_cdr3[sample_id] = seqs
        cohort_cdr3[cohort].update(seqs)
        for s in seqs:
            cdr3_to_samples[s].add(sample_id)
        for s, c in zip(clones["junction_aa"], clones["duplicate_count"]):
            if isinstance(s, str):
                cdr3_to_count[s] += int(c)
        for k in (3, 4, 5):
            for s in seqs:
                cohort_kmers[cohort][k].update(kmers(s, k))

    cohorts = sorted(cohort_cdr3.keys())

    # Pairwise overlap
    pair_rows = []
    for a, b in combinations(cohorts, 2):
        A, B = cohort_cdr3[a], cohort_cdr3[b]
        inter = A & B
        union = A | B
        jac = len(inter) / len(union) if union else 0.0
        row = {"cohort_a": a, "cohort_b": b,
               "n_shared_exact": len(inter),
               "jaccard": jac,
               "n_a": len(A), "n_b": len(B)}
        for k in (3, 4, 5):
            ka, kb = cohort_kmers[a][k], cohort_kmers[b][k]
            kinter = ka & kb
            kunion = ka | kb
            row[f"kmer{k}_overlap_pct"] = 100 * len(kinter) / len(kunion) if kunion else 0.0
            row[f"kmer{k}_shared"] = len(kinter)
        pair_rows.append(row)
    # also include diagonal pseudo entries (self)
    for c in cohorts:
        row = {"cohort_a": c, "cohort_b": c,
               "n_shared_exact": len(cohort_cdr3[c]),
               "jaccard": 1.0,
               "n_a": len(cohort_cdr3[c]), "n_b": len(cohort_cdr3[c])}
        for k in (3, 4, 5):
            row[f"kmer{k}_overlap_pct"] = 100.0
            row[f"kmer{k}_shared"] = len(cohort_kmers[c][k])
        pair_rows.append(row)
    pair_df = pd.DataFrame(pair_rows)
    OUT_PAIRS.parent.mkdir(parents=True, exist_ok=True)
    pair_df.to_csv(OUT_PAIRS, index=False)

    # Top shared CDR3-aa (those appearing in >=2 samples)
    rows = []
    for s, samples in cdr3_to_samples.items():
        if len(samples) >= 2:
            cohorts_s = sorted({sample_cohort[x] for x in samples})
            rows.append({
                "cdr3_aa": s,
                "n_samples": len(samples),
                "n_cohorts": len(cohorts_s),
                "cohorts": ";".join(cohorts_s),
                "samples": ";".join(sorted(samples)),
                "total_count": cdr3_to_count[s],
                "cdr3_len": len(s),
            })
    top_clones = (pd.DataFrame(rows)
                  .sort_values(["n_samples", "total_count"], ascending=[False, False])
                  .reset_index(drop=True)) if rows else pd.DataFrame()
    top_clones.to_csv(OUT_TOP_CLONES, index=False)

    # Top shared k-mers across cohorts (k=3,4,5)
    motif_rows = []
    for k in (3, 4, 5):
        # count cohorts each k-mer appears in
        kmer_cohort: dict[str, set[str]] = defaultdict(set)
        kmer_sample: dict[str, set[str]] = defaultdict(set)
        for sample_id, seqs in sample_cdr3.items():
            cohort = sample_cohort[sample_id]
            for s in seqs:
                for kk in kmers(s, k):
                    kmer_cohort[kk].add(cohort)
                    kmer_sample[kk].add(sample_id)
        # top motifs by n_cohorts then n_samples
        top_motifs = sorted(kmer_cohort.items(),
                            key=lambda kv: (len(kv[1]), len(kmer_sample[kv[0]])),
                            reverse=True)[:30]
        for motif, cohort_set in top_motifs:
            motif_rows.append({
                "motif_kmer": motif,
                "k": k,
                "n_cohorts": len(cohort_set),
                "n_samples": len(kmer_sample[motif]),
                "cohorts": ";".join(sorted(cohort_set)),
            })
    pd.DataFrame(motif_rows).to_csv(OUT_MOTIFS, index=False)

    # Cumulative shared count vs number of cohorts: how many CDR3 appear in >=N cohorts
    n_cohorts = len(cohorts)
    cum_rows = []
    for n in range(1, n_cohorts + 1):
        n_clones = sum(1 for s, samples in cdr3_to_samples.items()
                       if len({sample_cohort[x] for x in samples}) >= n)
        cum_rows.append({"n_cohorts_min": n, "n_clones_shared": n_clones})
    pd.DataFrame(cum_rows).to_csv(REPO / "data" / "metrics" / "cumulative_shared_by_cohort.csv",
                                  index=False)

    print(f"Wrote {OUT_PAIRS} ({len(pair_df)} rows, {len(cohorts)} cohorts)")
    print(f"Wrote {OUT_TOP_CLONES} ({len(top_clones)} multi-sample CDR3s)")
    print(f"Wrote {OUT_MOTIFS} ({len(motif_rows)} top motifs across k=3,4,5)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
