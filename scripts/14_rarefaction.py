#!/usr/bin/env python
"""Rarefaction curves per sample.

For each productive sample, expand the clone table into a read-pool (one row per
read using duplicate_count multiplicities), sample fractions 0.05, 0.10, ...,
1.00, count unique clones at each fraction (mean of 5 bootstraps for robustness).

Output: data/metrics/rarefaction.csv (sample_id, cohort, fraction, unique_clones)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from thoracictcr.repertoire.repertoire_io import load_trust4_clones  # noqa: E402

TRUST4_DIR = REPO / "data" / "repertoire" / "trust4" / "geo_pilot"
DIVERSITY_CSV = REPO / "data" / "metrics" / "atlas_diversity.csv"
OUT_CSV = REPO / "data" / "metrics" / "rarefaction.csv"

FRACTIONS = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
N_BOOTSTRAP = 5
SEED = 20260527


def main() -> int:
    rng = np.random.default_rng(SEED)
    div = pd.read_csv(DIVERSITY_CSV)
    sample_cohort = dict(zip(div.sample_id, div.cohort))

    rows = []
    for sample_id, cohort in sample_cohort.items():
        sdir = TRUST4_DIR / sample_id
        clones = load_trust4_clones(sdir, chain="TRB", productive_only=True)
        if clones is None or len(clones) == 0:
            continue
        # Build read pool
        clone_ids = np.arange(len(clones))
        counts = clones["duplicate_count"].astype(int).to_numpy()
        pool = np.repeat(clone_ids, counts)
        total_reads = len(pool)
        if total_reads == 0:
            continue
        for frac in FRACTIONS:
            n_pick = max(1, int(round(total_reads * frac)))
            uniq_counts = []
            for _ in range(N_BOOTSTRAP):
                sample = rng.choice(pool, size=n_pick, replace=False)
                uniq_counts.append(int(len(np.unique(sample))))
            rows.append({
                "sample_id": sample_id, "cohort": cohort,
                "fraction": frac,
                "n_reads_sampled": n_pick,
                "unique_clones": float(np.mean(uniq_counts)),
                "unique_clones_sd": float(np.std(uniq_counts)),
                "total_reads": total_reads,
            })
    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV} ({len(df)} rows; {df.sample_id.nunique()} samples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
