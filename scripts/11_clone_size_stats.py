#!/usr/bin/env python
"""Compute per-sample clone size distribution + CDR3 length distribution.

Outputs:
  - data/metrics/clone_size_distribution.csv
      (sample_id, cohort, rank, count, freq, log10_freq, cdr3_aa, cdr3_len)
  - data/metrics/clone_size_summary.csv
      per-sample summary including top1/top10/top100 fractions and # hyperexpanded
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
OUT_DIST = REPO / "data" / "metrics" / "clone_size_distribution.csv"
OUT_SUMMARY = REPO / "data" / "metrics" / "clone_size_summary.csv"

HYPEREXPANDED_FREQ = 0.05


def main() -> int:
    div = pd.read_csv(DIVERSITY_CSV)
    sample_cohort = dict(zip(div.sample_id, div.cohort))

    dist_rows = []
    summary_rows = []
    for sample_id, cohort in sample_cohort.items():
        sdir = TRUST4_DIR / sample_id
        clones = load_trust4_clones(sdir, chain="TRB", productive_only=True)
        if clones is None or len(clones) == 0:
            continue
        clones = clones.sort_values("duplicate_count", ascending=False).reset_index(drop=True)
        total = clones["duplicate_count"].sum()
        clones["rank"] = np.arange(1, len(clones) + 1)
        clones["freq"] = clones["duplicate_count"] / total
        clones["log10_freq"] = np.log10(clones["freq"])
        clones["cdr3_aa"] = clones["junction_aa"]
        clones["cdr3_len"] = clones["junction_aa"].str.len()
        clones["sample_id"] = sample_id
        clones["cohort"] = cohort
        dist_rows.append(clones[["sample_id", "cohort", "rank",
                                  "duplicate_count", "freq", "log10_freq",
                                  "cdr3_aa", "cdr3_len"]]
                         .rename(columns={"duplicate_count": "count"}))

        top1 = float(clones["freq"].iloc[0]) if len(clones) >= 1 else 0.0
        top10 = float(clones["freq"].iloc[:10].sum())
        top100 = float(clones["freq"].iloc[:100].sum())
        n_hyper = int((clones["freq"] > HYPEREXPANDED_FREQ).sum())
        summary_rows.append({
            "sample_id": sample_id, "cohort": cohort,
            "n_clones": len(clones), "total_reads": int(total),
            "top1_frac": top1, "top10_frac": top10, "top100_frac": top100,
            "n_hyperexpanded": n_hyper,
            "median_cdr3_len": float(clones["cdr3_len"].median()),
            "mean_cdr3_len": float(clones["cdr3_len"].mean()),
        })

    if not dist_rows:
        print("No data", file=sys.stderr)
        return 1
    dist = pd.concat(dist_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    OUT_DIST.parent.mkdir(parents=True, exist_ok=True)
    dist.to_csv(OUT_DIST, index=False)
    summary.to_csv(OUT_SUMMARY, index=False)
    print(f"Wrote {OUT_DIST} ({len(dist)} rows)")
    print(f"Wrote {OUT_SUMMARY} ({len(summary)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
