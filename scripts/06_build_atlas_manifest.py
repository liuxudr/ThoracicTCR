#!/usr/bin/env python3
"""Build expanded mini-atlas manifest from already-resolved GEO/SRA cohorts.

Strategy:
  - Reuse configs/geo_thoracic_manifest.csv produced by 01b_
  - Exclude GSE193258 (verified to produce 0 productive TRB clones — bad library prep)
  - Pick a stratified set across the 3 good cohorts:
      * GSE145370 (ESCC, 18 runs, mean 45 GB/sample!) -> take smallest 5
      * GSE135222 (NSCLC nivo, 27 runs, mean 3 GB) -> take smallest 12
      * GSE126044 (NSCLC nivo, 16 runs, mean 5 GB) -> take smallest 8
  - Total target: 25 samples, ~80-100 GB (disk-bounded)

Outputs:
  - configs/atlas_manifest.csv               : selected 30 samples
  - configs/atlas_download.urls.txt          : aria2 input with MD5
  - configs/atlas_manifest_summary.csv       : per-cohort summary
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import pandas as pd

MANIFEST_IN = REPO / "configs" / "geo_thoracic_manifest.csv"
OUT_DIR = REPO / "configs"

# Cohort selection rule: (cohort_id, n_to_select)
COHORT_QUOTA = {
    "GSE145370": 5,    # ESCC — strongest TCR signal, but huge per-sample size (mean 45 GB)
    "GSE135222": 12,   # NSCLC nivo — small samples (mean 3 GB), include more
    "GSE126044": 8,    # NSCLC nivo
}
EXCLUDE_COHORTS = {"GSE193258"}  # confirmed bad library prep for TCR

# Already-downloaded samples to avoid re-downloading
EXISTING_FASTQ = REPO / "data" / "raw" / "fastq"


def main() -> None:
    if not MANIFEST_IN.exists():
        print(f"[ERROR] Run 01b_build_geo_pilot_manifest.py first.")
        sys.exit(1)

    df = pd.read_csv(MANIFEST_IN)
    print(f"[INFO] Loaded {len(df)} runs from {MANIFEST_IN}")

    df = df[~df["queried_accession"].isin(EXCLUDE_COHORTS)]
    print(f"[INFO] After excluding {EXCLUDE_COHORTS}: {len(df)} runs")

    selected_frames = []
    for cohort, n in COHORT_QUOTA.items():
        sub = df[df["queried_accession"] == cohort].copy()
        sub = sub.dropna(subset=["fastq_r1_https", "fastq_r2_https"])
        sub = sub[sub["library_layout"].fillna("").str.upper() == "PAIRED"]
        sub = sub.sort_values("total_fastq_gb")
        sub = sub.drop_duplicates("sample_accession", keep="first")
        picked = sub.head(n)
        print(f"  {cohort}: picked {len(picked)}/{len(sub)} (target {n})")
        selected_frames.append(picked)

    atlas = pd.concat(selected_frames, ignore_index=True)
    atlas_csv = OUT_DIR / "atlas_manifest.csv"
    atlas.to_csv(atlas_csv, index=False)
    print(f"[OK]   atlas manifest: {atlas_csv}  (n={len(atlas)})")

    summary = (
        atlas.groupby("queried_accession")
        .agg(
            n_runs=("run_accession", "count"),
            total_gb=("total_fastq_gb", "sum"),
            mean_gb=("total_fastq_gb", "mean"),
        )
        .round(2)
    )
    print("\n[INFO] Per-cohort selection summary:")
    print(summary.to_string())
    summary.to_csv(OUT_DIR / "atlas_manifest_summary.csv")
    print(f"\n[INFO] Atlas total download estimate: {atlas['total_fastq_gb'].sum():.1f} GB")

    # Compute delta: which fastq pairs are NOT yet on disk
    existing_runs: set[str] = set()
    if EXISTING_FASTQ.exists():
        for p in EXISTING_FASTQ.glob("*_R1.fastq.gz"):
            existing_runs.add(p.name.replace("_R1.fastq.gz", ""))
    print(f"\n[INFO] Already on disk: {len(existing_runs)} runs")

    new_atlas = atlas[~atlas["run_accession"].isin(existing_runs)]
    delta_gb = new_atlas["total_fastq_gb"].sum()
    print(f"[INFO] New downloads needed: {len(new_atlas)} runs, {delta_gb:.1f} GB")

    # Generate aria2 URL list (with MD5)
    urls_path = OUT_DIR / "atlas_download.urls.txt"
    with urls_path.open("w") as f:
        for _, r in new_atlas.iterrows():
            for mate in ("r1", "r2"):
                url = r.get(f"fastq_{mate}_https", "")
                md5 = r.get(f"fastq_{mate}_md5", "")
                if not url:
                    continue
                fname = f"{r['run_accession']}_{mate.upper()}.fastq.gz"
                f.write(f"{url}\n  dir=data/raw/fastq\n  out={fname}\n")
                if md5:
                    f.write(f"  checksum=md5={md5}\n")
    n_url_lines = sum(1 for _ in urls_path.open()) if urls_path.exists() else 0
    print(f"[OK]   aria2 URL list (with MD5): {urls_path}  ({n_url_lines} lines)")


if __name__ == "__main__":
    main()
