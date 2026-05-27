#!/usr/bin/env python3
"""Build GEO/SRA pilot manifest for ThoracicTCR.

Uses ENA's filereport API (no toolkit, no API key). Produces:
  - configs/geo_thoracic_manifest.csv   : all queryable cohorts merged
  - configs/pilot_geo_10.csv            : 10-sample pilot subset (paired-end, smallest)
  - configs/pilot_geo_10.urls.txt       : aria2-compatible URL list for download

Default cohorts (open-access, NSCLC neoadjuvant immunotherapy with MPR labels):
  - GSE207422 / PRJNA856944 — neoadjuvant pembro/nivo, NSCLC (~21 patients)
  - GSE193258 / PRJNA793452 — neoadjuvant pembrolizumab, NSCLC (~25 patients)
  - GSE126044 / PRJNA517284 — nivolumab in NSCLC (~16, response labels)
  - GSE135222 / PRJNA532380 — nivolumab in NSCLC (~27, response labels)
  - GSE145370              — paired pre/post ESCC chemo (~10), Chinese cohort
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from thoracictcr.data.geo_manifest import build_geo_pilot_manifest, select_pilot

# GEO accessions (ENA resolves GSE -> SRP automatically)
PILOT_COHORTS = [
    "GSE207422",   # NSCLC neoadjuvant pembro/nivo (MPR labels)
    "GSE193258",   # NSCLC neoadjuvant pembro (MPR labels)
    "GSE126044",   # NSCLC nivolumab (response)
    "GSE135222",   # NSCLC nivolumab (response)
    "GSE145370",   # ESCC chemo paired (Chinese cohort)
]

PILOT_N = 10
OUT_DIR = REPO / "configs"
OUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    print(f"[INFO] Querying ENA for {len(PILOT_COHORTS)} cohorts: {PILOT_COHORTS}")
    df = build_geo_pilot_manifest(PILOT_COHORTS, rnaseq_only=True)

    if df.empty:
        print("[ERROR] ENA returned no RNA-Seq runs for any accession.")
        sys.exit(1)

    master_csv = OUT_DIR / "geo_thoracic_manifest.csv"
    df.to_csv(master_csv, index=False)
    print(f"[OK]   Master manifest: {master_csv}  (n_runs={len(df)})")

    print("\n[INFO] Per-cohort summary:")
    summary = (
        df.groupby("queried_accession")
        .agg(
            n_runs=("run_accession", "count"),
            n_samples=("sample_accession", "nunique"),
            total_gb=("total_fastq_gb", "sum"),
            mean_gb=("total_fastq_gb", "mean"),
            layouts=("library_layout", lambda s: ",".join(sorted(s.dropna().unique()))),
        )
        .round(2)
    )
    print(summary.to_string())
    summary.to_csv(OUT_DIR / "geo_thoracic_manifest_summary.csv")

    pilot = select_pilot(df, n=PILOT_N, seed=42, prefer_paired=True, max_gb_per_sample=20.0)
    pilot_csv = OUT_DIR / "pilot_geo_10.csv"
    pilot.to_csv(pilot_csv, index=False)
    print(f"\n[OK]   Pilot ({PILOT_N}) manifest: {pilot_csv}")

    show_cols = [
        "queried_accession", "run_accession", "sample_accession",
        "library_layout", "instrument_model", "total_fastq_gb",
    ]
    print(pilot[show_cols].to_string(index=False))
    print(f"\n[INFO] Pilot total download size: {pilot['total_fastq_gb'].sum():.2f} GB")

    # Generate aria2 URL list with MD5 checksum (one block per mate)
    urls_path = OUT_DIR / "pilot_geo_10.urls.txt"
    with urls_path.open("w") as f:
        for _, r in pilot.iterrows():
            for mate in ("r1", "r2"):
                url = r.get(f"fastq_{mate}_https", "")
                md5 = r.get(f"fastq_{mate}_md5", "")
                if not url:
                    continue
                fname = f"{r['run_accession']}_{mate.upper()}.fastq.gz"
                f.write(f"{url}\n  dir=data/raw/fastq\n  out={fname}\n")
                if md5:
                    f.write(f"  checksum=md5={md5}\n")
    print(f"[OK]   aria2 URL list (with MD5): {urls_path}")


if __name__ == "__main__":
    main()
