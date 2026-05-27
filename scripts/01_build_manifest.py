#!/usr/bin/env python3
"""Build master manifest for all 5 thoracic TCGA projects, plus a 10-sample THYM pilot."""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from thoracictcr.data.manifest import (
    build_thoracic_manifest,
    select_pilot,
    write_gdc_manifest,
)

PROJECTS = ["TCGA-LUAD", "TCGA-LUSC", "TCGA-ESCA", "TCGA-MESO", "TCGA-THYM"]

OUT_DIR = REPO / "configs"
OUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    print(f"[INFO] Building manifest for projects: {PROJECTS}")
    df = build_thoracic_manifest(projects=PROJECTS, access="open", size=10000)
    master_csv = OUT_DIR / "thoracic_manifest.csv"
    df.to_csv(master_csv, index=False)
    print(f"[OK]   Master manifest: {master_csv}  (rows={len(df)})")

    print("\n[INFO] Per-project case counts and total BAM size:")
    summary = (
        df.groupby("project")
        .agg(
            n_files=("file_id", "count"),
            n_cases=("case_id", "nunique"),
            total_gb=("file_size_gb", "sum"),
            mean_gb=("file_size_gb", "mean"),
        )
        .round(1)
    )
    print(summary.to_string())
    summary.to_csv(OUT_DIR / "thoracic_manifest_summary.csv")

    # Pilot: 10 THYM Primary Tumor BAMs
    pilot = select_pilot(df, project="TCGA-THYM", n=10, seed=42)
    pilot_csv = OUT_DIR / "pilot_thym_10.csv"
    pilot.to_csv(pilot_csv, index=False)
    print(f"\n[OK]   Pilot (10 THYM) manifest: {pilot_csv}")
    print(pilot[["case_id", "submitter_id", "file_name", "file_size_gb"]].to_string(index=False))

    # GDC-client compatible manifest
    gdc_manifest = OUT_DIR / "pilot_thym_10.gdc_manifest.txt"
    write_gdc_manifest(pilot, gdc_manifest)
    print(f"\n[OK]   GDC-client manifest: {gdc_manifest}")
    total_gb = pilot["file_size_gb"].sum()
    print(f"[INFO] Pilot total download size: {total_gb:.1f} GB")


if __name__ == "__main__":
    main()
