#!/usr/bin/env python3
"""Pilot QC + diversity for 10 TCGA-THYM TRUST4 outputs.

Validates the end-to-end pipeline:
  - TRUST4 produced parseable *_report.tsv
  - At least N productive TRB clones per sample
  - Diversity metrics computed without NaN/Inf
  - Top clones look reasonable (CDR3aa length 8-24, V/J gene names present)

Outputs:
  data/metrics/pilot_diversity.csv
  data/metrics/pilot_qc_report.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import pandas as pd
from rich.console import Console
from rich.table import Table

from thoracictcr.metrics.diversity import compute_diversity
from thoracictcr.repertoire.repertoire_io import load_trust4_clones, load_trust4_report

console = Console()

import os

# Resolve TRUST4 output dir from MODE env var (matches 03_trust4_pilot.sh)
# MODE=geo (default) -> geo_pilot/
# MODE=tcga         -> TCGA-THYM/
_MODE = os.environ.get("MODE", "geo").lower()
_TRUST4_SUBDIRS = {"geo": "geo_pilot", "tcga": "TCGA-THYM"}
TRUST4_ROOT = REPO / "data" / "repertoire" / "trust4" / _TRUST4_SUBDIRS.get(_MODE, "geo_pilot")
METRICS_DIR = REPO / "data" / "metrics"
METRICS_DIR.mkdir(parents=True, exist_ok=True)

MIN_PRODUCTIVE_TRB = 50  # pilot threshold; full pipeline can require more


def qc_sample(sample_dir: Path) -> tuple[dict, list[str]]:
    """Return (metrics_dict, list_of_qc_warnings)."""
    warnings: list[str] = []
    reports = list(sample_dir.glob("*_report.tsv"))
    if not reports:
        return {}, [f"NO_REPORT: {sample_dir.name}"]

    raw = load_trust4_report(reports[0])
    trb = raw[(raw["chain"] == "TRB") & raw["productive"]]
    n_trb = len(trb)

    if n_trb < MIN_PRODUCTIVE_TRB:
        warnings.append(f"LOW_TRB: {sample_dir.name} has only {n_trb} productive TRB clones")

    bad_len = trb[~trb["junction_aa"].str.len().between(8, 24)]
    if len(bad_len) > 0:
        warnings.append(
            f"BAD_CDR3_LEN: {sample_dir.name} {len(bad_len)} clones have CDR3 outside 8-24aa"
        )

    missing_v = trb["v_call"].isna().sum() + (trb["v_call"] == "").sum()
    if missing_v > 0:
        warnings.append(f"MISSING_V: {sample_dir.name} {missing_v} clones missing V gene")

    clones = load_trust4_clones(sample_dir, chain="TRB", productive_only=True, min_count=1)
    if clones is None or clones.empty:
        return {}, warnings + [f"EMPTY_AFTER_FILTER: {sample_dir.name}"]

    metrics = compute_diversity(clones)
    metrics["sample_id"] = sample_dir.name

    for k, v in metrics.items():
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            warnings.append(f"INVALID_METRIC: {sample_dir.name} {k}={v}")

    return metrics, warnings


def main() -> None:
    if not TRUST4_ROOT.exists():
        console.print(f"[red]TRUST4 output dir not found: {TRUST4_ROOT}[/red]")
        console.print("Run scripts/03_trust4_pilot.sh first.")
        sys.exit(1)

    sample_dirs = sorted([d for d in TRUST4_ROOT.iterdir() if d.is_dir()])
    console.print(f"[cyan]Found {len(sample_dirs)} pilot sample(s) under {TRUST4_ROOT}[/cyan]")

    rows: list[dict] = []
    all_warnings: list[str] = []

    for sd in sample_dirs:
        metrics, warns = qc_sample(sd)
        if metrics:
            rows.append(metrics)
        all_warnings.extend(warns)

    if not rows:
        console.print("[red]No valid TRUST4 outputs — pilot FAILED.[/red]")
        sys.exit(2)

    df = pd.DataFrame(rows)
    df = df[["sample_id"] + [c for c in df.columns if c != "sample_id"]]
    csv_path = METRICS_DIR / "pilot_diversity.csv"
    df.to_csv(csv_path, index=False)

    table = Table(title="Pilot diversity metrics (TRB)", show_lines=False)
    for col in ["sample_id", "n_clones", "shannon", "clonality", "top10_frac", "d50"]:
        table.add_column(col, justify="right" if col != "sample_id" else "left")
    for _, r in df.iterrows():
        table.add_row(
            r["sample_id"],
            str(int(r["n_clones"])),
            f"{r['shannon']:.3f}",
            f"{r['clonality']:.3f}",
            f"{r['top10_frac']:.3f}",
            str(int(r["d50"])),
        )
    console.print(table)

    qc_path = METRICS_DIR / "pilot_qc_report.txt"
    with qc_path.open("w") as f:
        f.write(f"Samples processed: {len(df)} / {len(sample_dirs)}\n")
        f.write(f"Warnings: {len(all_warnings)}\n\n")
        for w in all_warnings:
            f.write(w + "\n")
        f.write("\n=== Summary stats ===\n")
        f.write(
            df[["n_clones", "shannon", "clonality", "top10_frac", "d50"]]
            .describe()
            .round(3)
            .to_string()
        )

    console.print(f"\n[green]Diversity CSV: {csv_path}[/green]")
    console.print(f"[green]QC report:     {qc_path}[/green]")

    if all_warnings:
        console.print(f"\n[yellow]{len(all_warnings)} QC warnings (see report):[/yellow]")
        for w in all_warnings[:5]:
            console.print(f"  - {w}")
        if len(all_warnings) > 5:
            console.print(f"  ... and {len(all_warnings) - 5} more")
    else:
        console.print("\n[green]✓ All pilot samples passed QC[/green]")


if __name__ == "__main__":
    main()
