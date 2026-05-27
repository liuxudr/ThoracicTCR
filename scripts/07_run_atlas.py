#!/usr/bin/env python3
"""Orchestrate atlas construction after TRUST4 + (optional) HLA + (optional) DBs.

Pipeline:
    1. Load atlas_manifest.csv
    2. For each sample, load TRUST4 *_report.tsv -> filtered productive TRB clones
    3. Compute diversity metrics
    4. Optional: load arcasHLA outputs from data/hla/arcashla/
    5. Optional: run VDJdb / McPAS lookup on aggregate clones
    6. Save ThoracicTCRAtlas object to data/atlas/thoracic_tcr_atlas.joblib

Outputs:
    data/atlas/thoracic_tcr_atlas.joblib   - main atlas object
    data/metrics/atlas_diversity.csv        - per-sample diversity table
    data/atlas/atlas_summary.json           - small summary JSON
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import pandas as pd
from rich.console import Console

from thoracictcr.atlas.atlas_builder import build_atlas
from thoracictcr.repertoire.repertoire_io import load_trust4_clones

console = Console()

TRUST4_ROOT = REPO / "data" / "repertoire" / "trust4" / "geo_pilot"
HLA_ROOT    = REPO / "data" / "hla" / "arcashla"
VDJDB_DIR   = REPO / "data" / "external" / "vdjdb"
MCPAS_DIR   = REPO / "data" / "external" / "mcpas"
ATLAS_DIR   = REPO / "data" / "atlas"
METRICS_DIR = REPO / "data" / "metrics"
ATLAS_DIR.mkdir(parents=True, exist_ok=True)
METRICS_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    manifest = pd.read_csv(REPO / "configs" / "atlas_manifest.csv")
    manifest = manifest.rename(columns={"run_accession": "sample_id"})
    manifest["cohort"] = manifest["queried_accession"]
    console.print(f"[cyan]Atlas manifest: {len(manifest)} samples[/cyan]")

    atlas = build_atlas(
        trust4_root=TRUST4_ROOT,
        metadata=manifest,
        hla_root=HLA_ROOT if HLA_ROOT.exists() else None,
    )

    # Save diversity table (most useful intermediate)
    div_csv = METRICS_DIR / "atlas_diversity.csv"
    if not atlas.diversity.empty:
        # Merge with cohort/treatment from metadata
        div = atlas.diversity.merge(
            manifest[["sample_id", "cohort", "library_layout", "instrument_model"]],
            on="sample_id", how="left",
        )
        div.to_csv(div_csv, index=False)
        console.print(f"[green]Diversity CSV -> {div_csv}[/green]")
    else:
        console.print("[yellow]No samples with productive TRB clones![/yellow]")

    # Optional: VDJdb + McPAS annotation
    if VDJDB_DIR.exists() and list(VDJDB_DIR.glob("*.tsv")) + list(VDJDB_DIR.glob("*.txt")):
        try:
            from thoracictcr.antigen.vdjdb_lookup import VDJdbMatcher, load_vdjdb
            vdjdb = load_vdjdb(VDJDB_DIR, chain="TRB")
            matcher = VDJdbMatcher(vdjdb)
            ag = matcher.annotate_per_sample(atlas.repertoire, cdr3_col="junction_aa")
            if not ag.empty:
                ag_csv = ATLAS_DIR / "vdjdb_matches.csv"
                ag.to_csv(ag_csv, index=False)
                atlas.antigen_profiles = ag
                console.print(f"[green]VDJdb hits CSV -> {ag_csv} ({len(ag)} clone-hits)[/green]")
            else:
                console.print("[yellow]No VDJdb hits[/yellow]")
        except FileNotFoundError as e:
            console.print(f"[yellow]Skip VDJdb: {e}[/yellow]")
        except Exception as e:
            console.print(f"[yellow]VDJdb lookup failed: {e}[/yellow]")
    else:
        console.print("[yellow]Skip VDJdb (run scripts/05_download_public_tcr_dbs.sh)[/yellow]")

    # Save atlas
    atlas_path = ATLAS_DIR / "thoracic_tcr_atlas.joblib"
    atlas.save(atlas_path)

    # Write summary
    summary = atlas.summary()
    summary_path = ATLAS_DIR / "atlas_summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    console.print(f"[green]Atlas summary:[/green]")
    for k, v in summary.items():
        console.print(f"  {k}: {v}")
    console.print(f"[green]Summary JSON -> {summary_path}[/green]")


if __name__ == "__main__":
    main()
