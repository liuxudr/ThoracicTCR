"""``thoracictcr atlas ...`` — build the pan-thoracic TCR atlas object."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer
from rich.console import Console

from ._common import relpath
from ..atlas.atlas_builder import build_atlas

app = typer.Typer(
    name="atlas",
    help="Assemble per-sample TRUST4 / HLA / antigen annotations into a single atlas object.",
    no_args_is_help=True,
)
console = Console()


@app.command("build")
def build(
    manifest: Path = typer.Option(
        relpath("configs/atlas_manifest.csv"), "--manifest", "-m",
        help="Atlas manifest CSV (output of 'manifest atlas').",
    ),
    trust4_dir: Path = typer.Option(
        relpath("data/repertoire/trust4/geo_pilot"), "--trust4-dir", "-t",
    ),
    hla_dir: Path = typer.Option(
        relpath("data/hla/arcashla"), "--hla-dir",
        help="arcasHLA output dir (optional; skipped if missing).",
    ),
    vdjdb_dir: Path = typer.Option(
        relpath("data/external/vdjdb"), "--vdjdb-dir",
        help="VDJdb mirror dir (optional).",
    ),
    out_dir: Path = typer.Option(relpath("data/atlas"), "--out-dir", "-o"),
    metrics_dir: Path = typer.Option(
        relpath("data/metrics"), "--metrics-dir",
        help="Where to write per-sample atlas_diversity.csv (consumed by metrics commands).",
    ),
    chain: str = typer.Option("TRB", "--chain"),
) -> None:
    """Build the atlas joblib + diversity CSV + summary JSON.

    Optional steps (HLA, VDJdb antigen lookup) are silently skipped if their
    input directories are absent.
    """
    if not manifest.exists():
        console.print(f"[red]Manifest not found: {manifest}[/red]")
        raise typer.Exit(1)

    meta = pd.read_csv(manifest).rename(columns={"run_accession": "sample_id"})
    meta["cohort"] = meta.get("queried_accession", meta.get("cohort", ""))
    console.print(f"[cyan]Atlas manifest: {len(meta)} samples[/cyan]")

    atlas = build_atlas(
        trust4_root=trust4_dir,
        metadata=meta,
        hla_root=hla_dir if hla_dir.exists() else None,
        chain=chain,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    metrics_dir.mkdir(parents=True, exist_ok=True)

    div_csv = metrics_dir / "atlas_diversity.csv"
    if not atlas.diversity.empty:
        keep = ["sample_id"] + [c for c in ("cohort", "library_layout", "instrument_model")
                                 if c in meta.columns]
        div = atlas.diversity.merge(meta[keep], on="sample_id", how="left")
        div.to_csv(div_csv, index=False)
        console.print(f"[green]Diversity CSV -> {div_csv}[/green]")
    else:
        console.print("[yellow]No samples with productive clones — diversity CSV not written.[/yellow]")

    # Optional VDJdb annotation
    if vdjdb_dir.exists() and (list(vdjdb_dir.glob("*.tsv")) + list(vdjdb_dir.glob("*.txt"))):
        try:
            from ..antigen.vdjdb_lookup import VDJdbMatcher, load_vdjdb
            vdjdb = load_vdjdb(vdjdb_dir, chain=chain)
            matcher = VDJdbMatcher(vdjdb)
            ag = matcher.annotate_per_sample(atlas.repertoire, cdr3_col="junction_aa")
            if not ag.empty:
                ag_path = out_dir / "vdjdb_matches.csv"
                ag.to_csv(ag_path, index=False)
                atlas.antigen_profiles = ag
                console.print(f"[green]VDJdb hits -> {ag_path} ({len(ag)} clone-hits)[/green]")
            else:
                console.print("[yellow]No VDJdb hits[/yellow]")
        except Exception as e:
            console.print(f"[yellow]VDJdb lookup skipped: {e}[/yellow]")
    else:
        console.print("[yellow]VDJdb dir empty — skipping antigen lookup[/yellow]")

    atlas_path = out_dir / "thoracic_tcr_atlas.joblib"
    atlas.save(atlas_path)

    summary = atlas.summary()
    summary_path = out_dir / "atlas_summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    console.print("[green]Atlas summary:[/green]")
    for k, v in summary.items():
        console.print(f"  {k}: {v}")
    console.print(f"[green]Summary JSON -> {summary_path}[/green]")


@app.command("info")
def info(
    atlas_path: Path = typer.Option(
        relpath("data/atlas/thoracic_tcr_atlas.joblib"), "--atlas", "-a"
    ),
) -> None:
    """Load a previously saved atlas joblib and print its summary."""
    from ..atlas.atlas_builder import ThoracicTCRAtlas
    if not atlas_path.exists():
        console.print(f"[red]Atlas not found: {atlas_path}[/red]")
        raise typer.Exit(1)
    a = ThoracicTCRAtlas.load(atlas_path)
    for k, v in a.summary().items():
        console.print(f"  {k}: {v}")
