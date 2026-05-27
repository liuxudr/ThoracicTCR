"""``thoracictcr manifest ...`` — build sample manifests for the pipeline.

Subcommands
-----------
* ``tcga``   GDC API → 5 TCGA thoracic projects RNA-Seq BAM manifest
* ``geo``    ENA filereport → GEO/SRA pilot manifest (open-access fastq)
* ``atlas``  Filter + stratified-select runs from ``geo`` output for an atlas
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer
from rich.console import Console

from ._common import ensure_parent, parse_quota, print_summary, relpath
from ..data.geo_manifest import build_geo_pilot_manifest, select_pilot as select_geo_pilot
from ..data.manifest import (
    build_thoracic_manifest,
    select_pilot as select_tcga_pilot,
    write_gdc_manifest,
)

app = typer.Typer(
    name="manifest",
    help="Build TCGA / GEO sample manifests + pilot/atlas subset selection.",
    no_args_is_help=True,
)
console = Console()

DEFAULT_TCGA_PROJECTS = ["TCGA-LUAD", "TCGA-LUSC", "TCGA-ESCA", "TCGA-MESO", "TCGA-THYM"]
DEFAULT_GEO_COHORTS = ["GSE207422", "GSE193258", "GSE126044", "GSE135222", "GSE145370"]
# Cohort selection rule for the atlas (validated empirically — see scripts/06_)
DEFAULT_ATLAS_QUOTA = {"GSE145370": 5, "GSE135222": 12, "GSE126044": 8}
DEFAULT_ATLAS_EXCLUDE = ["GSE193258"]  # bad library prep for TCR


def _write_aria2_urls(df: pd.DataFrame, out_path: Path, dl_dir: str) -> int:
    """Write an aria2c-compatible URL list with one block per fastq mate."""
    n_lines = 0
    with out_path.open("w") as f:
        for _, r in df.iterrows():
            for mate in ("r1", "r2"):
                url = r.get(f"fastq_{mate}_https", "")
                md5 = r.get(f"fastq_{mate}_md5", "")
                if not url:
                    continue
                fname = f"{r['run_accession']}_{mate.upper()}.fastq.gz"
                f.write(f"{url}\n  dir={dl_dir}\n  out={fname}\n")
                if md5:
                    f.write(f"  checksum=md5={md5}\n")
                n_lines += 1
    return n_lines


@app.command("tcga")
def tcga(
    projects: list[str] = typer.Option(
        DEFAULT_TCGA_PROJECTS, "--project", "-p", help="TCGA project ID; repeatable."
    ),
    access: str = typer.Option(
        "open", "--access", help="GDC access tier: 'open' or 'controlled' (dbGaP)."
    ),
    out_dir: Path = typer.Option(
        relpath("configs"), "--out-dir", "-o", help="Where to write manifest CSVs."
    ),
    pilot_project: str = typer.Option(
        "TCGA-THYM", "--pilot-project", help="Project used for the small pilot subset."
    ),
    pilot_n: int = typer.Option(10, "--pilot-n", help="Pilot subset size."),
    seed: int = typer.Option(42, "--seed", help="Deterministic sampling seed."),
    size: int = typer.Option(10000, "--gdc-size", help="GDC API page size."),
) -> None:
    """Query GDC and write ``thoracic_manifest.csv`` + a pilot subset.

    Note: TCGA RNA-Seq BAMs are dbGaP-controlled. For ``--access controlled``
    you must have ``gdc-user-token.txt`` configured.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    df = build_thoracic_manifest(projects=projects, access=access, size=size)
    master_csv = out_dir / "thoracic_manifest.csv"
    df.to_csv(master_csv, index=False)
    console.print(f"[green]Master manifest -> {master_csv}[/green] ({len(df)} rows)")

    if df.empty:
        return

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
    summary.to_csv(out_dir / "thoracic_manifest_summary.csv")
    print_summary(summary, "Per-project case counts and total BAM size")

    pilot = select_tcga_pilot(df, project=pilot_project, n=pilot_n, seed=seed)
    if pilot.empty:
        console.print("[yellow]Pilot selection empty (no rows matched).[/yellow]")
        return
    pilot_csv = out_dir / f"pilot_{pilot_project.lower().replace('-', '_')}_{pilot_n}.csv"
    pilot.to_csv(pilot_csv, index=False)
    gdc_path = out_dir / pilot_csv.name.replace(".csv", ".gdc_manifest.txt")
    write_gdc_manifest(pilot, gdc_path)
    console.print(f"[green]Pilot manifest -> {pilot_csv}[/green]")
    console.print(f"[green]GDC-client manifest -> {gdc_path}[/green]")
    console.print(f"[cyan]Pilot total download size: {pilot['file_size_gb'].sum():.1f} GB[/cyan]")


@app.command("geo")
def geo(
    accessions: list[str] = typer.Option(
        DEFAULT_GEO_COHORTS, "--accession", "-a", help="GSE / SRP / PRJNA accession; repeatable."
    ),
    rnaseq_only: bool = typer.Option(True, "--rnaseq-only/--all-libraries"),
    pilot_n: int = typer.Option(10, "--pilot-n"),
    seed: int = typer.Option(42, "--seed"),
    max_gb_per_sample: float = typer.Option(
        20.0, "--max-gb-per-sample",
        help="Skip runs larger than this when picking the pilot (set 0 to disable).",
    ),
    out_dir: Path = typer.Option(relpath("configs"), "--out-dir", "-o"),
    download_dir: str = typer.Option(
        "data/raw/fastq", "--download-dir",
        help="Relative path written into the aria2 URL list 'dir=' header.",
    ),
) -> None:
    """Build a GEO/SRA pilot manifest via ENA filereport (no SRA toolkit needed).

    Writes:
      * ``geo_thoracic_manifest.csv`` — every queryable run
      * ``pilot_geo_N.csv``           — N smallest paired runs (one per sample)
      * ``pilot_geo_N.urls.txt``      — aria2c input list with MD5 checksums
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    console.print(f"[cyan]Querying ENA for {len(accessions)} accessions: {accessions}[/cyan]")
    df = build_geo_pilot_manifest(accessions, rnaseq_only=rnaseq_only)
    if df.empty:
        console.print("[red]ENA returned no RNA-Seq runs for any accession.[/red]")
        raise typer.Exit(1)

    master_csv = out_dir / "geo_thoracic_manifest.csv"
    df.to_csv(master_csv, index=False)
    console.print(f"[green]Master manifest -> {master_csv}[/green] ({len(df)} runs)")

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
    summary.to_csv(out_dir / "geo_thoracic_manifest_summary.csv")
    print_summary(summary, "Per-cohort summary")

    pilot = select_geo_pilot(
        df,
        n=pilot_n,
        seed=seed,
        prefer_paired=True,
        max_gb_per_sample=None if max_gb_per_sample <= 0 else max_gb_per_sample,
    )
    pilot_csv = out_dir / f"pilot_geo_{pilot_n}.csv"
    pilot.to_csv(pilot_csv, index=False)
    console.print(f"[green]Pilot manifest -> {pilot_csv}[/green]")

    urls_path = out_dir / f"pilot_geo_{pilot_n}.urls.txt"
    n_lines = _write_aria2_urls(pilot, urls_path, dl_dir=download_dir)
    console.print(
        f"[green]aria2 URL list -> {urls_path}[/green] ({n_lines} mates, "
        f"{pilot['total_fastq_gb'].sum():.2f} GB total)"
    )


@app.command("atlas")
def atlas(
    manifest_in: Path = typer.Option(
        relpath("configs/geo_thoracic_manifest.csv"), "--manifest", "-m",
        help="Master manifest CSV produced by 'manifest geo'.",
    ),
    quota: list[str] = typer.Option(
        None, "--quota",
        help="Per-cohort selection quota, e.g. --quota GSE145370=5 (repeatable). "
             "Default mirrors scripts/06_.",
    ),
    exclude: list[str] = typer.Option(
        DEFAULT_ATLAS_EXCLUDE, "--exclude", "-x",
        help="Cohorts to drop (repeatable).",
    ),
    existing_fastq_dir: Path = typer.Option(
        relpath("data/raw/fastq"), "--existing-fastq-dir",
        help="Skip runs already present here when building the download list.",
    ),
    out_dir: Path = typer.Option(relpath("configs"), "--out-dir", "-o"),
    download_dir: str = typer.Option("data/raw/fastq", "--download-dir"),
) -> None:
    """Filter + stratified-sample the GEO manifest to build the atlas cohort.

    Picks the N smallest *paired* runs per cohort (so disk footprint is bounded),
    excludes cohorts with verified bad library prep, and emits an aria2 list of
    only the runs not yet on disk.
    """
    if not manifest_in.exists():
        console.print(f"[red]Manifest not found: {manifest_in}[/red] — run 'manifest geo' first.")
        raise typer.Exit(1)

    df = pd.read_csv(manifest_in)
    console.print(f"[cyan]Loaded {len(df)} runs from {manifest_in}[/cyan]")

    exclude_set = set(exclude)
    df = df[~df["queried_accession"].isin(exclude_set)]
    console.print(f"[cyan]After excluding {exclude_set}: {len(df)} runs[/cyan]")

    quota_map = parse_quota(quota, DEFAULT_ATLAS_QUOTA)

    selected = []
    for cohort, n in quota_map.items():
        sub = df[df["queried_accession"] == cohort].copy()
        sub = sub.dropna(subset=["fastq_r1_https", "fastq_r2_https"])
        sub = sub[sub["library_layout"].fillna("").str.upper() == "PAIRED"]
        sub = sub.sort_values("total_fastq_gb").drop_duplicates("sample_accession", keep="first")
        picked = sub.head(n)
        console.print(f"  {cohort}: picked {len(picked)}/{len(sub)} (target {n})")
        selected.append(picked)

    atlas_df = pd.concat(selected, ignore_index=True) if selected else pd.DataFrame()
    out_dir.mkdir(parents=True, exist_ok=True)
    atlas_csv = out_dir / "atlas_manifest.csv"
    atlas_df.to_csv(atlas_csv, index=False)
    console.print(f"[green]Atlas manifest -> {atlas_csv}[/green] (n={len(atlas_df)})")

    if atlas_df.empty:
        return

    summary = (
        atlas_df.groupby("queried_accession")
        .agg(
            n_runs=("run_accession", "count"),
            total_gb=("total_fastq_gb", "sum"),
            mean_gb=("total_fastq_gb", "mean"),
        )
        .round(2)
    )
    summary.to_csv(out_dir / "atlas_manifest_summary.csv")
    print_summary(summary, "Per-cohort atlas selection")
    console.print(f"[cyan]Atlas total: {atlas_df['total_fastq_gb'].sum():.1f} GB[/cyan]")

    # Compute delta: which fastq pairs are NOT yet on disk
    existing_runs: set[str] = set()
    if existing_fastq_dir.exists():
        for p in existing_fastq_dir.glob("*_R1.fastq.gz"):
            existing_runs.add(p.name.replace("_R1.fastq.gz", ""))
    console.print(f"[cyan]Already on disk: {len(existing_runs)} runs[/cyan]")

    new = atlas_df[~atlas_df["run_accession"].isin(existing_runs)]
    delta_gb = new["total_fastq_gb"].sum()
    urls_path = out_dir / "atlas_download.urls.txt"
    n_lines = _write_aria2_urls(new, urls_path, dl_dir=download_dir)
    console.print(
        f"[green]New download list -> {urls_path}[/green] "
        f"({len(new)} runs, {delta_gb:.1f} GB, {n_lines} mates)"
    )
