"""``thoracictcr download ...`` — fetch public datasets used by the pipeline.

Subcommands
-----------
* ``fastq``    Pilot/atlas fastq from the ENA HTTP mirror (aria2c if available).
* ``tcr-dbs``  Public TCR-antigen DBs (VDJdb / McPAS-TCR / IEDB).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import typer
from rich.console import Console

from ._common import relpath
from .trust4_app import _find_repo_root

app = typer.Typer(
    name="download",
    help="Fetch public datasets (fastq, VDJdb / McPAS / IEDB).",
    no_args_is_help=True,
)
console = Console()


@app.command("fastq")
def fastq(
    urls: Path = typer.Option(
        relpath("configs/pilot_geo_10.urls.txt"), "--urls", "-u",
        help="aria2c-compatible URL list (produced by `manifest geo` / `manifest atlas`).",
    ),
    out_dir: Path = typer.Option(
        relpath("data/raw/fastq"), "--out-dir", "-o",
        help="Where downloaded fastq.gz files are written.",
    ),
    n_parallel: int = typer.Option(4, "--n-parallel", "-j"),
    n_conn_per_host: int = typer.Option(4, "--n-conn-per-host", "-x"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print plan and exit."),
    skip_integrity: bool = typer.Option(
        False, "--skip-integrity",
        help="Skip the post-download gzip integrity check.",
    ),
) -> None:
    """Download fastq from ENA listed in a URL manifest (aria2c parallel, resumable).

    Falls back to ``wget`` if ``aria2c`` is not on PATH. After download, runs
    ``gzip -t`` on every file and aborts if any are corrupt.
    """
    if not urls.exists():
        console.print(
            f"[red]URL list not found: {urls}[/red] — "
            "run `thoracictcr manifest geo` (or `manifest atlas`) first."
        )
        raise typer.Exit(1)

    n_files = sum(1 for line in urls.read_text().splitlines() if line.startswith("https"))
    console.print(f"[cyan]URL list:[/cyan] {urls} ({n_files} files)")
    console.print(f"[cyan]Output:[/cyan]   {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    if dry_run:
        console.print("[yellow]DRY-RUN — first 4 URLs:[/yellow]")
        for i, line in enumerate(urls.read_text().splitlines()):
            if line.startswith("https") and i < 8:
                console.print(f"   {line}")
        return

    aria2c = shutil.which("aria2c")
    if aria2c:
        console.print(f"[cyan]Using aria2c (j={n_parallel}, x={n_conn_per_host}, check-integrity=true)[/cyan]")
        rc = subprocess.run([
            aria2c, "-i", str(urls),
            "-j", str(n_parallel), "-x", str(n_conn_per_host), "-s", str(n_conn_per_host),
            "-c", "--check-integrity=true",
            "--auto-file-renaming=false",
            "--console-log-level=warn", "--summary-interval=10",
        ]).returncode
        if rc != 0:
            raise typer.Exit(rc)
    else:
        console.print("[yellow]aria2c not found; falling back to wget (sequential).[/yellow]")
        for line in urls.read_text().splitlines():
            if line.startswith("https"):
                subprocess.run(["wget", "-c", "-P", str(out_dir), line], check=False)

    if not skip_integrity:
        corrupt = 0
        for f in sorted(out_dir.glob("*.fastq.gz")):
            if subprocess.run(["gzip", "-t", str(f)], capture_output=True).returncode != 0:
                console.print(f"[red]CORRUPT[/red]: {f.name}")
                corrupt += 1
        if corrupt:
            console.print(f"[red]{corrupt} corrupt files — delete them and re-run.[/red]")
            raise typer.Exit(1)
        console.print("[green]All gzip integrity OK[/green]")


@app.command("tcr-dbs")
def tcr_dbs(
    out_root: Path = typer.Option(
        relpath("data/external"), "--out-root",
        help="Output root; per-database subdirs are created inside.",
    ),
) -> None:
    """Download VDJdb (CC-BY), McPAS-TCR (academic) and IEDB receptor subset (~30 MB total)."""
    repo = _find_repo_root(Path.cwd())
    if repo is None:
        console.print(
            "[red]Cannot find the download driver — please install ThoracicTCR from a repo checkout.[/red]"
        )
        raise typer.Exit(1)
    script = repo / "scripts" / "05_download_public_tcr_dbs.sh"
    env = os.environ.copy()
    env["EXT"] = str(out_root)
    rc = subprocess.run(["bash", str(script)], env=env, cwd=repo).returncode
    if rc != 0:
        raise typer.Exit(rc)
