"""``thoracictcr trust4 ...`` — TCR repertoire reconstruction with TRUST4.

Two paths:
  * ``run-bam``  Loop over BAM files in pure Python.
  * ``run``      Run TRUST4 over a directory of fastq (or BAM) samples in parallel.

``run`` also handles the BAM (TCGA) mode when invoked with ``--mode tcga``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import typer
from rich.console import Console

from ._common import relpath

app = typer.Typer(
    name="trust4",
    help="TCR repertoire reconstruction via TRUST4.",
    no_args_is_help=True,
)
console = Console()

TRUST4_REFS = {
    "hg38_bcrtcr.fa":  "https://raw.githubusercontent.com/liulab-dfci/TRUST4/master/hg38_bcrtcr.fa",
    "human_IMGT+C.fa": "https://raw.githubusercontent.com/liulab-dfci/TRUST4/master/human_IMGT+C.fa",
}


@app.command("fetch-refs")
def fetch_refs(
    out_dir: Path = typer.Option(
        relpath("data/external/trust4_ref"), "--out-dir", "-o",
        help="Where to place hg38_bcrtcr.fa and human_IMGT+C.fa.",
    ),
    force: bool = typer.Option(False, "--force/--no-force"),
) -> None:
    """Download TRUST4 reference files (~1 MB total) from the official GitHub mirror."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, url in TRUST4_REFS.items():
        dest = out_dir / name
        if dest.exists() and not force:
            console.print(f"[yellow]exists, skipping[/yellow]: {dest}")
            continue
        console.print(f"[cyan]fetching[/cyan]: {url}")
        urllib.request.urlretrieve(url, dest)
        console.print(f"[green]wrote[/green]: {dest} ({dest.stat().st_size/1e6:.2f} MB)")


@app.command("check")
def check() -> None:
    """Print whether ``run-trust4`` and reference files are visible to this shell."""
    trust4_bin = shutil.which("run-trust4")
    console.print(f"run-trust4 binary: {trust4_bin or '[red]NOT FOUND[/red]'}")
    ref_dir = relpath("data/external/trust4_ref")
    for fname in TRUST4_REFS:
        p = ref_dir / fname
        console.print(f"  ref {fname}: {'OK' if p.exists() else '[red]MISSING[/red]'} ({p})")


@app.command("run-bam")
def run_bam(
    bam_dir: Path = typer.Option(..., "--bam-dir", help="Directory of BAMs to process."),
    output_dir: Path = typer.Option(..., "--output", "-o"),
    threads: int = typer.Option(8, "--threads", "-t"),
) -> None:
    """Run TRUST4 on every ``*.bam`` under ``--bam-dir`` (pure Python loop)."""
    from ..repertoire.trust4_runner import TRUST4Runner

    runner = TRUST4Runner(n_threads=threads)
    output_dir.mkdir(parents=True, exist_ok=True)
    bams = sorted(bam_dir.glob("*.bam"))
    if not bams:
        console.print(f"[red]No BAMs under {bam_dir}[/red]")
        raise typer.Exit(1)
    for bam in bams:
        sample_id = bam.stem
        out_sub = output_dir / sample_id
        out_sub.mkdir(exist_ok=True)
        runner.run_on_bam(bam, out_sub)
        console.print(f"[green]done[/green]: {sample_id}")


def _find_repo_root(start: Path) -> Path | None:
    """Walk up from ``start`` looking for a directory that contains both
    ``scripts/03_trust4_pilot.sh`` and ``configs/paths.yaml``."""
    for d in [start, *start.parents]:
        if (d / "scripts" / "03_trust4_pilot.sh").exists() and (d / "configs" / "paths.yaml").exists():
            return d
    return None


@app.command("run")
def run(
    mode: str = typer.Option("geo", "--mode", help="'geo' (fastq) or 'tcga' (BAM)."),
    fastq_dir: Path = typer.Option(None, "--fastq-dir", help="Override $FASTQ_ROOT (mode=geo)."),
    bam_dir: Path = typer.Option(None, "--bam-dir", help="Override $BAM_ROOT (mode=tcga)."),
    n_parallel: int = typer.Option(2, "--n-parallel"),
    threads_per_job: int = typer.Option(4, "--threads-per-job"),
) -> None:
    """Run TRUST4 across a directory of samples (fastq or BAM) in parallel.

    ``--mode geo`` expects paired fastq.gz under ``data/raw/fastq/`` (the
    default; override with ``--fastq-dir``). ``--mode tcga`` expects BAM files
    under ``data/raw/bam/`` (override with ``--bam-dir``).
    """
    repo = _find_repo_root(Path.cwd())
    if repo is None:
        console.print(
            "[red]Cannot locate the TRUST4 driver. "
            "Please install ThoracicTCR from a repository checkout.[/red]"
        )
        raise typer.Exit(1)

    script = repo / "scripts" / "03_trust4_pilot.sh"
    env = os.environ.copy()
    env.update({
        "MODE": mode,
        "N_PARALLEL": str(n_parallel),
        "THREADS_PER_JOB": str(threads_per_job),
    })
    if fastq_dir:
        env["FASTQ_ROOT"] = str(fastq_dir)
    if bam_dir:
        env["BAM_ROOT"] = str(bam_dir)

    console.print(f"[cyan]bash {script}  (MODE={mode}, N_PARALLEL={n_parallel})[/cyan]")
    proc = subprocess.run(["bash", str(script)], env=env, cwd=repo)
    if proc.returncode != 0:
        raise typer.Exit(proc.returncode)
