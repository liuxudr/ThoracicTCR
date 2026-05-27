"""``thoracictcr figures ...`` — render the SCI manuscript figures via Rscript.

All figures are real R scripts under ``R/fig*.R`` rendered with ggplot2 + ggsci
(NPG palette). This sub-app just enumerates them and calls ``Rscript`` so the
pipeline is invokable from the unified CLI.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import typer
from rich.console import Console

from ._common import relpath
from .trust4_app import _find_repo_root  # reuse repo-root locator

app = typer.Typer(
    name="figures",
    help="Render manuscript figures (R/ggplot2 + ggsci NPG palette).",
    no_args_is_help=True,
)
console = Console()

DEFAULT_FIGURES = [
    "fig1_landscape.R", "fig2_diversity_clonality.R", "fig3_vj_usage.R",
    "fig4_clonal_expansion.R", "fig5_antigen.R", "fig6_shared_motifs.R",
    "figS1_pipeline.R", "figS2_cdr3_length.R", "figS3_depth.R",
    "figS4_gse193258.R", "figS5_hill.R", "figS6_rarefaction.R",
    "figS7_vdjdb_detail.R", "figS8_attrition.R",
]


@app.command("list")
def list_figures() -> None:
    """List the figure scripts shipped with the repo."""
    repo = _find_repo_root(Path.cwd())
    r_dir = (repo / "R") if repo else relpath("R")
    if not r_dir.exists():
        console.print(f"[red]R/ dir not found at {r_dir}[/red]")
        raise typer.Exit(1)
    for f in sorted(r_dir.glob("fig*.R")) + sorted(r_dir.glob("tableS*.R")):
        console.print(f"  {f.name}")


@app.command("render")
def render(
    figures: list[str] = typer.Option(
        None, "--figure", "-f",
        help="Specific figure script(s) to render (e.g. fig1_landscape.R). "
             "Default: all main + supplementary figures.",
    ),
    out_dir: Path = typer.Option(
        relpath("paper/figures"), "--out-dir", "-o",
        help="Where R scripts should write PDFs/PNGs (passed as $THORACICTCR_FIG_DIR).",
    ),
) -> None:
    """Render figures via ``Rscript R/fig*.R``.

    Requires R + ggplot2 + ggsci + patchwork. Each R script reads its inputs
    from ``data/metrics/`` and ``data/atlas/`` (set up by upstream commands).
    """
    if not shutil.which("Rscript"):
        console.print("[red]Rscript not found on PATH — install R first.[/red]")
        raise typer.Exit(1)

    repo = _find_repo_root(Path.cwd())
    if repo is None:
        console.print(
            "[red]Cannot locate the figure scripts. "
            "Please install ThoracicTCR from a repository checkout.[/red]"
        )
        raise typer.Exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)
    (repo / "paper" / "tables").mkdir(parents=True, exist_ok=True)

    target = list(figures) if figures else DEFAULT_FIGURES
    env = os.environ.copy()
    env["THORACICTCR_REPO"] = str(repo)
    env["THORACICTCR_FIG_DIR"] = str(out_dir)

    n_ok, n_skip = 0, 0
    for fr in target:
        path = repo / "R" / fr
        if not path.exists():
            console.print(f"[yellow]skip[/yellow]: {fr} not found")
            n_skip += 1
            continue
        console.print(f"=== Rendering {fr} ===")
        rc = subprocess.run(["Rscript", str(path)], env=env, cwd=repo).returncode
        if rc == 0:
            n_ok += 1
        else:
            console.print(f"[red]failed[/red]: {fr} (rc={rc})")

    console.print(f"[green]Rendered {n_ok}/{len(target)} figure scripts ({n_skip} skipped)[/green]")
    pdfs = sorted(out_dir.glob("*.pdf"))
    if pdfs:
        console.print(f"[cyan]PDFs in {out_dir}:[/cyan]")
        for p in pdfs:
            console.print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")
