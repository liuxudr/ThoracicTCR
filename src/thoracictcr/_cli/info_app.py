"""``thoracictcr info ...`` — environment, paths and software-version reporting."""

from __future__ import annotations

import os
import platform
import subprocess
import sys
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from .. import __version__
from ..utils.paths_loader import load_paths
from .trust4_app import _find_repo_root

app = typer.Typer(
    name="info",
    help="Environment + paths + tool version reporting.",
    no_args_is_help=True,
)
console = Console()


@app.command("env")
def env() -> None:
    """Print Python / OS / key-package versions."""
    table = Table(title="thoracictcr environment")
    table.add_column("item"); table.add_column("value")
    table.add_row("thoracictcr", __version__)
    table.add_row("python", sys.version.split()[0])
    table.add_row("executable", sys.executable)
    table.add_row("platform", platform.platform())
    for pkg in ("pandas", "numpy", "scipy", "scikit-learn", "lifelines", "lightgbm",
                "biopython", "requests", "typer"):
        try:
            mod = __import__(pkg.replace("-", "_"))
            ver = getattr(mod, "__version__", "?")
        except ImportError:
            ver = "[red]not installed[/red]"
        table.add_row(f"py:{pkg}", ver)
    console.print(table)


@app.command("paths")
def paths(
    config: Path = typer.Option(None, "--config", help="Override configs/paths.yaml location."),
) -> None:
    """Print fully-resolved paths from configs/paths.yaml."""
    cfg = load_paths(str(config) if config else None)

    def walk(node, prefix=""):
        if isinstance(node, dict):
            for k, v in node.items():
                walk(v, f"{prefix}.{k}" if prefix else k)
        else:
            console.print(f"  [cyan]{prefix}[/cyan] = {node}")
    walk(cfg)


@app.command("versions")
def versions(
    out: Path = typer.Option(
        None, "--out", "-o",
        help="If given, also write a TSV in the format used by tableS5_software_versions.tsv.",
    ),
) -> None:
    """Probe each external tool (TRUST4 / arcasHLA / samtools / aria2c / R / Rscript) for its version."""
    rows: list[tuple[str, str]] = []

    def probe(cmd: list[str]) -> str:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            return (r.stdout or r.stderr).splitlines()[0] if (r.stdout or r.stderr) else "(no output)"
        except FileNotFoundError:
            return "[red]not found[/red]"
        except Exception as e:
            return f"[red]error: {e}[/red]"

    rows.append(("TRUST4",     probe(["run-trust4", "-h"])))
    rows.append(("arcasHLA",   probe(["arcasHLA", "version"])))
    rows.append(("samtools",   probe(["samtools", "--version"])))
    rows.append(("aria2c",     probe(["aria2c", "--version"])))
    rows.append(("Rscript",    probe(["Rscript", "--version"])))
    rows.append(("gdc-client", probe(["gdc-client", "--version"])))

    table = Table(title="External tool versions")
    table.add_column("tool"); table.add_column("first-line")
    for k, v in rows:
        table.add_row(k, v)
    console.print(table)

    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w") as f:
            f.write("tool\tversion\tsource\n")
            for k, v in rows:
                # strip rich markup before writing
                v_plain = v.replace("[red]", "").replace("[/red]", "")
                f.write(f"{k}\t{v_plain}\tprobed\n")
        console.print(f"[green]wrote[/green] {out}")
