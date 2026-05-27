"""Shared helpers for thoracictcr CLI sub-apps.

These keep the per-command code small by centralising:
  * working-directory-relative defaults (``relpath``)
  * tabular summary printing
  * a small "quota" parser used by ``manifest atlas``
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from rich.console import Console
from rich.table import Table

console = Console()


def relpath(rel: str) -> Path:
    """Resolve ``rel`` relative to the current working directory.

    The CLI is meant to be invoked from the repo root, so ``configs/``,
    ``data/``, ``paper/`` etc. all resolve naturally.
    """
    return Path.cwd() / rel


def ensure_parent(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def print_summary(df: pd.DataFrame, title: str = "") -> None:
    """Pretty-print a small pandas summary frame to the console."""
    if title:
        console.print(f"\n[bold cyan]{title}[/bold cyan]")
    if df.empty:
        console.print("[yellow](empty)[/yellow]")
        return
    table = Table(show_lines=False)
    table.add_column(df.index.name or "", justify="left")
    for col in df.columns:
        table.add_column(str(col), justify="right")
    for idx, row in df.iterrows():
        table.add_row(str(idx), *[_fmt(v) for v in row.tolist()])
    console.print(table)


def _fmt(v: object) -> str:
    if isinstance(v, float):
        return f"{v:.2f}"
    return str(v)


def parse_quota(spec: list[str] | None, default: dict[str, int]) -> dict[str, int]:
    """Parse ``--quota GSE145370=5 --quota GSE135222=12`` -> dict.

    If ``spec`` is empty/None, return ``default`` (a copy).
    """
    if not spec:
        return dict(default)
    out: dict[str, int] = {}
    for item in spec:
        if "=" not in item:
            raise ValueError(f"Bad --quota '{item}', expected COHORT=N")
        k, v = item.split("=", 1)
        out[k.strip()] = int(v)
    return out
