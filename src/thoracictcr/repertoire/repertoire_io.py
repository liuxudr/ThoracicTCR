"""Standardized I/O for TRUST4/MiXCR repertoire outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

TRUST4_REPORT_COLS = [
    "count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "C", "cid", "cid_full_length",
]


def load_trust4_report(report_tsv: Path) -> pd.DataFrame:
    """Load TRUST4 *_report.tsv to a normalized DataFrame.

    Returns columns:
        chain, v_call, j_call, junction_aa, junction_nt, duplicate_count, productive
    """
    df = pd.read_csv(report_tsv, sep="\t")
    # TRUST4 column naming varies between versions; normalize defensively
    rename_map = {
        "#count": "count",
        "frequency": "frequency",
        "CDR3nt": "junction_nt",
        "CDR3aa": "junction_aa",
        "V": "v_call",
        "J": "j_call",
        "C": "c_call",
    }
    df = df.rename(columns=rename_map)

    out = pd.DataFrame()
    out["v_call"] = df.get("v_call", "")
    out["j_call"] = df.get("j_call", "")
    out["junction_aa"] = df.get("junction_aa", "")
    out["junction_nt"] = df.get("junction_nt", "")
    out["duplicate_count"] = df.get("count", 0).astype(int)
    # Infer chain from V/J gene prefix (TRA/TRB/TRG/TRD)
    out["chain"] = (
        out["v_call"].astype(str).str.extract(r"^(TR[ABGD])")[0]
        .fillna(out["j_call"].astype(str).str.extract(r"^(TR[ABGD])")[0])
    )
    # Productive = no stop codon AND not partial
    out["productive"] = (
        out["junction_aa"].astype(str).str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
        & ~out["junction_aa"].str.contains(r"\*", na=False)
        & ~out["junction_aa"].str.contains("partial", case=False, na=False)
    )
    return out


def load_trust4_clones(
    sample_dir: Path,
    chain: str = "TRB",
    productive_only: bool = True,
    min_count: int = 1,
) -> pd.DataFrame | None:
    """Load and filter clones from a TRUST4 output directory.

    Args:
        sample_dir: Directory containing TRUST4 *_report.tsv
        chain: TRA / TRB / TRG / TRD
        productive_only: filter to productive sequences
        min_count: minimum duplicate_count

    Returns:
        DataFrame with one row per unique CDR3aa (counts aggregated), or None if no report.
    """
    reports = list(sample_dir.glob("*_report.tsv"))
    if not reports:
        return None
    df = load_trust4_report(reports[0])
    df = df[df["chain"] == chain]
    if productive_only:
        df = df[df["productive"]]
    df = df[df["duplicate_count"] >= min_count]
    # Aggregate by CDR3aa (same protein from different nt sequences)
    df = (
        df.groupby("junction_aa", as_index=False)
        .agg({"duplicate_count": "sum", "v_call": "first", "j_call": "first"})
    )
    df = df[df["junction_aa"].str.len().between(8, 24)]
    df = df.sort_values("duplicate_count", ascending=False).reset_index(drop=True)
    return df
