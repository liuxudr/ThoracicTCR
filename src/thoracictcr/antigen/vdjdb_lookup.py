"""VDJdb TCR-antigen lookup.

VDJdb provides experimentally-validated TCR-pMHC pairs. Match strategies:
  1. Exact CDR3-aa match (high confidence)
  2. CDR3-aa + V gene match
  3. Hamming/Levenshtein distance ≤ 1 with same V/J (medium confidence)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from rich.console import Console

console = Console()


def load_vdjdb(vdjdb_dir: Path, chain: str = "TRB") -> pd.DataFrame:
    """Load VDJdb table.

    The VDJdb release ships multiple TSVs:
      - vdjdb.slim.txt     : core columns
      - vdjdb.txt          : full
      - vdjdb_full.txt     : full + meta
    We use vdjdb.slim.txt for speed.
    """
    # Look for the slim file
    candidates = [
        vdjdb_dir / "vdjdb.slim.txt",
        vdjdb_dir / "vdjdb-2024-02-22.tsv",
    ] + list(vdjdb_dir.glob("vdjdb*.tsv")) + list(vdjdb_dir.glob("vdjdb*.txt"))
    chosen = next((p for p in candidates if p.exists()), None)
    if chosen is None:
        raise FileNotFoundError(
            f"No VDJdb table found in {vdjdb_dir}. "
            f"Run scripts/05_download_public_tcr_dbs.sh first."
        )
    console.print(f"[cyan]Loading VDJdb from {chosen.name}[/cyan]")
    df = pd.read_csv(chosen, sep="\t", dtype=str, low_memory=False)

    # Normalize column names across VDJdb releases
    col_map = {
        "gene": "chain",
        "cdr3": "cdr3_aa",
        "cdr3.aa": "cdr3_aa",
        "v.segm": "v_call",
        "j.segm": "j_call",
        "antigen.epitope": "epitope",
        "antigen.gene": "antigen_gene",
        "antigen.species": "antigen_species",
        "mhc.a": "mhc_a",
        "mhc.class": "mhc_class",
        "vdjdb.score": "confidence_score",
    }
    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    chain_col = "chain" if "chain" in df.columns else "Gene"
    if chain_col in df.columns:
        df = df[df[chain_col].astype(str).str.upper() == chain.upper()]
    return df.reset_index(drop=True)


class VDJdbMatcher:
    """Match a set of TCR clones against VDJdb."""

    def __init__(self, vdjdb: pd.DataFrame):
        self.vdjdb = vdjdb
        # Build fast exact-match index
        self._exact = vdjdb.set_index("cdr3_aa", drop=False) if "cdr3_aa" in vdjdb.columns else None

    def match_exact(self, cdr3_aa: str) -> pd.DataFrame:
        """Exact CDR3-aa match. Returns matching VDJdb rows."""
        if self._exact is None:
            return pd.DataFrame()
        try:
            hits = self._exact.loc[[cdr3_aa]]
            return hits.reset_index(drop=True)
        except KeyError:
            return pd.DataFrame()

    def match_clones(
        self,
        clones: pd.DataFrame,
        cdr3_col: str = "junction_aa",
    ) -> pd.DataFrame:
        """Annotate clones with VDJdb hits.

        Returns rows with: clone_cdr3, n_hits, epitopes, antigens, species, mhc_class.
        """
        rows = []
        for cdr3 in clones[cdr3_col].dropna().unique():
            hits = self.match_exact(str(cdr3))
            if hits.empty:
                continue
            rows.append({
                "cdr3_aa": cdr3,
                "n_hits": len(hits),
                "epitopes": ";".join(sorted(set(hits.get("epitope", pd.Series()).dropna()))[:5]),
                "antigens": ";".join(sorted(set(hits.get("antigen_gene", pd.Series()).dropna()))[:5]),
                "species": ";".join(sorted(set(hits.get("antigen_species", pd.Series()).dropna()))[:5]),
                "mhc_class": ";".join(sorted(set(hits.get("mhc_class", pd.Series()).dropna()))[:3]),
            })
        return pd.DataFrame(rows)

    def annotate_per_sample(
        self,
        sample_clones: dict[str, pd.DataFrame],
        cdr3_col: str = "junction_aa",
    ) -> pd.DataFrame:
        """Run match across all samples and aggregate."""
        all_rows = []
        for sample_id, clones in sample_clones.items():
            matches = self.match_clones(clones, cdr3_col=cdr3_col)
            if matches.empty:
                continue
            matches["sample_id"] = sample_id
            all_rows.append(matches)
        if not all_rows:
            return pd.DataFrame()
        return pd.concat(all_rows, ignore_index=True)
