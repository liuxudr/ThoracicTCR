"""McPAS-TCR lookup. McPAS holds TCRs implicated in autoimmunity / pathogen / cancer.

Schema columns (as of 2024 release):
  CDR3.alpha.aa, CDR3.beta.aa, TRBV, TRBJ, Pathology, Antigen.protein, Epitope.peptide,
  Species, Category, MHC, T.Cell.Type, ...
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from rich.console import Console

console = Console()


def load_mcpas(mcpas_dir: Path) -> pd.DataFrame:
    """Load McPAS-TCR table."""
    csv_paths = list(mcpas_dir.glob("*.csv"))
    if not csv_paths:
        raise FileNotFoundError(
            f"No McPAS CSV found in {mcpas_dir}. "
            f"Run scripts/05_download_public_tcr_dbs.sh first."
        )
    console.print(f"[cyan]Loading McPAS from {csv_paths[0].name}[/cyan]")
    df = pd.read_csv(csv_paths[0], encoding="latin-1", dtype=str, low_memory=False)
    # Normalize key columns
    rename = {
        "CDR3.beta.aa": "cdr3_beta_aa",
        "CDR3.alpha.aa": "cdr3_alpha_aa",
        "TRBV": "v_beta",
        "TRBJ": "j_beta",
        "Pathology": "pathology",
        "Antigen.protein": "antigen",
        "Epitope.peptide": "epitope",
        "Species": "species",
        "Category": "category",
        "MHC": "mhc",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    return df.reset_index(drop=True)


class McPASMatcher:
    """Match TCRs against McPAS-TCR DB; categorize hits as autoimmune / viral / tumor."""

    def __init__(self, mcpas: pd.DataFrame):
        self.mcpas = mcpas
        # Build CDR3-beta exact-match index
        self._beta = None
        if "cdr3_beta_aa" in mcpas.columns:
            beta = mcpas.dropna(subset=["cdr3_beta_aa"])
            self._beta = beta.set_index("cdr3_beta_aa", drop=False)

    def match_beta(self, cdr3_aa: str) -> pd.DataFrame:
        if self._beta is None:
            return pd.DataFrame()
        try:
            hits = self._beta.loc[[cdr3_aa]]
            return hits.reset_index(drop=True)
        except KeyError:
            return pd.DataFrame()

    def match_clones(
        self,
        clones: pd.DataFrame,
        cdr3_col: str = "junction_aa",
    ) -> pd.DataFrame:
        rows = []
        for cdr3 in clones[cdr3_col].dropna().unique():
            hits = self.match_beta(str(cdr3))
            if hits.empty:
                continue
            categories = hits.get("category", pd.Series()).dropna().unique().tolist()
            pathologies = hits.get("pathology", pd.Series()).dropna().unique().tolist()
            rows.append({
                "cdr3_aa": cdr3,
                "n_hits": len(hits),
                "categories": ";".join(sorted(categories)[:5]),
                "pathologies": ";".join(sorted(pathologies)[:5]),
                "is_autoimmune": any("autoimmune" in str(c).lower() for c in categories),
                "is_pathogen":   any("pathogen"   in str(c).lower() for c in categories),
                "is_cancer":     any("cancer"     in str(c).lower() or "tumor" in str(c).lower()
                                     for c in categories),
            })
        return pd.DataFrame(rows)

    def annotate_per_sample(
        self,
        sample_clones: dict[str, pd.DataFrame],
        cdr3_col: str = "junction_aa",
    ) -> pd.DataFrame:
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
