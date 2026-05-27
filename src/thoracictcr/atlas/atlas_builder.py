"""Pan-thoracic TCR atlas builder.

Aggregates per-sample TRUST4 outputs + diversity metrics + HLA + antigen annotations
into a single `ThoracicTCRAtlas` object serializable via joblib.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import joblib
import pandas as pd
from rich.console import Console

from ..metrics.diversity import compute_diversity
from ..repertoire.repertoire_io import load_trust4_clones

console = Console()


@dataclass
class ThoracicTCRAtlas:
    """Unified atlas data structure.

    Attributes:
        metadata: per-sample DataFrame (sample_id, cohort, cancer, treatment, ...)
        repertoire: dict[sample_id, DataFrame] of clones
        diversity: per-sample diversity metrics
        hla: dict[sample_id, dict[gene, list[allele]]]
        antigen_profiles: per-sample antigen annotation hits (VDJdb / McPAS)
        version: schema version string
    """

    metadata: pd.DataFrame = field(default_factory=pd.DataFrame)
    repertoire: dict[str, pd.DataFrame] = field(default_factory=dict)
    diversity: pd.DataFrame = field(default_factory=pd.DataFrame)
    hla: dict[str, dict] = field(default_factory=dict)
    antigen_profiles: pd.DataFrame = field(default_factory=pd.DataFrame)
    version: str = "0.1.0"

    def save(self, path: Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, path)
        console.print(f"[green]Atlas saved -> {path}[/green]")

    @classmethod
    def load(cls, path: Path) -> "ThoracicTCRAtlas":
        return joblib.load(path)

    def summary(self) -> dict:
        return {
            "n_samples": len(self.metadata),
            "n_cohorts": int(self.metadata.get("cohort", pd.Series()).nunique())
                if not self.metadata.empty else 0,
            "n_repertoires": len(self.repertoire),
            "total_clones": int(sum(len(c) for c in self.repertoire.values())),
            "n_with_hla": len(self.hla),
            "n_antigen_hits": len(self.antigen_profiles),
            "version": self.version,
        }


def build_atlas(
    trust4_root: Path,
    metadata: pd.DataFrame,
    hla_root: Path | None = None,
    vdjdb_matches_csv: Path | None = None,
    chain: str = "TRB",
) -> ThoracicTCRAtlas:
    """Assemble atlas from per-sample TRUST4 outputs + optional HLA + antigen.

    Args:
        trust4_root: dir containing <sample_id>/<sample_id>_report.tsv
        metadata: per-sample metadata DataFrame with at least 'sample_id' column
        hla_root: dir containing <sample_id>/*.genotype.json (optional)
        vdjdb_matches_csv: pre-computed VDJdb match CSV (optional)
        chain: TCR chain to extract (default TRB)
    """
    atlas = ThoracicTCRAtlas()

    sample_ids = metadata["sample_id"].tolist()
    console.print(f"[cyan]Building atlas across {len(sample_ids)} samples[/cyan]")

    diversity_rows = []
    for sid in sample_ids:
        sample_dir = trust4_root / sid
        clones = load_trust4_clones(sample_dir, chain=chain, productive_only=True)
        if clones is None or clones.empty:
            console.print(f"[yellow]  skip {sid}: no productive {chain} clones[/yellow]")
            continue
        atlas.repertoire[sid] = clones
        m = compute_diversity(clones)
        m["sample_id"] = sid
        diversity_rows.append(m)

    atlas.diversity = pd.DataFrame(diversity_rows)
    atlas.metadata = metadata.copy()

    # HLA (optional)
    if hla_root and hla_root.exists():
        for sid in sample_ids:
            for gj in (hla_root / sid).glob("*genotype.json") if (hla_root / sid).exists() else []:
                with gj.open() as f:
                    atlas.hla[sid] = json.load(f)
                break

    # Antigen annotation (optional)
    if vdjdb_matches_csv and vdjdb_matches_csv.exists():
        atlas.antigen_profiles = pd.read_csv(vdjdb_matches_csv)

    return atlas
