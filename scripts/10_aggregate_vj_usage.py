#!/usr/bin/env python
"""Aggregate TRBV and TRBJ gene usage across all productive samples in the atlas.

Reads each TRUST4 *_report.tsv, filters to productive TRB clones via
`repertoire_io.load_trust4_clones`, counts V-gene and J-gene occurrences
weighted by `duplicate_count`, and writes a long-format table.

Output: data/metrics/vj_usage_long.csv
Columns: sample_id, cohort, gene_type (V|J), gene_name, count, frequency
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from thoracictcr.repertoire.repertoire_io import load_trust4_clones  # noqa: E402

TRUST4_DIR = REPO / "data" / "repertoire" / "trust4" / "geo_pilot"
DIVERSITY_CSV = REPO / "data" / "metrics" / "atlas_diversity.csv"
OUT_CSV = REPO / "data" / "metrics" / "vj_usage_long.csv"


def strip_allele(gene_call: str) -> str:
    """TRBV5-1*01 -> TRBV5-1; keep first if multiple comma-separated."""
    if not isinstance(gene_call, str) or not gene_call or gene_call == ".":
        return ""
    first = gene_call.split(",")[0]
    return first.split("*")[0]


def main() -> int:
    div = pd.read_csv(DIVERSITY_CSV)
    sample_cohort = dict(zip(div.sample_id, div.cohort))

    rows = []
    for sample_id, cohort in sample_cohort.items():
        sdir = TRUST4_DIR / sample_id
        if not sdir.exists():
            continue
        clones = load_trust4_clones(sdir, chain="TRB", productive_only=True)
        if clones is None or len(clones) == 0:
            continue
        clones["v_gene"] = clones["v_call"].map(strip_allele)
        clones["j_gene"] = clones["j_call"].map(strip_allele)
        total_count = clones["duplicate_count"].sum()

        # V genes
        v_tab = (clones.groupby("v_gene")["duplicate_count"].sum()
                 .reset_index().rename(columns={"v_gene": "gene_name",
                                                "duplicate_count": "count"}))
        v_tab = v_tab[v_tab["gene_name"] != ""]
        v_tab["gene_type"] = "V"
        v_tab["frequency"] = v_tab["count"] / total_count

        # J genes
        j_tab = (clones.groupby("j_gene")["duplicate_count"].sum()
                 .reset_index().rename(columns={"j_gene": "gene_name",
                                                "duplicate_count": "count"}))
        j_tab = j_tab[j_tab["gene_name"] != ""]
        j_tab["gene_type"] = "J"
        j_tab["frequency"] = j_tab["count"] / total_count

        for sub in (v_tab, j_tab):
            sub["sample_id"] = sample_id
            sub["cohort"] = cohort
            rows.append(sub[["sample_id", "cohort", "gene_type", "gene_name", "count", "frequency"]])

    if not rows:
        print("No data produced", file=sys.stderr)
        return 1
    out = pd.concat(rows, ignore_index=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV} ({len(out)} rows; "
          f"{out.sample_id.nunique()} samples, "
          f"{out[out.gene_type=='V'].gene_name.nunique()} V-genes, "
          f"{out[out.gene_type=='J'].gene_name.nunique()} J-genes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
