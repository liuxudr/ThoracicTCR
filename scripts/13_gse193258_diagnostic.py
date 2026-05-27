#!/usr/bin/env python
"""GSE193258 exclusion diagnostic.

Parses cdr3.out files for a problematic sample (SRR17491517) vs a good sample
(SRR11094249) and quantifies V/J alignment completeness.

cdr3.out columns (no header; TRUST4 internal format):
  0:assemble_id  1:row  2:v_call  3:d_call  4:j_call  5:c_call
  6:cdr1_nt      7:cdr2_nt  8:cdr3_nt  9:?  10:?  11:?  12:?

Outputs:
  - data/metrics/gse193258_diagnostic_summary.csv
  - data/metrics/gse193258_diagnostic_per_assembly.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
TRUST4_DIR = REPO / "data" / "repertoire" / "trust4" / "geo_pilot"
OUT_SUM = REPO / "data" / "metrics" / "gse193258_diagnostic_summary.csv"
OUT_PER = REPO / "data" / "metrics" / "gse193258_diagnostic_per_assembly.csv"

SAMPLES = [
    ("SRR17491517", "GSE193258", "Problematic (excluded)"),
    ("SRR17491496", "GSE193258", "Problematic (excluded)"),
    ("SRR11094249", "GSE145370", "Good cohort"),
    ("SRR8526721", "GSE126044", "Good cohort"),
]


def parse_cdr3_out(path: Path) -> pd.DataFrame:
    rows = []
    if not path.exists():
        return pd.DataFrame()
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            v_call = parts[2]
            j_call = parts[4]
            cdr3_nt = parts[8] if len(parts) > 8 else ""
            has_v = v_call not in ("", "*", ".")
            has_j = j_call not in ("", "*", ".")
            if has_v and has_j:
                completeness = "V+J"
            elif has_v:
                completeness = "V-only"
            elif has_j:
                completeness = "J-only"
            else:
                completeness = "none"
            chain = ""
            for prefix in ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"):
                if v_call.startswith(prefix) or j_call.startswith(prefix):
                    chain = prefix
                    break
            rows.append({
                "assemble_id": parts[0],
                "v_call": v_call,
                "j_call": j_call,
                "completeness": completeness,
                "chain": chain,
                "cdr3_nt_len": len(cdr3_nt) if cdr3_nt and cdr3_nt != "*" else 0,
            })
    return pd.DataFrame(rows)


def main() -> int:
    per_rows = []
    sum_rows = []
    for sample_id, cohort, label in SAMPLES:
        path = TRUST4_DIR / sample_id / f"{sample_id}_cdr3.out"
        df = parse_cdr3_out(path)
        if df.empty:
            print(f"WARN: no data for {sample_id}", file=sys.stderr)
            continue
        df["sample_id"] = sample_id
        df["cohort"] = cohort
        df["label"] = label
        per_rows.append(df)

        n_total = len(df)
        counts = df["completeness"].value_counts().to_dict()
        sum_rows.append({
            "sample_id": sample_id, "cohort": cohort, "label": label,
            "n_assemblies": n_total,
            "n_V_plus_J": counts.get("V+J", 0),
            "n_V_only": counts.get("V-only", 0),
            "n_J_only": counts.get("J-only", 0),
            "n_none": counts.get("none", 0),
            "pct_V_plus_J": 100 * counts.get("V+J", 0) / n_total if n_total else 0.0,
            "median_cdr3_nt_len": float(df["cdr3_nt_len"].median()),
            "mean_cdr3_nt_len": float(df["cdr3_nt_len"].mean()),
        })

    if not per_rows:
        return 1
    per = pd.concat(per_rows, ignore_index=True)
    summary = pd.DataFrame(sum_rows)
    OUT_SUM.parent.mkdir(parents=True, exist_ok=True)
    per.to_csv(OUT_PER, index=False)
    summary.to_csv(OUT_SUM, index=False)
    print(f"Wrote {OUT_SUM}\n{summary.to_string(index=False)}")
    print(f"Wrote {OUT_PER} ({len(per)} assemblies)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
