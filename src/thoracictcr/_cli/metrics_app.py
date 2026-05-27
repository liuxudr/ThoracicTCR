"""``thoracictcr metrics ...`` — diversity, clonality and motif metrics from TRUST4 outputs.

Each subcommand reads a directory of per-sample TRUST4 outputs and writes CSV
tables under ``data/metrics/`` (or wherever ``--out`` points).
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import typer
from rich.console import Console
from rich.table import Table

from ._common import ensure_parent, relpath
from ..metrics.diversity import compute_diversity
from ..repertoire.repertoire_io import load_trust4_clones, load_trust4_report

app = typer.Typer(
    name="metrics",
    help="Per-sample / cross-cohort repertoire metrics computed from TRUST4 outputs.",
    no_args_is_help=True,
)
console = Console()


def _iter_sample_dirs(trust4_dir: Path) -> list[Path]:
    if not trust4_dir.exists():
        console.print(f"[red]TRUST4 dir not found: {trust4_dir}[/red]")
        raise typer.Exit(1)
    return sorted([d for d in trust4_dir.iterdir() if d.is_dir()])


def _load_sample_cohort(diversity_csv: Path) -> dict[str, str]:
    if not diversity_csv.exists():
        console.print(
            f"[red]diversity CSV not found: {diversity_csv}[/red] — "
            "run 'metrics diversity' or 'atlas build' first."
        )
        raise typer.Exit(1)
    div = pd.read_csv(diversity_csv)
    if "cohort" not in div.columns:
        div["cohort"] = ""
    return dict(zip(div.sample_id.astype(str), div.cohort.astype(str)))


@app.command("diversity")
def diversity(
    trust4_dir: Path = typer.Option(
        relpath("data/repertoire/trust4/geo_pilot"), "--trust4-dir", "-t",
        help="Directory containing per-sample TRUST4 subdirs.",
    ),
    out: Path = typer.Option(
        relpath("data/metrics/diversity.csv"), "--out", "-o",
        help="Output CSV.",
    ),
    chain: str = typer.Option("TRB", "--chain", help="TRA / TRB / TRG / TRD"),
    productive_only: bool = typer.Option(True, "--productive-only/--include-nonproductive"),
    min_count: int = typer.Option(1, "--min-count"),
) -> None:
    """Compute Shannon/Simpson/Hill/Gini/clonality + top-N for every sample."""
    rows = []
    for sample_dir in _iter_sample_dirs(trust4_dir):
        clones = load_trust4_clones(
            sample_dir, chain=chain, productive_only=productive_only, min_count=min_count
        )
        if clones is None or clones.empty:
            console.print(f"[yellow]skip {sample_dir.name}: empty[/yellow]")
            continue
        metrics = compute_diversity(clones)
        metrics["sample_id"] = sample_dir.name
        rows.append(metrics)

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df[["sample_id"] + [c for c in df.columns if c != "sample_id"]]
    df.to_csv(ensure_parent(out), index=False)
    console.print(f"[green]Wrote diversity for {len(df)} samples -> {out}[/green]")


@app.command("pilot-qc")
def pilot_qc(
    trust4_dir: Path = typer.Option(
        relpath("data/repertoire/trust4/geo_pilot"), "--trust4-dir", "-t"
    ),
    out_dir: Path = typer.Option(relpath("data/metrics"), "--out-dir", "-o"),
    min_productive_trb: int = typer.Option(50, "--min-productive-trb"),
) -> None:
    """End-to-end pilot QC: diversity + warnings (low TRB, bad CDR3, missing V…)."""
    sample_dirs = _iter_sample_dirs(trust4_dir)
    console.print(f"[cyan]Found {len(sample_dirs)} pilot sample(s) under {trust4_dir}[/cyan]")

    rows: list[dict] = []
    all_warnings: list[str] = []

    for sd in sample_dirs:
        reports = list(sd.glob("*_report.tsv"))
        if not reports:
            all_warnings.append(f"NO_REPORT: {sd.name}")
            continue
        raw = load_trust4_report(reports[0])
        trb = raw[(raw["chain"] == "TRB") & raw["productive"]]
        if len(trb) < min_productive_trb:
            all_warnings.append(
                f"LOW_TRB: {sd.name} has only {len(trb)} productive TRB clones"
            )
        bad_len = trb[~trb["junction_aa"].str.len().between(8, 24)]
        if len(bad_len):
            all_warnings.append(
                f"BAD_CDR3_LEN: {sd.name} {len(bad_len)} clones outside 8-24aa"
            )
        missing_v = trb["v_call"].isna().sum() + (trb["v_call"] == "").sum()
        if missing_v:
            all_warnings.append(f"MISSING_V: {sd.name} {missing_v} clones missing V gene")

        clones = load_trust4_clones(sd, chain="TRB", productive_only=True, min_count=1)
        if clones is None or clones.empty:
            all_warnings.append(f"EMPTY_AFTER_FILTER: {sd.name}")
            continue
        metrics = compute_diversity(clones)
        metrics["sample_id"] = sd.name
        for k, v in metrics.items():
            if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                all_warnings.append(f"INVALID_METRIC: {sd.name} {k}={v}")
        rows.append(metrics)

    if not rows:
        console.print("[red]No valid TRUST4 outputs — pilot FAILED.[/red]")
        raise typer.Exit(2)

    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df = df[["sample_id"] + [c for c in df.columns if c != "sample_id"]]
    csv_path = out_dir / "pilot_diversity.csv"
    df.to_csv(csv_path, index=False)

    table = Table(title="Pilot diversity metrics (TRB)")
    for col in ["sample_id", "n_clones", "shannon", "clonality", "top10_frac", "d50"]:
        table.add_column(col, justify="right" if col != "sample_id" else "left")
    for _, r in df.iterrows():
        table.add_row(
            r["sample_id"], str(int(r["n_clones"])),
            f"{r['shannon']:.3f}", f"{r['clonality']:.3f}",
            f"{r['top10_frac']:.3f}", str(int(r["d50"])),
        )
    console.print(table)

    qc_path = out_dir / "pilot_qc_report.txt"
    with qc_path.open("w") as f:
        f.write(f"Samples processed: {len(df)} / {len(sample_dirs)}\n")
        f.write(f"Warnings: {len(all_warnings)}\n\n")
        for w in all_warnings:
            f.write(w + "\n")
        f.write("\n=== Summary stats ===\n")
        f.write(
            df[["n_clones", "shannon", "clonality", "top10_frac", "d50"]]
            .describe()
            .round(3)
            .to_string()
        )

    console.print(f"[green]Diversity CSV -> {csv_path}[/green]")
    console.print(f"[green]QC report -> {qc_path}[/green]")
    if all_warnings:
        console.print(f"[yellow]{len(all_warnings)} warnings (see QC report).[/yellow]")
    else:
        console.print("[green]All pilot samples passed QC[/green]")


def _strip_allele(gene_call: object) -> str:
    if not isinstance(gene_call, str) or not gene_call or gene_call == ".":
        return ""
    return gene_call.split(",")[0].split("*")[0]


@app.command("vj-usage")
def vj_usage(
    trust4_dir: Path = typer.Option(relpath("data/repertoire/trust4/geo_pilot"), "--trust4-dir", "-t"),
    diversity_csv: Path = typer.Option(
        relpath("data/metrics/atlas_diversity.csv"), "--diversity-csv",
        help="Sample list + cohort labels.",
    ),
    out: Path = typer.Option(relpath("data/metrics/vj_usage_long.csv"), "--out", "-o"),
) -> None:
    """Aggregate per-sample TRBV / TRBJ usage (count-weighted) as a long table."""
    sample_cohort = _load_sample_cohort(diversity_csv)
    rows = []
    for sample_id, cohort in sample_cohort.items():
        sdir = trust4_dir / sample_id
        if not sdir.exists():
            continue
        clones = load_trust4_clones(sdir, chain="TRB", productive_only=True)
        if clones is None or clones.empty:
            continue
        clones = clones.copy()
        clones["v_gene"] = clones["v_call"].map(_strip_allele)
        clones["j_gene"] = clones["j_call"].map(_strip_allele)
        total = clones["duplicate_count"].sum()
        for gene_type, col in (("V", "v_gene"), ("J", "j_gene")):
            tab = (
                clones.groupby(col)["duplicate_count"].sum()
                .reset_index().rename(columns={col: "gene_name", "duplicate_count": "count"})
            )
            tab = tab[tab["gene_name"] != ""]
            tab["gene_type"] = gene_type
            tab["frequency"] = tab["count"] / total
            tab["sample_id"] = sample_id
            tab["cohort"] = cohort
            rows.append(tab[["sample_id", "cohort", "gene_type", "gene_name", "count", "frequency"]])

    if not rows:
        console.print("[red]No data produced.[/red]")
        raise typer.Exit(1)
    out_df = pd.concat(rows, ignore_index=True)
    out_df.to_csv(ensure_parent(out), index=False)
    console.print(
        f"[green]VJ usage -> {out}[/green] ({len(out_df)} rows; "
        f"{out_df.sample_id.nunique()} samples; "
        f"{out_df[out_df.gene_type=='V'].gene_name.nunique()} V / "
        f"{out_df[out_df.gene_type=='J'].gene_name.nunique()} J genes)"
    )


@app.command("clone-size")
def clone_size(
    trust4_dir: Path = typer.Option(relpath("data/repertoire/trust4/geo_pilot"), "--trust4-dir", "-t"),
    diversity_csv: Path = typer.Option(relpath("data/metrics/atlas_diversity.csv"), "--diversity-csv"),
    out_dir: Path = typer.Option(relpath("data/metrics"), "--out-dir", "-o"),
    hyperexpanded_freq: float = typer.Option(0.05, "--hyperexpanded-freq"),
) -> None:
    """Per-sample clone-size distribution + CDR3 length + top-N fractions."""
    sample_cohort = _load_sample_cohort(diversity_csv)
    out_dir.mkdir(parents=True, exist_ok=True)
    dist_rows, summary_rows = [], []
    for sample_id, cohort in sample_cohort.items():
        sdir = trust4_dir / sample_id
        clones = load_trust4_clones(sdir, chain="TRB", productive_only=True)
        if clones is None or clones.empty:
            continue
        clones = clones.sort_values("duplicate_count", ascending=False).reset_index(drop=True)
        total = int(clones["duplicate_count"].sum())
        clones["rank"] = np.arange(1, len(clones) + 1)
        clones["freq"] = clones["duplicate_count"] / total
        clones["log10_freq"] = np.log10(clones["freq"])
        clones["cdr3_aa"] = clones["junction_aa"]
        clones["cdr3_len"] = clones["junction_aa"].str.len()
        clones["sample_id"] = sample_id
        clones["cohort"] = cohort
        dist_rows.append(
            clones[["sample_id", "cohort", "rank", "duplicate_count",
                    "freq", "log10_freq", "cdr3_aa", "cdr3_len"]]
            .rename(columns={"duplicate_count": "count"})
        )
        summary_rows.append({
            "sample_id": sample_id, "cohort": cohort,
            "n_clones": len(clones), "total_reads": total,
            "top1_frac": float(clones["freq"].iloc[0]),
            "top10_frac": float(clones["freq"].iloc[:10].sum()),
            "top100_frac": float(clones["freq"].iloc[:100].sum()),
            "n_hyperexpanded": int((clones["freq"] > hyperexpanded_freq).sum()),
            "median_cdr3_len": float(clones["cdr3_len"].median()),
            "mean_cdr3_len": float(clones["cdr3_len"].mean()),
        })

    if not dist_rows:
        console.print("[red]No data.[/red]")
        raise typer.Exit(1)
    dist = pd.concat(dist_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    dist_path = out_dir / "clone_size_distribution.csv"
    summary_path = out_dir / "clone_size_summary.csv"
    dist.to_csv(dist_path, index=False)
    summary.to_csv(summary_path, index=False)
    console.print(f"[green]Distribution -> {dist_path}[/green] ({len(dist)} rows)")
    console.print(f"[green]Summary -> {summary_path}[/green] ({len(summary)} rows)")


def _kmers(seq: str, k: int) -> set[str]:
    if not isinstance(seq, str) or len(seq) < k:
        return set()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


@app.command("shared-clones")
def shared_clones(
    trust4_dir: Path = typer.Option(relpath("data/repertoire/trust4/geo_pilot"), "--trust4-dir", "-t"),
    diversity_csv: Path = typer.Option(relpath("data/metrics/atlas_diversity.csv"), "--diversity-csv"),
    out_dir: Path = typer.Option(relpath("data/metrics"), "--out-dir", "-o"),
    ks: list[int] = typer.Option([3, 4, 5], "--k", help="k-mer sizes (repeatable)."),
) -> None:
    """Cross-cohort CDR3 sharing + k-mer motif overlap."""
    sample_cohort = _load_sample_cohort(diversity_csv)
    out_dir.mkdir(parents=True, exist_ok=True)

    cohort_cdr3: dict[str, set[str]] = defaultdict(set)
    cohort_kmers: dict[str, dict[int, set[str]]] = defaultdict(lambda: {k: set() for k in ks})
    sample_cdr3: dict[str, set[str]] = {}
    cdr3_to_samples: dict[str, set[str]] = defaultdict(set)
    cdr3_to_count: dict[str, int] = defaultdict(int)

    for sample_id, cohort in sample_cohort.items():
        sdir = trust4_dir / sample_id
        clones = load_trust4_clones(sdir, chain="TRB", productive_only=True)
        if clones is None or clones.empty:
            continue
        seqs = set(clones["junction_aa"].dropna().astype(str).tolist())
        sample_cdr3[sample_id] = seqs
        cohort_cdr3[cohort].update(seqs)
        for s in seqs:
            cdr3_to_samples[s].add(sample_id)
        for s, c in zip(clones["junction_aa"], clones["duplicate_count"]):
            if isinstance(s, str):
                cdr3_to_count[s] += int(c)
        for k in ks:
            for s in seqs:
                cohort_kmers[cohort][k].update(_kmers(s, k))

    cohorts = sorted(cohort_cdr3.keys())
    pair_rows = []
    for a, b in combinations(cohorts, 2):
        A, B = cohort_cdr3[a], cohort_cdr3[b]
        inter = A & B
        union = A | B
        row = {
            "cohort_a": a, "cohort_b": b,
            "n_shared_exact": len(inter),
            "jaccard": len(inter) / len(union) if union else 0.0,
            "n_a": len(A), "n_b": len(B),
        }
        for k in ks:
            ka, kb = cohort_kmers[a][k], cohort_kmers[b][k]
            kinter, kunion = ka & kb, ka | kb
            row[f"kmer{k}_overlap_pct"] = 100 * len(kinter) / len(kunion) if kunion else 0.0
            row[f"kmer{k}_shared"] = len(kinter)
        pair_rows.append(row)
    for c in cohorts:
        row = {"cohort_a": c, "cohort_b": c, "n_shared_exact": len(cohort_cdr3[c]),
               "jaccard": 1.0, "n_a": len(cohort_cdr3[c]), "n_b": len(cohort_cdr3[c])}
        for k in ks:
            row[f"kmer{k}_overlap_pct"] = 100.0
            row[f"kmer{k}_shared"] = len(cohort_kmers[c][k])
        pair_rows.append(row)
    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(out_dir / "shared_clones.csv", index=False)

    rows = []
    for s, samples in cdr3_to_samples.items():
        if len(samples) >= 2:
            ch = sorted({sample_cohort[x] for x in samples})
            rows.append({
                "cdr3_aa": s, "n_samples": len(samples), "n_cohorts": len(ch),
                "cohorts": ";".join(ch), "samples": ";".join(sorted(samples)),
                "total_count": cdr3_to_count[s], "cdr3_len": len(s),
            })
    top_clones = (
        pd.DataFrame(rows)
        .sort_values(["n_samples", "total_count"], ascending=[False, False])
        .reset_index(drop=True)
        if rows else pd.DataFrame()
    )
    top_clones.to_csv(out_dir / "top_shared_clones.csv", index=False)

    motif_rows = []
    for k in ks:
        kmer_cohort: dict[str, set[str]] = defaultdict(set)
        kmer_sample: dict[str, set[str]] = defaultdict(set)
        for sample_id, seqs in sample_cdr3.items():
            cohort = sample_cohort[sample_id]
            for s in seqs:
                for kk in _kmers(s, k):
                    kmer_cohort[kk].add(cohort)
                    kmer_sample[kk].add(sample_id)
        top_motifs = sorted(
            kmer_cohort.items(),
            key=lambda kv: (len(kv[1]), len(kmer_sample[kv[0]])),
            reverse=True,
        )[:30]
        for motif, cs in top_motifs:
            motif_rows.append({
                "motif_kmer": motif, "k": k,
                "n_cohorts": len(cs), "n_samples": len(kmer_sample[motif]),
                "cohorts": ";".join(sorted(cs)),
            })
    pd.DataFrame(motif_rows).to_csv(out_dir / "top_shared_motifs.csv", index=False)

    n_cohorts = len(cohorts)
    cum = [{
        "n_cohorts_min": n,
        "n_clones_shared": sum(
            1 for _s, samples in cdr3_to_samples.items()
            if len({sample_cohort[x] for x in samples}) >= n
        ),
    } for n in range(1, n_cohorts + 1)]
    pd.DataFrame(cum).to_csv(out_dir / "cumulative_shared_by_cohort.csv", index=False)

    console.print(
        f"[green]pairs={len(pair_df)} / multi_sample_CDR3s={len(top_clones)} / "
        f"motifs={len(motif_rows)} written under {out_dir}[/green]"
    )


@app.command("rarefaction")
def rarefaction(
    trust4_dir: Path = typer.Option(relpath("data/repertoire/trust4/geo_pilot"), "--trust4-dir", "-t"),
    diversity_csv: Path = typer.Option(relpath("data/metrics/atlas_diversity.csv"), "--diversity-csv"),
    out: Path = typer.Option(relpath("data/metrics/rarefaction.csv"), "--out", "-o"),
    fractions: list[float] = typer.Option(
        [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00],
        "--fraction",
    ),
    n_bootstrap: int = typer.Option(5, "--n-bootstrap"),
    seed: int = typer.Option(20260527, "--seed"),
) -> None:
    """Rarefaction curves: unique-clones vs subsampled read fraction."""
    sample_cohort = _load_sample_cohort(diversity_csv)
    rng = np.random.default_rng(seed)
    rows = []
    for sample_id, cohort in sample_cohort.items():
        sdir = trust4_dir / sample_id
        clones = load_trust4_clones(sdir, chain="TRB", productive_only=True)
        if clones is None or clones.empty:
            continue
        clone_ids = np.arange(len(clones))
        counts = clones["duplicate_count"].astype(int).to_numpy()
        pool = np.repeat(clone_ids, counts)
        total_reads = len(pool)
        if total_reads == 0:
            continue
        for frac in fractions:
            n_pick = max(1, int(round(total_reads * frac)))
            uniq = [
                int(len(np.unique(rng.choice(pool, size=n_pick, replace=False))))
                for _ in range(n_bootstrap)
            ]
            rows.append({
                "sample_id": sample_id, "cohort": cohort, "fraction": frac,
                "n_reads_sampled": n_pick,
                "unique_clones": float(np.mean(uniq)),
                "unique_clones_sd": float(np.std(uniq)),
                "total_reads": total_reads,
            })
    df = pd.DataFrame(rows)
    df.to_csv(ensure_parent(out), index=False)
    console.print(f"[green]Rarefaction -> {out}[/green] ({len(df)} rows; {df.sample_id.nunique()} samples)")


@app.command("gse-diagnostic")
def gse_diagnostic(
    trust4_dir: Path = typer.Option(relpath("data/repertoire/trust4/geo_pilot"), "--trust4-dir", "-t"),
    out_dir: Path = typer.Option(relpath("data/metrics"), "--out-dir", "-o"),
    samples: list[str] = typer.Option(
        [
            "SRR17491517:GSE193258:Problematic (excluded)",
            "SRR17491496:GSE193258:Problematic (excluded)",
            "SRR11094249:GSE145370:Good cohort",
            "SRR8526721:GSE126044:Good cohort",
        ],
        "--sample",
        help="SRR:COHORT:LABEL  — repeatable.",
    ),
) -> None:
    """Diagnostic for cohort exclusion: parse ``*_cdr3.out`` and quantify V+J completeness."""
    out_dir.mkdir(parents=True, exist_ok=True)

    per_rows, sum_rows = [], []
    for spec in samples:
        parts = spec.split(":")
        if len(parts) != 3:
            console.print(f"[red]bad --sample '{spec}' (expected SRR:COHORT:LABEL)[/red]")
            continue
        sample_id, cohort, label = parts
        path = trust4_dir / sample_id / f"{sample_id}_cdr3.out"
        if not path.exists():
            console.print(f"[yellow]missing {path}[/yellow]", file=sys.stderr)
            continue
        rows = []
        with open(path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 9:
                    continue
                v_call, j_call, cdr3_nt = p[2], p[4], p[8]
                has_v = v_call not in ("", "*", ".")
                has_j = j_call not in ("", "*", ".")
                completeness = ("V+J" if has_v and has_j else
                                "V-only" if has_v else
                                "J-only" if has_j else "none")
                chain = ""
                for prefix in ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"):
                    if v_call.startswith(prefix) or j_call.startswith(prefix):
                        chain = prefix
                        break
                rows.append({
                    "assemble_id": p[0], "v_call": v_call, "j_call": j_call,
                    "completeness": completeness, "chain": chain,
                    "cdr3_nt_len": len(cdr3_nt) if cdr3_nt and cdr3_nt != "*" else 0,
                })
        if not rows:
            continue
        df = pd.DataFrame(rows)
        df["sample_id"], df["cohort"], df["label"] = sample_id, cohort, label
        per_rows.append(df)
        n = len(df)
        counts = df["completeness"].value_counts().to_dict()
        sum_rows.append({
            "sample_id": sample_id, "cohort": cohort, "label": label,
            "n_assemblies": n,
            "n_V_plus_J": counts.get("V+J", 0),
            "n_V_only": counts.get("V-only", 0),
            "n_J_only": counts.get("J-only", 0),
            "n_none": counts.get("none", 0),
            "pct_V_plus_J": 100 * counts.get("V+J", 0) / n,
            "median_cdr3_nt_len": float(df["cdr3_nt_len"].median()),
            "mean_cdr3_nt_len": float(df["cdr3_nt_len"].mean()),
        })

    if not per_rows:
        console.print("[red]No data parsed.[/red]")
        raise typer.Exit(1)

    per = pd.concat(per_rows, ignore_index=True)
    summary = pd.DataFrame(sum_rows)
    per.to_csv(out_dir / "gse193258_diagnostic_per_assembly.csv", index=False)
    summary.to_csv(out_dir / "gse193258_diagnostic_summary.csv", index=False)
    console.print(summary.to_string(index=False))
    console.print(f"[green]Wrote diagnostic CSVs under {out_dir}[/green]")
