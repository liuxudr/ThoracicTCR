"""arcasHLA wrapper for HLA typing from paired RNA-Seq fastq.

arcasHLA expects pre-aligned BAM by default, but supports raw fastq via:
    arcasHLA extract <bam>           # OR
    arcasHLA genotype <r1.fq.gz> <r2.fq.gz>

We use the fastq path directly since our pilot data is fastq, not BAM.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


class ArcasHLARunner:
    """Run arcasHLA on paired fastq to call HLA-I + HLA-II alleles."""

    def __init__(
        self,
        binary: str = "arcasHLA",
        genes: tuple[str, ...] = ("A", "B", "C", "DRB1", "DQB1", "DPB1"),
        n_threads: int = 4,
    ):
        self.binary = binary
        self.genes = genes
        self.n_threads = n_threads
        if shutil.which(binary) is None:
            console.print(
                f"[yellow]WARNING: '{binary}' not on PATH. "
                f"Install: mamba install -c bioconda arcas-hla[/yellow]"
            )

    def genotype_from_fastq(
        self,
        fq1: Path,
        fq2: Path,
        output_dir: Path,
        sample_id: str | None = None,
    ) -> Path:
        """Call HLA alleles from paired fastq. Returns path to *_genotype.json."""
        output_dir.mkdir(parents=True, exist_ok=True)
        sample_id = sample_id or fq1.stem.replace("_R1.fastq", "").replace(".gz", "")

        cmd = [
            self.binary, "genotype",
            str(fq1), str(fq2),
            "--genes", ",".join(self.genes),
            "--outdir", str(output_dir),
            "-t", str(self.n_threads),
            "-v",
        ]
        console.print(f"[cyan]arcasHLA ({sample_id}):[/cyan] {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            console.print(f"[red]arcasHLA failed for {sample_id}:[/red]\n{result.stderr[-500:]}")
            raise RuntimeError(f"arcasHLA failed for {sample_id}")

        # arcasHLA writes <prefix>.genotype.json
        candidates = list(output_dir.glob(f"*genotype.json"))
        if not candidates:
            raise RuntimeError(f"No genotype.json produced for {sample_id}")
        return candidates[0]

    def parse_genotype(self, json_path: Path) -> dict[str, list[str]]:
        """Return {gene: [allele1, allele2]} from arcasHLA genotype.json."""
        with json_path.open() as f:
            raw = json.load(f)
        return raw


def batch_run_arcashla(
    fastq_dir: Path,
    output_dir: Path,
    n_threads: int = 4,
    skip_existing: bool = True,
) -> dict[str, dict | str]:
    """Run arcasHLA on all paired fastq under fastq_dir."""
    runner = ArcasHLARunner(n_threads=n_threads)
    results: dict[str, dict | str] = {}
    pairs: list[tuple[Path, Path, str]] = []
    for r1 in sorted(fastq_dir.glob("*_R1.fastq.gz")):
        sample = r1.name.replace("_R1.fastq.gz", "")
        r2 = fastq_dir / f"{sample}_R2.fastq.gz"
        if r2.exists():
            pairs.append((r1, r2, sample))

    console.print(f"[cyan]Running arcasHLA on {len(pairs)} samples[/cyan]")
    for r1, r2, sample in pairs:
        out_sub = output_dir / sample
        existing = list(out_sub.glob("*genotype.json"))
        if skip_existing and existing:
            console.print(f"[dim]SKIP {sample} (already done)[/dim]")
            results[sample] = runner.parse_genotype(existing[0])
            continue
        try:
            gjson = runner.genotype_from_fastq(r1, r2, out_sub, sample)
            results[sample] = runner.parse_genotype(gjson)
        except RuntimeError as e:
            results[sample] = f"FAILED: {e}"
    return results
