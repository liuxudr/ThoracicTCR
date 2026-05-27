"""TRUST4 runner: reconstruct TCR/BCR clones from RNA-seq BAM or FASTQ."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Literal

from rich.console import Console

console = Console()


class TRUST4Runner:
    """Wrap TRUST4 (`run-trust4`) calls.

    Required reference files (set via conda env or download manually):
        - bcrtcr.fa  : IMGT BCR/TCR reference
        - IMGT+C.fa  : IMGT human alleles
    These are bundled in TRUST4's `hg38_bcrtcr.fa` etc when installed via bioconda.
    """

    def __init__(
        self,
        trust4_bin: str = "run-trust4",
        bcrtcr_fa: Path | None = None,
        imgt_fa: Path | None = None,
        n_threads: int = 8,
        genome: Literal["hg38", "hg19"] = "hg38",
    ):
        self.trust4_bin = trust4_bin
        self.bcrtcr_fa = bcrtcr_fa
        self.imgt_fa = imgt_fa
        self.n_threads = n_threads
        self.genome = genome
        if shutil.which(trust4_bin) is None:
            console.print(
                f"[yellow]WARNING: '{trust4_bin}' not on PATH. "
                f"Install via: mamba install -c bioconda trust4[/yellow]"
            )

    def _resolve_refs(self) -> tuple[Path, Path]:
        """Locate bcrtcr.fa and IMGT+C.fa from TRUST4 install if not set."""
        if self.bcrtcr_fa and self.imgt_fa:
            return self.bcrtcr_fa, self.imgt_fa
        # bioconda installs reference next to binary
        trust4_path = shutil.which(self.trust4_bin)
        if trust4_path is None:
            raise FileNotFoundError(
                "TRUST4 binary not found. Install via 'mamba install -c bioconda trust4'."
            )
        share_dir = Path(trust4_path).resolve().parent.parent / "share" / "trust4"
        # Try common naming
        candidates_bcr = list(share_dir.rglob("hg38_bcrtcr.fa")) + list(share_dir.rglob("bcrtcr.fa"))
        candidates_imgt = list(share_dir.rglob("human_IMGT+C.fa")) + list(share_dir.rglob("IMGT+C.fa"))
        if not candidates_bcr or not candidates_imgt:
            raise FileNotFoundError(
                f"Could not auto-locate TRUST4 reference files under {share_dir}. "
                "Pass bcrtcr_fa and imgt_fa explicitly."
            )
        return candidates_bcr[0], candidates_imgt[0]

    def run_on_bam(
        self,
        bam_path: Path,
        output_dir: Path,
        prefix: str | None = None,
    ) -> Path:
        """Run TRUST4 on a BAM file. Returns path to *_report.tsv."""
        bcrtcr, imgt = self._resolve_refs()
        output_dir.mkdir(parents=True, exist_ok=True)
        prefix = prefix or bam_path.stem
        cmd = [
            self.trust4_bin,
            "-b", str(bam_path),
            "-f", str(bcrtcr),
            "--ref", str(imgt),
            "-t", str(self.n_threads),
            "-o", prefix,
            "--od", str(output_dir),
        ]
        console.print(f"[cyan]TRUST4 ({prefix}):[/cyan] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        report = output_dir / f"{prefix}_report.tsv"
        if not report.exists():
            raise RuntimeError(f"TRUST4 finished but {report} missing")
        return report

    def run_on_fastq(
        self,
        fq1: Path,
        fq2: Path | None,
        output_dir: Path,
        prefix: str,
    ) -> Path:
        bcrtcr, imgt = self._resolve_refs()
        output_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            self.trust4_bin,
            "-1", str(fq1),
            "-f", str(bcrtcr),
            "--ref", str(imgt),
            "-t", str(self.n_threads),
            "-o", prefix,
            "--od", str(output_dir),
        ]
        if fq2:
            cmd.insert(3, "-2")
            cmd.insert(4, str(fq2))
        subprocess.run(cmd, check=True)
        return output_dir / f"{prefix}_report.tsv"
