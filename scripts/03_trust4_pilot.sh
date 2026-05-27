#!/usr/bin/env bash
# Run TRUST4 on pilot samples (BAM or paired FASTQ).
#
# Auto-detects input mode:
#   - $BAM_ROOT (TCGA path)    -> BAM input
#   - $FASTQ_ROOT (GEO path)   -> paired FASTQ input
#
# Usage:
#   conda activate thoracictcr
#
#   # GEO pilot (default, after running 01b + 02b):
#   bash scripts/03_trust4_pilot.sh
#
#   # TCGA pilot (after running 01 + 02 with dbGaP token):
#   MODE=tcga bash scripts/03_trust4_pilot.sh
#
# Env vars:
#   MODE             "geo" (default) | "tcga"
#   N_PARALLEL       2  (samples in parallel)
#   THREADS_PER_JOB  4  (TRUST4 threads per sample)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE="${MODE:-geo}"
N_PARALLEL="${N_PARALLEL:-2}"
THREADS_PER_JOB="${THREADS_PER_JOB:-4}"

case "$MODE" in
  geo)
    INPUT_ROOT="${FASTQ_ROOT:-${REPO_ROOT}/data/raw/fastq}"
    OUT_ROOT="${REPO_ROOT}/data/repertoire/trust4/geo_pilot"
    INPUT_KIND="fastq"
    ;;
  tcga)
    INPUT_ROOT="${BAM_ROOT:-${REPO_ROOT}/data/raw/bam/TCGA-THYM}"
    OUT_ROOT="${REPO_ROOT}/data/repertoire/trust4/TCGA-THYM"
    INPUT_KIND="bam"
    ;;
  *)
    echo "[ERROR] Unknown MODE=$MODE (use 'geo' or 'tcga')"
    exit 1
    ;;
esac
LOG_DIR="${REPO_ROOT}/logs/trust4_pilot"

mkdir -p "$OUT_ROOT" "$LOG_DIR"

if ! command -v run-trust4 >/dev/null 2>&1; then
  echo "[ERROR] run-trust4 not found on PATH."
  echo "        Run: conda activate thoracictcr"
  exit 1
fi

# Locate TRUST4 reference files. Search order:
#   1. $BCRTCR_FA / $IMGT_FA env vars (highest priority)
#   2. $REPO_ROOT/data/external/trust4_ref/    (downloaded from TRUST4 github)
#   3. $(which run-trust4)/../share/trust4/    (bioconda — empty in some builds)
FALLBACK_REF_DIR="$REPO_ROOT/data/external/trust4_ref"
TRUST4_SHARE="$(dirname "$(which run-trust4)")/../share/trust4"
: "${BCRTCR_FA:=}"
: "${IMGT_FA:=}"
[[ -z "$BCRTCR_FA" && -f "$FALLBACK_REF_DIR/hg38_bcrtcr.fa" ]] && BCRTCR_FA="$FALLBACK_REF_DIR/hg38_bcrtcr.fa"
[[ -z "$IMGT_FA"   && -f "$FALLBACK_REF_DIR/human_IMGT+C.fa" ]] && IMGT_FA="$FALLBACK_REF_DIR/human_IMGT+C.fa"
[[ -z "$BCRTCR_FA" ]] && BCRTCR_FA="$(find "$TRUST4_SHARE" -name 'hg38_bcrtcr.fa' -o -name 'bcrtcr.fa' 2>/dev/null | head -1)"
[[ -z "$IMGT_FA"   ]] && IMGT_FA="$(find "$TRUST4_SHARE" -name 'human_IMGT+C.fa' -o -name 'IMGT+C.fa' 2>/dev/null | head -1)"

if [[ -z "$BCRTCR_FA" || -z "$IMGT_FA" ]]; then
  echo "[ERROR] Could not locate TRUST4 reference files."
  echo "        Tried:"
  echo "          \$BCRTCR_FA / \$IMGT_FA env vars"
  echo "          $FALLBACK_REF_DIR/{hg38_bcrtcr.fa,human_IMGT+C.fa}"
  echo "          $TRUST4_SHARE"
  echo "        To fetch references, run:"
  echo "          mkdir -p $FALLBACK_REF_DIR && cd $FALLBACK_REF_DIR &&"
  echo "          curl -fsSL -O https://raw.githubusercontent.com/liulab-dfci/TRUST4/master/hg38_bcrtcr.fa &&"
  echo "          curl -fsSL -O https://raw.githubusercontent.com/liulab-dfci/TRUST4/master/human_IMGT+C.fa"
  exit 1
fi
echo "[INFO] MODE: $MODE ($INPUT_KIND)"
echo "[INFO] Input root: $INPUT_ROOT"
echo "[INFO] Output:     $OUT_ROOT"
echo "[INFO] BCR/TCR ref: $BCRTCR_FA"
echo "[INFO] IMGT+C ref:  $IMGT_FA"

# Worker functions
run_one_bam() {
  local bam="$1"
  local sample_id; sample_id="$(basename "$(dirname "$bam")")"
  local out_dir="$OUT_ROOT/$sample_id"
  local log="$LOG_DIR/${sample_id}.log"

  if [[ -f "$out_dir/${sample_id}_report.tsv" ]]; then
    echo "[SKIP] $sample_id (already done)"
    return 0
  fi
  mkdir -p "$out_dir"
  echo "[RUN]  $sample_id (bam)"
  run-trust4 \
    -b "$bam" \
    -f "$BCRTCR_FA" \
    --ref "$IMGT_FA" \
    -t "$THREADS_PER_JOB" \
    -o "$sample_id" \
    --od "$out_dir" \
      > "$log" 2>&1 \
  && echo "[OK]   $sample_id" \
  || echo "[FAIL] $sample_id (see $log)"
}

run_one_fastq_pair() {
  local r1="$1"
  local sample_id; sample_id="$(basename "$r1" | sed -E 's/_R?1\.fastq(\.gz)?$//; s/_1\.fastq(\.gz)?$//')"
  local r2="${r1/_R1./_R2.}"
  [[ ! -f "$r2" ]] && r2="${r1/_1./_2.}"
  if [[ ! -f "$r2" ]]; then
    echo "[FAIL] $sample_id (cannot find R2 mate for $r1)"
    return 1
  fi
  local out_dir="$OUT_ROOT/$sample_id"
  local log="$LOG_DIR/${sample_id}.log"
  if [[ -f "$out_dir/${sample_id}_report.tsv" ]]; then
    echo "[SKIP] $sample_id (already done)"
    return 0
  fi
  mkdir -p "$out_dir"
  echo "[RUN]  $sample_id (fastq paired)"
  run-trust4 \
    -1 "$r1" -2 "$r2" \
    -f "$BCRTCR_FA" \
    --ref "$IMGT_FA" \
    -t "$THREADS_PER_JOB" \
    -o "$sample_id" \
    --od "$out_dir" \
      > "$log" 2>&1 \
  && echo "[OK]   $sample_id" \
  || echo "[FAIL] $sample_id (see $log)"
}

export -f run_one_bam run_one_fastq_pair
export OUT_ROOT LOG_DIR BCRTCR_FA IMGT_FA THREADS_PER_JOB

# Discover inputs and dispatch
if [[ "$INPUT_KIND" == "bam" ]]; then
  mapfile -t INPUTS < <(find "$INPUT_ROOT" -name "*.bam" -not -name "*.bai" | sort)
  [[ ${#INPUTS[@]} -eq 0 ]] && { echo "[ERROR] No BAM files under $INPUT_ROOT"; exit 1; }
  echo "[INFO] Found ${#INPUTS[@]} BAM(s)"
  printf '%s\n' "${INPUTS[@]}" | parallel -j "$N_PARALLEL" run_one_bam {}
else
  # Find R1 mates (handles both `_R1.fastq.gz` and `_1.fastq.gz` naming)
  mapfile -t INPUTS < <(find "$INPUT_ROOT" \( -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" \) | sort -u)
  [[ ${#INPUTS[@]} -eq 0 ]] && { echo "[ERROR] No paired fastq under $INPUT_ROOT"; exit 1; }
  echo "[INFO] Found ${#INPUTS[@]} fastq R1 mate(s)"
  printf '%s\n' "${INPUTS[@]}" | parallel -j "$N_PARALLEL" run_one_fastq_pair {}
fi

echo ""
echo "[INFO] TRUST4 outputs:"
find "$OUT_ROOT" -name "*_report.tsv" | sort
echo ""
echo "[INFO] Clone counts per sample:"
for f in $(find "$OUT_ROOT" -name "*_report.tsv" | sort); do
  n=$(awk 'NR>1 {c++} END {print c+0}' "$f")
  echo "  $(basename "$(dirname "$f")"): $n unique clones"
done
