#!/usr/bin/env bash
# Run arcasHLA on all paired fastq for HLA-I + HLA-II typing.
#
# Inputs:  data/raw/fastq/SRR*_R{1,2}.fastq.gz
# Outputs: data/hla/arcashla/<sample>/<sample>.genotype.json
#
# Time:   ~10-30 min per sample on 4 threads (HLA-only K-mer matching)
# Memory: ~8 GB RAM per active worker

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FQ_DIR="$REPO_ROOT/data/raw/fastq"
OUT_ROOT="$REPO_ROOT/data/hla/arcashla"
LOG_DIR="$REPO_ROOT/logs/arcashla"
N_PARALLEL="${N_PARALLEL:-2}"
THREADS_PER_JOB="${THREADS_PER_JOB:-4}"
GENES="${GENES:-A,B,C,DRB1,DQB1,DPB1}"

mkdir -p "$OUT_ROOT" "$LOG_DIR"

if ! command -v arcasHLA >/dev/null 2>&1; then
  echo "[ERROR] arcasHLA not found. conda activate thoracictcr"
  exit 1
fi

run_one() {
  local r1="$1"
  local sample; sample="$(basename "$r1" _R1.fastq.gz)"
  local r2="$(dirname "$r1")/${sample}_R2.fastq.gz"
  [[ ! -f "$r2" ]] && { echo "[FAIL] $sample (no R2)"; return 1; }
  local out_dir="$OUT_ROOT/$sample"
  local log="$LOG_DIR/${sample}.log"

  if [[ -f "$out_dir/${sample}.genotype.json" ]] \
     || ls "$out_dir"/*.genotype.json >/dev/null 2>&1; then
    echo "[SKIP] $sample"
    return 0
  fi
  mkdir -p "$out_dir"
  echo "[RUN]  $sample"
  arcasHLA genotype "$r1" "$r2" \
    --genes "$GENES" \
    --outdir "$out_dir" \
    -t "$THREADS_PER_JOB" \
    -v \
      > "$log" 2>&1 \
  && echo "[OK]   $sample" \
  || echo "[FAIL] $sample (see $log)"
}
export -f run_one
export OUT_ROOT LOG_DIR THREADS_PER_JOB GENES

mapfile -t R1S < <(find "$FQ_DIR" -name "*_R1.fastq.gz" | sort)
if [[ ${#R1S[@]} -eq 0 ]]; then
  echo "[ERROR] No fastq pairs under $FQ_DIR"
  exit 1
fi
echo "[INFO] Found ${#R1S[@]} samples; running $N_PARALLEL parallel × $THREADS_PER_JOB threads each"

printf '%s\n' "${R1S[@]}" | parallel -j "$N_PARALLEL" run_one {}

echo ""
echo "[INFO] HLA results:"
find "$OUT_ROOT" -name "*.genotype.json" | sort | while read -r gj; do
  s=$(basename "$(dirname "$gj")")
  echo "  $s: $(jq -c . "$gj" 2>/dev/null | cut -c1-180)"
done
