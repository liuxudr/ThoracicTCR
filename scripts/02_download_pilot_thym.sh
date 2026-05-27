#!/usr/bin/env bash
# Download 10 TCGA-THYM RNA-seq BAMs using gdc-client.
# Expects scripts/01_build_manifest.py to have produced configs/pilot_thym_10.gdc_manifest.txt
#
# Usage:
#   bash scripts/02_download_pilot_thym.sh           # standard
#   DRY_RUN=1 bash scripts/02_download_pilot_thym.sh # show what would download
#
# Note: TCGA "open" access BAMs are publicly downloadable WITHOUT a dbGaP token.
# For controlled-access fastq, set GDC_TOKEN env var to path of token file.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="${REPO_ROOT}/configs/pilot_thym_10.gdc_manifest.txt"
OUT_DIR="${REPO_ROOT}/data/raw/bam/TCGA-THYM"
N_JOBS="${N_JOBS:-4}"

if [[ ! -f "$MANIFEST" ]]; then
  echo "[ERROR] Manifest not found: $MANIFEST"
  echo "        Run: python scripts/01_build_manifest.py"
  exit 1
fi

mkdir -p "$OUT_DIR"

echo "[INFO] Manifest: $MANIFEST"
echo "[INFO] Output:   $OUT_DIR"
echo "[INFO] Total bytes in manifest:"
awk -F'\t' 'NR>1 {sum+=$4} END {printf "       %.2f GB\n", sum/1e9}' "$MANIFEST"

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY-RUN] Would run:"
  echo "  gdc-client download -m $MANIFEST -d $OUT_DIR -n $N_JOBS"
  exit 0
fi

GDC_ARGS=(-m "$MANIFEST" -d "$OUT_DIR" -n "$N_JOBS")
if [[ -n "${GDC_TOKEN:-}" && -f "$GDC_TOKEN" ]]; then
  GDC_ARGS+=(-t "$GDC_TOKEN")
  echo "[INFO] Using GDC token: $GDC_TOKEN"
fi

echo "[INFO] Starting gdc-client download with $N_JOBS parallel jobs..."
gdc-client download "${GDC_ARGS[@]}"

echo ""
echo "[OK] Downloaded files:"
find "$OUT_DIR" -name "*.bam" | head
echo ""
echo "[INFO] Total downloaded size:"
du -sh "$OUT_DIR"
