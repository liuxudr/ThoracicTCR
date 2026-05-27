#!/usr/bin/env bash
# Download 10-sample GEO/SRA pilot fastq via aria2 from ENA mirror.
# Reads URL list from configs/pilot_geo_10.urls.txt (produced by 01b_build_geo_pilot_manifest.py).
#
# Usage:
#   bash scripts/02b_download_pilot_geo_fastq.sh           # download
#   DRY_RUN=1 bash scripts/02b_download_pilot_geo_fastq.sh # preview total size + first URLs

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
URL_LIST="${REPO_ROOT}/configs/pilot_geo_10.urls.txt"
OUT_DIR="${REPO_ROOT}/data/raw/fastq"
LOG_DIR="${REPO_ROOT}/logs/geo_download"
N_PARALLEL="${N_PARALLEL:-4}"
N_CONN_PER_HOST="${N_CONN_PER_HOST:-4}"

if [[ ! -f "$URL_LIST" ]]; then
  echo "[ERROR] URL list not found: $URL_LIST"
  echo "        Run: python scripts/01b_build_geo_pilot_manifest.py"
  exit 1
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"

N_FILES=$(grep -c "^https" "$URL_LIST" || true)
echo "[INFO] URL list: $URL_LIST  (n_files=$N_FILES)"
echo "[INFO] Output:   $OUT_DIR"

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY-RUN] First 4 URLs:"
  grep "^https" "$URL_LIST" | head -4 | sed 's/^/   /'
  echo "[DRY-RUN] Would run:"
  echo "  aria2c -i $URL_LIST -j $N_PARALLEL -x $N_CONN_PER_HOST -c --auto-file-renaming=false"
  exit 0
fi

# Prefer aria2 (parallel + resumable); fall back to wget if absent.
if command -v aria2c >/dev/null 2>&1; then
  echo "[INFO] Using aria2c (j=$N_PARALLEL, x=$N_CONN_PER_HOST, --check-integrity=true)"
  aria2c \
    -i "$URL_LIST" \
    -j "$N_PARALLEL" \
    -x "$N_CONN_PER_HOST" \
    -s "$N_CONN_PER_HOST" \
    -c \
    --check-integrity=true \
    --auto-file-renaming=false \
    --console-log-level=warn \
    --summary-interval=10 \
    --log="$LOG_DIR/aria2.log"

  echo ""
  echo "[INFO] Post-download gzip integrity check"
  N_CORRUPT=0
  while IFS= read -r url; do
    # Parse the next 'out=' line to get filename
    :
  done < /dev/null
  for f in "$OUT_DIR"/*.fastq.gz; do
    if ! gzip -t "$f" 2>/dev/null; then
      echo "  ✗ CORRUPT: $(basename "$f")"
      N_CORRUPT=$((N_CORRUPT + 1))
    fi
  done
  if [[ $N_CORRUPT -gt 0 ]]; then
    echo "[ERROR] $N_CORRUPT corrupt files detected. Delete them and re-run this script."
    exit 1
  fi
  echo "  ✓ All gzip integrity OK"
else
  echo "[WARN] aria2c not found; falling back to wget (slower, sequential per URL)"
  grep "^https" "$URL_LIST" | while read -r url; do
    fname=$(basename "$url")
    wget -c -P "$OUT_DIR" "$url" 2>>"$LOG_DIR/wget.log"
  done
fi

echo ""
echo "[INFO] Downloaded files:"
find "$OUT_DIR" -name "*.fastq.gz" | sort | head
echo ""
echo "[INFO] Total downloaded size:"
du -sh "$OUT_DIR"
