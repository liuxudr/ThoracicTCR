#!/usr/bin/env bash
# Render all SCI figures via Rscript.
# Reads:  data/metrics/atlas_diversity.csv + data/atlas/vdjdb_matches.csv
# Writes: paper/figures/fig*.{pdf,png}

set -euo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found. Install R first."
  exit 1
fi

export THORACICTCR_REPO="$REPO_ROOT"
mkdir -p "$REPO_ROOT/paper/figures" "$REPO_ROOT/paper/tables"

FIGURES=(
  fig1_landscape.R
  fig2_diversity_clonality.R
  fig3_vj_usage.R
  fig4_clonal_expansion.R
  fig5_antigen.R
  fig6_shared_motifs.R
  figS1_pipeline.R
  figS2_cdr3_length.R
  figS3_depth.R
  figS4_gse193258.R
  figS5_hill.R
  figS6_rarefaction.R
  figS7_vdjdb_detail.R
  figS8_attrition.R
)

for fr in "${FIGURES[@]}"; do
  echo "=== Rendering $fr ==="
  if [[ ! -f "$REPO_ROOT/R/$fr" ]]; then
    echo "  [skip] $fr not found"
    continue
  fi
  Rscript "$REPO_ROOT/R/$fr" 2>&1 | tail -4
done

echo ""
echo "[INFO] Figures rendered:"
ls -la "$REPO_ROOT/paper/figures/"*.pdf 2>/dev/null
