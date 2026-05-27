#!/usr/bin/env bash
# Download public TCR-antigen databases used in Phase 5 antigen lookup.
#
# Targets:
#   - VDJdb        (CC-BY 4.0)     https://github.com/antigenomics/vdjdb-db/releases
#   - McPAS-TCR    (free academic) http://friedmanlab.weizmann.ac.il/McPAS-TCR/
#   - IEDB         (free)          https://www.iedb.org/database_export_v3.php
#
# Output: data/external/<db>/
# Total size ~30 MB. Safe to run anytime, no quota / no auth.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EXT="$REPO_ROOT/data/external"

mkdir -p "$EXT/vdjdb" "$EXT/mcpas" "$EXT/iedb"

###############################################################################
# 1. VDJdb — latest release tarball
###############################################################################
echo "[INFO] Fetching VDJdb release list..."
VDJDB_URL=$(curl -fsSL https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest \
            | grep -E '"browser_download_url".*vdjdb.*\.zip"' \
            | head -1 \
            | sed -E 's/.*"([^"]+)".*/\1/')

if [[ -n "$VDJDB_URL" ]]; then
  echo "[INFO] VDJdb URL: $VDJDB_URL"
  curl -fSL -o "$EXT/vdjdb/vdjdb.zip" "$VDJDB_URL"
  ( cd "$EXT/vdjdb" && unzip -o vdjdb.zip > /dev/null && rm -f vdjdb.zip )
  echo "  ✓ VDJdb extracted to $EXT/vdjdb"
  ls "$EXT/vdjdb" | head -5 | sed 's/^/    /'
else
  echo "  ✗ Could not resolve VDJdb release URL"
fi

###############################################################################
# 2. McPAS-TCR — single CSV
###############################################################################
MCPAS_URL="http://friedmanlab.weizmann.ac.il/McPAS-TCR/session/c5cf60bc63f9c3e1c1fdf8b81f56c2dc/download/downloadDB?w="
# Stable mirror (when site is slow):
# MCPAS_FALLBACK="https://github.com/.../McPAS-TCR.csv"  # community mirrors exist

if curl -fSL -o "$EXT/mcpas/McPAS-TCR.csv" "$MCPAS_URL"; then
  echo "  ✓ McPAS-TCR -> $EXT/mcpas/McPAS-TCR.csv"
  awk 'NR<=3' "$EXT/mcpas/McPAS-TCR.csv" | cut -c1-100 | sed 's/^/    /'
else
  echo "  ! McPAS-TCR download failed (site may be temporarily down)."
  echo "    Retry later: bash $0 mcpas"
fi

###############################################################################
# 3. IEDB — receptor + epitope subsets
###############################################################################
IEDB_BASE="https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"
if curl -fSL -o "$EXT/iedb/receptor_full_v3.zip" "$IEDB_BASE"; then
  ( cd "$EXT/iedb" && unzip -o receptor_full_v3.zip > /dev/null && rm -f receptor_full_v3.zip )
  echo "  ✓ IEDB TCR receptors -> $EXT/iedb"
  ls "$EXT/iedb" | head -5 | sed 's/^/    /'
else
  echo "  ! IEDB receptor download failed (rate-limit?). Skip."
fi

echo ""
echo "[INFO] Public TCR DB mirror sizes:"
du -sh "$EXT"/* 2>/dev/null
