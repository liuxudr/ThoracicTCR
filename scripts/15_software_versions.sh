#!/usr/bin/env bash
# Generate paper/tables/tableS5_software_versions.tsv by probing each tool
REPO=/home/bio1/workdata1/WZ/XXG2/ideaB
OUT=$REPO/paper/tables/tableS5_software_versions.tsv
PY=/home/bio1/biosoft/miniconda3/envs/thoracictcr/bin/python
RBIN=/usr/bin/Rscript
mkdir -p "$REPO/paper/tables"

# Always source the conda env so tool wrappers resolve
source /home/bio1/biosoft/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate thoracictcr 2>/dev/null || true

emit () {
  # tool<TAB>version<TAB>source
  printf "%s\t%s\t%s\n" "$1" "$2" "$3" >> "$OUT"
}

safe_run () {
  out=$("$@" 2>&1 | head -1 || true)
  printf "%s" "$out"
}

> "$OUT"
printf "tool\tversion\tsource\n" > "$OUT"

# arcasHLA does not implement --version; try `arcasHLA version` then fall back
arcas_ver="$(arcasHLA version 2>&1 | head -1 || true)"
if [ -z "$arcas_ver" ] || echo "$arcas_ver" | grep -qi usage; then
  arcas_ver="$(conda list -n thoracictcr arcas-hla 2>/dev/null | awk '/^arcas-hla/ {print $2; exit}')"
fi
emit "arcasHLA"   "${arcas_ver:-installed (version not reported)}"  "conda thoracictcr"
emit "TRUST4"     "$(safe_run run-trust4 -h | head -1)"          "conda thoracictcr"
emit "samtools"   "$(safe_run samtools --version)"               "conda thoracictcr"
emit "gdc-client" "$(safe_run /home/bio1/biosoft/gdc-client --version)" "local install"
emit "aria2c"     "$(safe_run aria2c --version)"                  "conda thoracictcr"
emit "Python"     "$(safe_run python --version)"                  "conda thoracictcr"

# pip packages
for pkg in pandas numpy scipy scikit-learn requests lifelines lightgbm biopython; do
  ver=$($PY -m pip show "$pkg" 2>/dev/null | awk -F': ' '/^Version:/ {print $2}')
  emit "py:$pkg" "${ver:-not-installed}" "pip (conda thoracictcr)"
done

# R packages
$RBIN -e '
pkgs <- c("ggplot2","ggsci","patchwork","dplyr","tidyr","readr","knitr",
          "purrr","scales")
out <- file("/home/bio1/workdata1/WZ/XXG2/ideaB/paper/tables/_rver.tsv", "w")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error=function(e) "missing")
  cat(sprintf("r:%s\t%s\tR (system Rscript)\n", p, v), file=out)
}
close(out)
' 2>/dev/null
if [ -f $REPO/paper/tables/_rver.tsv ]; then
  cat $REPO/paper/tables/_rver.tsv >> "$OUT"
  rm $REPO/paper/tables/_rver.tsv
fi

# OS & R version
emit "OS-kernel"  "$(uname -r)" "host"
emit "R"          "$($RBIN --version 2>&1 | head -1)" "system"

# Also write a kable-style .txt version
$RBIN -e '
df <- read.delim("'"$OUT"'", check.names=FALSE)
lines <- c("Supplementary Table S5. Software versions used in the ThoracicTCR pipeline.",
           "",
           knitr::kable(df, format="simple"))
writeLines(lines, "'"${OUT%.tsv}.txt"'")
' 2>/dev/null

echo "Wrote $OUT"
cat "$OUT"
