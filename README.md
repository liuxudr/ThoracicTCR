# ThoracicTCR

A unified command-line toolkit for pan-thoracic T-cell receptor repertoire analysis from bulk RNA-Seq.

Reconstructs paired TCR repertoires from fastq or BAM, computes diversity / clonality / V-J / clone-size / motif metrics, builds a cross-cohort atlas, and renders publication-ready figures — all from one `thoracictcr` command.

Cancers supported: **LUAD, LUSC, ESCA, MESO, THYM** (≈ 1,400 samples).

## Installation

Requires Python ≥ 3.10 and an external TRUST4 install (provided by conda).

```bash
# 1. Create the conda environment (Python + R + TRUST4 + arcasHLA + aria2 …)
conda env create -f environment.yml
conda activate thoracictcr

# 2. Install the CLI
pip install -e ".[dev]"

# 3. Fetch TRUST4 reference files (~1 MB)
thoracictcr trust4 fetch-refs
```

Verify the install:

```bash
thoracictcr version
thoracictcr info env       # Python + key package versions
thoracictcr info paths     # data paths resolved from configs/paths.yaml
thoracictcr trust4 check   # check that run-trust4 + refs are visible
```

## Command map

```
thoracictcr
├── version                                        # print version
├── info        env | paths | versions             # diagnostics
├── manifest    tcga | geo | atlas                 # build sample manifests
├── download    fastq | tcr-dbs                    # fetch fastq / VDJdb / McPAS / IEDB
├── trust4      fetch-refs | check | run-bam | run # TCR reconstruction
├── metrics     diversity | pilot-qc | vj-usage |  # repertoire metrics
│               clone-size | shared-clones |
│               rarefaction | gse-diagnostic
├── atlas       build | info                       # cross-cohort atlas
└── figures     list | render                      # publication figures (R)
```

Every node accepts `--help` for full options.

## Quickstart — open-access pilot (10 samples, ≈ 54 GB)

The default pilot uses four open-access GEO/SRA cohorts of NSCLC/ESCC bulk RNA-Seq
with anti-PD-1 response labels (GSE126044, GSE135222, GSE145370, GSE193258).
No controlled-access credentials required.

```bash
conda activate thoracictcr

# 1. Build the pilot manifest from ENA (no download, just metadata)
thoracictcr manifest geo
#   -> configs/pilot_geo_10.csv
#   -> configs/pilot_geo_10.urls.txt

# 2. Download pilot fastq (aria2c, parallel, resumable; preview first)
thoracictcr download fastq --dry-run
thoracictcr download fastq

# 3. Reconstruct TCR repertoires with TRUST4
thoracictcr trust4 run --mode geo

# 4. End-to-end QC: diversity table + warnings report
thoracictcr metrics pilot-qc
```

Step 4 prints `All pilot samples passed QC` if every sample produced parseable
productive TRB clones — that confirms the pipeline runs end-to-end on your machine.

## Scaling to the full atlas

```bash
# 5. Pick a stratified atlas cohort and emit the delta download list
thoracictcr manifest atlas

# 6. Download the additional fastq + run TRUST4 on the full atlas
thoracictcr download fastq --urls configs/atlas_download.urls.txt
thoracictcr trust4 run --mode geo

# 7. Assemble the atlas object + per-sample diversity table
thoracictcr atlas build

# 8. Compute repertoire metrics
thoracictcr metrics vj-usage
thoracictcr metrics clone-size
thoracictcr metrics shared-clones
thoracictcr metrics rarefaction

# 9. (Optional) Antigen lookup against public TCR databases
thoracictcr download tcr-dbs
thoracictcr atlas build      # re-run with VDJdb annotations

# 10. Render manuscript figures (requires R + ggplot2 + ggsci + patchwork)
thoracictcr figures render
```

All intermediate tables are written under `data/metrics/` and `data/atlas/`;
figures go to `paper/figures/`. Filenames and seeds are deterministic so the
pipeline is reproducible from raw fastq to final PDF.

## Controlled-access TCGA path

TCGA RNA-Seq BAMs are dbGaP-controlled; you need an approved project and a GDC
download token before running this path.

```bash
# Once $GDC_TOKEN is set
thoracictcr manifest tcga --access controlled
gdc-client download -m configs/pilot_tcga_thym_10.gdc_manifest.txt -t $GDC_TOKEN
thoracictcr trust4 run --mode tcga
thoracictcr metrics pilot-qc --trust4-dir data/repertoire/trust4/TCGA-THYM
```

## Configuration

All data paths resolve through [`configs/paths.yaml`](configs/paths.yaml).
Override per-machine by setting `THORACICTCR_PATHS_CONFIG` to point at a copy:

```bash
export THORACICTCR_PATHS_CONFIG=/path/to/my/paths.yaml
thoracictcr info paths    # confirm the resolved values
```

You can also override any path inline via the command's `--out-dir`,
`--trust4-dir`, `--manifest`, … flags — see `--help` on each command.

## Repository layout

```
src/thoracictcr/   Python package (CLI + library functions)
R/                 Figure scripts (ggplot2 + ggsci NPG palette)
configs/           YAML configuration (paths, cohorts, panels)
scripts/           External-tool drivers wrapped by the CLI
data/              (gitignored) raw inputs and pipeline outputs
paper/             Manuscript figures and tables
```

## Data availability

- **Open-access inputs** — GEO/SRA accessions GSE126044, GSE135222, GSE145370, GSE193258 (NSCLC / ESCC bulk RNA-Seq); VDJdb, McPAS-TCR, IEDB.
- **Controlled inputs** — TCGA RNA-Seq BAMs require dbGaP authorization for the relevant project. Not redistributable.
- **Pipeline outputs** — every diversity / VJ / clone-size / motif table is reproducible from the raw fastq via the commands above.

## Citing

If you use ThoracicTCR in your research, please cite the entry in [`CITATION.cff`](CITATION.cff).

## License

MIT — see [`LICENSE`](LICENSE).
