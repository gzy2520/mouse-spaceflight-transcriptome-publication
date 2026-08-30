# Mouse spaceflight transcriptome: publication figures and tables

This repository is the compact, manuscript-facing release of the mouse spaceflight transcriptome analysis. It reproduces the figures and tables selected in `24_teacher_selected_figure_release_20260829` from frozen publication inputs. The earlier exploratory analyses, caches, repeated exports and superseded figure versions are intentionally excluded.

The checked-in files under `results/` are the approved reference outputs. A reproduction run always writes to a separate directory and never overwrites them.

The Fig. S5 PNG is freshly rendered during a run. Its companion PDF is copied byte-for-byte from the approved vector export because Cairo writes the current creation time into otherwise identical PDF content; keeping that file immutable makes the release archive and its checksum reproducible.

## Repository contents

```text
data/publication_input/    Frozen, figure-level input tables
R/figures/                 R scripts for the R-rendered figures
python/                    Matplotlib script retained for four heatmaps
workflow/run_pipeline.R    One-command reproduction entry point
tests/                     File, hash, stable-ID and analysis-scope checks
results/                   Approved figures and manuscript tables
provenance/                Input/output hashes and figure-source mapping
```

The release contains 8 main-figure PNGs, 8 supplementary PNGs, the PDF version of Fig. S5, and 27 CSV tables. The expected 44-file set and SHA-256 values are recorded in `provenance/reference_sha256.csv`.

## Reproduce the release

R 4.4.3 and Python 3.13 were used for the approved render. The exact R package versions are recorded in `renv.lock`; the Python packages are pinned in `requirements.txt`.

```bash
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")'
Rscript -e 'renv::restore()'
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
PYTHON=.venv/bin/python Rscript workflow/run_pipeline.R
```

The command creates `reproduced_results/` and then checks every file against the approved release. To validate an existing directory without rendering again:

```bash
Rscript tests/validate_publication_contract.R results
```

The validator checks the exact file list and SHA-256 values, as well as the main stable-ID and scope contracts used by the manuscript tables. A successful run ends with `PUBLICATION_CONTRACT_PASS`.

## Analysis conventions retained in this release

- Ensembl gene IDs are the analysis keys. Gene symbols are retained only as display labels.
- The tissue analyses cover 26 mouse tissues. The seven-pathway source-sample table contains 761 source sample columns across seven pathways.
- Fig. 4b and Fig. 5b use the median of mission-level median log2FC values; their displayed nominal tissue-meta P values use signed Stouffer aggregation.
- The YARN/qsmooth heatmaps retain the approved normalized matrix values. Their dendrogram order is computed from the corresponding frozen component matrix.
- Neil2 (`ENSMUSG00000035121`) is excluded from the final 45-component SSB display and is retained only in the dedicated exclusion audit table.
- A positive GSEA NES indicates enrichment toward the spaceflight-ranked end of the gene list; it is not described as pathway activation.

This is a direct publication-output reproduction package, not a replacement for the archived raw-data-to-analysis workflow. `provenance/figure_sources.csv` records the script and primary frozen input used for every figure.
