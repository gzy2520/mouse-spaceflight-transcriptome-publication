# Mouse spaceflight transcriptome: publication figures and tables

This repository is the compact, manuscript-facing release of the mouse spaceflight transcriptome analysis. It reproduces the figures and tables selected in `24_teacher_selected_figure_release_20260829` from frozen publication inputs. The earlier exploratory analyses, caches, repeated exports and superseded figure versions are intentionally excluded.

The checked-in files under `results/` are the approved reference outputs. A reproduction run always writes to a separate directory and never overwrites them.

The Fig. S5 PNG is freshly rendered during a run. Its companion PDF is copied byte-for-byte from the approved vector export because Cairo writes the current creation time into otherwise identical PDF content; keeping that file immutable makes the release archive and its checksum reproducible.

## Repository contents

```text
data/publication_input/    Frozen, figure-level input tables
R/data/                    Source-refresh and metadata-audit scripts
R/figures/                 R scripts for the R-rendered figures
python/                    Matplotlib script retained for four heatmaps
workflow/run_pipeline.R    One-command reproduction entry point
tests/                     File, hash, stable-ID and analysis-scope checks
results/                   Approved figures and manuscript tables
provenance/                Input/output hashes and figure-source mapping
```

The release contains 9 main-figure PNGs (including the mouse-composition overview Fig. 6), one main-figure PDF vector export, 8 supplementary PNGs, the PDF version of Fig. S5, and 28 CSV tables. The expected 47-file set and SHA-256 values are recorded in `provenance/reference_sha256.csv`.

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

## Figure 6 display convention

`Main/Fig_6_mouse_metadata_bubble.png` (and its PDF companion) is a sample-design
overview. The x-axis follows the 26-tissue order in the cross-tissue GO table.
Each mission occupies two rows (female, then male); the mission name is printed
between the paired rows. Age is encoded by a fixed, ordered discrete palette with
no continuous colour interpolation, while point area is the number of unique mice
in that mission–tissue–sex–age stratum. Flight and Ground samples are both included;
their separate counts remain in `results/tables/Fig_6_mouse_metadata_grouped.csv`.

The Fig. 6 metadata was refreshed directly from all 48 NASA OSDR ISA metadata
ZIPs on 2026-08-30. The frozen publication input contains 1,588 complete ISA
source/sample rows and a one-to-one audited join for all 761 selected analysis
sample columns. OSD-162 is assigned to the current dataset-API mission
`SpaceX-8`; its legacy ISA value `SpaceX-3` is retained only in the conflict
audit. Age ranges are shown only when the downloaded ISA cell is itself a range.
Exact source values such as `32 week` are not converted to months.

To refresh the source snapshot from NASA and rebuild the Fig. 6 input tables:

```bash
Rscript R/data/01_refresh_osdr_mouse_metadata.R
```

This network-dependent refresh is intentionally separate from the offline
publication-output reproduction command. URLs, download time, sizes and SHA-256
values are recorded in
`data/publication_input/mouse_metadata/02_osdr_download_manifest.csv`.

## Conservative visual refresh (Style A)

The original `results/` directory is the approved reference set. A visual-only refresh can be rendered separately with:

```bash
Rscript workflow/run_style_A.R results_style_A_rebuilt
```

Style A keeps the same plotting inputs, filters, ordering, membership counts and numeric matrices while applying a consistent publication palette, typography, spacing and high-resolution export. It writes the same 44-file set to the chosen output directory, never overwrites `results/`, and ends with `STYLE_A_CONTRACT_PASS`.

## Narrative figure variant (Style C)

The C-story variant is an independent, higher-level manuscript layout. It is
rendered with:

```bash
Rscript workflow/run_style_C.R results_style_C_rebuilt
```

It keeps the frozen analysis inputs and publication tables unchanged, but uses
functional-block facets for Fig. 1, a single lower-triangle display for Fig. 2,
UpSet geometry for the four- and five-set membership panels, dot matrices for
sparse meta-log2FC displays, and pathway-faceted small multiples for Fig. S5.
The output contains 16 PNG figures, matching PDF figures, and the 27 copied
publication tables. The complete figure mapping and display-only changes are
recorded in `provenance/style_C_figure_manifest.csv`.

Validate an existing C-story directory with:

```bash
Rscript tests/validate_style_C_contract.R results_style_C_rebuilt
```

`results/` remains the approved reference set; Style C is a candidate visual
variant until the manuscript team selects it.

## Figure 1 display convention

The checked-in `Main/Fig_1.png` follows the teacher-requested heatmap presentation: the frozen 26 × 15 matrix is unchanged, in-cell NES numbers and the decorative top bar are omitted, and a fine continuous blue–white–red scale carries the value detail. Exact NES values remain in the publication tables. Figure 2 is the explicit exception and keeps its Spearman numbers.

## Analysis conventions retained in this release

- Ensembl gene IDs are the analysis keys. Gene symbols are retained only as display labels.
- The tissue analyses cover 26 mouse tissues. The seven-pathway source-sample table contains 761 source sample columns across seven pathways.
- Fig. 4b and Fig. 5b use the median of mission-level median log2FC values; their displayed nominal tissue-meta P values use signed Stouffer aggregation.
- The YARN/qsmooth heatmaps retain the approved normalized matrix values. Their dendrogram order is computed from the corresponding frozen component matrix.
- Neil2 (`ENSMUSG00000035121`) is excluded from the final 45-component SSB display and is retained only in the dedicated exclusion audit table.
- A positive GSEA NES indicates enrichment toward the spaceflight-ranked end of the gene list; it is not described as pathway activation.

This is a direct publication-output reproduction package, not a replacement for the archived raw-data-to-analysis workflow. `provenance/figure_sources.csv` records the script and primary frozen input used for every figure.
