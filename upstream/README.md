# Reproduce from NASA processed data

The manuscript starts from NASA's published differential-expression CSV files,
which also contain sample normalized-expression columns. This workflow rebuilds
the stable-ID registry, GO membership, GSEA, pathway counts, expression
normalization, meta summaries, overlap groups and sample correlations from those
source files. It does not realign FASTQ or refit NASA's differential-expression
models. Ensembl IDs are analysis keys; the seed is 25.

## Run

Use R 4.4.3 and restore the repository's `renv.lock`, then install Python packages
from `requirements.txt`. Run from the repository root:

```sh
Rscript -e 'renv::restore(prompt=FALSE)'
publication_r_library=$(Rscript -e 'cat(renv::paths$library())')
.venv/bin/python upstream/run.py init --work reproduced_results/from_sources
.venv/bin/python upstream/run.py acquire --work reproduced_results/from_sources --download
.venv/bin/python upstream/run.py run --work reproduced_results/from_sources --r-library "$publication_r_library"
.venv/bin/python upstream/run.py render --work reproduced_results/from_sources --r-library "$publication_r_library"
```

The uncompressed processed inputs total approximately 4.64 GB.

`init` needs a new directory. `acquire` downloads the manifest's exact inputs;
`--local-source /path/to/original-data-layout` can instead reuse an existing copy
whose hashes match. `plan` lists the stages. `run --from-stage NAME --until NAME`
resumes a selected range after an interrupted run. Existing GSEA caches belong
to this newly computed workspace. A failed step names its log under `logs/`.

`render` exports the computed inputs into `candidate_publication/`, copies only
the presentation code and externally supplied statistics, and renders into
`reproduced_figures/`. It does not fill missing computed inputs from the frozen
publication tables. The original repository and approved `release/` are untouched.

## Fixed inputs

- `nasa_processed_files.csv` fixes 49 processed source paths and SHA-256 values.
  Nine historical names are absent from the current NASA files API. Their exact
  source bytes are preserved in the downloadable GitHub release asset described
  by `historical_snapshot.json`. The downloader verifies its SHA-256, then each
  source file. The remaining files are obtained from NASA.
- `annotation_snapshots.csv` identifies the five bundled historical GO/MGI files
  in `data/annotation_snapshot_20260822.tar.gz`.
- `msigdb_snapshot.csv` identifies the bundled MSigDB 2026.1.Mm Ensembl membership.
  The registry uses that fixed resource rather than silently refreshing MSigDB.
- `design_inputs/` contains the reviewed study/sample/contrast scope, component
  selection, tissue identities and mapping adjudications. These are design inputs,
  not substituted expression or differential-expression results.
- `approved_label_overrides.csv` preserves the final Eye `missions` display label.
  Its explicit before/after rule changes no GO statistic or numerical selection.
- The final two-way ANOVA workbook is supplied software output. Its bytes are
  preserved; recreating the software model is outside this code's scope.

## Computation

Stages 01–06 build annotation, the stable-ID registry and enrichment summaries.
Stages 07–14 build the seven-pathway membership/counts, component meta matrices
and YARN qsmooth expression. Stages 15–17b generate DDR/common-direction groups
and UpSet tables. Stages 18–20 reconstruct metadata and the sample-level
Spearman/Fisher-z analysis. All source scripts and required helpers are packaged
under `sources/`; `source_manifest.csv` records original and packaged hashes.

Targeted GO inference retains its original `fgseaMultilevel` calculation and
P-value precision. The global descriptive summaries consume NES only. Their
unused adaptive P-value refinement is omitted, and global P-value columns are
not exported. Initial enrichment-score/permutation normalization and directional
fallback are retained. A direct comparison across 5,254 GO-BP sets found identical
pathways and leading-edge genes, with NES differences below 6e-15 (CSV round-trip
precision). Final figure and workbook comparisons remain the decisive check.

`export` records source-to-candidate hashes. Gzip headers, paths, download times
and floating-point serialization can differ without changing analytical values;
it does not claim that every candidate is byte-identical. Missing computed inputs
are errors. Current reconstruction evidence is summarized in
`../provenance/cleanup_report.md`.
