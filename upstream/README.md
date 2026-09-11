# Acquisition and upstream analysis

This directory reconstructs the manuscript's numerical inputs from NASA's
published processed differential-expression CSVs, which also contain the sample
normalized-expression columns. It does not download FASTQ files, realign reads,
or refit the original differential-expression models. The separate publication
renderer then produces the selected manuscript figures and supplementary tables.

Analysis joins and gene membership use stable Ensembl IDs. Reviewed component
maps retain MGI/Entrez identifiers and symbols for display. Existing arithmetic,
contrast orientation, filtering and categorical orders are preserved; the random
seed is 25 where the original analysis uses randomness.

## Inputs and provenance

- `nasa_processed_files.csv`: 49 exact processed files, including the contrasts
  used by the study registry; each has an OSDR files-API endpoint, workspace path,
  size and SHA-256. The publication scope is narrower than the complete registry.
- `annotation_snapshots.csv`: five historical annotation files. The exact copies
  are bundled in `../data/annotation_snapshot_20260822.tar.gz`. Acquisition checks
  its member inventory and each file's SHA-256. Current GO/MGI downloads are not
  substitutes for this historical snapshot.
- `msigdb_snapshot.csv`: the exact MSigDB 2026.1.Mm mouse Ensembl membership
  snapshot. Acquisition looks for its basename under `../data/` or for the
  specified `local_source_path` under `--local-source`. Without that snapshot,
  stage 03 retains the original msigdbr route; do not assume a current database
  download reproduces the approved membership.
- `design_inputs/` and `design_manifest.csv`: 13 curated study, contrast, sample,
  tissue-order and reviewed component-mapping inputs, including the teacher's
  component workbook. These are fixed design inputs, not regenerated analysis
  results, even where their historical paths contain `03_analysis_results`.
- `sources/` and `source_manifest.csv`: original computation scripts and the
  vendored seven-pathway catalog/matrix implementation. Original and packaged
  hashes record portability edits separately from unchanged sources.
- `approved_artifact_map.csv`: publication-relative destinations, expected
  approved hashes and possible upstream producers. The existing `status` column
  describes the historical source search, not a successful fresh reconstruction.
  Fresh export results are recorded separately in `logs/export_comparison.csv`.

## Environment

The publication `renv.lock` covers the selected renderer, not the complete
upstream Bioconductor environment. Upstream execution also requires `fgsea`,
`Biobase`, `yarn` (and its qsmooth dependencies), `readxl`, `ggrepel`, `jsonlite`
and `ragg`, in addition to the renderer's `data.table`, `ggplot2`, `patchwork`,
`scales` and `digest`. `msigdbr` is needed for the database-fetch fallback.
Base/recommended R packages such as `grid` and `grDevices` come with R.
Use a dedicated R library and pass its absolute path with `--r-library`;
`run.py` does not install packages. Record the actual installed versions before
claiming a reproducible upstream environment.

Use the repository Python environment with `../requirements.txt`; upstream
scripts require SciPy as well as NumPy, pandas and matplotlib. The current
requirements include SciPy 1.18.1. R scripts also use `shasum` and Cairo graphics;
matching fonts and graphics-library versions matter for identical image pixels.

## Isolated execution

Run these commands from the repository root. Choose a new work directory;
`init` refuses an existing directory. Never use the original source checkout as
`--work`. The `--local-source` directory is read only to the acquisition command.

```sh
python3 upstream/run.py plan
python3 upstream/run.py init --work reproduced_results/upstream_new
python3 upstream/run.py acquire --work reproduced_results/upstream_new --download
.venv/bin/python upstream/run.py run --work reproduced_results/upstream_new \
  --python "$PWD/.venv/bin/python" --r-library /absolute/path/to/upstream-R-library
python3 upstream/run.py export --work reproduced_results/upstream_new
```

For an archival local acquisition, replace the acquire command with:

```sh
python3 upstream/run.py acquire --work reproduced_results/upstream_new \
  --local-source /absolute/path/to/archival-source
```

`--download` may also be combined with `--local-source`. Available local files
are checksum-verified before copying; NASA files are downloaded when needed.
The five annotation snapshots come from the bundled archive. The MSigDB snapshot
must also be available to guarantee the frozen database route. Stages 18–19 use
NASA metadata independently of processed-data acquisition: stage 18 refreshes
live dataset/ISA metadata and therefore still needs network access in a local
acquisition run.

`--from-stage` and `--until` select an inclusive stage interval, using names from
`plan`. A resumed interval requires the earlier products to exist. Logs and
completion markers are written under the isolated workspace's `logs/`; a marker
means that command exited successfully, not that its outputs match the manuscript.
`init` copies packaged sources once: later package fixes are not automatically
propagated into an already initialized run directory. Coordinate source updates
with the runner before resuming an active local or server run.

## Computation chain

| Stage | Inputs and products needed downstream |
| --- | --- |
| 01–03 | Historical GO/MGI/Ensembl snapshots and curated decisions → concrete GO membership, validated stable-ID registry, full gene statistics and frozen MSigDB membership. |
| 04 | Registry statistics and gene sets → global GSEA, concrete GO unit/tissue statistics and NES matrices. |
| 05–06 | Stage 04 → GO correlations and ordered matrices; global descriptive term selection used by Hallmark Fig. S1. |
| 07–08 | Processed data, sample registry and ontology catalog → 26-tissue RNA meta tables and seven-pathway membership. The two catalog runs intentionally use different approved relation/rule contracts. |
| 09 | Stage 07 sample scope and stage 08 catalog → per-source-sample seven-pathway counts and repertoire tables. |
| 10–12 | Stage 07 RNA meta tables and reviewed component maps → DSB/all-component/SSB meta-log2FC tables. |
| 13–14 | Processed expression, sample scope and reviewed components → YARN/qsmooth values, tissue order and the selected 74-component output after Neil2 exclusion. Stage 13 already writes the `18_essential_components_yarn_qsmooth_spearman_reordered_20260827` directory required by stage 14. |
| 15–16 | Concrete GO membership and stable-ID statistics → two-GO DDR membership and four common-direction group tables. Helper modules are packaged beside their callers. |
| 17 | Sample-level processed expression, stage 07 tissue mapping and stage 04 NES → dataset-presence membership, exact regions and set summaries. |
| UpSet, after 17 | `plot_teacher_selected_go_upsets_20260829.R` → the two `intersection_counts_*_tissues.csv` tables and `upset_render_audit.csv` from stage 17 membership. These export-mapped products need this additional original script. |
| 18–19 | Stage 09 sample scope plus live OSDR ISA/dataset metadata → mouse metadata and Fig. 6 complete audit. |
| 20 | Stage 18 metadata, stage 04 GO membership/statistics and stage 05 legacy correlation bridge → within-tissue sample-log2FC Spearman matrices and equal-weight Fisher-z pooling. |

The vendored UpSet script was added during the dependency audit; its runner
registration is being coordinated separately. Check `plan` for its presence
before running the full sequence. Its default stage-17 input is
`03_analysis_results/19_ddr_go0006974_dataset_level_venn_20260826` and its output
is `03_analysis_results/24_teacher_selected_figure_release_20260829/go_upsets`.

## Export and publication boundary

`export` writes a candidate under the work directory, never into the approved
publication inputs. Byte-identical candidates are reported as such. Changed or
missing candidates require review; gzip headers, run paths, metadata download
times and fresh numerical results must be distinguished in a semantic comparison.
Do not update approved hashes merely to make a candidate pass.

The exporter does not itself invoke the final publication renderer or install a
candidate as its input. The separate renderer in `../workflow/run_final_release.R`
currently reads the repository's approved numerical inputs. Candidate installation
and final figure/table comparison therefore remain an explicit integration step;
running that renderer alone does not establish acquisition-to-figure equivalence.
The externally supplied significance workbook is retained verbatim and is not
recomputed by this upstream chain.

## Inspection status: 2026-09-10

The independent stage 04–20 inspection checked packaged helper scripts, principal
input/output bridges and acquisition coverage. All 49 distinct processed paths
in the dataset manifest, all 48 paths in the sample-expression audit and all three
paths in the four-thymus audit are covered by the acquisition manifest. Packaged
R/Python sources parse without execution. The missing original UpSet producer
has been added unchanged, with its source hash recorded.

The inspection did not execute stages 04–20 or the full pipeline. At handoff,
the runner reported stages 01–02 complete and stage 03 in progress. Separate
publication-render checks reported 51 matching PNGs and seven cell-identical
workbooks; those checks do not validate a fresh upstream reconstruction.

Two export mappings need attention from the runner/integration owner: the GO
`tissue_statistics` entry and the Fig. 6 complete metadata audit have empty
`source_candidates` in the inspected map. The former is produced by stage 04;
the latter by stage 19 under `publication_release/results/tables/`.
A read-only comparison of the original local GO source and current publication
input found 390 rows with identical headers and all cells identical except 15
Eye `missions` labels: `SpaceX-4;SpaceX-8+SpaceX-9` in the original and
`SpaceX-4;SpaceX-8` in the publication input. A documented label bridge or explicit
review is needed; this observation does not authorize replacing the approved
baseline hash. The inspection left `run.py`, the artifact map and approved
publication outputs unchanged.
