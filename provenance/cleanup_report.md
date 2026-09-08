# Publication cleanup — 2026-09-08

Baseline: `0cfd9d7` (744 tracked files, 387.00 MiB). The final publication
entry point is `workflow/run_final_release.R`; selected outputs live only in
`release/`. Git history is retained, not rewritten or purged.

Final tracked package: **321 files, 106.07 MiB**, a **72.6% size reduction**
from the baseline (excluding Git objects, the local Python environment and
the recovery archive). Scripts/helpers/tests decreased from 35 to 19 files.
The clean-directory build and subsequent independent validation both passed:
`FINAL_RELEASE_BUILD_PASS` and `FINAL_RELEASE_CONTRACT_PASS`.
The complete output SHA-256 manifest and all gallery links were verified.

## Removed from the working release

- Alternative Style A/C render pipelines, their validators and result copies.
- Superseded base, log2FC-only, integrated and acceptance result directories.
- Old acceptance ZIP and pathway ZIP; obsolete logarithmic boxplot output.
- Duplicate panel exports under three alternative filenames and multiple folders.
- Original figure renderers replaced by the selected final renderers.
- Legacy NES-correlation Fig. 2 rendering (replaced by sample-log2FC correlations).
- Unreachable Style A/C branches in the selected Fig. S5 renderer, unused
  Hallmark helper functions and the unused Chinese tissue-label dictionary.
- Hard-coded local input fallback and automatic writes into historical results.

Frozen publication inputs, all 29 frozen source tables, source metadata audits,
environment locks, manuscript draft and metadata-refresh utilities are retained.
No upstream raw-data directory was modified.

## Reproducibility and numerical checks

- A fresh build renders 16 main panels, 9 supplementary panels and 26 individual
  tissue panels: 51 PNG/PDF pairs, with a local HTML review gallery.
- The original 27 manuscript tables are checked against their frozen SHA-256.
- Linear boxplot selection, gene ANOVA and summary grouping use Ensembl IDs.
  Symbols remain display labels; the summary now also exports EnsemblID.
- Both pathway and gene ANOVAs, tissue summaries (shared columns) and sample
  expression values match the pre-cleanup linear output at tolerance 1e-12.
- Tissue clustering now keys its matrix by EnsemblID. All tree merges, leaf
  order and heights agree with the approved baseline at tolerance 1e-12.
- The validator independently recomputes each tissue's Spearman matrix from
  frozen sample-term log2FC values and compares it at tolerance 1e-10.
- Stable-ID membership, GO scope, mouse metadata, palettes, dotplot scope,
  Fig. S5 scope/seed and the complete PNG/PDF file set are checked.
- Input and output checksums, runtime and source commit are exported. Output
  manifest paths are relative and resolve within the delivered release.

The figure gallery and main contact sheet support visual review. These are
reproductions from frozen analysis-level inputs; this cleanup does not claim
to rerun raw reads, enrichment, normalization or upstream meta-analysis.
The existing nominal ANOVA model and sample-unit assumptions are preserved.

## Recovery

Removed material is recoverable from Git at baseline `0cfd9d7`. A local recovery
copy (including untracked folder metadata) was also moved outside this worktree
to `../../99_archive_legacy/publication-cleanup-20260908.9mqI4K/`.
The archive is excluded from the publication package. Git checkout size is
reduced; Git historical object storage is intentionally unchanged.

## Remaining author decisions

The original methods Word document remains a draft; reconcile its historical
script references and sample/animal wording before submission. License and
citation decisions remain with the authors (see PUBLICATION_CHECKLIST.md).
