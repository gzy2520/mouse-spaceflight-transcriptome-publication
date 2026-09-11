# Publication cleanup, 2026-09-10

The baseline is commit 67bf6ab. All 102 desktop PNG/PDF files matched its release
files by SHA-256 before editing. The selected renderer set was rebuilt into a
new directory after removing the unused one-way ANOVA computations and the
implicit writes to duplicate result directories. All 51 PNGs are pixel-identical
to the approved desktop reference; all 51 PDFs were generated.

The portable workbook exporter reproduces all sheet names, dimensions and cells
of Tables S1–S3 and S5–S8. The supplied `2way_anova_7_pathways.xlsx` is copied
byte-for-byte. The unused old S4 and its one-way ANOVA are not claimed as the
paper's statistics. No new statistical model has been substituted.

Removed: obsolete linear-expression renderer, duplicate log-expression results,
historical document editing scripts, repeated frozen-output validator, and
unused one-way ANOVA output files. Approved figure bytes were preserved.
The source checkout's unrelated staged/untracked changes were not modified.
Removed tracked content remains available in Git history.

The drawing entry point and workbook exporter write only into the requested new
destination. Upstream reconstruction is separate; its scope and run status are
documented under `upstream/`. A successful frozen-input figure rebuild is not
represented as evidence of an upstream reanalysis.
