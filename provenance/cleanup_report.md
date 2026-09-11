# Publication cleanup and source reconstruction, 2026-09-11

Baseline: commit 67bf6ab. Before editing, all 102 approved desktop PNG/PDF
files matched the release by SHA-256. Approved figure bytes remain unchanged.

## Completed reconstruction

The isolated acquisition and analysis workflow was executed using all 49
hash-verified NASA processed inputs and the fixed GO/MGI and MSigDB resources.
All computation stages completed across isolated local/server directories.
Rendering used newly computed inputs, with no fallback to frozen numerical
publication results. The source boundary is NASA processed differential-expression
and normalized-expression CSVs, not FASTQ alignment or refitting NASA models.

The final comparison against the desktop manuscript package found:

- 51 of 51 PNGs pixel-identical; all 51 PDFs generated. PDF byte identity is
  not asserted because creation metadata can differ.
- Seven analysis workbooks (S1–S3, S5–S8): all sheet names, dimensions and
  97,577 cells checked. Six workbooks have exact cell equality. S5 has 27
  floating-point differences, maximum absolute difference 2.9976021664879227e-15;
  all other values and text agree. Comparison tolerance: relative 1e-12,
  absolute 1e-14.
- The supplied seven-sheet `2way_anova_7_pathways.xlsx` was byte-identical
  during reconstruction validation. Subsequently only its internal theme names
  were normalized; all other archive members remain byte-identical.
  It remains external software output; neither an inferred replacement model
  nor the unused old S4/one-way ANOVA is represented as the paper's statistics.
- Of 57 exported numerical/provenance artifacts, 49 are byte-identical.
  Remaining differences were checked: floating-point serialization, gzip
  headers, source hashes, local/server paths and refreshed metadata timestamps.
  Sample information agrees. The final Eye mission label is preserved through
  the explicit label override without changing numerical statistics.

Acquisition checks verified the manifest hashes of all locally reused sources,
current NASA API availability of 40 exact filenames, and one complete current
NASA download against its original hash. The nine unavailable historical names
are preserved in the GitHub release asset `reproduction-inputs-20260911`;
its remotely reported SHA-256 and size match the local archive. This is not a
claim that all 49 files were freshly downloaded from NASA during validation.

## Scope of cleanup

Removed the obsolete linear-expression renderer, duplicate log-expression
results, historical document editing scripts, repeated frozen-output validator,
unused one-way ANOVA calculations/exports and non-manuscript global GSEA plots.
Added the previously missing acquisition and upstream analysis sources, pinned
resources and portable figure/table entry points. Essential input checks remain.
Global descriptive NES skips unused adaptive P-value refinement; targeted GO
inference retains the original calculation. A 5,254-pathway comparison confirmed
NES equivalence to floating-point precision and identical leading-edge genes.

Frozen-input and full-source builds both passed. The detailed final comparison
is saved in `provenance/reproduction_validation.json`. Validation scripts and
intermediate outputs stay outside the versioned publication source package.
The original working checkout and unrelated exploration were not modified.
Removed tracked content remains recoverable from Git history.
