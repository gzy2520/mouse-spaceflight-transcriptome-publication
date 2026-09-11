# Mouse spaceflight transcriptome

The approved manuscript figures are in `release/Main/` and `release/Suppl/`.
This repository preserves the selected log2-expression version. Fig. 1a is
prepared separately. Gene selection and analysis use stable Ensembl identifiers;
gene symbols are display labels. Random seed is 25.

## Rebuild figures and supplementary tables

The publication renderer starts from the versioned numerical inputs. It writes
only to a new output directory, including the workbook exports.

```sh
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv"); renv::restore(prompt=FALSE); renv::activate()'
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
PYTHON="$PWD/.venv/bin/python" Rscript workflow/run_final_release.R reproduced_results/final
```

Python 3.12 or later is required; Python 3.14.7 was used for verification.
R 4.4.3 and the versions in `renv.lock` were used for the approved figures.
Use the same fonts and graphics libraries for pixel-identical rendering; PDF
metadata such as creation dates can differ. Python dependencies are pinned in
`requirements.txt`.

## Outputs and statistics

- `Main/`: 16 panels, each PNG and PDF. Fig. 4b–d and 5b–e show log2 expression;
  Fig. 4a/5a are tissue expression heatmaps and Fig. 4e/5f are meta-log2FC dotplots.
- `Suppl/`: 9 panels and 26 individual tissue correlation panels, each PNG/PDF.
  Fig. S1 is Hallmark enrichment. Fig. S6 is sample-log2FC Spearman correlation,
  with equal-weight Fisher-z pooling, not a meta-log2FC effect-size plot.
- `Supplementary_Tables/`: Tables S1–S3 and S5–S8, plus the supplied seven-sheet
  `2way_anova_7_pathways.xlsx` used for significance.
- `tables/log/`: regenerated expression summaries and sample values.
- `provenance/`: figure/source mapping, sample scope, tissue order and runtime.

The old `Table_S4_Repair_Genes_Expression_and_ANOVA.xlsx` was not used for the
paper's significance analysis. The supplied `2way_anova_7_pathways.xlsx` is an
external software result, preserved verbatim under `data/external_statistics/`.
The renderer neither recalculates it nor substitutes the obsolete one-way ANOVA.
The external result alone does not encode a reproducible software analysis setup.

For the complete source-data workflow, follow [`upstream/README.md`](upstream/README.md).

The frozen `results/tables/` CSVs are numerical inputs to the publication package.
Copying them is not an upstream analysis rerun. See `upstream/README.md` for the
acquisition and analysis source chain and its validation status.

The complete source reconstruction passed: 51 PNGs match the approved figures
pixel-for-pixel; all seven analysis workbooks agree within floating-point
precision. See [validation evidence](provenance/cleanup_report.md).

## Repository layout

`R/` and `python/` contain the selected renderers and metadata preparation.
`workflow/` contains the publication entry point, workbook exporter and gallery.
`data/publication_input/` contains the exact publication inputs.
`release/` preserves the approved outputs. Historical alternatives and removed
scripts remain recoverable from Git history. The retained Word document is an
editorial manuscript draft, not executable analysis documentation.

The publication input scope is 761 sample columns across 48 accessions and
26 tissues. Sample columns are not automatically unique animals. Fig. 2 uses
360 flight sample observations; tissues with fewer than five flight samples are
shown individually but excluded from the cross-tissue correlation pooling.
