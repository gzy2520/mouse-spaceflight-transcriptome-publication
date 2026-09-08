# Mouse spaceflight transcriptome: publication release

Open `release/index.html` for the clickable figure review gallery.
The single final output is `release/`. Figures are freshly rendered from frozen
publication inputs, with the latest manuscript numbering: Fig. 4b–d and 5b–e
are linear expression boxplots; Fig. 4e and 5f are meta-log2FC dotplots.
Fig. 1a is prepared separately and is intentionally absent.

## Rebuild

R 4.4.3 was used. R dependencies are recorded in `renv.lock`; Python dependencies
are pinned in `requirements.txt`.

```sh
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv"); renv::restore()'
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
PYTHON="$PWD/.venv/bin/python" Rscript workflow/run_final_release.R reproduced_results/check
Rscript tests/validate_final_release.R reproduced_results/check full
```

The destination must not exist. The build does not overwrite earlier results.
This is a frozen-input publication reproduction package: it regenerates figures
and the seven-pathway linear-expression summaries and ANOVAs. It does not rerun
raw-read processing, differential expression, GO enrichment or YARN normalization.

## Final outputs

- `release/Main/`: 16 main panels, each in PNG and PDF.
- `release/Suppl/`: 9 supplementary panels, each in PNG and PDF.
- `release/Suppl/Fig_S6_per_tissue/`: 26 tissue correlation panels, each in PNG and PDF.
- `release/tables/`: 27 unchanged manuscript tables.
- `release/tables/linear/`: four regenerated linear-expression statistical tables.
- `release/provenance/`: input/output SHA-256 manifests, plotting audits,
  sample metadata audit, runtime and source revision.

No descriptive-name duplicates, alternative style exports or obsolete log-scale
boxplots are included. Historical versions remain recoverable through Git.
See `provenance/cleanup_report.md` for the cleanup scope and verification.

## Source layout

`data/publication_input/` contains frozen numerical inputs and source audits.
`results/tables/` contains 29 frozen source tables: 27 manuscript tables and
two metadata audits. These are retained inputs, not disposable output.
`R/final_figures/`, the two retained scripts in `R/figures/` and the two Python
renderers implement the selected figures. The retained Hallmark style helper
supports its approved appearance. `R/data/` contains optional network-dependent
metadata-refresh scripts and is separate from the offline release build.

## Analysis conventions

Ensembl IDs are used for gene selection, grouping and clustering. Symbols are
display labels. Random seed is 25. The metadata scope is 761 sample columns
across 48 accessions and 26 tissues; sample columns must not automatically be
interpreted as unique animals.

Fig. 1c displays the frozen 26 × 15 mission-equal mean NES matrix. Positive NES
means enrichment toward the flight-ranked end, not pathway activation.
Fig. 2 and S6 use ground-referenced sample log2FC, within-tissue Spearman
correlations and equal-weight Fisher-z pooling across eligible tissues.

Fig. 4a and 5a both display YARN qsmooth log2 expression with average-linkage
tissue trees based on 1 − Spearman rho. Fig. 4e and 5f display tissue meta-log2FC
(median of mission-level medians) with nominal signed-Stouffer P < 0.05.
Neil2 (ENSMUSG00000035121) is excluded from the 45-component SSB display.

Linear boxplots use 2^YARNNormalizedLog2 − 1. The central line is the mean;
box edges are Q1/Q3, with the existing 1.5-IQR whisker convention. Gene and
pathway one-way ANOVAs describe tissue heterogeneity. Their nominal P values
are retained; the cleanup does not change the statistical model.

The existing methods Word document is retained as a manuscript draft; its
historical script references and sample-unit wording need editorial reconciliation
before submission. Author-selected licensing and citation metadata are tracked
in `provenance/PUBLICATION_CHECKLIST.md`.
