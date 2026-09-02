# Revised Fig. 2: sample-level log2FC GO-term relationships

This dated release is an isolated revision of Fig. 2. The approved
`final_result/` directory is unchanged. Fig. 1c and every other figure are not
part of this revision.

## Analysis contract

- Observation unit: Flight sample (mouse) within each tissue.
- Input effect: `log2(flight expression) - mean(log2(matched Ground expression))`,
  using the same OSD dataset and tissue for the Ground reference. Source
  matrices already contain a positive offset; no extra pseudocount is added.
- Term value: arithmetic mean of member-gene log2FC values. All analysis keys
  are Ensembl IDs; symbols are not used for matching or correlation.
- Correlation: first 15 × 15 Spearman matrices within each tissue, then an
  equal-weight Fisher-z average across tissues with at least five Flight
  samples. Twenty-three tissues enter the overall matrix; all 26 remain in the
  per-tissue outputs.
- Display: the existing Fig. 2 blue–white–red palette (`#3B6FB6`, `#FFFFFF`,
  `#C65A5A`) is retained. The fourth label is “Intrinsic apoptotic signaling”
  and the seventh is “Telomeric region”; only display wording changes.

## Contents

`Main/Fig_2.png` is the revised lower-triangle overall matrix.
`Suppl/Fig_S6_log2fc_per_tissue_merged.png` is the 26-tissue small-multiple
figure, and `Suppl/Fig_S6_per_tissue/` contains one full-size PNG per tissue.
The frozen calculation tables are in
`../data/publication_input/go/log2fc_spearman_20260902/`.

## Reproduce

From this repository root, with Python packages in `requirements.txt`:

```bash
bash workflow/run_log2fc_per_tissue_pathway_spearman_20260902.sh
```

With no argument the command writes a fresh
`reproduced_results/log2fc_spearman_20260902/` directory, leaving this
checked-in snapshot untouched. Pass a relative or absolute output directory to
choose another destination; an existing `Main/Fig_2.png` is never overwritten.
The renderer creates all single-tissue PNG/PDF files in a temporary directory,
then assembles the compact `Main/` and `Suppl/` tree. The table-level
validation summary and figure manifest are copied to `provenance/`.

For compact bottom labels (GO IDs only on the horizontal axes and merged term
key), pass `go_id` as the second argument:

```bash
bash workflow/run_log2fc_per_tissue_pathway_spearman_20260902.sh \
  reproduced_results/log2fc_spearman_go_id_20260902 go_id
```

The full raw-matrix-to-table builder is retained in the main project at
`07_scripts/06_cross_tissue_integration/build_log2fc_per_tissue_pathway_spearman_20260902.py`;
this compact publication repository freezes its resulting tables so a reader
can reproduce the released figures without hidden local paths.
