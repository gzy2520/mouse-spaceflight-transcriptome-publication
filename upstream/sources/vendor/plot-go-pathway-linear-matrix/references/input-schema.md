# Input schema

## RNA table

CSV or TSV with one row per group and stable gene identifier.

Required:

- `EnsemblID`: Ensembl gene ID. Version suffixes such as `.12` are removed.

Optional:

- `Group`: biological comparison or sample group. Defaults to `all`.
- `DisplayLabel`: gene symbol or short display name; never used for matching.
- `log2FoldChange`, VST, log-CPM, z-score, or another numeric column selected with `--effect-column`.
- `padj` or another adjusted-p-value column selected with `--padj-column`.

Each `Group + EnsemblID` pair must be unique. Aggregate biological replicates before running, using a design appropriate for the RNA experiment.

## Direct GO annotation table

CSV or TSV containing:

- `EnsemblID`: same species and namespace as the RNA table.
- `GO_ID`: direct GO identifier such as `GO:0006284`.
- `GO_Name`: optional display/audit field.

Provide all direct annotations for the RNA genes, not only terms previously
observed in another protein or gene set. The plotting script matches them to
the bundled frozen seven-pathway catalog. Canonical and GO alternative IDs
are supported; ontology expansion is performed when the catalog is built,
not separately for each input dataset.

## RNA value selection

| Goal | Effect column |
|---|---|
| Differential expression | log2FoldChange |
| Relative expression across samples | DESeq2 VST or edgeR log-CPM |
| Cross-sample pattern | gene-wise z-score derived from normalized expression |

If adjusted P values are supplied, cells with `padj` above the cutoff are rendered as unassigned RNA direction while retaining GO membership.

## Outputs

All matching and ordering use Ensembl ID. A gene is plotted only if it has at least one direct GO term assigned to the seven pathways. Terms assigned only to `Others` remain in the audit tables.

The summary figure reports pathway membership. The accompanying
`pathway_summary.csv` additionally records positive, negative, and
zero/non-significant RNA counts. RNA direction and magnitude are shown in the
linear matrix, not used to reorder genes.

Every run also copies the exact catalog to
`seven_pathway_go_catalog_used.csv` and records its ontology release and
SHA-256 hashes in `run_metadata.csv`.
