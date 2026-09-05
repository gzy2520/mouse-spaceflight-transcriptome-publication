#!/usr/bin/env bash

# Build the dated, integrated figure release without modifying final_result/.
set -euo pipefail

root_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
output_arg="${1:-final_result_integrated_20260905}"
if [[ "$output_arg" = /* ]]; then
  output_dir="$output_arg"
else
  output_dir="$root_dir/$output_arg"
fi

if [[ -e "$output_dir" ]]; then
  printf 'Refusing to overwrite an existing integration output: %s\n' "$output_dir" >&2
  exit 2
fi

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/publication-integrated-release.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

# Start from the frozen package so all untouched figures, tables and audits are
# carried forward byte-for-byte.
ditto "$root_dir/final_result" "$output_dir"

mkdir -p "$tmp_dir/fig1b" "$tmp_dir/expression"
PROJECT_ROOT="$root_dir" PUBLICATION_OUTPUT_DIR="$tmp_dir/fig1b" \
  Rscript --vanilla "$root_dir/R/final_figures/02_fig1b_metadata_bubble.R"
PROJECT_ROOT="$root_dir" PUBLICATION_OUTPUT_DIR="$tmp_dir/expression" FIGURE_FONT=Arial \
  Rscript --vanilla "$root_dir/R/final_figures/05_expression_tree_heatmaps.R"
PYTHON="${PYTHON:-python}" "$root_dir/workflow/run_log2fc_per_tissue_pathway_spearman_20260902.sh" \
  "$tmp_dir/log2fc" name_only

mkdir -p "$output_dir/Suppl/Fig_S6_per_tissue"
cp "$tmp_dir/fig1b/Main/Fig_1b.png" "$output_dir/Main/Fig_1b.png"
cp "$tmp_dir/fig1b/Main/Fig_1b.pdf" "$output_dir/Main/Fig_1b.pdf"
cp "$tmp_dir/fig1b/provenance/Fig_1b_mouse_metadata_grouped.csv" \
  "$output_dir/provenance/Fig_1b_mouse_metadata_grouped.csv"
cp "$tmp_dir/fig1b/provenance/figure_palette_audit.csv" \
  "$output_dir/provenance/Fig_1b_palette_audit_20260905.csv"

for figure_id in Fig_4a Fig_5a; do
  cp "$tmp_dir/expression/Main/$figure_id.png" "$output_dir/Main/$figure_id.png"
  cp "$tmp_dir/expression/Main/$figure_id.pdf" "$output_dir/Main/$figure_id.pdf"
done
cp "$tmp_dir/expression/provenance/expression_tree_structure_audit.csv" \
  "$output_dir/provenance/expression_tree_structure_audit.csv"
cp "$tmp_dir/expression/provenance/figure_palette_audit.csv" \
  "$output_dir/provenance/expression_tree_palette_audit_20260905.csv"

for suffix in png pdf; do
  cp "$tmp_dir/log2fc/Main/Fig_2.$suffix" "$output_dir/Main/Fig_2.$suffix"
  cp "$tmp_dir/log2fc/Suppl/Fig_S6_log2fc_per_tissue_merged.$suffix" \
    "$output_dir/Suppl/Fig_S6_log2fc_per_tissue_merged.$suffix"
  cp "$tmp_dir/log2fc/Suppl/Fig_S6_per_tissue/"*".$suffix" \
    "$output_dir/Suppl/Fig_S6_per_tissue/"
done
cp "$tmp_dir/log2fc/provenance/figure_manifest.csv" \
  "$output_dir/provenance/log2fc_spearman_figure_manifest_20260905.csv"
cp "$tmp_dir/log2fc/provenance/12_validation_summary.csv" \
  "$output_dir/provenance/log2fc_spearman_validation_20260905.csv"
cp "$tmp_dir/log2fc/provenance/renderer_README.md" \
  "$output_dir/provenance/log2fc_spearman_renderer_20260905.md"
cp "$root_dir/README_integrated_final_result_20260905.md" "$output_dir/README.md"

# Contract checks: untouched tables, retained qsmooth tree geometry, frozen
# Fig. 1b reference palette and complete GO-name Spearman output.
diff -qr "$root_dir/final_result/tables" "$output_dir/tables"
cmp "$root_dir/final_result/provenance/expression_tree_structure_audit.csv" \
  "$output_dir/provenance/expression_tree_structure_audit.csv"
expected_fig1b_palette='#FDE333,#D4E02D,#A6DA42,#73D25B,#25C771,#00BA82,#00AC8E,#009B95,#008A98,#007796,#006290,#1E4D85,#3C3777,#471D67,#4B0055'
observed_fig1b_palette=$(awk -F, 'NR > 1 {print $3}' "$tmp_dir/fig1b/provenance/figure_palette_audit.csv" | paste -sd, -)
[[ "$observed_fig1b_palette" == "$expected_fig1b_palette" ]]
[[ $(find "$output_dir/Suppl/Fig_S6_per_tissue" -maxdepth 1 -name '*.png' -type f | wc -l | tr -d ' ') == 26 ]]
[[ $(find "$output_dir/Suppl/Fig_S6_per_tissue" -maxdepth 1 -name '*.pdf' -type f | wc -l | tr -d ' ') == 26 ]]

git -C "$root_dir" rev-parse HEAD > "$output_dir/provenance/integration_git_revision.txt"
shasum -a 256 \
  "$output_dir/Main/Fig_1b.png" "$output_dir/Main/Fig_2.png" \
  "$output_dir/Main/Fig_4a.png" "$output_dir/Main/Fig_5a.png" \
  "$output_dir/Suppl/Fig_S6_log2fc_per_tissue_merged.png" \
  "$output_dir/tables/"* | sort > "$output_dir/provenance/integrated_replacement_and_table_sha256.txt"

printf 'INTEGRATED_FINAL_RESULT_PASS: %s\n' "$output_dir"
