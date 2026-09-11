#!/usr/bin/env bash

set -euo pipefail

root_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
input_dir="$root_dir/data/publication_input/go/log2fc_spearman_20260902"
output_arg="${1:-reproduced_results/log2fc_spearman_20260902}"
if [[ "$output_arg" = /* ]]; then
  output_dir="$output_arg"
else
  output_dir="$root_dir/$output_arg"
fi
label_mode="${2:-name_go_id}"
if [[ "$label_mode" != "name_go_id" && "$label_mode" != "go_id" && "$label_mode" != "name_only" ]]; then
  printf 'Label mode must be name_go_id, go_id or name_only: %s\n' "$label_mode" >&2
  exit 2
fi
renderer="$root_dir/python/plot_log2fc_per_tissue_pathway_spearman_20260902.py"
python_bin="${PYTHON:-python3}"

if [[ ! -d "$input_dir" || ! -f "$renderer" ]]; then
  printf 'Missing frozen inputs or renderer.\n' >&2
  exit 1
fi
if [[ -e "$output_dir/Main/Fig_2.png" ]]; then
  printf 'Refusing to overwrite an existing dated output: %s\n' "$output_dir" >&2
  exit 1
fi

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/log2fc-spearman-release.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

mkdir -p "$output_dir/Main" "$output_dir/Suppl/Fig_S6_per_tissue" "$output_dir/provenance"
cd "$root_dir"
"$python_bin" "$renderer" --root "$root_dir" --input-dir "$input_dir" --out-dir "$tmp_dir" --label-mode "$label_mode"

for suffix in png pdf; do
  cp "$tmp_dir/Fig_2_log2fc_fisher_z_overall.$suffix" "$output_dir/Main/Fig_2.$suffix"
  cp "$tmp_dir/Fig_S6_log2fc_fisher_z_per_tissue_merged.$suffix" \
    "$output_dir/Suppl/Fig_S6_log2fc_per_tissue_merged.$suffix"
  cp "$tmp_dir/per_tissue/"*".$suffix" "$output_dir/Suppl/Fig_S6_per_tissue/"
done
cp "$tmp_dir/figure_manifest.csv" "$output_dir/provenance/figure_manifest.csv"
cp "$tmp_dir/README.md" "$output_dir/provenance/renderer_README.md"

printf 'LOG2FC_RELEASE_PASS: %s\n' "$output_dir"
