#!/usr/bin/env python3
"""Render sample-level log2FC GO-term Spearman figures.

The numeric tables are produced by
``build_log2fc_per_tissue_pathway_spearman_20260902.py``.  This script is a
display-only renderer: it never re-computes expression values or correlations.
It writes one full-size heatmap for every tissue, a 26-tissue small-multiple
figure, and the overall 15 x 15 matrix used for the revised Fig. 2.

The blue-white-red palette is intentionally the approved Fig. 2 palette.  GO
term keys remain stable Ensembl-backed analysis identifiers; display labels are
only the frozen GO names with the previously requested wording cleanup:
``Intrinsic apoptotic signaling`` and ``Telomeric region``.
"""

from __future__ import annotations

import argparse
import os
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize


EXPECTED_N_TERMS = 15
EXPECTED_N_TISSUES = 26
PAIR_COUNT = EXPECTED_N_TERMS * (EXPECTED_N_TERMS - 1) // 2
PALETTE = {
    "low": "#3B6FB6",
    "mid": "#FFFFFF",
    "high": "#C65A5A",
    "ink": "#263238",
    "muted": "#52606D",
    "grid": "#FFFFFF",
}
FIG2_CMAP = LinearSegmentedColormap.from_list(
    "approved_fig2_blue_white_red",
    [PALETTE["low"], PALETTE["mid"], PALETTE["high"]],
    N=257,
)


def configure_style() -> None:
    """Use a print-safe sans-serif style with deliberately large labels."""

    plt.rcParams.update(
        {
            "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 13,
            "axes.titlesize": 18,
            "axes.labelsize": 14,
            "axes.linewidth": 1.1,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def path_for_manifest(path: Path, root: Path) -> str:
    """Use a repository-relative path, hiding ephemeral render locations."""

    try:
        return str(path.relative_to(root))
    except ValueError:
        return f"<temporary-render-output>/{path.name}"


def read_term_metadata(path: Path) -> tuple[pd.DataFrame, list[str]]:
    meta = pd.read_csv(path)
    required = {"term_key", "term_name", "go_id", "term_order"}
    _assert(required.issubset(meta.columns), f"Term metadata lacks {required - set(meta.columns)}")
    meta = meta.sort_values("term_order").reset_index(drop=True)
    _assert(len(meta) == EXPECTED_N_TERMS, f"Expected {EXPECTED_N_TERMS} terms, found {len(meta)}")
    _assert(meta.term_key.is_unique, "Term metadata contains duplicate term_key values")
    return meta, meta.term_key.tolist()


def display_name(row: pd.Series, width: int = 17) -> str:
    """Return compact, readable display wording without changing the key.

    The full source term name remains in ``10_term_metadata.csv``.  These
    shorter labels are used only on dense axes so that enlarged type does not
    collide; the two requested wording edits are explicit in the mapping.
    """

    compact_names = {
        "DNA_damage_response_GO0006974": "DNA damage response",
        "DNA_repair_GO0006281": "DNA repair",
        "DNA_damage_signal_transduction_GO0042770": "DNA-damage signal transduction",
        "intrinsic_apoptotic_signaling_GO0008630": "Intrinsic apoptotic signaling",
        "DNA_damage_tolerance_GO0006301": "DNA damage tolerance",
        "telomere_maintenance_GO0043247": "Telomere maintenance",
        "telomere_region_GO0000781": "Telomeric region",
        "cell_cycle_GO0007049": "Cell cycle",
        "DNA_replication_GO0006260": "DNA replication",
        "DNA_templated_transcription_GO0006351": "DNA-templated transcription",
        "translation_GO0006412": "Translation",
        "mitochondrion_GO0005739": "Mitochondrion",
        "mechanosensation_GO0050954": "Mechanical sensing",
        "cytoskeleton_GO0005856": "Cytoskeleton",
        "cell_death_GO0008219": "Cell death",
    }
    name = compact_names.get(str(row.term_key), str(row.term_name))
    wrapped = textwrap.fill(name, width=width, break_long_words=False)
    return f"{wrapped}\n{row.go_id}"


def read_overall(path: Path, term_keys: list[str]) -> pd.DataFrame:
    table = pd.read_csv(path)
    _assert(table.columns[0] == "term_key", "Overall matrix must start with term_key")
    table = table.set_index("term_key")
    _assert(table.index.is_unique, "Overall matrix has duplicate term keys")
    _assert(set(table.index) == set(term_keys), "Overall matrix term set differs from metadata")
    table = table.loc[term_keys, term_keys].astype(float)
    values = table.to_numpy()
    _assert(np.isfinite(values).all(), "Overall matrix contains non-finite values")
    _assert(np.max(np.abs(values - values.T)) < 1e-10, "Overall matrix is not symmetric")
    _assert(np.allclose(np.diag(values), 1.0), "Overall matrix diagonal is not one")
    _assert(values.min() >= -1.0000001 and values.max() <= 1.0000001, "Overall rho outside [-1, 1]")
    return table


def read_tissue_matrices(
    path: Path, term_keys: list[str]
) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    pairs = pd.read_csv(path)
    required = {
        "analysis_tissue",
        "n_flight_samples",
        "pathway_1_key",
        "pathway_2_key",
        "spearman_rho",
    }
    _assert(required.issubset(pairs.columns), f"Pair table lacks {required - set(pairs.columns)}")
    _assert(
        len(pairs) == EXPECTED_N_TISSUES * PAIR_COUNT,
        f"Expected {EXPECTED_N_TISSUES} x {PAIR_COUNT} pairs, found {len(pairs)}",
    )
    _assert(
        not pairs.duplicated(["analysis_tissue", "pathway_1_key", "pathway_2_key"]).any(),
        "Tissue pair table contains duplicate pairs",
    )
    _assert(set(pairs.pathway_1_key).issubset(term_keys), "Unknown pathway_1_key in pair table")
    _assert(set(pairs.pathway_2_key).issubset(term_keys), "Unknown pathway_2_key in pair table")

    index = {key: i for i, key in enumerate(term_keys)}
    matrices: dict[str, np.ndarray] = {}
    sample_counts: dict[str, int] = {}
    for tissue, chunk in pairs.groupby("analysis_tissue", sort=True):
        _assert(len(chunk) == PAIR_COUNT, f"{tissue}: expected {PAIR_COUNT} pairs, found {len(chunk)}")
        counts = chunk.n_flight_samples.dropna().astype(int).unique()
        _assert(len(counts) == 1, f"{tissue}: inconsistent n_flight_samples")
        matrix = np.full((len(term_keys), len(term_keys)), np.nan, dtype=float)
        np.fill_diagonal(matrix, 1.0)
        for row in chunk.itertuples(index=False):
            i = index[row.pathway_1_key]
            j = index[row.pathway_2_key]
            rho = float(row.spearman_rho) if pd.notna(row.spearman_rho) else np.nan
            _assert(not np.isfinite(rho) or -1.0000001 <= rho <= 1.0000001, f"{tissue}: rho outside [-1, 1]")
            matrix[i, j] = rho
            matrix[j, i] = rho
        matrices[tissue] = matrix
        sample_counts[tissue] = int(counts[0])
    _assert(len(matrices) == EXPECTED_N_TISSUES, f"Expected {EXPECTED_N_TISSUES} tissues, found {len(matrices)}")
    return matrices, sample_counts


def read_cluster_order(path: Path, term_keys: list[str]) -> list[str]:
    order = pd.read_csv(path)
    _assert({"cluster_position", "term_key"}.issubset(order.columns), "Cluster order schema changed")
    order = order.sort_values("cluster_position")
    ordered = order.term_key.tolist()
    _assert(len(ordered) == len(term_keys) and set(ordered) == set(term_keys), "Invalid cluster order")
    return ordered


def ordered_matrix(matrix: np.ndarray, current_keys: list[str], order: list[str]) -> np.ndarray:
    indices = [current_keys.index(key) for key in order]
    return matrix[np.ix_(indices, indices)]


def add_heatmap_text(ax: plt.Axes, matrix: np.ndarray, fontsize: float = 8.4) -> None:
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if not np.isfinite(value):
                continue
            colour = "white" if abs(value) >= 0.58 else PALETTE["ink"]
            ax.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=fontsize,
                color=colour,
                fontweight="normal",
                clip_on=True,
            )


def draw_single(
    matrix: np.ndarray,
    labels: list[str],
    title: str,
    subtitle: str,
    output_stem: Path,
    mask_upper: bool = False,
    x_labels: list[str] | None = None,
    x_label_rows: int = 3,
) -> list[Path]:
    fig, ax = plt.subplots(figsize=(14.5, 14.5))
    display_matrix = matrix.copy()
    if mask_upper:
        display_matrix[np.triu_indices_from(display_matrix, k=1)] = np.nan
    masked = np.ma.masked_invalid(display_matrix)
    im = ax.imshow(masked, cmap=FIG2_CMAP, norm=Normalize(-1, 1), interpolation="nearest")
    n = matrix.shape[0]
    ticks = np.arange(n)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    # Keep the y-axis compact with readable single-line term names
    y_labels = [label.rsplit("\n", 1)[0].replace("\n", " ") for label in labels]
    ax.set_yticklabels(y_labels, fontsize=10.4, linespacing=0.88)
    ax.tick_params(axis="both", length=0, pad=8)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color=PALETTE["grid"], linewidth=0.55)
    ax.tick_params(which="minor", bottom=False, left=False)
    if x_labels is None:
        x_labels = labels
    _assert(len(x_labels) == n, "x-axis label count does not match matrix size")
    clean_x_labels = [lbl.replace("\n", " ") if "\nGO:" not in lbl else lbl for lbl in x_labels]
    ax.set_xticklabels(
        clean_x_labels,
        rotation=55,
        ha="right",
        rotation_mode="anchor",
        fontsize=10.8,
        linespacing=0.88,
        color=PALETTE["ink"],
    )
    ax.tick_params(axis="x", length=0, pad=6)
    add_heatmap_text(ax, display_matrix)
    fig.suptitle(title, x=0.27, y=0.978, ha="left", fontweight="bold", fontsize=19, color=PALETTE["ink"])
    fig.text(0.27, 0.938, subtitle, ha="left", va="bottom", fontsize=12.2, color=PALETTE["muted"])
    cbar = fig.colorbar(im, ax=ax, fraction=0.038, pad=0.04, ticks=[-1, -0.5, 0, 0.5, 1])
    cbar.set_label("Spearman rho", fontsize=13.5, fontweight="bold", labelpad=12)
    cbar.ax.tick_params(labelsize=11, length=3)
    fig.text(
        0.01,
        0.008,
        "Sample-level log2FC; cells with undefined rho are left blank. Stable GO IDs are retained in frozen metadata.",
        ha="left",
        va="bottom",
        fontsize=10.2,
        color=PALETTE["muted"],
    )
    fig.subplots_adjust(left=0.27, right=0.90, top=0.88, bottom=0.24)
    outputs = []
    for suffix, kwargs in ((".png", {"dpi": 400}), (".pdf", {})):
        path = output_stem.with_suffix(suffix)
        fig.savefig(path, bbox_inches="tight", pad_inches=0.10, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def draw_merged(
    matrices: dict[str, np.ndarray],
    sample_counts: dict[str, int],
    order: list[str],
    term_labels_by_key: dict[str, str],
    output_stem: Path,
    legend_columns: int = 5,
) -> list[Path]:
    tissues = sorted(matrices)
    n_panels = len(tissues)
    ncols, nrows = 5, int(np.ceil(n_panels / 5))
    fig, axes = plt.subplots(nrows, ncols, figsize=(19.0, 3.55 * nrows), squeeze=False)
    axes_flat = axes.ravel()
    last_im = None
    for ax, tissue in zip(axes_flat, tissues):
        matrix = ordered_matrix(matrices[tissue], list(term_labels_by_key), order)
        last_im = ax.imshow(
            np.ma.masked_invalid(matrix), cmap=FIG2_CMAP, norm=Normalize(-1, 1), interpolation="nearest"
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(
            f"{tissue}\n$n$ = {sample_counts[tissue]}",
            fontsize=11.2,
            fontweight="bold",
            color=PALETTE["ink"],
            pad=5,
        )
        ax.set_xticks(np.arange(-0.5, EXPECTED_N_TERMS, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, EXPECTED_N_TERMS, 1), minor=True)
        ax.grid(which="minor", color=PALETTE["grid"], linewidth=0.22)
        ax.tick_params(which="minor", bottom=False, left=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    fig.suptitle(
        "Per-tissue GO-term relationships from sample-level log2FC",
        x=0.50,
        y=0.995,
        fontsize=21,
        fontweight="bold",
        color=PALETTE["ink"],
    )
    fig.text(
        0.50,
        0.978,
        "Each panel is a tissue-level Spearman matrix; n is the number of Flight samples contributing all 15 terms."
        " Term order follows the overall average-linkage order.",
        ha="center",
        va="top",
        fontsize=12.5,
        color=PALETTE["muted"],
    )
    legend_lines = []
    for pos, key in enumerate(order, start=1):
        legend_lines.append(f"{pos}. {term_labels_by_key[key].replace(chr(10), ' ')}")
    legend_text = "   ".join(legend_lines)
    _assert(legend_columns >= 1, "legend_columns must be positive")
    legend_rows = [
        "   ".join(legend_lines[i : i + legend_columns])
        for i in range(0, len(legend_lines), legend_columns)
    ]
    fig.text(
        0.50,
        0.016,
        "Term key (cluster order):\n" + "\n".join(legend_rows),
        ha="center",
        va="bottom",
        fontsize=9.5 if legend_columns >= 8 else 8.8,
        color=PALETTE["muted"],
        wrap=True,
    )
    if last_im is not None:
        cax = fig.add_axes([0.965, 0.24, 0.012, 0.54])
        cbar = fig.colorbar(last_im, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1])
        cbar.set_label("Spearman rho", fontsize=12.5, fontweight="bold", labelpad=9)
        cbar.ax.tick_params(labelsize=10)
    fig.subplots_adjust(
        left=0.025,
        right=0.94,
        top=0.94,
        bottom=0.11 if legend_columns >= 8 else 0.14,
        wspace=0.12,
        hspace=0.52,
    )
    outputs = []
    for suffix, kwargs in ((".png", {"dpi": 350}), (".pdf", {})):
        path = output_stem.with_suffix(suffix)
        fig.savefig(path, bbox_inches="tight", pad_inches=0.10, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Main project root")
    parser.add_argument(
        "--input-dir",
        default=None,
        help="Dated builder output (default: 03_analysis_results/06_cross_tissue_integration/10_log2fc_per_tissue_pathway_spearman_20260902)",
    )
    parser.add_argument("--out-dir", default=None, help="Figure output directory (default: <input-dir>/figures)")
    parser.add_argument(
        "--label-mode",
        choices=("name_go_id", "go_id", "name_only"),
        default="name_go_id",
        help="Axis/term-key display mode: compact name plus GO ID, GO IDs only, or compact names only",
    )
    args = parser.parse_args(argv)
    root = Path(args.root).resolve()
    input_dir = Path(args.input_dir).resolve() if args.input_dir else root / "03_analysis_results/06_cross_tissue_integration/10_log2fc_per_tissue_pathway_spearman_20260902"
    figure_dir = Path(args.out_dir).resolve() if args.out_dir else input_dir / "figures"
    single_dir = figure_dir / "per_tissue"
    figure_dir.mkdir(parents=True, exist_ok=True)
    single_dir.mkdir(parents=True, exist_ok=True)

    configure_style()
    meta, term_keys = read_term_metadata(input_dir / "10_term_metadata.csv")
    overall = read_overall(input_dir / "04_overall_fisher_z_spearman_rho_matrix.csv", term_keys)
    matrices, sample_counts = read_tissue_matrices(
        input_dir / "03_tissue_pairwise_spearman_log2fc_long.csv", term_keys
    )
    order = read_cluster_order(input_dir / "07_cluster_order_overall_log2fc.csv", term_keys)
    labels_by_key = {row.term_key: display_name(row) for _, row in meta.iterrows()}
    ordered_labels = [labels_by_key[key] for key in order]
    names_by_key = {key: label.rsplit("\n", 1)[0].replace("\n", " ") for key, label in labels_by_key.items()}
    meta_by_key = meta.set_index("term_key")
    x_labels_by_key = {
        key: (
            str(meta_by_key.loc[key, "go_id"])
            if args.label_mode == "go_id"
            else names_by_key[key]
            if args.label_mode == "name_only"
            else labels_by_key[key]
        )
        for key in term_keys
    }
    ordered_x_labels = [x_labels_by_key[key] for key in order]

    inclusion = pd.read_csv(input_dir / "06_fisher_z_inclusion_audit.csv")
    _assert(set(inclusion.analysis_tissue) == set(matrices), "Inclusion audit tissue set differs from pair table")
    included_n = int(inclusion.included_in_overall.sum())
    _assert(included_n >= 1, "No tissues are included in the overall Fisher-z matrix")

    manifest_rows: list[dict[str, object]] = []
    overall_ordered = overall.loc[order, order].to_numpy()
    overall_outputs = draw_single(
        overall_ordered,
        ordered_labels,
        "Cross-tissue relationships among GO response terms",
        f"Sample-level log2FC → tissue-wise Spearman → equal-weight Fisher-z; {included_n} tissues included (n ≥ 5)",
        figure_dir / "Fig_2_log2fc_fisher_z_overall",
        mask_upper=True,
        x_labels=ordered_x_labels,
        x_label_rows=2 if args.label_mode == "go_id" else 3,
    )
    for path in overall_outputs:
        manifest_rows.append(
            {
                "figure_id": "Fig_2_log2fc_fisher_z_overall",
                "kind": "overall",
                "analysis_tissue": "ALL_INCLUDED_TISSUES",
                "n_flight_samples": int(sum(sample_counts.values())),
                "path": path_for_manifest(path, root),
                "palette": "#3B6FB6 | #FFFFFF | #C65A5A",
                "label_mode": args.label_mode,
            }
        )

    for tissue in sorted(matrices):
        tissue_ordered = ordered_matrix(matrices[tissue], term_keys, order)
        tissue_outputs = draw_single(
            tissue_ordered,
            ordered_labels,
            f"GO-term relationships — {tissue}",
            f"Flight samples n = {sample_counts[tissue]}; tissue-level Spearman rho from sample-level log2FC",
            single_dir / f"Spearman_log2FC_{re.sub(r'[^A-Za-z0-9]+', '_', tissue).strip('_')}",
            x_labels=ordered_x_labels,
            x_label_rows=2 if args.label_mode == "go_id" else 3,
        )
        for path in tissue_outputs:
            manifest_rows.append(
                {
                    "figure_id": "Fig_S6_per_tissue_spearman_log2fc",
                    "kind": "per_tissue",
                    "analysis_tissue": tissue,
                    "n_flight_samples": sample_counts[tissue],
                    "path": path_for_manifest(path, root),
                    "palette": "#3B6FB6 | #FFFFFF | #C65A5A",
                    "label_mode": args.label_mode,
                }
            )

    merged_outputs = draw_merged(
        matrices,
        sample_counts,
        order,
        x_labels_by_key,
        figure_dir / "Fig_S6_log2fc_fisher_z_per_tissue_merged",
        legend_columns=8 if args.label_mode == "go_id" else 3 if args.label_mode == "name_only" else 5,
    )
    for path in merged_outputs:
        manifest_rows.append(
            {
                "figure_id": "Fig_S6_log2fc_fisher_z_per_tissue_merged",
                "kind": "merged_per_tissue",
                "analysis_tissue": "ALL_26_TISSUES",
                "n_flight_samples": int(sum(sample_counts.values())),
                "path": path_for_manifest(path, root),
                "palette": "#3B6FB6 | #FFFFFF | #C65A5A",
                "label_mode": args.label_mode,
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(figure_dir / "figure_manifest.csv", index=False)
    label_description = {
        "go_id": "横轴和合并图 term key 只显示 GO ID；纵轴保留简短 term 名称。",
        "name_only": "横轴、纵轴和合并图 term key 均显示相同的简短 GO 名称；GO ID 保留在冻结元数据中。",
        "name_go_id": "横轴显示简短 term 名称与 GO ID；纵轴显示简短 term 名称。",
    }[args.label_mode]
    readme = f"""# log2FC 逐组织 Spearman 图（2026-09-02）

本目录只读取同目录的冻结计算表，不重新计算表达量或相关性。输入方法是：每个 Flight
样本相对于同一 OSD 数据集、同一组织 Ground log2 均值的样本级 log2FC；先在每个组织中以
小鼠样本为观测单位计算 15 条 GO term 的 Spearman rho，再对 {included_n} 个 n ≥ 5
组织的 rho 做等权 Fisher-z 平均并反变换，得到 revised Fig. 2。n < 5 的组织仍保留在
逐组织图中，但不进入整体汇总。

输出：

- `Fig_2_log2fc_fisher_z_overall.png/pdf`：整体跨组织矩阵（蓝—白—红配色保持原 Fig. 2）。
- `per_tissue/Spearman_log2FC_*.png/pdf`：26 个组织的完整单图，含相关系数数值和所选 term 标签。
- `Fig_S6_log2fc_fisher_z_per_tissue_merged.png/pdf`：26 个组织的合并小 multiples 图。
- `figure_manifest.csv`：每个文件、组织、样本数和配色的审计清单。

第 4 个 term 的显示名为 `Intrinsic apoptotic signaling`，第 7 个为 `Telomeric region`；
这只改显示文字，不改变 term_key 或数值。合并图中的 term 编号按整体 average-linkage
顺序，完整名称和 GO ID 见 `10_term_metadata.csv` 与各组织单图。
{label_description}
"""
    (figure_dir / "README.md").write_text(readme, encoding="utf-8")
    print(
        f"LOG2FC_PER_TISSUE_FIGURES_PASS: overall=2; per_tissue={len(matrices)}x2; merged=2; "
        f"included_tissues={included_n}; output={figure_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
