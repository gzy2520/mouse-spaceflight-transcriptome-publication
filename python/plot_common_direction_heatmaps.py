#!/usr/bin/env python3
"""Render the four approved common-direction log2FC heatmaps."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-manuscript-release-24")

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("PROJECT_ROOT", Path.cwd())).resolve()
OUT = Path(os.environ["PUBLICATION_OUTPUT_DIR"]).resolve()
INPUT = ROOT / "data/publication_input/common_direction"
STYLE_C = os.environ.get("FIGURE_STYLE", "approved") == "C"
STYLE_A = os.environ.get("FIGURE_STYLE", "approved") == "A"

SPECS = {
    "A_kidney": {
        "title": "A | Kidney",
        "cutoff": 0.05,
        "output": "Suppl/Fig_S2b.png",
        "labels": {
            "OSD-163": "OSD-163\nSpaceX-8",
            "OSD-253": "OSD-253\nSpaceX-15",
            "OSD-462": "OSD-462\nSpaceX-21",
            "OSD-513": "OSD-513\nSpaceX-21",
            "OSD-771": "OSD-771\nSpaceX-18",
        },
    },
    "B_thymus": {
        "title": "B | Four thymus comparisons",
        "cutoff": 0.05,
        "output": "Main/Fig_3b.png",
        "labels": {
            "OSD289_MHU1_thymus": "OSD-289 MHU-1\nuG - GC",
            "OSD289_MHU2_thymus": "OSD-289 MHU-2\nuG - vivarium",
            "OSD421_thymus": "OSD-421\nFlight - Ground",
            "OSD515_thymus": "OSD-515\nFlight - Ground",
        },
    },
    "C_NES_lt_minus1_excl_Lung_Thymus": {
        "title": "C | NES < -1, excluding Lung and Thymus",
        "cutoff": 0.10,
        "output": "Suppl/Fig_S4b.png",
    },
    "D_NES_gt1_excl_Kidney": {
        "title": "D | NES > 1, excluding Kidney",
        "cutoff": 0.10,
        "output": "Suppl/Fig_S3b.png",
    },
}


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 9 if not STYLE_A else 9.4,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.0,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def signed_cmap() -> LinearSegmentedColormap:
    cmap = LinearSegmentedColormap.from_list(
        "negative_blue_positive_red",
        ["#3775BA", "#F7F7F7", "#B64342"] if not STYLE_A else ["#3B6FB6", "#F7F7F5", "#C65A5A"],
    )
    cmap.set_bad("#D9D9D9")
    return cmap


def row_labels(selected: pd.DataFrame) -> list[str]:
    symbols = selected["display_symbol"].where(
        selected["display_symbol"].ne(""), selected["ensembl_id"]
    )
    arrows = selected["common_direction"].map(
        {"common_up": "↑", "common_down": "↓"}
    ).fillna("")
    if STYLE_A:
        return [
            f"{symbol} {arrow}  q={q:.2g}  min|FC|={score:.2f}"
            for symbol, arrow, q, score in zip(
                symbols, arrows, selected["final_fdr"], selected["ranking_score"]
            )
        ]
    return [
        f"{symbol} {arrow} | q={q:.2g} | minimum |log2FC|={score:.2f}"
        for symbol, arrow, q, score in zip(
            symbols, arrows, selected["final_fdr"], selected["ranking_score"]
        )
    ]


def row_labels_c(selected: pd.DataFrame) -> list[str]:
    symbols = selected["display_symbol"].where(
        selected["display_symbol"].ne(""), selected["ensembl_id"]
    )
    arrows = selected["common_direction"].map(
        {"common_up": "↑", "common_down": "↓"}
    ).fillna("")
    return [f"{symbol} {arrow}" for symbol, arrow in zip(symbols, arrows)]


def render_one(group_id: str, columns: pd.DataFrame, norm: Normalize) -> None:
    spec = SPECS[group_id]
    source = INPUT / f"{group_id}_DDR_common_direction_top30_min_abs_log2fc_direction_ordered.csv"
    selected = pd.read_csv(source, encoding="utf-8-sig")
    unit_order = tuple(columns["column_id"])
    if "labels" in spec:
        labels = spec["labels"]
    else:
        labels = {
            row.column_id: (
                f"{row.column_id[len('Heart / '):] if row.column_id.startswith('Heart / ') else row.column_id}"
                f"\nNES={float(row.DDR_NES):.2f}"
            )
            for row in columns.itertuples(index=False)
        }
    values = selected[[f"log2fc__{unit}" for unit in unit_order]].to_numpy(float)
    n_rows, n_cols = values.shape
    if STYLE_C:
        render_one_c(group_id, selected, columns, values, norm, labels)
        return
    fig, ax = plt.subplots(
        figsize=(max(10.5, 1.25 * n_cols + 4.5), max(7.2, 3.0 + 0.30 * n_rows))
    )
    cmap = signed_cmap()
    ax.imshow(
        np.ma.masked_invalid(values), aspect="auto", interpolation="none",
        cmap=cmap, norm=norm
    )
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([labels[unit] for unit in unit_order], rotation=35, ha="right")
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels(selected), fontsize=6.2 if n_rows > 30 else 7.2)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", length=1.8, pad=1.5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    up_n = int(selected["common_direction"].eq("common_up").sum())
    down_n = int(selected["common_direction"].eq("common_down").sum())
    if STYLE_A:
        title = (
            f"{spec['title']}\nq < {spec['cutoff']:.2f} · pooled Top 30 · "
            f"direction ordered (up {up_n}, down {down_n})"
        )
    else:
        title = (
            f"{spec['title']}\nq < {spec['cutoff']:.2f}; same pooled Top 30 ordered by common direction "
            f"(up {up_n}, down {down_n}), then minimum |log2FC|"
        )
    ax.set_title(title, loc="left", weight="bold", pad=8)
    for row in range(n_rows):
        for column in range(n_cols):
            value = values[row, column]
            if np.isfinite(value):
                ax.text(
                    column, row, "0.00" if abs(value) < 0.005 else f"{value:.2f}",
                    ha="center", va="center", fontsize=5.7 if n_rows > 30 else 6.2,
                    color="white" if abs(value) > 0.65 * norm.vmax else "#222222",
                )
    fig.tight_layout(pad=1.0, rect=(0.01, 0.02, 0.86, 0.96))
    colorbar_axis = fig.add_axes([0.89, 0.24, 0.018, 0.52])
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=colorbar_axis)
    colorbar.set_label("Unified log2FC: spaceflight - ground/control")
    fig.savefig(
        OUT / spec["output"], dpi=600 if STYLE_A else 300, bbox_inches="tight", pad_inches=0.08,
        facecolor="white",
        # Keep the approved archive metadata stable across matplotlib patch versions.
        metadata={"Software": "Matplotlib version3.10.8, https://matplotlib.org/"},
    )
    plt.close(fig)


def render_one_c(
    group_id: str,
    selected: pd.DataFrame,
    columns: pd.DataFrame,
    values: np.ndarray,
    norm: Normalize,
    labels: dict[str, str],
) -> None:
    """C-story heatmap: concise gene labels and a quieter cell layer."""
    n_rows, n_cols = values.shape
    fig, ax = plt.subplots(
        figsize=(max(9.5, 1.12 * n_cols + 3.8), max(6.8, 2.8 + 0.24 * n_rows))
    )
    cmap = signed_cmap()
    ax.imshow(
        np.ma.masked_invalid(values), aspect="auto", interpolation="none",
        cmap=cmap, norm=norm
    )
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([labels[unit] for unit in columns["column_id"]], rotation=38, ha="right")
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels_c(selected), fontsize=7.1 if n_rows > 30 else 7.6)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.45)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", length=1.5, pad=2.0)
    direction_colours = {"common_up": "#C65A5A", "common_down": "#3B6FB6"}
    for tick, direction in zip(ax.get_yticklabels(), selected["common_direction"]):
        tick.set_color(direction_colours.get(direction, "#263442"))
    spec = SPECS[group_id]
    up_n = int(selected["common_direction"].eq("common_up").sum())
    down_n = int(selected["common_direction"].eq("common_down").sum())
    ax.set_title(
        f"{spec['title']}\nTop 30 pooled genes; common direction (up {up_n}, down {down_n})",
        loc="left", weight="bold", pad=8,
    )
    fig.tight_layout(pad=1.0, rect=(0.01, 0.02, 0.86, 0.96))
    colorbar_axis = fig.add_axes([0.89, 0.24, 0.018, 0.52])
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=colorbar_axis)
    colorbar.set_label("Unified log2FC: spaceflight − ground/control")
    stem = OUT / spec["output"]
    fig.savefig(stem, dpi=600, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig)


def main() -> None:
    apply_style()
    columns = pd.read_csv(INPUT / "06_group_column_order.csv", encoding="utf-8-sig")
    scale = pd.read_csv(INPUT / "07_figure_colour_scale_audit.csv", encoding="utf-8-sig")
    norm = Normalize(vmin=float(scale.loc[0, "vmin"]), vmax=float(scale.loc[0, "vmax"]))
    for group_id in SPECS:
        group_columns = columns.loc[columns["group_id"].eq(group_id)].sort_values("column_position")
        if group_columns.empty:
            raise ValueError(f"Missing column order for {group_id}")
        render_one(group_id, group_columns, norm)


if __name__ == "__main__":
    main()
