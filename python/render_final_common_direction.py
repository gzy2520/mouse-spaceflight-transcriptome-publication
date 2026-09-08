#!/usr/bin/env python3
"""Render final common-direction log2FC heatmaps from frozen matrices."""

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

SELECTED_CMAP = {"low": "#3775BA", "mid": "#F7F7F7", "high": "#B64342"}


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            # Keep the selected supplementary heatmap typography unchanged;
            # Fig. 3b receives its larger sizes explicitly in render_final_matrix().
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.0,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def signed_cmap() -> LinearSegmentedColormap:
    palette = SELECTED_CMAP
    cmap = LinearSegmentedColormap.from_list(
        "negative_blue_positive_red",
        [palette["low"], palette["mid"], palette["high"]],
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
    render_final_matrix(group_id, selected, columns, values, norm, labels)


def render_final_matrix(
    group_id: str,
    selected: pd.DataFrame,
    columns: pd.DataFrame,
    values: np.ndarray,
    norm: Normalize,
    labels: dict[str, str],
) -> None:
    """Final heatmap with concise display labels and paired PNG/PDF export."""
    n_rows, n_cols = values.shape
    is_final_fig3b = group_id == "B_thymus"
    fig, ax = plt.subplots(
        figsize=(
            13.2 if is_final_fig3b else max(9.5, 1.12 * n_cols + 3.8),
            12.2 if is_final_fig3b else max(6.8, 2.8 + 0.24 * n_rows),
        )
    )
    cmap = signed_cmap()
    ax.imshow(
        np.ma.masked_invalid(values), aspect="auto", interpolation="none",
        cmap=cmap, norm=norm
    )
    ax.set_xticks(np.arange(n_cols))
    if is_final_fig3b:
        final_labels = {
            "OSD289_MHU1_thymus": "OSD-289 MHU-1\nFlight-Ground",
            "OSD289_MHU2_thymus": "OSD-289 MHU-2\nFlight-Ground (vivarium)",
            "OSD421_thymus": "OSD-421\nFlight-Ground",
            "OSD515_thymus": "OSD-515\nFlight-Ground",
        }
        displayed_labels = [final_labels[unit] for unit in columns["column_id"]]
        ax.set_xticklabels(displayed_labels, rotation=0, ha="center", linespacing=1.05)
    else:
        displayed_labels = [labels[unit] for unit in columns["column_id"]]
        ax.set_xticklabels(displayed_labels, rotation=38, ha="right")
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(
        row_labels(selected),
        fontsize=(12.8 if is_final_fig3b else (7.1 if n_rows > 30 else 7.6)),
    )
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.45)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="x", labelsize=14.5 if is_final_fig3b else 8, pad=10 if is_final_fig3b else 3.5)
    ax.tick_params(axis="y", length=1.5, pad=2.0)
    direction_colours = {"common_up": "#C65A5A", "common_down": "#3B6FB6"}
    for tick, direction in zip(ax.get_yticklabels(), selected["common_direction"]):
        tick.set_color(direction_colours.get(direction, "#263442"))
    spec = SPECS[group_id]
    up_n = int(selected["common_direction"].eq("common_up").sum())
    down_n = int(selected["common_direction"].eq("common_down").sum())
    ax.set_title(
        f"{spec['title']}\nTop 30 pooled genes; common direction (up {up_n}, down {down_n})",
        loc="left", weight="bold", pad=10 if is_final_fig3b else 8,
        fontsize=21 if is_final_fig3b else None,
    )
    fig.tight_layout(pad=1.0, rect=(0.01, 0.03 if is_final_fig3b else 0.02, 0.86, 0.96))
    colorbar_axis = fig.add_axes([0.89, 0.22 if is_final_fig3b else 0.24, 0.020, 0.56 if is_final_fig3b else 0.52])
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=colorbar_axis)
    colorbar.set_label(
        "Unified log2FC: Flight-Ground" if is_final_fig3b
        else "Unified log2FC: spaceflight − ground/control",
        fontsize=15 if is_final_fig3b else None,
    )
    if is_final_fig3b:
        colorbar.ax.tick_params(labelsize=13)
    stem = OUT / spec["output"]
    fig.savefig(stem, dpi=600, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig)


def main() -> None:
    apply_style()
    columns = pd.read_csv(INPUT / "06_group_column_order.csv", encoding="utf-8-sig")
    scale = pd.read_csv(INPUT / "07_figure_colour_scale_audit.csv", encoding="utf-8-sig")
    norm = Normalize(vmin=float(scale.loc[0, "vmin"]), vmax=float(scale.loc[0, "vmax"]))
    label_audit: list[dict[str, str]] = []
    for group_id in SPECS:
        group_columns = columns.loc[columns["group_id"].eq(group_id)].sort_values("column_position")
        if group_columns.empty:
            raise ValueError(f"Missing column order for {group_id}")
        render_one(group_id, group_columns, norm)
        if group_id == "B_thymus":
            final_labels = [
                "OSD-289 MHU-1 / Flight-Ground",
                "OSD-289 MHU-2 / Flight-Ground (vivarium)",
                "OSD-421 / Flight-Ground",
                "OSD-515 / Flight-Ground",
            ]
            for position, (column_id, display_label) in enumerate(
                zip(group_columns["column_id"], final_labels), start=1
            ):
                label_audit.append(
                    {"position": str(position), "column_id": column_id, "display_label": display_label}
                )
    provenance = OUT / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(label_audit).to_csv(
        provenance / "Fig_3b_column_label_audit.csv", index=False
    )
    palette_path = provenance / "figure_palette_audit.csv"
    current = pd.read_csv(palette_path) if palette_path.exists() else pd.DataFrame(
        columns=["figure", "role", "hex"]
    )
    figure_ids = [Path(spec["output"]).stem for spec in SPECS.values()]
    current = current.loc[~current["figure"].isin(figure_ids)]
    palette_rows = [
        {"figure": figure_id, "role": role, "hex": value.upper()}
        for figure_id in figure_ids
        for role, value in SELECTED_CMAP.items()
    ]
    pd.concat([current, pd.DataFrame(palette_rows)], ignore_index=True).sort_values(
        ["figure", "role"]
    ).to_csv(palette_path, index=False)


if __name__ == "__main__":
    main()
