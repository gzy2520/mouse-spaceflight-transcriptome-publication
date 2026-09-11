#!/usr/bin/env python3
"""Build the four-group DDR common-direction log2FC heatmap release.

This release implements the latest selection contract:

1. Construct four prespecified groups (Kidney, four-thymus comparisons,
   NES < -1 excluding Lung/Thymus, and NES > 1 excluding Kidney).
2. Apply a group-specific final BH-FDR cutoff (0.05 for Kidney/thymus and
   0.10 for the two NES-defined tissue groups).
3. After significance filtering, retain only genes whose log2FC has one sign
   in every column of the group.
4. Rank the resulting common-up/common-down genes by either mean |log2FC| or
   minimum |log2FC|, take one global pooled Top 30, and make both a pooled
   figure and a direction-ordered single-panel display of that same Top 30.

Thus each group has four visible figures: mean-pooled, mean-split,
minimum-pooled, and minimum-split. PNG and PDF are written for every figure.
Ensembl Gene ID is the analysis key; Symbol is display-only.
"""

from __future__ import annotations

import hashlib
import importlib.util
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-ddr-four-groups-log2fc")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BUILDER_PATH = ROOT / "07_scripts/06_cross_tissue_integration/build_ddr_common_direction_tissue_groups_current.py"

RUN_DATE = os.environ.get("MOUSE_GO_DDR_LOG2FC_RUN_DATE", "20260826").strip()
if not RUN_DATE.isdigit() or len(RUN_DATE) != 8:
    raise ValueError("MOUSE_GO_DDR_LOG2FC_RUN_DATE must be an YYYYMMDD string")

OUT = ROOT / (
    "03_analysis_results/10_full_stable_id_rerun_20260823/"
    f"05e_ddr_common_direction_four_groups_log2fc_{RUN_DATE}"
)
OUT = Path(os.environ.get("MOUSE_GO_DDR_LOG2FC_OUTPUT_DIR", str(OUT)))
if not OUT.is_absolute():
    OUT = ROOT / OUT
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

DISPLAY_N = 30
RANK_METHODS = {
    "mean_abs_log2fc": {
        "column": "mean_abs_log2fc",
        "label": "mean |log2FC|",
    },
    "min_abs_log2fc": {
        "column": "min_abs_log2fc",
        "label": "minimum |log2FC|",
    },
}

GROUP_ORDER = (
    "A_kidney",
    "B_thymus",
    "C_NES_lt_minus1_excl_Lung_Thymus",
    "D_NES_gt1_excl_Kidney",
)
GROUP_TITLES = {
    "A_kidney": "A | Kidney",
    "B_thymus": "B | Four thymus comparisons",
    "C_NES_lt_minus1_excl_Lung_Thymus": "C | NES < -1, excluding Lung and Thymus",
    "D_NES_gt1_excl_Kidney": "D | NES > 1, excluding Kidney",
}
GROUP_CUTOFFS = {
    "A_kidney": 0.05,
    "B_thymus": 0.05,
    "C_NES_lt_minus1_excl_Lung_Thymus": 0.10,
    "D_NES_gt1_excl_Kidney": 0.10,
}

# C/D order is inherited from the displayed 26-tissue DDR-NES ordering used
# in the preceding release.  It is not re-sorted by this script.
GROUP_TISSUES = {
    "C_NES_lt_minus1_excl_Lung_Thymus": (
        "Heart / Heart right ventricle",
        "Retina",
        "Spleen-distal",
        "Bone marrow",
        "Femoral skin",
        "Cecum",
    ),
    "D_NES_gt1_excl_Kidney": (
        "Quadriceps femoris",
        "Tibialis anterior",
        "Optic nerve",
        "Spleen",
        "Femoral lateral skin",
        "Left lobe of the liver",
        "Colon",
        "Cerebellum",
        "Heart",
        "Eye",
    ),
}

THYMUS_ORDER = (
    "OSD289_MHU1_thymus",
    "OSD289_MHU2_thymus",
    "OSD421_thymus",
    "OSD515_thymus",
)
KIDNEY_ORDER = (
    "OSD-163",
    "OSD-253",
    "OSD-462",
    "OSD-513",
    "OSD-771",
)


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "ddr_common_direction_group_builder", BUILDER_PATH
    )
    if spec is None or spec.loader is None:
        raise ImportError(BUILDER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def signed_cmap() -> LinearSegmentedColormap:
    cmap = LinearSegmentedColormap.from_list(
        "negative_blue_positive_red",
        ["#3775BA", "#F7F7F7", "#B64342"],
    )
    cmap.set_bad("#D9D9D9")
    return cmap


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.0,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def normalize_ensembl(values: pd.Series) -> pd.Series:
    return (
        values.astype(str)
        .str.strip()
        .str.upper()
        .str.replace(r"\.\d+$", "", regex=True)
    )


def add_rank_metrics(common: pd.DataFrame, unit_order: tuple[str, ...]) -> pd.DataFrame:
    result = common.copy()
    value_columns = [f"log2fc__{unit}" for unit in unit_order]
    absolute_values = result[value_columns].abs()
    result["mean_abs_log2fc"] = absolute_values.mean(axis=1)
    result["min_abs_log2fc"] = absolute_values.min(axis=1)
    result["mean_log2fc"] = result[value_columns].mean(axis=1)
    result["min_log2fc"] = result[value_columns].min(axis=1)
    result["max_log2fc"] = result[value_columns].max(axis=1)
    result["common_direction"] = np.where(
        result["effect_pattern"].eq("all_positive"), "common_up", "common_down"
    )
    return result


def select_ranked(
    common: pd.DataFrame,
    rank_column: str,
    direction: str | None,
    n: int = DISPLAY_N,
) -> pd.DataFrame:
    part = common if direction is None else common.loc[
        common["common_direction"].eq(direction)
    ]
    part = part.sort_values(
        [rank_column, "final_fdr", "final_integrated_p", "ensembl_id"],
        ascending=[False, True, True, True],
        kind="mergesort",
    ).head(n).copy()
    if direction is None:
        part["pooled_rank"] = np.arange(1, len(part) + 1)
        part["rank_scope"] = "pooled common-up and common-down"
    else:
        part["rank_within_direction"] = np.arange(1, len(part) + 1)
        part["rank_scope"] = f"common {direction} only"
    part["ranking_score"] = part[rank_column]
    return part.reset_index(drop=True)


def order_directional_display(selected: pd.DataFrame) -> pd.DataFrame:
    """Order one pooled Top-30 table for a single, direction-grouped figure.

    Selection has already happened globally in ``select_ranked``.  This
    function only changes the row order for display: common-up rows first,
    followed by common-down rows, with the requested absolute-effect score
    decreasing within each direction.  It never adds genes to the pooled
    Top-30 set.
    """
    if selected.empty:
        result = selected.copy()
        result["direction_display_order"] = pd.Series(dtype="int64")
        result["rank_scope"] = "direction-ordered display of pooled Top 30"
        return result
    result = selected.copy()
    result["direction_display_order"] = result["common_direction"].map(
        {"common_up": 0, "common_down": 1}
    )
    result = result.sort_values(
        ["direction_display_order", "ranking_score", "final_fdr", "final_integrated_p", "ensembl_id"],
        ascending=[True, False, True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    result["display_order"] = np.arange(1, len(result) + 1)
    result["rank_scope"] = "direction-ordered display of pooled Top 30"
    return result


def _row_labels(
    selected: pd.DataFrame,
    rank_label: str,
    fontsize: float,
) -> list[str]:
    symbols = selected["display_symbol"].where(
        selected["display_symbol"].ne(""), selected["ensembl_id"]
    )
    arrows = selected["common_direction"].map(
        {"common_up": "↑", "common_down": "↓"}
    ).fillna("")
    return [
        f"{symbol} {arrow} | q={q:.2g} | {rank_label}={score:.2f}"
        for symbol, arrow, q, score in zip(
            symbols,
            arrows,
            selected["final_fdr"],
            selected["ranking_score"],
        )
    ]


def draw_heatmap_axis(
    ax: plt.Axes,
    selected: pd.DataFrame,
    unit_order: tuple[str, ...],
    labels: dict[str, str],
    rank_label: str,
    title: str,
    norm: Normalize,
    cmap: LinearSegmentedColormap,
) -> None:
    if selected.empty:
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No genes pass the group-specific cutoff and common-direction rule",
            ha="center",
            va="center",
            fontsize=10,
        )
        ax.set_title(title, loc="left", weight="bold", pad=11)
        return

    values = selected[[f"log2fc__{unit}" for unit in unit_order]].to_numpy(float)
    image = ax.imshow(
        np.ma.masked_invalid(values),
        aspect="auto",
        interpolation="none",
        cmap=cmap,
        norm=norm,
    )
    del image
    n_rows, n_cols = values.shape
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(
        [labels[unit] for unit in unit_order],
        rotation=35,
        ha="right",
    )
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(
        _row_labels(
            selected,
            rank_label,
            fontsize=6.2 if n_rows > 30 else 7.2,
        ),
        fontsize=6.2 if n_rows > 30 else 7.2,
    )
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", length=1.8, pad=1.5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, loc="left", weight="bold", pad=8)
    for row in range(n_rows):
        for column in range(n_cols):
            value = values[row, column]
            if np.isfinite(value):
                ax.text(
                    column,
                    row,
                    "0.00" if abs(value) < 0.005 else f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=5.7 if n_rows > 30 else 6.2,
                    color="white" if abs(value) > 0.65 * norm.vmax else "#222222",
                )


def save_figure(
    fig: plt.Figure,
    stem: str,
    *,
    apply_tight_layout: bool = True,
) -> tuple[Path, Path]:
    if apply_tight_layout:
        fig.tight_layout(pad=1.0)
    outputs = []
    for suffix in ("png", "pdf"):
        path = FIGURES / f"{stem}.{suffix}"
        fig.savefig(
            path,
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.08,
            facecolor="white",
        )
        outputs.append(path)
    plt.close(fig)
    return outputs[0], outputs[1]


def plot_pooled(
    selected: pd.DataFrame,
    unit_order: tuple[str, ...],
    labels: dict[str, str],
    group_title: str,
    cutoff: float,
    rank_label: str,
    rank_method: str,
    stem: str,
    norm: Normalize,
    cmap: LinearSegmentedColormap,
) -> tuple[Path, Path]:
    n_rows = len(selected)
    fig, ax = plt.subplots(
        figsize=(
            max(10.5, 1.25 * len(unit_order) + 4.5),
            max(7.2, 3.0 + 0.30 * n_rows),
        )
    )
    draw_heatmap_axis(
        ax,
        selected,
        unit_order,
        labels,
        rank_label,
        (
            f"{group_title}\n"
            f"q < {cutoff:.2f}; common up/down pooled Top 30 by {rank_label}"
        ),
        norm,
        cmap,
    )
    colorbar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        ax=ax,
        fraction=0.035,
        pad=0.02,
    )
    colorbar.set_label("Unified log2FC: spaceflight - ground/control")
    return save_figure(fig, stem)


def plot_direction_ordered(
    selected: pd.DataFrame,
    unit_order: tuple[str, ...],
    labels: dict[str, str],
    group_title: str,
    cutoff: float,
    rank_label: str,
    rank_method: str,
    stem: str,
    norm: Normalize,
    cmap: LinearSegmentedColormap,
) -> tuple[Path, Path]:
    n_rows = len(selected)
    fig, ax = plt.subplots(
        figsize=(
            max(10.5, 1.25 * len(unit_order) + 4.5),
            max(7.2, 3.0 + 0.30 * n_rows),
        )
    )
    up_n = int(selected["common_direction"].eq("common_up").sum())
    down_n = int(selected["common_direction"].eq("common_down").sum())
    draw_heatmap_axis(
        ax,
        selected,
        unit_order,
        labels,
        rank_label,
        (
            f"{group_title}\n"
            f"q < {cutoff:.2f}; same pooled Top 30 ordered by common direction "
            f"(up {up_n}, down {down_n}), then {rank_label}"
        ),
        norm,
        cmap,
    )
    # Reserve a right-hand strip for the shared colorbar before exporting.
    fig.tight_layout(pad=1.0, rect=(0.01, 0.02, 0.86, 0.96))
    colorbar_axis = fig.add_axes([0.89, 0.24, 0.018, 0.52])
    colorbar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        cax=colorbar_axis,
    )
    colorbar.set_label("Unified log2FC: spaceflight - ground/control")
    return save_figure(fig, stem, apply_tight_layout=False)


def build_group_inputs(builder, support, tissue_frames, tissue_nes, unit_nes, unit_frames):
    thymus_frames, thymus_audit = support.load_four_thymus()
    thymus_groups = {
        mission: tuple(group["dataset_id"])
        for mission, group in thymus_audit.groupby("mission_cluster", sort=False)
    }
    kidney_frames, kidney_audit = builder.build_kidney_accession_frames(
        unit_frames,
        unit_nes,
        KIDNEY_ORDER,
    )
    nes_lookup = tissue_nes.set_index("analysis_tissue")[builder.DDR_NES_COLUMN]
    groups = {
        "A_kidney": (
            kidney_frames,
            KIDNEY_ORDER,
            {accession: (accession,) for accession in KIDNEY_ORDER},
            builder.KIDNEY_LABELS,
        ),
        "B_thymus": (thymus_frames, THYMUS_ORDER, thymus_groups, builder.THYMUS_LABELS),
    }
    for group_id, tissues in GROUP_TISSUES.items():
        frames = {tissue: tissue_frames[tissue] for tissue in tissues}
        labels = {
            tissue: (
                f"{tissue[len('Heart / '):] if tissue.startswith('Heart / ') else tissue}"
                f"\nNES={float(nes_lookup[tissue]):.2f}"
            )
            for tissue in tissues
        }
        groups[group_id] = (
            frames,
            tissues,
            {tissue: (tissue,) for tissue in tissues},
            labels,
        )
    return groups, thymus_audit, kidney_audit


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style()

    builder = load_builder()
    support = builder.load_base_module()
    support.GO_SET_ID_TYPE = "ensembl_gene_id"
    # build_final_results uses this module-level value only to form its
    # descriptive selected table.  We apply the group-specific cutoffs below.
    support.FINAL_FDR_CUTOFF = 1.0

    go_sets, membership = builder.load_go_sets()
    membership.to_csv(
        TABLES / "00_two_GO_union_membership_used.csv",
        index=False,
        encoding="utf-8-sig",
    )

    tissue_nes = pd.read_csv(builder.NES_MATRIX, encoding="utf-8-sig")
    unit_nes = pd.read_csv(builder.UNIT_NES, encoding="utf-8-sig")
    unit_nes = unit_nes.drop_duplicates("analysis_unit_id", keep="first")
    unit_frames = builder.load_unit_frames(set(unit_nes["analysis_unit_id"]))
    tissue_frames, tissue_audit = builder.build_tissue_frames(
        unit_frames,
        unit_nes,
        tuple(tissue_nes["analysis_tissue"]),
    )
    nes_lookup = tissue_nes.set_index("analysis_tissue")[builder.DDR_NES_COLUMN]
    expected_c = {
        tissue
        for tissue, value in nes_lookup.items()
        if float(value) < -1
        and tissue not in {"Lung", "Thymus"}
    }
    expected_d = {
        tissue
        for tissue, value in nes_lookup.items()
        if float(value) > 1 and tissue != "Kidney"
    }
    if expected_c != set(GROUP_TISSUES["C_NES_lt_minus1_excl_Lung_Thymus"]):
        raise ValueError(
            "C tissue list does not match the NES < -1 rule after Lung/Thymus exclusion"
        )
    if expected_d != set(GROUP_TISSUES["D_NES_gt1_excl_Kidney"]):
        raise ValueError(
            "D tissue list does not match the NES > 1 rule after Kidney exclusion"
        )
    tissue_audit = tissue_audit.merge(
        tissue_nes[["analysis_tissue", builder.DDR_NES_COLUMN]],
        on="analysis_tissue",
        how="left",
        validate="one_to_one",
    ).rename(columns={builder.DDR_NES_COLUMN: "DDR_NES"})
    tissue_audit.to_csv(
        TABLES / "01_tissue_meta_aggregation_audit.csv",
        index=False,
        encoding="utf-8-sig",
    )

    group_inputs, thymus_audit, kidney_audit = build_group_inputs(
        builder, support, tissue_frames, tissue_nes, unit_nes, unit_frames
    )
    thymus_audit.to_csv(
        TABLES / "02_four_thymus_input_contrast_audit_used.csv",
        index=False,
        encoding="utf-8-sig",
    )
    kidney_audit.to_csv(
        TABLES / "03_kidney_accession_aggregation_audit.csv",
        index=False,
        encoding="utf-8-sig",
    )

    all_common: dict[str, pd.DataFrame] = {}
    all_results: dict[str, pd.DataFrame] = {}
    summary_rows: list[dict[str, object]] = []
    display_rows: list[dict[str, object]] = []
    figure_requests: list[dict[str, object]] = []
    column_rows: list[dict[str, object]] = []

    nes_lookup = tissue_nes.set_index("analysis_tissue")[builder.DDR_NES_COLUMN]
    for group_id in GROUP_ORDER:
        frames, unit_order, mission_groups, labels = group_inputs[group_id]
        result, _unused_selected, algorithm_summary = support.build_final_results(
            frames, unit_order, mission_groups, go_sets
        )
        cutoff = GROUP_CUTOFFS[group_id]
        significant = result.loc[result["final_fdr"].lt(cutoff)].copy()
        common = significant.loc[
            significant["effect_pattern"].isin(["all_positive", "all_negative"])
        ].copy()
        common = add_rank_metrics(common, unit_order)
        all_results[group_id] = result
        all_common[group_id] = common

        cutoff_tag = f"q{int(round(cutoff * 100)):03d}"
        result.to_csv(
            TABLES / f"{group_id}_DDR_final_test_all.csv",
            index=False,
            encoding="utf-8-sig",
        )
        significant.to_csv(
            TABLES / f"{group_id}_DDR_significant_{cutoff_tag}.csv",
            index=False,
            encoding="utf-8-sig",
        )
        common.to_csv(
            TABLES / f"{group_id}_DDR_common_direction_{cutoff_tag}.csv",
            index=False,
            encoding="utf-8-sig",
        )

        summary_base = {
            "group_id": group_id,
            "group_title": GROUP_TITLES[group_id],
            "significance_cutoff": cutoff,
            "n_columns": len(unit_order),
            "column_ids": ";".join(unit_order),
            "common_complete_background": int(algorithm_summary["common_complete_background"]),
            "two_GO_candidate_family": int(len(result)),
            "significant_n": int(len(significant)),
            "common_direction_n": int(len(common)),
            "common_up_n": int(common["common_direction"].eq("common_up").sum()),
            "common_down_n": int(common["common_direction"].eq("common_down").sum()),
        }
        for rank_method, info in RANK_METHODS.items():
            rank_column = info["column"]
            pooled = select_ranked(common, rank_column, None)
            # The direction-ordered figure is a reordering of the *same*
            # global pooled Top 30.  Do not rank/select the two directions
            # independently: that would allow the display to contain genes
            # that are absent from the pooled heatmap.
            up = pooled.loc[pooled["common_direction"].eq("common_up")].copy()
            down = pooled.loc[pooled["common_direction"].eq("common_down")].copy()
            for table in (up, down):
                table["rank_scope"] = "directional audit subset of pooled Top 30"
                table["rank_within_direction"] = np.arange(1, len(table) + 1)
                table.reset_index(drop=True, inplace=True)
            direction_ordered = order_directional_display(pooled)
            if set(direction_ordered["ensembl_id"]) != set(pooled["ensembl_id"]):
                raise AssertionError("direction-ordered display changed the pooled Top 30 set")
            pooled_path = TABLES / (
                f"{group_id}_DDR_common_direction_top30_{rank_method}_pooled.csv"
            )
            direction_ordered_path = TABLES / (
                f"{group_id}_DDR_common_direction_top30_{rank_method}_direction_ordered.csv"
            )
            separate_up_path = TABLES / (
                f"{group_id}_DDR_common_direction_top30_{rank_method}_up.csv"
            )
            separate_down_path = TABLES / (
                f"{group_id}_DDR_common_direction_top30_{rank_method}_down.csv"
            )
            pooled.to_csv(pooled_path, index=False, encoding="utf-8-sig")
            direction_ordered.to_csv(direction_ordered_path, index=False, encoding="utf-8-sig")
            up.to_csv(separate_up_path, index=False, encoding="utf-8-sig")
            down.to_csv(separate_down_path, index=False, encoding="utf-8-sig")

            for display_mode, table in (("pooled", pooled), ("direction_ordered", direction_ordered)):
                for rank, record in enumerate(table.to_dict("records"), start=1):
                    display_rows.append(
                        {
                            "group_id": group_id,
                            "rank_method": rank_method,
                            "display_mode": display_mode,
                            "display_rank": rank,
                            "ensembl_id": record["ensembl_id"],
                            "display_symbol": record["display_symbol"],
                            "common_direction": record["common_direction"],
                            "final_fdr": record["final_fdr"],
                            "final_integrated_p": record["final_integrated_p"],
                            "ranking_score": record["ranking_score"],
                            "ranking_score_name": rank_method,
                        }
                    )
            summary_rows.append(
                {
                    **summary_base,
                    "rank_method": rank_method,
                    "rank_definition": info["label"],
                    "pooled_top30_n": int(len(pooled)),
                    "pooled_top30_up_n": int(pooled["common_direction"].eq("common_up").sum()),
                    "pooled_top30_down_n": int(pooled["common_direction"].eq("common_down").sum()),
                    "direction_ordered_n": int(len(direction_ordered)),
                    "direction_ordered_up_n": int(direction_ordered["common_direction"].eq("common_up").sum()),
                    "direction_ordered_down_n": int(direction_ordered["common_direction"].eq("common_down").sum()),
                    "selection_order": "group-specific q filter -> common direction -> descending |log2FC| score -> global Top 30; direction-ordered figure then displays common up before common down",
                }
            )
            figure_requests.extend(
                [
                    {
                        "group_id": group_id,
                        "rank_method": rank_method,
                        "rank_label": info["label"],
                        "unit_order": unit_order,
                        "labels": labels,
                        "cutoff": cutoff,
                        "pooled": pooled,
                        "direction_ordered": direction_ordered,
                    }
                ]
            )

        for position, unit_id in enumerate(unit_order, start=1):
            if group_id == "A_kidney":
                analysis_tissue = "Kidney"
                nes = float(nes_lookup["Kidney"])
                n_units = int(
                    kidney_audit.loc[
                        kidney_audit["kidney_accession"].eq(unit_id), "n_analysis_units"
                    ].iloc[0]
                )
                n_missions = 1
            elif group_id == "B_thymus":
                analysis_tissue = "Thymus"
                nes = float(nes_lookup["Thymus"])
                row = thymus_audit.loc[thymus_audit["dataset_id"].eq(unit_id)].iloc[0]
                n_units = 1
                n_missions = 1
                del row
            else:
                analysis_tissue = unit_id
                nes = float(nes_lookup[unit_id])
                audit = tissue_audit.loc[
                    tissue_audit["analysis_tissue"].eq(unit_id)
                ].iloc[0]
                n_units = int(audit["n_analysis_units"])
                n_missions = int(audit["n_missions"])
            column_rows.append(
                {
                    "group_id": group_id,
                    "column_position": position,
                    "column_id": unit_id,
                    "analysis_tissue": analysis_tissue,
                    "DDR_NES": nes,
                    "n_analysis_units": n_units,
                    "n_missions": n_missions,
                    "group_definition": GROUP_TITLES[group_id],
                }
            )

    all_values = []
    for request in figure_requests:
        unit_order = request["unit_order"]
        for table in (request["pooled"], request["direction_ordered"]):
            if len(table):
                all_values.append(
                    table[[f"log2fc__{unit}" for unit in unit_order]].to_numpy(float).ravel()
                )
    limit = max(1.0, float(np.nanpercentile(np.abs(np.concatenate(all_values)), 98)))
    norm = Normalize(vmin=-limit, vmax=limit)
    cmap = signed_cmap()
    pd.DataFrame(
        [
            {
                "colour_scale": "shared across all four groups and both rank methods",
                "vmin": -limit,
                "vmax": limit,
                "calculation": "max(1.0, 98th percentile of absolute displayed log2FC values)",
            }
        ]
    ).to_csv(TABLES / "07_figure_colour_scale_audit.csv", index=False, encoding="utf-8-sig")

    manifest_rows = []
    for request in figure_requests:
        group_id = request["group_id"]
        rank_method = request["rank_method"]
        rank_label = request["rank_label"]
        pooled_stem = f"{group_id}_DDR_common_direction_{rank_method}_pooled_top30"
        direction_ordered_stem = f"{group_id}_DDR_common_direction_{rank_method}_direction_ordered_top30"
        pooled_png, pooled_pdf = plot_pooled(
            request["pooled"],
            request["unit_order"],
            request["labels"],
            GROUP_TITLES[group_id],
            request["cutoff"],
            rank_label,
            rank_method,
            pooled_stem,
            norm,
            cmap,
        )
        direction_png, direction_pdf = plot_direction_ordered(
            request["direction_ordered"],
            request["unit_order"],
            request["labels"],
            GROUP_TITLES[group_id],
            request["cutoff"],
            rank_label,
            rank_method,
            direction_ordered_stem,
            norm,
            cmap,
        )
        manifest_rows.extend(
            [
                {
                    "group_id": group_id,
                    "rank_method": rank_method,
                    "display_mode": "pooled",
                    "n_rows": len(request["pooled"]),
                    "png": str(pooled_png.relative_to(OUT)),
                    "pdf": str(pooled_pdf.relative_to(OUT)),
                },
                {
                    "group_id": group_id,
                    "rank_method": rank_method,
                    "display_mode": "direction_ordered_single_panel",
                    "n_rows": len(request["direction_ordered"]),
                    "n_rows_up": int(request["direction_ordered"]["common_direction"].eq("common_up").sum()),
                    "n_rows_down": int(request["direction_ordered"]["common_direction"].eq("common_down").sum()),
                    "png": str(direction_png.relative_to(OUT)),
                    "pdf": str(direction_pdf.relative_to(OUT)),
                },
            ]
        )

    pd.DataFrame(summary_rows).to_csv(
        TABLES / "04_group_rank_summary.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(display_rows).to_csv(
        TABLES / "05_all_displayed_genes_membership.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(column_rows).to_csv(
        TABLES / "06_group_column_order.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(manifest_rows).to_csv(
        TABLES / "08_figure_manifest.csv", index=False, encoding="utf-8-sig"
    )

    provenance_rows = [
        {
            "item": "gene_statistics",
            "path": str(builder.LONG_STATS),
            "sha256": sha256(builder.LONG_STATS),
            "note": "stable Ensembl long table",
        },
        {
            "item": "tissue_DDR_NES",
            "path": str(builder.NES_MATRIX),
            "sha256": sha256(builder.NES_MATRIX),
            "note": "explicit GO tissue mission-equal NES matrix",
        },
        {
            "item": "analysis_unit_tissue_mapping",
            "path": str(builder.UNIT_NES),
            "sha256": sha256(builder.UNIT_NES),
            "note": "analysis-unit tissue and mission mapping",
        },
        {
            "item": "GO_membership",
            "path": str(builder.GO_MEMBERSHIP),
            "sha256": sha256(builder.GO_MEMBERSHIP),
            "note": "two-GO Ensembl membership used by the stable-ID rerun",
        },
    ]
    pd.DataFrame(provenance_rows).to_csv(
        TABLES / "09_input_provenance.csv", index=False, encoding="utf-8-sig"
    )

    group_lines = []
    for group_id in GROUP_ORDER:
        cutoff = GROUP_CUTOFFS[group_id]
        summary = pd.DataFrame(summary_rows).loc[
            lambda frame: frame["group_id"].eq(group_id)
        ].iloc[0]
        columns = "; ".join(
            tissue[len("Heart / "):] if tissue.startswith("Heart / ") else tissue
            for tissue in group_inputs[group_id][1]
        )
        group_lines.append(
            f"- **{GROUP_TITLES[group_id]}**：q < {cutoff:.2f}；列（固定顺序）：{columns}；"
            f"共同方向候选 {int(summary['common_direction_n'])} 个（上调 {int(summary['common_up_n'])}、下调 {int(summary['common_down_n'])}）。"
        )
    readme = f"""# 四组小鼠 DDR 共同方向热图（|log2FC| mean/min）

## 本版规则

按组先应用不同的最终 BH-q 阈值，再筛选共同方向，最后按 `|log2FC|` 排序：

1. 在预先指定的两个 GO term（GO:0006974 与 GO:0006281）并集、且各列都有完整数据的 Ensembl Gene ID 背景中，计算每组的最终 integrated signed-Z、双侧 P 和 BH-q。
2. A/B 组保留最终 q < 0.05；C/D 组保留最终 q < 0.10。
3. 在显著基因中，只保留该组所有列的 log2FC 均大于 0（common up）或均小于 0（common down）的基因。
4. 对共同方向基因计算两个排序分数：跨列 `mean(|log2FC|)` 和跨列 `min(|log2FC|)`；分数从大到小排序，q、最终 P、Ensembl ID 仅作为并列时的次级排序。

每组输出四张可见图（每张图同时提供 PNG/PDF）：

- `mean_abs_log2fc_pooled_top30`：上下调合并后全局 Top 30；
- `mean_abs_log2fc_direction_ordered_top30`：仍为一张完整热图，使用同一批全局 Top 30，先排共同上调、再排共同下调，方向内按分数排序；
- `min_abs_log2fc_pooled_top30`：上下调合并后全局 Top 30；
- `min_abs_log2fc_direction_ordered_top30`：仍为一张完整热图，使用同一批全局 Top 30，先排共同上调、再排共同下调，方向内按分数排序。

## 四组定义

{chr(10).join(group_lines)}

- C 的 NES 条件为 `< -1`，并去除 Lung、Thymus。
- D 的 NES 条件为 `> 1`，并去除 Kidney。
- C/D 组织列沿用此前 26 组织 DDR-NES 图的固定显示顺序，不在本脚本中重新排序；NES 只用于定义分组和列顺序，不用于基因显著性或活性判定。

## 聚合和 ID

- A（Kidney）：OSD-163、OSD-253、OSD-462、OSD-513、OSD-771 五个 accession；同一 accession 内的子分析等权合并，最终五列。
- B（Thymus）：OSD-289 MHU-1、OSD-289 MHU-2、OSD-421、OSD-515 四个比较；沿用四胸腺的 mission 聚类和等权整合。
- C/D：同一组织内先在每个 mission 对分析单位取均值，再对 mission 等权，每个组织形成一列。
- 全部分析和 GO 连接使用 Ensembl Gene ID；Symbol 仅用于图中显示，缺失 Symbol 时显示 Ensembl ID。
- GO 注释沿用 MGI MOUSE-mod GAF（date-generated 2026-08-04）、GO basic release 2026-07-26，以及 `is_a`、`part_of`、`regulates`、`positively_regulates`、`negatively_regulates` 五关系闭包；小鼠 taxon 10090，排除 NOT。

## 输出

- `figures/`：4 组 × 2 种排序 × 2 种展示，共 16 张 PNG 和 16 张 PDF。
- `tables/04_group_rank_summary.csv`：按组、排序方法汇总阈值、候选数和 Top 30 数量。
- `tables/05_all_displayed_genes_membership.csv`：所有图中基因的 Ensembl、Symbol、方向、q 和排序分数。
- `tables/*_up.csv`、`tables/*_down.csv`：全局 Top 30 的方向审计子集，不代表分别重新取 Top 30。
- `tables/06_group_column_order.csv`：组内列顺序、组织 NES 和聚合审计。
- `tables/08_figure_manifest.csv`：每张图的组、排序方式、展示方式和文件名。
- `tables/09_input_provenance.csv`：输入路径和 SHA256。

上一版七组 pooled 图目录 `05b_ddr_common_direction_tissue_groups_20260826` 与 `05c_ddr_common_direction_tissue_groups_abs_integrated_z_20260826` 保留不变；本目录是本次四组、双 `|log2FC|` 排序规则的独立版本。
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
