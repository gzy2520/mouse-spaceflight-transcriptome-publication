#!/usr/bin/env python3
"""Build pooled common-direction DDR heatmaps for seven tissue groups.

The existing four-thymus figure is intentionally left untouched.  This script
adds a versioned extension in which common-up and common-down genes are pooled
and the best 30 rows are displayed.  The ranking mode is controlled by
``MOUSE_GO_DDR_GROUPED_RANK_MODE``: ``q`` (the default) ranks by final BH q,
whereas ``absolute_integrated_z`` ignores q and ranks by the absolute
gene-level integrated signed Z.

The four-thymus panel keeps one column per original comparison.  The other
panels use one column per tissue, with within-tissue analysis units aggregated
by mission (mean signed Z and mean log2FC) and missions given equal weight.
All joins and GO membership decisions use Ensembl Gene IDs; symbols are display
annotations only.
"""

from __future__ import annotations

import importlib.util
import hashlib
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-ddr-common-tissue-groups")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE_PATH = ROOT / "07_scripts/06_cross_tissue_integration/_ddr_heatmap_support.py"

RANK_MODE = os.environ.get("MOUSE_GO_DDR_GROUPED_RANK_MODE", "q").strip().lower()
if RANK_MODE not in {"q", "absolute_integrated_z"}:
    raise ValueError(
        "MOUSE_GO_DDR_GROUPED_RANK_MODE must be 'q' or 'absolute_integrated_z'"
    )
RUN_DATE = os.environ.get("MOUSE_GO_DDR_GROUPED_RUN_DATE", "20260826").strip()
if not RUN_DATE.isdigit() or len(RUN_DATE) != 8:
    raise ValueError("MOUSE_GO_DDR_GROUPED_RUN_DATE must be an YYYYMMDD string")
DEFAULT_OUT = (
    ROOT
    / "03_analysis_results/10_full_stable_id_rerun_20260823"
    / (
        f"05b_ddr_common_direction_tissue_groups_{RUN_DATE}"
        if RANK_MODE == "q"
        else f"05c_ddr_common_direction_tissue_groups_abs_integrated_z_{RUN_DATE}"
    )
)
Q_OUTPUT_DIR = (
    ROOT
    / "03_analysis_results/10_full_stable_id_rerun_20260823"
    / f"05b_ddr_common_direction_tissue_groups_{RUN_DATE}"
)
ABS_OUTPUT_DIR = (
    ROOT
    / "03_analysis_results/10_full_stable_id_rerun_20260823"
    / f"05c_ddr_common_direction_tissue_groups_abs_integrated_z_{RUN_DATE}"
)
OUT = Path(os.environ.get("MOUSE_GO_DDR_GROUPED_OUTPUT_DIR", str(DEFAULT_OUT)))
if not OUT.is_absolute():
    OUT = ROOT / OUT
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

LONG_STATS = ROOT / (
    "03_analysis_results/10_full_stable_id_rerun_20260823/"
    "01_stable_id_registry/03_full_gene_statistics_ensembl_long.csv.gz"
)
NES_MATRIX = ROOT / (
    "03_analysis_results/10_full_stable_id_rerun_20260823/"
    "03_mouse_go_concrete_terms_20260824/tables/"
    "10_tissue_mission_equal_NES_pathway_matrix.csv"
)
UNIT_NES = ROOT / (
    "03_analysis_results/10_full_stable_id_rerun_20260823/"
    "03_mouse_go_concrete_terms_20260824/tables/"
    "11_analysis_unit_NES_pathway_matrix.csv"
)
GO_MEMBERSHIP = ROOT / (
    "03_analysis_results/10_full_stable_id_rerun_20260823/"
    "05_direct_ddr_final_fdr/tables/00_two_GO_union_membership.csv"
)

DDR_NES_COLUMN = "DNA_damage_response_GO0006974"
FINAL_FDR_CUTOFF = 0.10
DISPLAY_N = 30
RANK_FILE_TAG = "pooled_top30" if RANK_MODE == "q" else "pooled_top30_abs_integrated_z"

THYMUS_ORDER = (
    "OSD289_MHU1_thymus",
    "OSD289_MHU2_thymus",
    "OSD421_thymus",
    "OSD515_thymus",
)
THYMUS_LABELS = {
    "OSD289_MHU1_thymus": "OSD-289 MHU-1\nuG - GC",
    "OSD289_MHU2_thymus": "OSD-289 MHU-2\nuG - vivarium",
    "OSD421_thymus": "OSD-421\nFlight - Ground",
    "OSD515_thymus": "OSD-515\nFlight - Ground",
}

# The order is frozen from the displayed 26-tissue DDR NES heatmap, not from a
# re-sorting performed by this script.  This prevents accidental redefinition
# of the teacher-requested Lung-to-Thymus and Quadriceps-to-Eye intervals.
GROUP_TISSUES = {
    "G2_low_NES_Lung_to_Thymus": (
        "Lung",
        "Heart / Heart right ventricle",
        "Retina",
        "Spleen-distal",
        "Bone marrow",
        "Femoral skin",
        "Cecum",
        "Thymus",
    ),
    "G3_low_NES_excluding_Thymus": (
        "Lung",
        "Heart / Heart right ventricle",
        "Retina",
        "Spleen-distal",
        "Bone marrow",
        "Femoral skin",
        "Cecum",
    ),
    "G4_high_NES_Quadriceps_to_Eye": (
        "Quadriceps femoris",
        "Tibialis anterior",
        "Optic nerve",
        "Spleen",
        "Femoral lateral skin",
        "Left lobe of the liver",
        "Colon",
        "Cerebellum",
        "Heart",
        "Kidney",
        "Eye",
    ),
    "G5_high_NES_top4": (
        "Quadriceps femoris",
        "Tibialis anterior",
        "Optic nerve",
        "Spleen",
    ),
    "G7_low_NES_four_tissues": (
        "Spleen-distal",
        "Bone marrow",
        "Femoral skin",
        "Cecum",
    ),
}

GROUP_TITLES = {
    "G1_four_thymus": "Group 1 | Four thymus comparisons",
    "G2_low_NES_Lung_to_Thymus": "Group 2 | DDR NES interval: Lung → Thymus",
    "G3_low_NES_excluding_Thymus": "Group 3 | Group 2 minus Thymus",
    "G4_high_NES_Quadriceps_to_Eye": "Group 4 | DDR NES interval: Quadriceps femoris → Eye",
    "G5_high_NES_top4": "Group 5 | Four highest-DDR-NES tissues",
    "G6_kidney_five_datasets": "Group 6 | Kidney five-dataset representative",
    "G7_low_NES_four_tissues": "Group 7 | Four low-DDR-NES tissues",
}

KIDNEY_ACCESSION_ORDER = (
    "OSD-163",
    "OSD-253",
    "OSD-462",
    "OSD-513",
    "OSD-771",
)
KIDNEY_LABELS = {
    "OSD-163": "OSD-163\nSpaceX-8",
    "OSD-253": "OSD-253\nSpaceX-15",
    "OSD-462": "OSD-462\nSpaceX-21",
    "OSD-513": "OSD-513\nSpaceX-21",
    "OSD-771": "OSD-771\nSpaceX-18",
}


def load_base_module():
    spec = importlib.util.spec_from_file_location("ddr_heatmap_support_grouped", BASE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(BASE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_ensembl(values: pd.Series) -> pd.Series:
    return (
        values.astype(str)
        .str.strip()
        .str.upper()
        .str.replace(r"\.\d+$", "", regex=True)
    )


def normalize_symbol(values: pd.Series) -> pd.Series:
    return values.fillna("").astype(str).str.strip().str.upper()


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


def signed_cmap() -> LinearSegmentedColormap:
    cmap = LinearSegmentedColormap.from_list(
        "negative_blue_positive_red",
        ["#3775BA", "#F7F7F7", "#B64342"],
    )
    cmap.set_bad("#D9D9D9")
    return cmap


def load_go_sets() -> tuple[dict[str, set[str]], pd.DataFrame]:
    membership = pd.read_csv(GO_MEMBERSHIP, encoding="utf-8-sig")
    required = {
        "ensembl_id",
        "in_DNA_damage_response_GO0006974",
        "in_DNA_repair_GO0006281",
    }
    missing = required - set(membership.columns)
    if missing:
        raise ValueError(f"GO membership is missing columns: {sorted(missing)}")
    membership["ensembl_id"] = normalize_ensembl(membership["ensembl_id"])

    def as_bool(values: pd.Series) -> pd.Series:
        if pd.api.types.is_bool_dtype(values):
            return values.fillna(False)
        return (
            values.astype(str)
            .str.strip()
            .str.lower()
            .isin({"true", "1", "yes", "y"})
        )

    membership["_damage_flag"] = as_bool(
        membership["in_DNA_damage_response_GO0006974"]
    )
    membership["_repair_flag"] = as_bool(membership["in_DNA_repair_GO0006281"])
    damage = set(
        membership.loc[
            membership["_damage_flag"],
            "ensembl_id",
        ]
    )
    repair = set(
        membership.loc[membership["_repair_flag"], "ensembl_id"]
    )
    if not damage or not repair:
        raise ValueError("Empty DDR or DNA-repair Ensembl gene set")
    return {
        "DNA Damage Response (GO:0006974)": damage,
        "DNA Repair (GO:0006281)": repair,
    }, membership


def deduplicate_unit_frame(frame: pd.DataFrame, unit_id: str) -> pd.DataFrame:
    frame = frame.loc[
        ~frame["ensembl_id"].isin(["", "NAN"])
        & frame["log2fc"].notna()
        & frame["wald_z"].notna()
    ].copy()
    duplicated = frame["ensembl_id"].duplicated(keep=False)
    if duplicated.any():
        checks = frame.loc[duplicated].groupby("ensembl_id").agg(
            n_log2fc=("log2fc", "nunique"),
            n_z=("wald_z", "nunique"),
            n_fdr=("sample_fdr", "nunique"),
        )
        if (checks[["n_log2fc", "n_z", "n_fdr"]] > 1).any(axis=None):
            raise ValueError(f"Non-identical duplicate Ensembl rows in {unit_id}")
        frame = frame.drop_duplicates("ensembl_id", keep="first")
    return frame.set_index("ensembl_id").sort_index()


def load_unit_frames(unit_ids: set[str]) -> dict[str, pd.DataFrame]:
    usecols = [
        "ensembl_id",
        "gene_symbol",
        "analysis_unit_id",
        "log2fc_spaceflight_vs_ground",
        "stat_spaceflight_vs_ground",
        "adjusted_p_value",
    ]
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(LONG_STATS, usecols=usecols, chunksize=200_000):
        part = chunk.loc[chunk["analysis_unit_id"].isin(unit_ids)]
        if not part.empty:
            pieces.append(part)
    if not pieces:
        raise ValueError("No stable-ID gene statistics matched the tissue units")
    long = pd.concat(pieces, ignore_index=True)
    frames: dict[str, pd.DataFrame] = {}
    for unit_id in sorted(unit_ids):
        part = long.loc[long["analysis_unit_id"].eq(unit_id)].copy()
        if part.empty:
            raise ValueError(f"Missing gene statistics for {unit_id}")
        frame = pd.DataFrame(
            {
                "ensembl_id": normalize_ensembl(part["ensembl_id"]),
                "gene_symbol": normalize_symbol(part["gene_symbol"]),
                "log2fc": pd.to_numeric(
                    part["log2fc_spaceflight_vs_ground"], errors="coerce"
                ),
                "wald_z": pd.to_numeric(
                    part["stat_spaceflight_vs_ground"], errors="coerce"
                ),
                "sample_fdr": pd.to_numeric(
                    part["adjusted_p_value"], errors="coerce"
                ),
            }
        )
        frames[unit_id] = deduplicate_unit_frame(frame, unit_id)
    return frames


def resolve_symbol(values: list[str]) -> str:
    valid = [value for value in values if value and value != "NAN"]
    if not valid:
        return ""
    counts = pd.Series(valid).value_counts()
    best_count = counts.max()
    return sorted(counts.index[counts.eq(best_count)].tolist())[0]


def build_tissue_frames(
    unit_frames: dict[str, pd.DataFrame],
    unit_meta: pd.DataFrame,
    tissue_order: tuple[str, ...],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Aggregate each tissue with equal mission weight.

    A gene is retained in a tissue only when it is present with finite log2FC
    and signed Z in every underlying analysis unit for that tissue.  This is
    deliberately conservative and makes the later cross-tissue common
    background auditable.
    """
    frames: dict[str, pd.DataFrame] = {}
    audit_rows: list[dict[str, object]] = []
    for tissue in tissue_order:
        meta = unit_meta.loc[unit_meta["analysis_tissue"].eq(tissue)].copy()
        if meta.empty:
            raise ValueError(f"No analysis units for tissue {tissue}")
        units = meta["analysis_unit_id"].tolist()
        missions = meta["mission_cluster"].drop_duplicates().tolist()
        common = set.intersection(*(set(unit_frames[u].index) for u in units))
        if not common:
            raise ValueError(f"No complete gene background for tissue {tissue}")
        common_ids = sorted(common)
        mission_log2fc: dict[str, pd.Series] = {}
        mission_z: dict[str, pd.Series] = {}
        mission_fdr: dict[str, pd.Series] = {}
        symbols: dict[str, list[str]] = defaultdict(list)
        for unit_id in units:
            frame = unit_frames[unit_id].loc[common_ids]
            for gene_id, symbol in frame["gene_symbol"].items():
                symbols[gene_id].append(str(symbol))
        for mission in missions:
            mission_units = meta.loc[
                meta["mission_cluster"].eq(mission), "analysis_unit_id"
            ].tolist()
            log2fc_values = pd.concat(
                [unit_frames[u].loc[common_ids, "log2fc"] for u in mission_units],
                axis=1,
            )
            z_values = pd.concat(
                [unit_frames[u].loc[common_ids, "wald_z"] for u in mission_units],
                axis=1,
            )
            fdr_values = pd.concat(
                [unit_frames[u].loc[common_ids, "sample_fdr"] for u in mission_units],
                axis=1,
            )
            mission_log2fc[mission] = log2fc_values.mean(axis=1)
            mission_z[mission] = z_values.mean(axis=1)
            mission_fdr[mission] = fdr_values.median(axis=1)
        mission_log2fc_wide = pd.DataFrame(mission_log2fc, index=common_ids)
        mission_z_wide = pd.DataFrame(mission_z, index=common_ids)
        mission_fdr_wide = pd.DataFrame(mission_fdr, index=common_ids)
        tissue_frame = pd.DataFrame(index=common_ids)
        tissue_frame["gene_symbol"] = [resolve_symbol(symbols[g]) for g in common_ids]
        tissue_frame["log2fc"] = mission_log2fc_wide.mean(axis=1)
        tissue_frame["wald_z"] = mission_z_wide.sum(axis=1) / np.sqrt(len(missions))
        tissue_frame["sample_fdr"] = mission_fdr_wide.mean(axis=1)
        frames[tissue] = tissue_frame
        audit_rows.append(
            {
                "analysis_tissue": tissue,
                "n_analysis_units": len(units),
                "analysis_units": ";".join(units),
                "n_missions": len(missions),
                "missions": ";".join(missions),
                "complete_gene_background": len(common_ids),
                "tissue_log2fc_rule": "mean of mission-level mean log2FC",
                "tissue_signed_z_rule": "sum mission-level mean signed Z / sqrt(n missions)",
                "mission_weight_rule": "equal weight",
            }
        )
    return frames, pd.DataFrame(audit_rows)


def build_kidney_accession_frames(
    unit_frames: dict[str, pd.DataFrame],
    unit_meta: pd.DataFrame,
    accession_order: tuple[str, ...],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Collapse Kidney analysis units into the five requested accessions.

    OSD-253 and OSD-771 each have documented subtype/age analysis units.  They
    are averaged within accession so that the standalone Kidney panel has one
    column per OSDR dataset rather than giving one accession extra visual or
    statistical weight.  All units for an accession belong to one mission in
    the current manifest.
    """
    frames: dict[str, pd.DataFrame] = {}
    audit_rows: list[dict[str, object]] = []
    kidney_meta = unit_meta.loc[unit_meta["analysis_tissue"].eq("Kidney")].copy()
    for accession in accession_order:
        meta = kidney_meta.loc[kidney_meta["accession"].eq(accession)].copy()
        if meta.empty:
            raise ValueError(f"No Kidney analysis units for accession {accession}")
        units = meta["analysis_unit_id"].tolist()
        missions = meta["mission_cluster"].drop_duplicates().tolist()
        if len(missions) != 1:
            raise ValueError(
                f"Kidney accession {accession} spans multiple missions: {missions}"
            )
        common = set.intersection(*(set(unit_frames[u].index) for u in units))
        if not common:
            raise ValueError(f"No complete gene background for Kidney {accession}")
        common_ids = sorted(common)
        log2fc_values = pd.concat(
            [unit_frames[u].loc[common_ids, "log2fc"] for u in units], axis=1
        )
        z_values = pd.concat(
            [unit_frames[u].loc[common_ids, "wald_z"] for u in units], axis=1
        )
        fdr_values = pd.concat(
            [unit_frames[u].loc[common_ids, "sample_fdr"] for u in units], axis=1
        )
        symbols: dict[str, list[str]] = defaultdict(list)
        for unit_id in units:
            for gene_id, symbol in unit_frames[unit_id].loc[
                common_ids, "gene_symbol"
            ].items():
                symbols[gene_id].append(str(symbol))
        frame = pd.DataFrame(index=common_ids)
        frame["gene_symbol"] = [resolve_symbol(symbols[g]) for g in common_ids]
        frame["log2fc"] = log2fc_values.mean(axis=1)
        frame["wald_z"] = z_values.mean(axis=1)
        frame["sample_fdr"] = fdr_values.median(axis=1)
        frames[accession] = frame
        audit_rows.append(
            {
                "kidney_accession": accession,
                "analysis_tissue": "Kidney",
                "analysis_units": ";".join(units),
                "n_analysis_units": len(units),
                "mission_cluster": missions[0],
                "complete_gene_background": len(common_ids),
                "aggregation_rule": "mean log2FC and mean signed Z across accession units; median descriptive sample FDR",
                "unit_weight_rule": "equal weight within accession",
            }
        )
    return frames, pd.DataFrame(audit_rows)


def select_top_pooled_common_direction(
    selected: pd.DataFrame, n: int = DISPLAY_N
) -> pd.DataFrame:
    common = selected.loc[
        selected["effect_pattern"].isin(["all_positive", "all_negative"])
    ].copy()
    common["common_direction"] = np.where(
        common["effect_pattern"].eq("all_positive"), "common_up", "common_down"
    )
    if RANK_MODE == "q":
        common["ranking_score"] = common["final_fdr"]
        common["ranking_score_name"] = "final_fdr"
        common = common.sort_values(
            ["final_fdr", "final_integrated_p", "ensembl_id"],
            ascending=[True, True, True],
            kind="mergesort",
        ).head(n)
        selection_rule = (
            "common up/down pooled; final FDR ascending, then final P and Ensembl ID; "
            f"top {n}"
        )
    else:
        # There is no gene-level GSEA NES in the direct DDR tables.  The
        # available gene-level signed integration score is the Stouffer
        # statistic, which is retained with its sign for common-direction
        # filtering and ranked by absolute magnitude only.  q is not used for
        # filtering or ranking in this mode.
        common["ranking_score"] = common["final_integrated_z"].abs()
        common["ranking_score_name"] = "abs_final_integrated_z"
        common = common.sort_values(
            ["ranking_score", "ensembl_id"],
            ascending=[False, True],
            kind="mergesort",
        ).head(n)
        selection_rule = (
            "common up/down pooled; no q filter; absolute final integrated Z "
            f"descending, then Ensembl ID; top {n}"
        )
    common["pooled_common_rank"] = np.arange(1, len(common) + 1)
    common["heatmap_selection_rule"] = selection_rule
    return common.reset_index(drop=True)


def plot_heatmap(
    selected: pd.DataFrame,
    unit_order: tuple[str, ...],
    labels: dict[str, str],
    title: str,
    stem: str,
    shared_limit: float,
) -> None:
    values = selected[[f"log2fc__{unit}" for unit in unit_order]].to_numpy(float)
    n_rows = len(selected)
    n_cols = len(unit_order)
    fig, ax = plt.subplots(
        figsize=(max(10.5, 1.25 * n_cols + 4.5), max(7.2, 3.0 + 0.30 * n_rows))
    )
    image = ax.imshow(
        np.ma.masked_invalid(values),
        aspect="auto",
        interpolation="none",
        cmap=signed_cmap(),
        vmin=-shared_limit,
        vmax=shared_limit,
    )
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([labels[unit] for unit in unit_order], rotation=35, ha="right")
    ax.set_yticks(np.arange(n_rows))
    gene_names = selected["display_symbol"].where(
        selected["display_symbol"].ne(""), selected["ensembl_id"]
    )
    if RANK_MODE == "q":
        row_labels = [
            f"{symbol} | q={q:.2g}"
            for symbol, q in zip(gene_names, selected["final_fdr"])
        ]
    else:
        row_labels = [
            f"{symbol} | |Z|={score:.2f}"
            for symbol, score in zip(gene_names, selected["ranking_score"])
        ]
    ax.set_yticklabels(
        row_labels,
        fontsize=7.4,
    )
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", length=1.8, pad=1.5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, loc="left", weight="bold", pad=11)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    colorbar.set_label("Unified log2FC: spaceflight - ground/control")

    def format_value(value: float) -> str:
        # Avoid displaying a tiny negative value as the visually confusing
        # string "-0.00" while retaining the unrounded value in the CSV.
        return "0.00" if abs(value) < 0.005 else f"{value:.2f}"

    for row in range(n_rows):
        for column in range(n_cols):
            value = values[row, column]
            if np.isfinite(value):
                ax.text(
                    column,
                    row,
                    format_value(value),
                    ha="center",
                    va="center",
                    fontsize=6.0,
                    color="white" if abs(value) > 0.65 * shared_limit else "#222222",
                )
    fig.tight_layout(pad=1.0)
    for suffix in ("png", "pdf"):
        fig.savefig(
            FIGURES / f"{stem}.{suffix}",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.08,
            facecolor="white",
        )
    plt.close(fig)


def group_summary_row(
    group_id: str,
    result: pd.DataFrame,
    selected: pd.DataFrame,
    display: pd.DataFrame,
    unit_order: tuple[str, ...],
    common_background: int,
) -> dict[str, object]:
    return {
        "group_id": group_id,
        "group_title": GROUP_TITLES[group_id],
        "rank_mode": RANK_MODE,
        "n_columns": len(unit_order),
        "column_ids": ";".join(unit_order),
        "common_background": int(common_background),
        "two_GO_candidate_family": int(len(result)),
        "final_FDR_lt_0.10": int(len(selected)),
        "final_FDR_lt_0.10_common_up": int(
            selected["effect_pattern"].eq("all_positive").sum()
        ),
        "final_FDR_lt_0.10_common_down": int(
            selected["effect_pattern"].eq("all_negative").sum()
        ),
        "ranking_pool_n": int(len(result) if RANK_MODE == "absolute_integrated_z" else len(selected)),
        "ranking_pool_common_up": int(
            (result if RANK_MODE == "absolute_integrated_z" else selected)["effect_pattern"].eq("all_positive").sum()
        ),
        "ranking_pool_common_down": int(
            (result if RANK_MODE == "absolute_integrated_z" else selected)["effect_pattern"].eq("all_negative").sum()
        ),
        "pooled_common_top30_n": int(len(display)),
        "pooled_common_top30_up": int(display["common_direction"].eq("common_up").sum()),
        "pooled_common_top30_down": int(
            display["common_direction"].eq("common_down").sum()
        ),
        "min_display_q": float(display["final_fdr"].min()) if len(display) else np.nan,
        "max_display_q": float(display["final_fdr"].max()) if len(display) else np.nan,
        "min_display_abs_integrated_z": float(display["ranking_score"].min()) if len(display) and RANK_MODE == "absolute_integrated_z" else np.nan,
        "max_display_abs_integrated_z": float(display["ranking_score"].max()) if len(display) and RANK_MODE == "absolute_integrated_z" else np.nan,
    }


def main() -> None:
    for path in (LONG_STATS, NES_MATRIX, UNIT_NES, GO_MEMBERSHIP):
        if not path.exists():
            raise FileNotFoundError(path)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style()

    base = load_base_module()
    base.GO_SET_ID_TYPE = "ensembl_gene_id"
    go_sets, membership = load_go_sets()
    membership.drop(columns=["_damage_flag", "_repair_flag"], errors="ignore").to_csv(
        TABLES / "00_two_GO_union_membership_used.csv",
        index=False,
        encoding="utf-8-sig",
    )

    tissue_nes = pd.read_csv(NES_MATRIX, encoding="utf-8-sig")
    unit_nes = pd.read_csv(UNIT_NES, encoding="utf-8-sig")
    required_nes = {"analysis_tissue", DDR_NES_COLUMN}
    if required_nes - set(tissue_nes.columns):
        raise ValueError("Tissue NES table lacks the DDR column")
    if {"analysis_unit_id", "analysis_tissue", "mission_cluster"} - set(unit_nes.columns):
        raise ValueError("Analysis-unit NES table lacks tissue/mission mapping")
    tissue_nes = tissue_nes.drop_duplicates("analysis_tissue", keep="first")
    unit_nes = unit_nes.drop_duplicates("analysis_unit_id", keep="first")
    unit_ids = set(unit_nes["analysis_unit_id"])
    unit_frames = load_unit_frames(unit_ids)
    tissue_frames, tissue_audit = build_tissue_frames(
        unit_frames,
        unit_nes,
        tuple(tissue_nes["analysis_tissue"]),
    )
    tissue_audit = tissue_audit.merge(
        tissue_nes[["analysis_tissue", DDR_NES_COLUMN]],
        on="analysis_tissue",
        how="left",
        validate="one_to_one",
    ).rename(columns={DDR_NES_COLUMN: "DDR_NES"})
    tissue_audit.to_csv(TABLES / "01_tissue_meta_aggregation_audit.csv", index=False, encoding="utf-8-sig")
    kidney_frames, kidney_audit = build_kidney_accession_frames(
        unit_frames,
        unit_nes,
        KIDNEY_ACCESSION_ORDER,
    )
    kidney_audit.to_csv(
        TABLES / "02_kidney_accession_aggregation_audit.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # The first panel is the established four-comparison thymus input.
    thymus_frames, thymus_audit = base.load_four_thymus()
    thymus_groups = {
        mission: tuple(group["dataset_id"])
        for mission, group in thymus_audit.groupby("mission_cluster", sort=False)
    }
    group_inputs: dict[str, tuple[dict[str, pd.DataFrame], tuple[str, ...], dict[str, tuple[str, ...]]]] = {
        "G1_four_thymus": (thymus_frames, THYMUS_ORDER, thymus_groups),
    }
    for group_id, tissues in GROUP_TISSUES.items():
        missing = [t for t in tissues if t not in tissue_frames]
        if missing:
            raise ValueError(f"{group_id} missing tissue frames: {missing}")
        frames = {t: tissue_frames[t] for t in tissues}
        # Each tissue column is already mission-equal; one synthetic mission per
        # tissue keeps build_final_results' Stouffer implementation exact.
        groups = {t: (t,) for t in tissues}
        group_inputs[group_id] = (frames, tissues, groups)
    group_inputs["G6_kidney_five_datasets"] = (
        kidney_frames,
        KIDNEY_ACCESSION_ORDER,
        {accession: (accession,) for accession in KIDNEY_ACCESSION_ORDER},
    )
    # Keep the report and figure order aligned with the numbered definitions;
    # G7 is stored in GROUP_TISSUES but is intentionally rendered after G6.
    group_order = (
        "G1_four_thymus",
        "G2_low_NES_Lung_to_Thymus",
        "G3_low_NES_excluding_Thymus",
        "G4_high_NES_Quadriceps_to_Eye",
        "G5_high_NES_top4",
        "G6_kidney_five_datasets",
        "G7_low_NES_four_tissues",
    )
    group_inputs = {group_id: group_inputs[group_id] for group_id in group_order}

    all_displays: dict[str, pd.DataFrame] = {}
    all_results: dict[str, pd.DataFrame] = {}
    all_selected: dict[str, pd.DataFrame] = {}
    group_summaries: list[dict[str, object]] = []
    audit_membership: list[dict[str, object]] = []

    for group_id, (frames, unit_order, mission_groups) in group_inputs.items():
        result, selected, summary = base.build_final_results(
            frames, unit_order, mission_groups, go_sets
        )
        ranking_pool = result if RANK_MODE == "absolute_integrated_z" else selected
        display = select_top_pooled_common_direction(ranking_pool, DISPLAY_N)
        all_results[group_id] = result
        all_selected[group_id] = selected
        all_displays[group_id] = display
        result.to_csv(TABLES / f"{group_id}_DDR_final_test_all.csv", index=False, encoding="utf-8-sig")
        selected.to_csv(TABLES / f"{group_id}_DDR_final_FDR01_selected.csv", index=False, encoding="utf-8-sig")
        display.to_csv(
            TABLES / f"{group_id}_DDR_common_direction_{RANK_FILE_TAG}.csv",
            index=False,
            encoding="utf-8-sig",
        )
        row = group_summary_row(
            group_id,
            result,
            selected,
            display,
            unit_order,
            int(summary["common_complete_background"]),
        )
        row.update({f"algorithm_{key}": value for key, value in summary.items()})
        group_summaries.append(row)

        for rank, record in enumerate(display.to_dict("records"), start=1):
            audit_membership.append(
                {
                    "group_id": group_id,
                    "pooled_common_rank": rank,
                    "ensembl_id": record["ensembl_id"],
                    "display_symbol": record["display_symbol"],
                    "common_direction": record["common_direction"],
                    "final_fdr": record["final_fdr"],
                    "final_integrated_p": record["final_integrated_p"],
                    "ranking_score": record["ranking_score"],
                    "ranking_score_name": record["ranking_score_name"],
                    "effect_pattern": record["effect_pattern"],
                }
            )

    # A common colour limit is computed after the first pass so all seven panels
    # use the same log2FC scale.  This is preferable for cross-group visual
    # comparison and is recorded in the figure provenance below.
    pooled_values = []
    for display, (_, unit_order, _) in zip(
        all_displays.values(), group_inputs.values()
    ):
        if len(display):
            pooled_values.append(display[[f"log2fc__{u}" for u in unit_order]].to_numpy(float).ravel())
    if pooled_values:
        shared_limit = max(1.0, float(np.nanpercentile(np.abs(np.concatenate(pooled_values)), 98)))
    else:
        shared_limit = 1.0
    pd.DataFrame(
        [{
            "colour_scale": "shared across all seven pooled Top-30 panels",
            "vmin": -shared_limit,
            "vmax": shared_limit,
            "calculation": "max(1.0, 98th percentile of absolute displayed log2FC values)",
        }]
    ).to_csv(TABLES / "07_figure_colour_scale_audit.csv", index=False, encoding="utf-8-sig")
    for group_id, display in all_displays.items():
        _, unit_order, _ = group_inputs[group_id]
        if group_id == "G1_four_thymus":
            labels = THYMUS_LABELS
            subtitle = (
                "All 4 comparison log2FC values share sign; pooled top 30 by "
                + ("final q (q<0.10)" if RANK_MODE == "q" else "absolute integrated Z")
            )
        elif group_id == "G6_kidney_five_datasets":
            labels = KIDNEY_LABELS
            subtitle = (
                "Five Kidney datasets; all columns share sign; pooled top 30 by "
                + ("final q (q<0.10)" if RANK_MODE == "q" else "absolute integrated Z")
            )
        else:
            nes_lookup = tissue_nes.set_index("analysis_tissue")[DDR_NES_COLUMN]
            labels = {t: f"{t}\nNES={float(nes_lookup[t]):.2f}" for t in unit_order}
            subtitle = (
                "Tissue-level mission-equal log2FC; all columns share sign; "
                + "pooled top 30 by "
                + ("final q (q<0.10)" if RANK_MODE == "q" else "absolute integrated Z")
            )
        plot_heatmap(
            display,
            unit_order,
            labels,
            f"{GROUP_TITLES[group_id]}\n{subtitle}",
            f"{group_id}_DDR_common_direction_{RANK_FILE_TAG}",
            shared_limit=shared_limit,
        )

    pd.DataFrame(group_summaries).to_csv(
        TABLES / "02_group_summary.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(audit_membership).to_csv(
        TABLES / "03_all_displayed_genes_membership.csv",
        index=False,
        encoding="utf-8-sig",
    )
    group_map_rows = []
    nes_lookup = tissue_nes.set_index("analysis_tissue")[DDR_NES_COLUMN]
    for group_id, tissues in GROUP_TISSUES.items():
        for position, tissue in enumerate(tissues, start=1):
            meta = tissue_audit.loc[tissue_audit["analysis_tissue"].eq(tissue)].iloc[0]
            group_map_rows.append(
                {
                    "group_id": group_id,
                    "column_position": position,
                    "column_id": tissue,
                    "analysis_tissue": tissue,
                    "DDR_NES": float(nes_lookup[tissue]),
                    "n_analysis_units": int(meta["n_analysis_units"]),
                    "n_missions": int(meta["n_missions"]),
                    "group_definition": GROUP_TITLES[group_id],
                }
            )
    kidney_nes = float(
        tissue_nes.loc[tissue_nes["analysis_tissue"].eq("Kidney"), DDR_NES_COLUMN].iloc[0]
    )
    for position, accession in enumerate(KIDNEY_ACCESSION_ORDER, start=1):
        meta = kidney_audit.loc[kidney_audit["kidney_accession"].eq(accession)].iloc[0]
        group_map_rows.append(
            {
                "group_id": "G6_kidney_five_datasets",
                "column_position": position,
                "column_id": accession,
                "analysis_tissue": "Kidney",
                "DDR_NES": kidney_nes,
                "n_analysis_units": int(meta["n_analysis_units"]),
                "n_missions": 1,
                "group_definition": GROUP_TITLES["G6_kidney_five_datasets"],
            }
        )
    pd.DataFrame(group_map_rows).to_csv(
        TABLES / "04_tissue_group_column_order.csv", index=False, encoding="utf-8-sig"
    )
    thymus_audit.to_csv(TABLES / "05_four_thymus_input_contrast_audit_used.csv", index=False, encoding="utf-8-sig")

    comparison_rows = []
    for group_id in group_inputs:
        q_path = Q_OUTPUT_DIR / "tables" / f"{group_id}_DDR_common_direction_pooled_top30.csv"
        abs_path = ABS_OUTPUT_DIR / "tables" / f"{group_id}_DDR_common_direction_pooled_top30_abs_integrated_z.csv"
        if not q_path.exists() or not abs_path.exists():
            continue
        q_table = pd.read_csv(q_path)
        abs_table = pd.read_csv(abs_path)
        q_ids = set(q_table["ensembl_id"])
        abs_ids = set(abs_table["ensembl_id"])
        q_dir = dict(zip(q_table["ensembl_id"], q_table["common_direction"]))
        abs_dir = dict(zip(abs_table["ensembl_id"], abs_table["common_direction"]))
        comparison_rows.append(
            {
                "group_id": group_id,
                "q_version_n": len(q_ids),
                "absolute_integrated_z_version_n": len(abs_ids),
                "overlap_n": len(q_ids & abs_ids),
                "q_only_n": len(q_ids - abs_ids),
                "absolute_integrated_z_only_n": len(abs_ids - q_ids),
                "q_version_up_n": sum(value == "common_up" for value in q_dir.values()),
                "q_version_down_n": sum(value == "common_down" for value in q_dir.values()),
                "absolute_integrated_z_version_up_n": sum(value == "common_up" for value in abs_dir.values()),
                "absolute_integrated_z_version_down_n": sum(value == "common_down" for value in abs_dir.values()),
                "overlap_direction_changes_n": sum(
                    q_dir[gene] != abs_dir[gene] for gene in q_ids & abs_ids
                ),
            }
        )
    pd.DataFrame(comparison_rows).to_csv(
        TABLES / "08_q_vs_absolute_integrated_z_selection_comparison.csv",
        index=False,
        encoding="utf-8-sig",
    )

    overlap = len(go_sets["DNA Damage Response (GO:0006974)"] & go_sets["DNA Repair (GO:0006281)"])
    if RANK_MODE == "q":
        readme_title = "七组小鼠 DDR 共同方向热图（按 q pooled Top 30）"
        selection_description = (
            "每组先保留最终 BH q < 0.10 的基因，再要求所有该组列的 log2FC 全部大于 0（common up）或全部小于 0（common down）；两类基因合并后按最终 q 从小到大排序，q 相同依次按最终双侧 P 值和 Ensembl Gene ID 排序，取前 30 个。"
        )
        ranking_description = "本版排序和入选使用最终 q 值；Group 2/3/4/5/7 的组织 NES 只用于冻结分组和列顺序。"
        output_note = "`tables/03_all_displayed_genes_membership.csv`：逐基因 Ensembl、Symbol、方向、q 和排序。"
    else:
        readme_title = "七组小鼠 DDR 共同方向热图（按绝对 integrated Z pooled Top 30）"
        selection_description = (
            "完全不使用 q 值进行筛选或排序。每组直接在完整的 GO 候选基因家族中要求所有该组列的 log2FC 全部大于 0（common up）或全部小于 0（common down）；两类基因合并后按基因层面的 |final_integrated_z| 从大到小排序，绝对值相同按 Ensembl Gene ID 排序，取前 30 个。"
        )
        ranking_description = "本版的基因排序分数是 |final_integrated_z|（各列 signed Z 的等权 Stouffer 整合绝对值），不是基因层面的 GSEA NES；q 仅保留在审计表中，不参与筛选、排序或作图标签。跨组织 DDR NES 仍只用于冻结分组和列顺序。"
        output_note = "`tables/03_all_displayed_genes_membership.csv`：逐基因 Ensembl、Symbol、方向、|final_integrated_z| 排序分数及审计用 q。"
    readme = f"""# {readme_title}

## 本版的关键修改

共同上调和共同下调不再分别取 15 个。{selection_description} 因此每组上下调数量可以不相等；如果某组符合条件的基因少于 30 个，则全部保留。

## 数据和 ID

- 基因统计：`{LONG_STATS.relative_to(ROOT)}`，Ensembl Gene ID 作为唯一分析键。
- GO 注释：MGI MOUSE-mod GAF，date-generated 2026-08-04；GO basic ontology release 2026-07-26；关系为 `is_a`、`part_of`、`regulates`、`positively_regulates`、`negatively_regulates`；小鼠 taxon 10090；排除 NOT；不按 evidence、日期、参考文献或注释机构额外过滤。
- 候选基因家族：GO:0006974 DNA damage response 与 GO:0006281 DNA repair 的并集；当前表中 DDR {len(go_sets['DNA Damage Response (GO:0006974)'])} 个、DNA repair {len(go_sets['DNA Repair (GO:0006281)'])} 个，重叠 {overlap} 个，合计 {len(go_sets['DNA Damage Response (GO:0006974)'] | go_sets['DNA Repair (GO:0006281)'])} 个 Ensembl Gene ID。
- Symbol 仅作图中显示；缺失 Symbol 时显示 Ensembl ID。

## 七组定义

1. Group 1：当前四个胸腺比较（OSD-289 MHU-1、OSD-289 MHU-2、OSD-421、OSD-515）。
2. Group 2：跨组织 DDR NES 图中明确的 `Lung → Thymus` 区间，共 8 个组织。
3. Group 3：Group 2 去掉 Thymus，共 7 个组织。
4. Group 4：跨组织 DDR NES 图中明确的 `Quadriceps femoris → Eye` 区间，共 11 个组织。
5. Group 5：Group 4 中 DDR NES 最高的 4 个组织：Quadriceps femoris、Tibialis anterior、Optic nerve、Spleen，使列数与四胸腺图一致。
6. Group 6：Kidney 的五个 OSDR accession（OSD-163、OSD-253、OSD-462、OSD-513、OSD-771）。OSD-253 的品系子分析、OSD-771 的年龄子分析分别在 accession 内等权合并，因此最终为五列。
7. Group 7：四个低 DDR-NES 组织：Spleen-distal、Bone marrow、Femoral skin、Cecum，列顺序固定。

Group 2/3/4/5/7 的列是组织层面值：同一组织内先在每个 mission 中对分析单位取均值，再对 mission 等权；组织 log2FC 为 mission-level log2FC 均值，组织 signed Z 为 mission-level signed Z 的 Stouffer 合并。随后各组织列等权形成该组的 gene-level integrated Z。Group 1 保留原四胸腺比较列和原 mission 聚类算法。{ranking_description}

跨组织 DDR NES 只用于冻结 Group 2/3/4/5/7 的分组和列顺序，不作为基因 q 值或通路活性判定；NES 仍表示富集方向，不等同于通路激活/抑制。

## 输出

- `figures/G1_four_thymus_DDR_common_direction_{RANK_FILE_TAG}.png/pdf`
- `figures/G2_low_NES_Lung_to_Thymus_DDR_common_direction_{RANK_FILE_TAG}.png/pdf`
- `figures/G3_low_NES_excluding_Thymus_DDR_common_direction_{RANK_FILE_TAG}.png/pdf`
- `figures/G4_high_NES_Quadriceps_to_Eye_DDR_common_direction_{RANK_FILE_TAG}.png/pdf`
- `figures/G5_high_NES_top4_DDR_common_direction_{RANK_FILE_TAG}.png/pdf`
- `figures/G6_kidney_five_datasets_DDR_common_direction_{RANK_FILE_TAG}.png/pdf`
- `figures/G7_low_NES_four_tissues_DDR_common_direction_{RANK_FILE_TAG}.png/pdf`
- `tables/02_group_summary.csv`：每组背景、q 审计统计、排序池和 pooled Top 30 上下调数量。
- {output_note}
- `tables/04_tissue_group_column_order.csv`：组织 NES、列顺序和组织聚合审计。
- `tables/02_kidney_accession_aggregation_audit.csv`：Kidney 五个 accession 的子分析合并审计。
- `tables/08_q_vs_absolute_integrated_z_selection_comparison.csv`：与上一版 q 排序结果的重叠和差异。

旧版 `02c_four_thymus_DDR_common_up_down_final_FDR01_top15_per_direction` 和 20260825 的六组结果均未覆盖；本目录是加入 Group 7 后的独立版本。
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    pd.DataFrame(
        [
            {"item": "gene_statistics", "path": str(LONG_STATS), "sha256": sha256(LONG_STATS), "note": "stable Ensembl long table"},
            {"item": "tissue_DDR_NES", "path": str(NES_MATRIX), "sha256": sha256(NES_MATRIX), "note": "explicit GO tissue mission-equal NES matrix"},
            {"item": "analysis_unit_tissue_mapping", "path": str(UNIT_NES), "sha256": sha256(UNIT_NES), "note": "analysis-unit tissue and mission mapping"},
            {"item": "GO_membership", "path": str(GO_MEMBERSHIP), "sha256": sha256(GO_MEMBERSHIP), "note": "two-GO union membership used by prior stable-ID DDR rerun"},
        ]
    ).to_csv(TABLES / "06_input_provenance.csv", index=False, encoding="utf-8-sig")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
