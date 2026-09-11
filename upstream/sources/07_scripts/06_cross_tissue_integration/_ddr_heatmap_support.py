#!/usr/bin/env python3
"""Build direct DDR log2FC heatmaps using final, not per-comparison, FDR.

For each tissue group:
1. Restrict the common analyzable gene background to the prespecified union of
   DNA damage response (GO:0006974) and DNA repair (GO:0006281).
2. Use every comparison's signed test statistic without filtering its FDR.
3. Give each mission one vote by averaging signed Z values within mission.
4. Apply BH correction to the final integrated P values within the prespecified
   two-GO candidate family.
5. Plot all genes with final FDR < 0.10; no top-N or direction-consistency rule.
"""

from __future__ import annotations

from collections import Counter
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-direct-ddr-final-fdr01")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[2]
GO_GMT = (
    ROOT / "07_scripts/config/rigorous_v2_gene_sets/GO_Biological_Process_2026.gmt"
)
def env_path(name: str, default: Path) -> Path:
    value = Path(os.environ.get(name, str(default)))
    return value if value.is_absolute() else ROOT / value


STAGE01 = env_path(
    "MOUSE_GO_STABLE_STAGE01_DIR",
    ROOT / "03_analysis_results/01_main_data_selection/current",
)
STAGE04 = env_path(
    "MOUSE_GO_STAGE04_DIR",
    ROOT / "03_analysis_results/04_dna_damage_repair_focus/current",
)
V4 = env_path(
    "MOUSE_GO_FOUR_THYMUS_DIR",
    ROOT / "03_analysis_results/06_cross_tissue_integration/01_four_thymus",
)
OUT = ROOT / "03_analysis_results/06_cross_tissue_integration/04_mouse_go_ddr"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

LONG_STATS = env_path(
    "MOUSE_GO_STABLE_LONG_STATS",
    STAGE01 / "03_full_gene_statistics_long.csv.gz",
)
EXACT_DATASET_RESULTS = env_path(
    "MOUSE_GO_EXACT_DATASET_RESULTS",
    STAGE04 / "01_exact_go_dataset_results.csv",
)
THYMUS_AUDIT_SOURCE = V4 / "tables/00_four_thymus_input_contrast_audit.csv"

FINAL_FDR_CUTOFF = 0.10
TARGET_GO_TERMS = {
    "DNA Damage Response (GO:0006974)",
    "DNA Repair (GO:0006281)",
}
GO_SET_ID_TYPE = "gene_symbol"

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

LIVER_TISSUES = {"Liver", "Left lobe of the liver"}
LIVER_ORDER = (
    "OSD-47",
    "OSD-48",
    "OSD-137",
    "OSD-242",
    "OSD-686",
    "OSD-379",
    "OSD-379__age_32w",
    "OSD-463",
)
LIVER_LABELS = {
    "OSD-47": "OSD-47\nLiver | SX-4",
    "OSD-48": "OSD-48\nLiver | SX-4",
    "OSD-137": "OSD-137\nLiver | SX-8",
    "OSD-242": "OSD-242\nLiver | SX-12",
    "OSD-686": "OSD-686\nLiver | SX-12",
    "OSD-379": "OSD-379 10-12w\nLeft lobe | SX-16",
    "OSD-379__age_32w": "OSD-379 32w\nLeft lobe | SX-16",
    "OSD-463": "OSD-463\nLeft lobe | SX-21",
}


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


def save_figure(fig: plt.Figure, stem: str) -> None:
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


def bh(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    if not len(p):
        return np.array([], dtype=float)
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(len(p), dtype=float)
    result[order] = np.clip(adjusted, 0, 1)
    return result


def normalize_ensembl(values: pd.Series) -> pd.Series:
    return (
        values.astype(str)
        .str.strip()
        .str.upper()
        .str.replace(r"\.\d+$", "", regex=True)
    )


def normalize_symbol(values: pd.Series) -> pd.Series:
    return values.fillna("").astype(str).str.strip().str.upper()


def parse_target_go() -> tuple[dict[str, set[str]], pd.DataFrame]:
    go_sets: dict[str, set[str]] = {}
    with GO_GMT.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] in TARGET_GO_TERMS:
                go_sets[fields[0]] = {
                    gene.strip().upper() for gene in fields[2:] if gene.strip()
                }
    if set(go_sets) != TARGET_GO_TERMS:
        raise ValueError("The two prespecified GO terms were not both found")

    damage = go_sets["DNA Damage Response (GO:0006974)"]
    repair = go_sets["DNA Repair (GO:0006281)"]
    membership = pd.DataFrame({"gene_symbol": sorted(damage | repair)})
    membership["in_DNA_damage_response_GO0006974"] = membership[
        "gene_symbol"
    ].isin(damage)
    membership["in_DNA_repair_GO0006281"] = membership["gene_symbol"].isin(
        repair
    )
    membership["GO_membership"] = np.select(
        [
            membership["in_DNA_damage_response_GO0006974"]
            & membership["in_DNA_repair_GO0006281"],
            membership["in_DNA_damage_response_GO0006974"],
        ],
        ["damage+repair", "damage"],
        default="repair",
    )
    return go_sets, membership


def deduplicate_frame(frame: pd.DataFrame, unit_id: str) -> pd.DataFrame:
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


def load_four_thymus() -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    audit = pd.read_csv(THYMUS_AUDIT_SOURCE)
    audit = (
        audit.loc[audit["dataset_id"].isin(THYMUS_ORDER)]
        .set_index("dataset_id")
        .loc[list(THYMUS_ORDER)]
        .reset_index()
    )
    labels = (
        audit["raw_left_label"].astype(str)
        + " "
        + audit["raw_right_label"].astype(str)
    )
    if labels.str.contains(r"\bAG\b", regex=True).any():
        raise AssertionError("An AG comparison entered the four-thymus analysis")

    frames: dict[str, pd.DataFrame] = {}
    for record in audit.to_dict("records"):
        unit_id = str(record["dataset_id"])
        raw = pd.read_csv(ROOT / str(record["input_relative_path"]), low_memory=False)
        log2fc_column = str(record["raw_log2fc_column"])
        stat_column = log2fc_column.replace("Log2fc_", "Stat_", 1)
        fdr_column = log2fc_column.replace("Log2fc_", "Adj.p.value_", 1)
        missing = [
            column
            for column in (
                "ENSEMBL",
                "SYMBOL",
                log2fc_column,
                stat_column,
                fdr_column,
            )
            if column not in raw.columns
        ]
        if missing:
            raise KeyError(f"{unit_id} missing columns: {missing}")
        multiplier = int(record["orientation_multiplier"])
        frame = pd.DataFrame(
            {
                "ensembl_id": normalize_ensembl(raw["ENSEMBL"]),
                "gene_symbol": normalize_symbol(raw["SYMBOL"]),
                "log2fc": pd.to_numeric(raw[log2fc_column], errors="coerce")
                * multiplier,
                "wald_z": pd.to_numeric(raw[stat_column], errors="coerce")
                * multiplier,
                "sample_fdr": pd.to_numeric(raw[fdr_column], errors="coerce"),
            }
        )
        frames[unit_id] = deduplicate_frame(frame, unit_id)
    audit["sample_FDR_filter"] = "none"
    audit["final_selection"] = "two-GO candidate-family final FDR < 0.10"
    return frames, audit


def liver_input_audit() -> pd.DataFrame:
    columns = [
        "analysis_unit_id",
        "accession",
        "mission_cluster",
        "biological_experiment_id",
        "animal_cohort_id",
        "duration_days",
        "raw_left_label",
        "raw_right_label",
        "orientation_multiplier",
        "unified_contrast_label",
    ]
    exact_header = pd.read_csv(EXACT_DATASET_RESULTS, nrows=0)
    tissue_column = (
        "analysis_tissue"
        if "analysis_tissue" in exact_header.columns
        else "tissue"
    )
    exact = pd.read_csv(
        EXACT_DATASET_RESULTS,
        usecols=columns + [tissue_column],
    ).drop_duplicates()
    if tissue_column != "analysis_tissue":
        exact = exact.rename(columns={tissue_column: "analysis_tissue"})
    audit = (
        exact.loc[
            exact["analysis_unit_id"].isin(LIVER_ORDER)
        ]
        .set_index("analysis_unit_id")
        .loc[list(LIVER_ORDER)]
        .reset_index()
    )
    if not audit["unified_contrast_label"].eq("Space Flight vs Ground").all():
        raise ValueError("Liver effects are not uniformly spaceflight minus ground")
    audit["sample_FDR_filter"] = "none"
    audit["final_selection"] = "two-GO candidate-family final FDR < 0.10"
    return audit


def load_integrated_liver() -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    usecols = [
        "ensembl_id",
        "gene_symbol",
        "analysis_unit_id",
        "log2fc_spaceflight_vs_ground",
        "stat_spaceflight_vs_ground",
        "adjusted_p_value",
    ]
    chunks = []
    for chunk in pd.read_csv(
        LONG_STATS,
        usecols=usecols,
        chunksize=200_000,
        low_memory=False,
    ):
        selected = chunk.loc[chunk["analysis_unit_id"].isin(LIVER_ORDER)]
        if not selected.empty:
            chunks.append(selected)
    if not chunks:
        raise ValueError("No integrated-liver gene statistics found")
    long = pd.concat(chunks, ignore_index=True)

    frames = {}
    for unit_id in LIVER_ORDER:
        part = long.loc[long["analysis_unit_id"].eq(unit_id)]
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
        frames[unit_id] = deduplicate_frame(frame, unit_id)
    return frames, liver_input_audit()


def resolve_symbols(values: list[str]) -> tuple[str, bool, str]:
    valid = [value for value in values if value and value != "NAN"]
    if not valid:
        return "", False, ""
    counts = Counter(valid)
    best_count = max(counts.values())
    best = sorted(symbol for symbol, count in counts.items() if count == best_count)[0]
    unique = sorted(counts)
    return best, len(unique) > 1, ";".join(unique)


def build_final_results(
    frames: dict[str, pd.DataFrame],
    unit_order: tuple[str, ...],
    mission_groups: dict[str, tuple[str, ...]],
    go_sets: dict[str, set[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    common = sorted(set.intersection(*(set(frames[unit].index) for unit in unit_order)))
    table = pd.DataFrame({"ensembl_id": common}).set_index("ensembl_id")
    symbol_columns = []
    for unit_id in unit_order:
        frame = frames[unit_id].reindex(common)
        for column in ("gene_symbol", "log2fc", "wald_z", "sample_fdr"):
            table[f"{column}__{unit_id}"] = frame[column]
        symbol_columns.append(f"gene_symbol__{unit_id}")

    resolved = table[symbol_columns].apply(
        lambda row: resolve_symbols(row.fillna("").astype(str).tolist()),
        axis=1,
    )
    table.insert(0, "display_symbol", [value[0] for value in resolved])
    table.insert(1, "symbol_annotation_conflict", [value[1] for value in resolved])
    table.insert(2, "all_observed_symbols", [value[2] for value in resolved])

    damage = go_sets["DNA Damage Response (GO:0006974)"]
    repair = go_sets["DNA Repair (GO:0006281)"]
    if GO_SET_ID_TYPE == "ensembl_gene_id":
        table.insert(
            3,
            "in_DNA_damage_response_GO0006974",
            table.index.isin(damage),
        )
        table.insert(
            4,
            "in_DNA_repair_GO0006281",
            table.index.isin(repair),
        )
    elif GO_SET_ID_TYPE == "gene_symbol":
        observed = table["all_observed_symbols"].map(
            lambda value: set(value.split(";")) if value else set()
        )
        table.insert(
            3,
            "in_DNA_damage_response_GO0006974",
            observed.map(lambda symbols: bool(symbols & damage)),
        )
        table.insert(
            4,
            "in_DNA_repair_GO0006281",
            observed.map(lambda symbols: bool(symbols & repair)),
        )
    else:
        raise ValueError(f"Unsupported GO_SET_ID_TYPE: {GO_SET_ID_TYPE}")
    table.insert(
        5,
        "GO_membership",
        np.select(
            [
                table["in_DNA_damage_response_GO0006974"]
                & table["in_DNA_repair_GO0006281"],
                table["in_DNA_damage_response_GO0006974"],
            ],
            ["damage+repair", "damage"],
            default="repair",
        ),
    )
    candidate = table.loc[
        table["in_DNA_damage_response_GO0006974"]
        | table["in_DNA_repair_GO0006281"]
    ].copy()
    if candidate.empty:
        raise ValueError("No two-GO genes in the complete common background")

    mission_z = pd.DataFrame(index=candidate.index)
    for mission, members in mission_groups.items():
        mission_z[mission] = candidate[
            [f"wald_z__{unit_id}" for unit_id in members]
        ].mean(axis=1)
    candidate["final_integrated_z"] = mission_z.sum(axis=1) / np.sqrt(
        mission_z.shape[1]
    )
    candidate["final_integrated_p"] = 2 * norm.sf(
        np.abs(candidate["final_integrated_z"])
    )
    candidate["final_fdr"] = bh(candidate["final_integrated_p"].to_numpy())
    log2fc = candidate[
        [f"log2fc__{unit_id}" for unit_id in unit_order]
    ]
    candidate["n_positive_effects"] = log2fc.gt(0).sum(axis=1)
    candidate["n_negative_effects"] = log2fc.lt(0).sum(axis=1)
    candidate["effect_pattern"] = np.select(
        [
            candidate["n_positive_effects"].eq(len(unit_order)),
            candidate["n_negative_effects"].eq(len(unit_order)),
        ],
        ["all_positive", "all_negative"],
        default="mixed",
    )
    candidate["sample_FDR_filter_applied"] = False
    candidate["final_test_family"] = (
        "complete common Ensembl IDs annotated to GO:0006974 or GO:0006281"
    )
    candidate = (
        candidate.reset_index()
        .sort_values(
            ["final_integrated_z", "final_fdr", "ensembl_id"],
            ascending=[False, True, True],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )
    selected = candidate.loc[candidate["final_fdr"].lt(FINAL_FDR_CUTOFF)].copy()
    if selected.empty:
        raise ValueError("No two-GO genes passed final FDR < 0.10")
    summary = {
        "common_complete_background": len(common),
        "two_GO_candidate_family": len(candidate),
        "final_FDR01": len(selected),
        "final_FDR01_positive_z": int(selected["final_integrated_z"].gt(0).sum()),
        "final_FDR01_negative_z": int(selected["final_integrated_z"].lt(0).sum()),
        "final_FDR01_mixed_effect_pattern": int(
            selected["effect_pattern"].eq("mixed").sum()
        ),
        "final_FDR01_common_up": int(
            selected["effect_pattern"].eq("all_positive").sum()
        ),
        "final_FDR01_common_down": int(
            selected["effect_pattern"].eq("all_negative").sum()
        ),
        "mission_clusters": len(mission_groups),
    }
    return candidate, selected, summary


def select_top_per_final_direction(
    selected: pd.DataFrame,
    n_per_direction: int = 15,
) -> pd.DataFrame:
    """Select display rows only; the complete final-FDR table remains unchanged."""
    display = selected.copy()
    display["final_direction"] = np.where(
        display["final_integrated_z"].gt(0), "up", "down"
    )
    parts = []
    for direction in ("up", "down"):
        part = (
            display.loc[display["final_direction"].eq(direction)]
            .sort_values(
                ["final_fdr", "final_integrated_p", "ensembl_id"],
                ascending=[True, True, True],
                kind="mergesort",
            )
            .head(n_per_direction)
        )
        parts.append(part)
    result = pd.concat(parts, ignore_index=True)
    if result.empty:
        raise ValueError("No genes available for top-per-direction display")
    result["heatmap_selection_rule"] = (
        f"up to {n_per_direction} genes per final integrated-Z direction, "
        "ranked by final FDR then final P"
    )
    return result


def select_top_common_direction(
    selected: pd.DataFrame,
    n_per_direction: int = 15,
) -> pd.DataFrame:
    """Require every comparison's log2FC to have the same sign."""
    definitions = (
        ("all_positive", "common_up"),
        ("all_negative", "common_down"),
    )
    parts = []
    for effect_pattern, direction in definitions:
        part = (
            selected.loc[selected["effect_pattern"].eq(effect_pattern)]
            .sort_values(
                ["final_fdr", "final_integrated_p", "ensembl_id"],
                ascending=[True, True, True],
                kind="mergesort",
            )
            .head(n_per_direction)
            .copy()
        )
        part["common_direction"] = direction
        parts.append(part)
    result = pd.concat(parts, ignore_index=True)
    result["heatmap_selection_rule"] = (
        f"all comparisons share sign; up to {n_per_direction} per common "
        "direction ranked by final FDR then final P"
    )
    return result


def plot_direct_heatmap(
    selected: pd.DataFrame,
    unit_order: tuple[str, ...],
    labels: dict[str, str],
    title: str,
    stem: str,
) -> None:
    values = selected[
        [f"log2fc__{unit_id}" for unit_id in unit_order]
    ].to_numpy(dtype=float)
    limit = max(1.0, float(np.nanpercentile(np.abs(values), 98)))
    n_rows = len(selected)
    fig, ax = plt.subplots(
        figsize=(
            max(10.5, 1.35 * len(unit_order) + 4.5),
            max(7.0, min(28.0, 3.1 + 0.28 * n_rows)),
        )
    )
    image = ax.imshow(
        np.ma.masked_invalid(values),
        aspect="auto",
        interpolation="none",
        cmap=signed_cmap(),
        vmin=-limit,
        vmax=limit,
    )
    ax.set_xticks(np.arange(len(unit_order)))
    ax.set_xticklabels([labels[unit] for unit in unit_order], rotation=0)
    ax.set_yticks(np.arange(n_rows))
    gene_names = selected["display_symbol"].where(
        selected["display_symbol"].ne(""), selected["ensembl_id"]
    )
    ax.set_yticklabels(
        [
            f"{symbol} | q={q:.2g}"
            for symbol, q in zip(gene_names, selected["final_fdr"])
        ],
        fontsize=6.5 if n_rows > 60 else 7.5,
    )
    ax.set_xticks(np.arange(-0.5, len(unit_order), 1), minor=True)
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
    if n_rows <= 40:
        for row in range(n_rows):
            for column in range(len(unit_order)):
                value = values[row, column]
                if np.isfinite(value):
                    ax.text(
                        column,
                        row,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=6.2,
                        color="white" if abs(value) > 0.65 * limit else "#222222",
                    )
    save_figure(fig, stem)


def plot_no_common_direction_audit(
    candidate: pd.DataFrame,
    stem: str,
) -> None:
    patterns = ["all_positive", "all_negative"]
    labels = ["8/8 common up", "8/8 common down"]
    colors = ["#B64342", "#3775BA"]
    subsets = [candidate.loc[candidate["effect_pattern"].eq(value)] for value in patterns]
    minima = [
        float(subset["final_fdr"].min()) if not subset.empty else np.nan
        for subset in subsets
    ]
    counts = [len(subset) for subset in subsets]

    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    bars = ax.bar(labels, minima, color=colors, width=0.58)
    ax.axhline(
        FINAL_FDR_CUTOFF,
        color="#24313F",
        linewidth=1.6,
        linestyle="--",
        label="Required final FDR = 0.10",
    )
    for bar, minimum, count in zip(bars, minima, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            minimum + 0.008,
            f"minimum q = {minimum:.3f}\n{count} direction-consistent candidates",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    ax.set_ylim(0, max(0.24, max(minima) + 0.06))
    ax.set_ylabel("Minimum final FDR within direction-consistent DDR genes")
    ax.set_xlabel("")
    ax.set_title(
        "Integrated liver | no common-direction DDR gene passes final FDR < 0.10\n"
        "Single-comparison FDR is not restricted; all 8 log2FC values must share sign",
        loc="left",
        weight="bold",
    )
    ax.legend(loc="upper left")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.6)
    save_figure(fig, stem)


def write_summary(
    go_membership: pd.DataFrame,
    thymus_summary: dict[str, int],
    liver_summary: dict[str, int],
) -> None:
    summary = pd.DataFrame(
        [
            {"analysis": "four_thymus", **thymus_summary},
            {"analysis": "integrated_liver", **liver_summary},
        ]
    )
    summary["final_FDR_cutoff"] = FINAL_FDR_CUTOFF
    summary["sample_FDR_filter"] = "none"
    summary["BH_family"] = "two prespecified GO terms after common-background mapping"
    summary.to_csv(
        TABLES / "07_final_FDR01_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    overlap = int(
        (
            go_membership["in_DNA_damage_response_GO0006974"]
            & go_membership["in_DNA_repair_GO0006281"]
        ).sum()
    )
    membership_unit = (
        "Ensembl Gene ID"
        if GO_SET_ID_TYPE == "ensembl_gene_id"
        else "基因符号"
    )
    text = f"""# DDR直接log2FC热图：仅限制最终FDR

## 固定候选范围

- DNA damage response (GO:0006974): {int(go_membership['in_DNA_damage_response_GO0006974'].sum())}个{membership_unit}。
- DNA repair (GO:0006281): {int(go_membership['in_DNA_repair_GO0006281'].sum())}个{membership_unit}。
- 两者并集: {len(go_membership)}个；重叠: {overlap}个。

## 正确的阈值层级

单个比较的FDR不限制。每个组织先在所有比较均有有限log2FC和signed-Z的公共
Ensembl背景中，与两个GO term的并集取交集；随后同一任务内signed-Z先平均，
使每个任务贡献一票。仅对这个预先指定的DDR候选家族的最终双侧P值做BH校正，
以最终FDR<0.10筛选热图行。

不设单比较FDR或log2FC阈值，不要求所有比较方向一致，也不使用top-N。热图
直接显示所有最终FDR<0.10基因在各比较中的统一方向log2FC，行按最终整合Z值
从高到低排列。

完整热图继续保留；另新增汇报版胸腺热图，按最终整合Z的正负方向分别选择
最终FDR最低的Top 15。Top 15仅控制展示行数，不改变完整显著基因结果。

## 共同上下调主筛选

共同上调要求所有比较log2FC均大于0，共同下调要求所有比较log2FC均小于0；
随后还必须满足最终FDR<0.10。四胸腺满足条件的共同上调为
{thymus_summary['final_FDR01_common_up']}个、共同下调为
{thymus_summary['final_FDR01_common_down']}个，汇报图分别展示最终FDR最低的
Top 15。

整合肝脏在该规则下共同上调和共同下调均为0，因此不绘制空热图，改为审计图
展示最接近阈值的结果。{liver_summary['final_FDR01']}基因直接热图继续保留，
但不能称为共同上下调图。

## 四胸腺

- 公共完整背景: {thymus_summary['common_complete_background']:,}个Ensembl ID。
- DDR候选检验家族: {thymus_summary['two_GO_candidate_family']}个。
- 最终FDR<0.10: {thymus_summary['final_FDR01']}个，其中最终Z为正
  {thymus_summary['final_FDR01_positive_z']}个、为负
  {thymus_summary['final_FDR01_negative_z']}个。
- MHU-1 AG未纳入，只保留MHU-1 uG-GC。

## 整合肝脏

- Liver与Left lobe of the liver合为肝脏展示，8个比较仍分别显示。
- 公共完整背景: {liver_summary['common_complete_background']:,}个Ensembl ID。
- DDR候选检验家族: {liver_summary['two_GO_candidate_family']}个。
- 最终FDR<0.10: {liver_summary['final_FDR01']}个，其中最终Z为正
  {liver_summary['final_FDR01_positive_z']}个、为负
  {liver_summary['final_FDR01_negative_z']}个。

## 解释边界

- 红色表示统一方向下正log2FC，蓝色表示负log2FC。
- 单列原始FDR保留在结果表中用于审计，但没有参与筛选。
- 最终FDR<0.10是本轮指定的探索性标准，不能写成FDR<0.05的严格证据。
- 热图保留跨比较方向不一致，不再强制共同方向。
"""
    (OUT / "00_DDR直接热图_最终FDR01说明.md").write_text(text, encoding="utf-8")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style()

    go_sets, go_membership = parse_target_go()
    go_membership.to_csv(
        TABLES / "00_two_GO_union_membership.csv",
        index=False,
        encoding="utf-8-sig",
    )

    thymus_frames, thymus_audit = load_four_thymus()
    thymus_groups = {
        mission: tuple(group["dataset_id"])
        for mission, group in thymus_audit.groupby("mission_cluster", sort=False)
    }
    thymus_all, thymus_selected, thymus_summary = build_final_results(
        thymus_frames, THYMUS_ORDER, thymus_groups, go_sets
    )
    thymus_top15 = select_top_per_final_direction(thymus_selected)
    thymus_common_top15 = select_top_common_direction(thymus_selected)
    thymus_audit.to_csv(
        TABLES / "01_four_thymus_input_contrast_audit.csv",
        index=False,
        encoding="utf-8-sig",
    )
    thymus_all.to_csv(
        TABLES / "02_four_thymus_DDR_final_test_all.csv",
        index=False,
        encoding="utf-8-sig",
    )
    thymus_selected.to_csv(
        TABLES / "03_four_thymus_DDR_final_FDR01_heatmap_genes.csv",
        index=False,
        encoding="utf-8-sig",
    )
    thymus_top15.to_csv(
        TABLES
        / "03b_four_thymus_DDR_final_FDR01_heatmap_top15_per_direction.csv",
        index=False,
        encoding="utf-8-sig",
    )
    thymus_common_top15.to_csv(
        TABLES
        / "03c_four_thymus_DDR_common_direction_final_FDR01_top15.csv",
        index=False,
        encoding="utf-8-sig",
    )
    plot_direct_heatmap(
        thymus_selected,
        THYMUS_ORDER,
        THYMUS_LABELS,
        (
            "Four thymus comparisons | direct DDR log2FC heatmap "
            "(final FDR < 0.10)\n"
            "No per-comparison FDR filter; no direction-consistency or top-N filter"
        ),
        "02_four_thymus_DDR_direct_log2FC_final_FDR01_heatmap",
    )
    plot_direct_heatmap(
        thymus_top15,
        THYMUS_ORDER,
        THYMUS_LABELS,
        (
            "Four thymus comparisons | direct DDR log2FC heatmap "
            "(final FDR < 0.10)\n"
            "Top 15 per final direction by final FDR; "
            "no per-comparison FDR filter"
        ),
        "02b_four_thymus_DDR_direct_log2FC_final_FDR01_top15_per_direction",
    )
    plot_direct_heatmap(
        thymus_common_top15,
        THYMUS_ORDER,
        THYMUS_LABELS,
        (
            "Four thymus comparisons | common DDR up/down "
            "(final FDR < 0.10)\n"
            "All 4 log2FC values share sign; Top 15 per common direction by final FDR"
        ),
        "02c_four_thymus_DDR_common_up_down_final_FDR01_top15_per_direction",
    )

    liver_frames, liver_audit = load_integrated_liver()
    liver_groups = {
        mission: tuple(group["analysis_unit_id"])
        for mission, group in liver_audit.groupby("mission_cluster", sort=False)
    }
    liver_all, liver_selected, liver_summary = build_final_results(
        liver_frames, LIVER_ORDER, liver_groups, go_sets
    )
    liver_common = select_top_common_direction(
        liver_selected,
        n_per_direction=max(1, len(liver_selected)),
    )
    liver_audit.to_csv(
        TABLES / "04_integrated_liver_input_contrast_audit.csv",
        index=False,
        encoding="utf-8-sig",
    )
    liver_all.to_csv(
        TABLES / "05_integrated_liver_DDR_final_test_all.csv",
        index=False,
        encoding="utf-8-sig",
    )
    liver_selected.to_csv(
        TABLES / "06_integrated_liver_DDR_final_FDR01_heatmap_genes.csv",
        index=False,
        encoding="utf-8-sig",
    )
    liver_common.to_csv(
        TABLES
        / "06b_integrated_liver_DDR_common_direction_final_FDR01.csv",
        index=False,
        encoding="utf-8-sig",
    )
    plot_direct_heatmap(
        liver_selected,
        LIVER_ORDER,
        LIVER_LABELS,
        (
            "Integrated liver comparisons | direct DDR log2FC heatmap "
            "(final FDR < 0.10)\n"
            "Liver + left liver lobe; no per-comparison FDR, direction, or top-N filter"
        ),
        "03_integrated_liver_DDR_direct_log2FC_final_FDR01_heatmap",
    )
    plot_no_common_direction_audit(
        liver_all,
        "03b_integrated_liver_DDR_common_direction_final_FDR01_no_genes",
    )

    write_summary(go_membership, thymus_summary, liver_summary)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
