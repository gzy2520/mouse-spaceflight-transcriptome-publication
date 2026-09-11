#!/usr/bin/env python3
"""Build dataset-level GO:0006974 membership tables for Venn/UpSet figures.

The set unit in this deliverable is an OSDR accession (dataset).  If one
accession contains several curated analysis strata (for example strain or
age), all of its included flight/control sample columns are pooled before
the accession-level set is made.  A tissue set is then the union of the
accession sets assigned to that tissue.

The teacher-requested Kidney five-dataset and Thymus four-comparison views
are retained explicitly.  The latter uses the four curated comparison units
(OSD-289 MHU-1, OSD-289 MHU-2, OSD-421 and OSD-515), while all tissue-level
sets remain accession-level.

Analysis key: Ensembl Gene ID.  Symbol is display/audit only.
Presence rule: source value - 1 > 0 in at least one included flight/control
sample.  This is a repertoire-presence rule, not a differential-expression
significance rule.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RUN_DATE = os.environ.get("DDR_VENN_RUN_DATE", "20260826").strip()
if not re.fullmatch(r"\d{8}", RUN_DATE):
    raise ValueError("DDR_VENN_RUN_DATE must be YYYYMMDD")

OUTPUT_DEFAULT = ROOT / "03_analysis_results" / f"20_ddr_go0006974_four_primary_figures_{RUN_DATE}"
OUT = Path(os.environ.get("DDR_VENN_OUTPUT_DIR", str(OUTPUT_DEFAULT)))
if not OUT.is_absolute():
    OUT = ROOT / OUT
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

UNIT_AUDIT = ROOT / "03_analysis_results" / "06_cross_tissue_integration" / "05_proliferation_state" / "tables" / "02_analysis_unit_expression_and_group_audit.csv"
TISSUE_MAPPING = ROOT / "03_analysis_results" / "10_full_stable_id_rerun_20260823" / "04_seven_pathway_response_26_tissues_20260824" / "prepared_inputs" / "00_normalized_original_material_26_tissue_mapping.csv"
TISSUE_NES = ROOT / "03_analysis_results" / "10_full_stable_id_rerun_20260823" / "03_mouse_go_concrete_terms_20260824" / "tables" / "04_tissue_statistics_concrete_terms_and_context.csv"
GO_MEMBERSHIP = ROOT / "08_GO_annotation" / "mgi_mod_gaf_five_relation_concrete_terms_20260824" / "final_manual_curated" / "DNA_damage_response_GO0006974_with_gene_ids.tsv.gz"

EXCLUDED_UP_TISSUES = {"Kidney"}
EXCLUDED_DOWN_TISSUES = {"Lung", "Thymus"}
KIDNEY_ACCESSIONS = ("OSD-163", "OSD-253", "OSD-462", "OSD-513", "OSD-771")
THYMUS_UNITS = ("OSD-289", "OSD-289__MHU2", "OSD-421", "OSD-515")

TISSUE_ZH = {
    "Adrenal gland": "肾上腺", "Bone marrow": "骨髓", "Cecum": "盲肠",
    "Cerebellum": "小脑", "Colon": "结肠", "Dorsal skin": "背部皮肤",
    "Extensor digitorum longus": "趾长伸肌", "Eye": "眼",
    "Femoral lateral skin": "股外侧皮肤", "Femoral skin": "股部皮肤",
    "Gastrocnemius": "腓肠肌", "Heart": "心脏",
    "Heart / Heart right ventricle": "心脏 / 右心室", "Kidney": "肾脏",
    "Left lobe of the liver": "肝左叶", "Liver": "肝脏", "Lung": "肺",
    "Mammary gland": "乳腺", "Optic nerve": "视神经",
    "Quadriceps femoris": "股四头肌", "Retina": "视网膜", "Soleus": "比目鱼肌",
    "Spleen": "脾脏", "Spleen-distal": "脾脏（远端）", "Thymus": "胸腺",
    "Tibialis anterior": "胫骨前肌",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_ensembl(values: pd.Series) -> pd.Series:
    return values.astype("string").str.strip().str.upper().str.replace(r"\.\d+$", "", regex=True)


def split_samples(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    return [x.strip() for x in text.split(";") if x.strip()]


def parse_bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype("string").str.strip().str.lower().isin(["true", "t", "1", "yes"])


def canonical_symbol(values: pd.Series) -> str:
    vals = sorted({str(x).strip() for x in values.dropna() if str(x).strip()})
    canonical = [x for x in vals if "/" not in x]
    return (canonical or vals or [""])[0]


def sanitize(text: str) -> str:
    value = re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_")
    return value or "set"


def load_go_set() -> pd.DataFrame:
    go = pd.read_csv(GO_MEMBERSHIP, sep="\t", compression="infer", low_memory=False)
    required = {"ensembl_gene_id", "SYMBOL", "included_in_gene_set"}
    missing = required - set(go.columns)
    if missing:
        raise ValueError(f"GO membership missing columns: {sorted(missing)}")
    go = go.loc[parse_bool_series(go["included_in_gene_set"])].copy()
    go["ensembl_id"] = norm_ensembl(go["ensembl_gene_id"])
    go = go.loc[go["ensembl_id"].str.fullmatch(r"ENSMUSG\d+").fillna(False)].copy()
    id_map = (
        go.groupby("ensembl_id", sort=True, as_index=False)
        .agg(
            display_symbol=("SYMBOL", canonical_symbol),
            all_source_symbols=("SYMBOL", lambda x: ";".join(sorted({str(v).strip() for v in x.dropna() if str(v).strip()}))),
            mgi_id=("mgi_id", lambda x: ";".join(sorted({str(v).strip() for v in x.dropna() if str(v).strip()}))),
            n_go_annotation_rows=("ensembl_id", "size"),
        )
    )
    if id_map.empty or id_map["ensembl_id"].duplicated().any():
        raise ValueError("GO:0006974 Ensembl map is empty or non-unique")
    if id_map["ensembl_id"].nunique() != 893:
        raise ValueError(f"Expected frozen GO:0006974 set size 893; observed {id_map['ensembl_id'].nunique()}")
    return id_map


def load_metadata() -> tuple[pd.DataFrame, pd.DataFrame]:
    units = pd.read_csv(UNIT_AUDIT, low_memory=False)
    mapping = pd.read_csv(TISSUE_MAPPING, low_memory=False)
    nes = pd.read_csv(TISSUE_NES, low_memory=False)
    required = {"analysis_unit_id", "accession", "input_relative_path", "flight_samples", "control_samples"}
    if required - set(units.columns):
        raise ValueError(f"Expression audit missing: {sorted(required - set(units.columns))}")
    if {"accession", "analysis_tissue"} - set(mapping.columns):
        raise ValueError("26-tissue mapping must contain accession and analysis_tissue")
    if mapping["accession"].duplicated().any():
        raise ValueError("26-tissue mapping must be one row per accession")
    units = units.merge(mapping[["accession", "analysis_tissue", "original_material_type"]], on="accession", how="left", validate="many_to_one")
    if units["analysis_tissue"].isna().any():
        raise ValueError("Some analysis units lack a 26-tissue mapping")
    if units["analysis_unit_id"].duplicated().any():
        raise ValueError("Expression audit has duplicate analysis_unit_id")
    nes = nes.loc[nes["go_id"].eq("GO:0006974")].copy()
    if nes.empty or nes["analysis_tissue"].duplicated().any():
        raise ValueError("GO:0006974 tissue NES is not one-row-per-tissue")
    nes["mean_mission_nes"] = pd.to_numeric(nes["mean_mission_nes"], errors="coerce")
    if nes["mean_mission_nes"].isna().any():
        raise ValueError("GO:0006974 tissue NES contains NA")
    if nes["analysis_tissue"].nunique() != 26:
        raise ValueError(f"Expected 26 tissues; observed {nes['analysis_tissue'].nunique()}")
    return units, nes


def unique_sample_lists(group: pd.DataFrame) -> dict[str, list[str]]:
    flight: list[str] = []
    control: list[str] = []
    for value in group["flight_samples"]:
        flight.extend(split_samples(value))
    for value in group["control_samples"]:
        control.extend(split_samples(value))
    flight = list(dict.fromkeys(flight))
    control = list(dict.fromkeys(control))
    overlap = set(flight) & set(control)
    if overlap:
        raise ValueError(f"Sample appears in both flight and control groups: {sorted(overlap)[:5]}")
    return {"flight": flight, "control": control, "included": list(dict.fromkeys(flight + control))}


def build_dataset_metadata(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for accession, group in units.groupby("accession", sort=True):
        tissues = sorted(group["analysis_tissue"].dropna().unique())
        paths = sorted(group["input_relative_path"].dropna().unique())
        if len(tissues) != 1 or len(paths) != 1:
            raise ValueError(f"Accession {accession} has inconsistent tissue or input path")
        samples = unique_sample_lists(group)
        if not samples["included"]:
            raise ValueError(f"No included samples for accession {accession}")
        rows.append({
            "dataset_id": str(accession),
            "accession": str(accession),
            "analysis_tissue": tissues[0],
            "input_relative_path": paths[0],
            "n_analysis_units": int(len(group)),
            "analysis_unit_ids": ";".join(group["analysis_unit_id"].astype(str)),
            "n_flight_samples": len(samples["flight"]),
            "n_control_samples": len(samples["control"]),
            "flight_samples": ";".join(samples["flight"]),
            "control_samples": ";".join(samples["control"]),
        })
    datasets = pd.DataFrame(rows).sort_values("dataset_id").reset_index(drop=True)
    if len(datasets) != units["accession"].nunique():
        raise AssertionError("Dataset metadata count mismatch")
    return datasets


def source_frame(path: Path, sample_names: list[str], target_ids: set[str], id_order: list[str]) -> tuple[pd.DataFrame, int]:
    header = pd.read_csv(path, nrows=0)
    candidates = {"ENSEMBL", "ENSEMBL_ID", "GENE_ID", "GENEID"}
    id_col = next((c for c in header.columns if str(c).strip().upper() in candidates), None)
    if id_col is None:
        raise ValueError(f"No Ensembl-like column found in {path}")
    missing = sorted(set(sample_names) - set(header.columns))
    if missing:
        raise ValueError(f"Missing sample columns in {path}: {missing[:5]}")
    frame = pd.read_csv(path, usecols=[id_col] + sample_names, dtype={s: "float32" for s in sample_names}, low_memory=False)
    frame = frame.rename(columns={id_col: "ensembl_id"})
    frame["ensembl_id"] = norm_ensembl(frame["ensembl_id"])
    frame = frame.loc[frame["ensembl_id"].isin(target_ids)].copy()
    duplicate_rows = int(frame["ensembl_id"].duplicated().sum())
    if duplicate_rows:
        frame = frame.groupby("ensembl_id", as_index=False, sort=False)[sample_names].max()
    frame = frame.set_index("ensembl_id").reindex(id_order)
    return frame, duplicate_rows


def presence_for_values(values: np.ndarray) -> dict[str, np.ndarray]:
    valid = np.isfinite(values)
    adjusted = values - 1.0
    positive = valid & (adjusted > 0)
    return {
        "valid": valid,
        "positive": positive,
        "any_present": positive.any(axis=1),
    }


def read_presence(datasets: pd.DataFrame, units: pd.DataFrame, id_map: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    target_ids = set(id_map["ensembl_id"])
    id_order = id_map["ensembl_id"].tolist()
    unit_rows: list[pd.DataFrame] = []
    dataset_rows: list[pd.DataFrame] = []
    unit_audits: list[dict[str, object]] = []
    dataset_audits: list[dict[str, object]] = []

    for source_path, path_units in units.groupby("input_relative_path", sort=True):
        path = ROOT / str(source_path)
        if not path.exists():
            raise FileNotFoundError(path)
        dataset_part = datasets.loc[datasets["input_relative_path"].eq(source_path)].copy()
        all_samples: list[str] = []
        for _, row in dataset_part.iterrows():
            all_samples.extend(split_samples(row["flight_samples"]))
            all_samples.extend(split_samples(row["control_samples"]))
        all_samples = list(dict.fromkeys(all_samples))
        frame, duplicate_rows = source_frame(path, all_samples, target_ids, id_order)

        for _, row in dataset_part.iterrows():
            flight = split_samples(row["flight_samples"])
            control = split_samples(row["control_samples"])
            included = list(dict.fromkeys(flight + control))
            values = frame[included].to_numpy(dtype=float)
            info = presence_for_values(values)
            valid, positive = info["valid"], info["positive"]
            source_min = np.full(len(id_order), np.nan)
            source_max = np.full(len(id_order), np.nan)
            has_values = valid.any(axis=1)
            if has_values.any():
                masked = np.where(valid[has_values], values[has_values], np.nan)
                source_min[has_values] = np.nanmin(masked, axis=1)
                source_max[has_values] = np.nanmax(masked, axis=1)
            dataset_rows.append(pd.DataFrame({
                "dataset_id": row["dataset_id"], "accession": row["accession"], "analysis_tissue": row["analysis_tissue"],
                "n_analysis_units": row["n_analysis_units"], "n_included_samples": len(included),
                "n_flight_samples": len(flight), "n_control_samples": len(control),
                "n_valid_values": valid.sum(axis=1).astype(int), "n_positive_after_subtract1": positive.sum(axis=1).astype(int),
                "present_any_included_sample": info["any_present"].astype(int),
                "source_value_min": source_min, "source_value_max": source_max,
                "ensembl_id": id_order,
            }))
            dataset_audits.append({
                "dataset_id": row["dataset_id"], "accession": row["accession"], "analysis_tissue": row["analysis_tissue"],
                "input_relative_path": source_path, "source_sha256": sha256(path), "n_analysis_units": row["n_analysis_units"],
                "n_included_samples": len(included), "n_flight_samples": len(flight), "n_control_samples": len(control),
                "n_target_ensembl_ids": len(id_order), "n_gene_ids_found_in_source": int(frame.notna().any(axis=1).sum()),
                "n_present_any_included_sample": int(info["any_present"].sum()), "n_absent_after_subtract1": int((~info["any_present"]).sum()),
                "n_source_values_exactly_one": int(np.isclose(values, 1.0, rtol=0, atol=1e-7).sum()),
                "n_source_values_le_one": int((valid & ((values - 1.0) <= 0)).sum()),
                "n_nonfinite_source_values": int((~valid).sum()), "n_duplicate_target_rows_collapsed": duplicate_rows,
                "presence_rule": "any included flight/control sample with source value - 1 > 0",
            })

        for _, row in path_units.iterrows():
            flight, control = split_samples(row["flight_samples"]), split_samples(row["control_samples"])
            included = list(dict.fromkeys(flight + control))
            values = frame[included].to_numpy(dtype=float)
            info = presence_for_values(values)
            valid, positive = info["valid"], info["positive"]
            unit_rows.append(pd.DataFrame({
                "analysis_unit_id": row["analysis_unit_id"], "dataset_id": row["accession"], "accession": row["accession"],
                "analysis_tissue": row["analysis_tissue"], "n_included_samples": len(included),
                "n_flight_samples": len(flight), "n_control_samples": len(control),
                "n_valid_values": valid.sum(axis=1).astype(int), "n_positive_after_subtract1": positive.sum(axis=1).astype(int),
                "present_any_included_sample": info["any_present"].astype(int), "ensembl_id": id_order,
            }))
            unit_audits.append({
                "analysis_unit_id": row["analysis_unit_id"], "dataset_id": row["accession"], "accession": row["accession"],
                "analysis_tissue": row["analysis_tissue"], "input_relative_path": source_path, "source_sha256": sha256(path),
                "n_included_samples": len(included), "n_flight_samples": len(flight), "n_control_samples": len(control),
                "n_target_ensembl_ids": len(id_order), "n_present_any_included_sample": int(info["any_present"].sum()),
                "n_absent_after_subtract1": int((~info["any_present"]).sum()), "n_duplicate_target_rows_collapsed": duplicate_rows,
                "presence_rule": "any included flight/control sample with source value - 1 > 0",
            })

    dataset_presence = pd.concat(dataset_rows, ignore_index=True)
    unit_presence = pd.concat(unit_rows, ignore_index=True)
    if len(dataset_presence) != len(datasets) * len(id_order):
        raise AssertionError("Dataset presence row count mismatch")
    if len(unit_presence) != len(units) * len(id_order):
        raise AssertionError("Unit presence row count mismatch")
    return dataset_presence, unit_presence, pd.DataFrame(dataset_audits), pd.DataFrame(unit_audits)


def matrix_outputs(family: str, labels: list[str], members: OrderedDict[str, set[str]], id_map: pd.DataFrame) -> None:
    columns = OrderedDict((label, f"set_{i+1:02d}_{sanitize(label)}") for i, label in enumerate(labels))
    result = id_map[["ensembl_id", "display_symbol", "all_source_symbols", "mgi_id"]].copy()
    for label, col in columns.items():
        result[col] = result["ensembl_id"].isin(members[label]).astype(int)
    result.insert(0, "family", family)
    result.to_csv(TABLES / f"04_membership_matrix_{family}.csv", index=False, encoding="utf-8-sig")
    long_parts = []
    for label, col in columns.items():
        part = result.loc[result[col].eq(1), ["ensembl_id", "display_symbol", "all_source_symbols", "mgi_id"]].copy()
        part.insert(0, "set_name", label)
        part.insert(0, "family", family)
        long_parts.append(part)
    pd.concat(long_parts, ignore_index=True).to_csv(TABLES / f"05_membership_long_{family}.csv", index=False, encoding="utf-8-sig")
    exact_parts = []
    flag_cols = list(columns.values())
    grouped = result.groupby(flag_cols, sort=False, dropna=False)
    for bits, group in grouped:
        if not isinstance(bits, tuple):
            bits = (bits,)
        included = [label for label, bit in zip(labels, bits) if int(bit) == 1]
        exact_parts.append({
            "family": family, "set_membership_pattern": "".join(str(int(x)) for x in bits),
            "exact_region": " & ".join(included) if included else "none", "n_genes": len(group),
            "ensembl_ids": ";".join(group["ensembl_id"].astype(str)),
        })
    exact = pd.DataFrame(exact_parts).sort_values(["n_genes", "set_membership_pattern"], ascending=[False, True])
    exact.to_csv(TABLES / f"06_exact_region_counts_{family}.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame({"set_name": labels, "matrix_column": list(columns.values()), "n_genes": [len(members[label]) for label in labels]}).to_csv(TABLES / f"07_set_sizes_{family}.csv", index=False, encoding="utf-8-sig")


def build_families(dataset_presence: pd.DataFrame, unit_presence: pd.DataFrame, datasets: pd.DataFrame, units: pd.DataFrame, nes: pd.DataFrame) -> tuple[dict[str, OrderedDict[str, set[str]]], pd.DataFrame, pd.DataFrame]:
    dataset_sets = dataset_presence.loc[dataset_presence["present_any_included_sample"].eq(1)].groupby("dataset_id")["ensembl_id"].apply(set).to_dict()
    tissue_sets = {}
    for tissue, group in datasets.groupby("analysis_tissue", sort=True):
        tissue_sets[tissue] = set().union(*(dataset_sets.get(acc, set()) for acc in group["dataset_id"]))
    nes = nes.copy()
    nes["up_selection"] = nes["mean_mission_nes"].gt(1) & ~nes["analysis_tissue"].isin(EXCLUDED_UP_TISSUES)
    nes["down_selection"] = nes["mean_mission_nes"].lt(-1) & ~nes["analysis_tissue"].isin(EXCLUDED_DOWN_TISSUES)
    nes["analysis_tissue_zh"] = nes["analysis_tissue"].map(TISSUE_ZH)
    if nes["analysis_tissue_zh"].isna().any():
        raise ValueError("Missing Chinese tissue label")
    selected_up = nes.loc[nes["up_selection"]].sort_values(["mean_mission_nes", "analysis_tissue"], ascending=[False, True])
    selected_down = nes.loc[nes["down_selection"]].sort_values(["mean_mission_nes", "analysis_tissue"], ascending=[True, True])
    if len(selected_up) != 10 or len(selected_down) != 6:
        raise AssertionError(f"Expected 10 up and 6 down tissues; observed {len(selected_up)} and {len(selected_down)}")
    families: dict[str, OrderedDict[str, set[str]]] = {
        "up_tissues": OrderedDict((r.analysis_tissue, tissue_sets[r.analysis_tissue]) for r in selected_up.itertuples()),
        "down_tissues": OrderedDict((r.analysis_tissue, tissue_sets[r.analysis_tissue]) for r in selected_down.itertuples()),
        "kidney_five_datasets": OrderedDict((acc, dataset_sets[acc]) for acc in KIDNEY_ACCESSIONS),
    }
    unit_sets = unit_presence.loc[unit_presence["present_any_included_sample"].eq(1)].groupby("analysis_unit_id")["ensembl_id"].apply(set).to_dict()
    families["thymus_four_comparisons"] = OrderedDict([
        ("OSD-289 MHU-1", unit_sets["OSD-289"]), ("OSD-289 MHU-2", unit_sets["OSD-289__MHU2"]),
        ("OSD-421", unit_sets["OSD-421"]), ("OSD-515", unit_sets["OSD-515"]),
    ])
    # Every tissue with at least two OSDR accessions gets a dataset-to-dataset Venn.
    for tissue, group in datasets.groupby("analysis_tissue", sort=True):
        accessions = list(group.sort_values("dataset_id")["dataset_id"])
        if len(accessions) >= 2:
            families[f"tissue_{sanitize(tissue)}_datasets"] = OrderedDict((acc, dataset_sets[acc]) for acc in accessions)
    rows = []
    for family, sets in families.items():
        family_type = "tissue_internal_datasets" if family.startswith("tissue_") else family
        tissue_name = family[len("tissue_"):-len("_datasets")].replace("_", " ") if family.startswith("tissue_") else ""
        for label, genes in sets.items():
            rows.append({"family": family, "family_type": family_type, "tissue": tissue_name, "set_name": label, "n_genes": len(genes), "gene_ids": ";".join(sorted(genes))})
    return families, nes, pd.DataFrame(rows)


def write_all_tissue_tables(tissue_sets: dict[str, set[str]], dataset_sets: dict[str, set[str]], datasets: pd.DataFrame, id_map: pd.DataFrame) -> None:
    labels = sorted(tissue_sets)
    tissue_members = OrderedDict((label, tissue_sets[label]) for label in labels)
    matrix_outputs("all_tissues_dataset_union", labels, tissue_members, id_map)
    pd.DataFrame({"analysis_tissue": labels, "analysis_tissue_zh": [TISSUE_ZH[x] for x in labels], "n_datasets": [int((datasets["analysis_tissue"] == x).sum()) for x in labels], "n_genes": [len(tissue_sets[x]) for x in labels], "accessions": [";".join(datasets.loc[datasets["analysis_tissue"].eq(x), "dataset_id"]) for x in labels]}).to_csv(TABLES / "12_tissue_dataset_union_summary.csv", index=False, encoding="utf-8-sig")
    rows = [{"dataset_id": acc, "analysis_tissue": datasets.loc[datasets["dataset_id"].eq(acc), "analysis_tissue"].iloc[0], "n_genes": len(genes), "ensembl_ids": ";".join(sorted(genes))} for acc, genes in sorted(dataset_sets.items())]
    pd.DataFrame(rows).to_csv(TABLES / "13_dataset_set_summary.csv", index=False, encoding="utf-8-sig")


def write_sample_selection(families: dict[str, OrderedDict[str, set[str]]], datasets: pd.DataFrame, units: pd.DataFrame) -> None:
    rows = []
    for family, sets in families.items():
        for label in sets:
            if family in {"up_tissues", "down_tissues"}:
                part = datasets.loc[datasets["analysis_tissue"].eq(label)]
            elif family == "kidney_five_datasets":
                part = datasets.loc[datasets["dataset_id"].eq(label)]
            elif family == "thymus_four_comparisons":
                unit_id = {"OSD-289 MHU-1": "OSD-289", "OSD-289 MHU-2": "OSD-289__MHU2", "OSD-421": "OSD-421", "OSD-515": "OSD-515"}[label]
                part = units.loc[units["analysis_unit_id"].eq(unit_id)]
            else:
                tissue = family[len("tissue_"):-len("_datasets")].replace("_", " ")
                part = datasets.loc[datasets["dataset_id"].eq(label) & datasets["analysis_tissue"].str.replace(r"[^A-Za-z0-9]+", "_", regex=True).str.strip("_").eq(sanitize(tissue))]
                if part.empty:
                    part = datasets.loc[datasets["dataset_id"].eq(label)]
            for _, row in part.iterrows():
                rows.append({"family": family, "set_name": label, "dataset_id": row.get("dataset_id", row.get("accession")), "analysis_unit_id": row.get("analysis_unit_ids", row.get("analysis_unit_id", "")), "accession": row.get("accession", ""), "analysis_tissue": row.get("analysis_tissue", ""), "n_flight_samples": row.get("n_flight_samples", ""), "n_control_samples": row.get("n_control_samples", ""), "flight_samples": row.get("flight_samples", ""), "control_samples": row.get("control_samples", "")})
    pd.DataFrame(rows).to_csv(TABLES / "08_input_dataset_selection.csv", index=False, encoding="utf-8-sig")


def write_provenance(id_map: pd.DataFrame, datasets: pd.DataFrame, units: pd.DataFrame, nes: pd.DataFrame) -> None:
    paths = [GO_MEMBERSHIP, UNIT_AUDIT, TISSUE_MAPPING, TISSUE_NES] + [ROOT / str(x) for x in datasets["input_relative_path"].drop_duplicates()]
    records = [{"path": str(p.relative_to(ROOT)), "exists": p.exists(), "sha256": sha256(p) if p.exists() else "", "size_bytes": p.stat().st_size if p.exists() else np.nan} for p in paths]
    pd.DataFrame(records).drop_duplicates("path").sort_values("path").to_csv(TABLES / "09_input_provenance_sha256.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([
        {"metric": "go_id", "value": "GO:0006974"}, {"metric": "go_set_unique_ensembl_ids", "value": int(id_map["ensembl_id"].nunique())},
        {"metric": "n_osdr_datasets", "value": int(datasets["dataset_id"].nunique())}, {"metric": "n_analysis_units", "value": int(units["analysis_unit_id"].nunique())},
        {"metric": "n_tissues", "value": int(nes["analysis_tissue"].nunique())}, {"metric": "expression_presence_rule", "value": "any included sample with source value - 1 > 0"},
        {"metric": "analysis_key", "value": "Ensembl Gene ID; Symbol display-only"}, {"metric": "random_seed", "value": 25},
    ]).to_csv(TABLES / "00_run_summary.csv", index=False, encoding="utf-8-sig")


def write_readme(families: dict[str, OrderedDict[str, set[str]]], datasets: pd.DataFrame, units: pd.DataFrame, id_map: pd.DataFrame) -> None:
    internal = [k for k in families if k.startswith("tissue_")]
    lines = [
        "# GO:0006974 数据集级 DDR 集合与交集图（2026-08-26）", "",
        "本版本按老师最新口径重建。分析键为 Ensembl Gene ID，Symbol 仅作为显示和审计标签。",
        "", "## 集合定义", "",
        f"- GO 基因集：GO:0006974（DNA damage response），五关系后代闭包的冻结结果，共 {id_map['ensembl_id'].nunique()} 个唯一小鼠 Ensembl ID。",
        "- OSDR 数据集集合：同一 accession 内的年龄、品系等分析分层先合并；该 accession 的全部纳入 flight/control 样本合并后计算集合。",
        "- 表达存在规则：对 GeneLab 源值先减去 1，只有调整后值 > 0 才算存在；任一纳入样本满足即纳入该数据集集合。值等于 1 或缺失不会计入。",
        "- 组织集合：该组织所属 accession 集合的并集。NES>1（排除 Kidney）和 NES<-1（排除 Lung、Thymus）只用于选择组织分组，不改变表达集合定义。",
        f"- 当前数据规模：{datasets['dataset_id'].nunique()} 个 OSDR accession、{units['analysis_unit_id'].nunique()} 个分析单元、26 个组织。",
        "", "## 图形解释", "",
        "- 2–7 个集合使用普通非面积拟合 Venn 图；每个区域的数字来自精确 Ensembl membership pattern。圆/多边形位置是示意，不代表区域面积。",
        "- 10 个上调组织无法用可读的普通 Venn 展示，因此使用精确 UpSet 交集矩阵；其柱高和点阵直接对应精确交集计数，不是面积拟合。",
        "- 胸腺四比较按老师指定保留 OSD-289 MHU-1、OSD-289 MHU-2、OSD-421、OSD-515；这是比较层视图，胸腺组织并集仍按 accession 级别计算。",
        "", "## 输出", "",
        "- `figures/`：本次老师要求的四张主图——上调组织、下调组织、Kidney 五数据集、Thymus 四比较。组织内部 dataset 家族仅保留在 tables 中供审计，不作为本次主图。",
        "- `tables/04–07_*`：membership 矩阵、长表、精确区域计数和集合大小；`tables/06_exact_region_counts_*.csv` 是区域数字的权威来源。",
        "- `tables/02_dataset_presence_audit.csv`：每个 accession 的样本数、缺失/等于 1/减 1 后为正的审计。",
        "- `tables/12_tissue_dataset_union_summary.csv`：26 个组织的 accession 组成和并集大小。",
        "",
        "本版本保留旧版输出，不覆盖旧结果；旧版将 analysis unit 直接作为集合，不能代替本数据集级版本。",
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    id_map = load_go_set()
    units, nes = load_metadata()
    datasets = build_dataset_metadata(units)
    dataset_presence, unit_presence, dataset_audit, unit_audit = read_presence(datasets, units, id_map)
    dataset_presence = dataset_presence.merge(id_map, on="ensembl_id", how="left", validate="many_to_one")
    unit_presence = unit_presence.merge(id_map, on="ensembl_id", how="left", validate="many_to_one")
    dataset_presence.to_csv(TABLES / "03_dataset_gene_presence_after_subtract1.csv.gz", index=False, compression="gzip")
    unit_presence.to_csv(TABLES / "04_analysis_unit_gene_presence_after_subtract1.csv.gz", index=False, compression="gzip")
    dataset_audit.to_csv(TABLES / "02_dataset_presence_audit.csv", index=False, encoding="utf-8-sig")
    unit_audit.to_csv(TABLES / "02b_analysis_unit_presence_audit.csv", index=False, encoding="utf-8-sig")

    dataset_sets = dataset_presence.loc[dataset_presence["present_any_included_sample"].eq(1)].groupby("dataset_id")["ensembl_id"].apply(set).to_dict()
    tissue_sets = {t: set().union(*(dataset_sets[a] for a in g["dataset_id"])) for t, g in datasets.groupby("analysis_tissue", sort=True)}
    write_all_tissue_tables(tissue_sets, dataset_sets, datasets, id_map)
    families, nes, family_summary = build_families(dataset_presence, unit_presence, datasets, units, nes)
    nes.sort_values("analysis_tissue").to_csv(TABLES / "01_go0006974_tissue_nes_selection.csv", index=False, encoding="utf-8-sig")
    family_summary.to_csv(TABLES / "10_family_set_summary.csv", index=False, encoding="utf-8-sig")
    for family, sets in families.items():
        matrix_outputs(family, list(sets), sets, id_map)
    write_sample_selection(families, datasets, units)
    write_provenance(id_map, datasets, units, nes)
    write_readme(families, datasets, units, id_map)
    manifest = {
        "run_date": RUN_DATE, "output": str(OUT), "go_id": "GO:0006974", "go_set_unique_ensembl": int(id_map["ensembl_id"].nunique()),
        "n_osdr_datasets": int(datasets["dataset_id"].nunique()), "n_analysis_units": int(units["analysis_unit_id"].nunique()), "n_tissues": int(nes["analysis_tissue"].nunique()),
        "presence_rule": "any included flight/control sample with (source value - 1) > 0", "dataset_rule": "all strata within accession pooled before dataset set construction",
        "tissue_rule": "union of accession sets assigned to tissue", "analysis_key": "Ensembl Gene ID", "symbol_role": "display and audit only",
        "venn_geometry": "schematic non-area-proportional geometry; exact region counts in CSV", "families": {k: {"n_sets": len(v), "set_names": list(v)} for k, v in families.items()},
    }
    (OUT / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: {label: len(genes) for label, genes in v.items()} for k, v in families.items()}, ensure_ascii=False, indent=2))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
