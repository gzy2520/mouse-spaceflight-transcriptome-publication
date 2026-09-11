#!/usr/bin/env python3
"""Sample-level log2FC pathway relationships: per-tissue Spearman, Fisher-z overall.

Teacher-requested method (2026-09-02) for the publication Fig_2 pipeline:
1. Input is per-sample log2FC (spaceflight sample vs same-tissue, same-dataset
   ground mean), NOT mission-equal NES.
2. First compute the 15 x 15 Spearman relationship matrix WITHIN each tissue,
   with the flight sample (mouse) as the observation unit.
3. Then combine the per-tissue correlations across tissues into one overall
   15 x 15 matrix via equal-weight Fisher-z averaging.

Feasibility note (audited 2026-09-02): 21 of 26 tissues have fewer than 3 OSD
datasets, so the dataset/mission grain cannot support within-tissue
correlations. The sample grain is the only feasible within-tissue unit
(flight samples per tissue: 2-55). Tissues with fewer than 5 flight samples
(Colon 4, Heart 3, Femoral skin 2) remain in the per-tissue tables and
figures but are excluded from the overall Fisher-z combination.

All analysis keys are Ensembl Gene IDs; symbols are display-only.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

# The selected source matrices are normalized expression values with a
# positive (+1) source offset already present.  The established project
# sample-wise effect definition therefore uses log2 values directly and does
# not add another pseudocount.
PSEUDOCOUNT = 0.0
MIN_FLIGHT_SAMPLES_FOR_OVERALL = 5
Z_CAP = 5.0
EXPECTED_N_TERMS = 15
EXPECTED_N_TISSUES = 26
EXPECTED_N_FLIGHT_SAMPLES = 360


def normalize_ensembl(series: pd.Series) -> pd.Series:
    """Uppercase and strip trailing version numbers, matching stage01."""
    cleaned = series.astype(str).str.strip().str.upper()
    return cleaned.str.replace(r"\.[0-9]+$", "", regex=True)


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compute_sample_term_log2fc(
    expression: np.ndarray,
    term_member_masks: dict[str, np.ndarray],
    ground_index: np.ndarray,
    flight_index: np.ndarray,
    pseudocount: float = PSEUDOCOUNT,
) -> pd.DataFrame:
    """Per flight sample, per term: mean member-gene log2FC vs ground log2 mean.

    expression: genes x samples matrix (same normalized scale per dataset).
    term_member_masks: term_key -> boolean gene membership (same gene order).
    Returns: DataFrame with flight samples as rows and term keys as columns.
    """
    if ground_index.size == 0 or flight_index.size == 0:
        raise ValueError("Need at least one ground and one flight sample")
    if (expression <= 0).any():
        raise ValueError("Sample-wise log2FC requires strictly positive source expression values")
    ground_log_ref = np.log2(expression[:, ground_index] + pseudocount).mean(axis=1)
    sample_log2fc = np.log2(expression[:, flight_index] + pseudocount) - ground_log_ref[:, None]
    data = {}
    for term_key, mask in term_member_masks.items():
        if not mask.any():
            raise ValueError(f"Term {term_key} has no member genes")
        data[term_key] = sample_log2fc[mask, :].mean(axis=0)
    result = pd.DataFrame(data)
    return result


def tissue_spearman_matrix(sample_term: pd.DataFrame) -> np.ndarray:
    """15 x 15 Spearman matrix of term values across samples (rows)."""
    return sample_term.corr(method="spearman").to_numpy(dtype=float)


def fisher_z_aggregate(
    per_tissue_rho: dict[str, np.ndarray],
    term_keys: list[str],
    included_tissues: list[str],
    z_cap: float = Z_CAP,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict]]:
    """Equal-weight Fisher-z combination of per-tissue Spearman matrices.

    Returns (overall_rho long-freeze frame, per_tissue_z long frame, audit).
    |z| above z_cap (including rho = +/-1) is capped and flagged.
    """
    z_long_rows = []
    audit = {}
    overall = pd.DataFrame(np.nan, index=term_keys, columns=term_keys)
    for i, key_i in enumerate(term_keys):
        for j, key_j in enumerate(term_keys):
            if j <= i:
                continue
            z_values = []
            for tissue in included_tissues:
                rho = per_tissue_rho[tissue][i, j]
                if not np.isfinite(rho):
                    continue
                z_raw = np.arctanh(np.clip(rho, -1.0 + 1e-12, 1.0 - 1e-12))
                z = float(np.clip(z_raw, -z_cap, z_cap))
                capped = abs(z_raw) > z_cap
                z_long_rows.append(
                    {
                        "analysis_tissue": tissue,
                        "pathway_1_key": key_i,
                        "pathway_2_key": key_j,
                        "fisher_z": z,
                        "z_capped": bool(capped),
                    }
                )
                z_values.append(z)
            overall.loc[key_i, key_j] = np.tanh(float(np.mean(z_values))) if z_values else np.nan
    z_long = pd.DataFrame(z_long_rows)
    for tissue in included_tissues:
        rho_matrix = per_tissue_rho[tissue]
        pair_slice = z_long[
            (z_long["analysis_tissue"] == tissue)
        ]
        audit[tissue] = {
            "n_pairs_na": int(np.sum(~np.isfinite(rho_matrix[np.triu_indices_from(rho_matrix, 1)]))),
            "n_pairs_capped": int(pair_slice["z_capped"].sum()) if len(pair_slice) else 0,
        }
    overall = overall.fillna(overall.T)
    return overall, z_long, audit


def average_linkage_order(rho_matrix: np.ndarray, term_keys: list[str]) -> list[str]:
    """Average-linkage (UPGMA) order on distance = (1 - rho) / 2, as in R hclust."""
    rho = np.where(np.isfinite(rho_matrix), rho_matrix, 0.0)
    rho = (rho + rho.T) / 2.0
    np.fill_diagonal(rho, 1.0)
    distance = (1.0 - rho) / 2.0
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=False)
    merged = linkage(condensed, method="average")
    order = [term_keys[idx] for idx in leaves_list(merged)]
    return order


def load_inputs(root: str) -> dict:
    publication_root = os.path.join(root, "publication_release")
    metadata_path = os.path.join(
        publication_root, "data/publication_input/mouse_metadata/03_mouse_level_metadata.csv"
    )
    term_stats_path = os.path.join(
        publication_root, "data/publication_input/go/04_tissue_statistics_concrete_terms_and_context.csv"
    )
    membership_path = os.path.join(
        root,
        "03_analysis_results/10_full_stable_id_rerun_20260823/"
        "03_mouse_go_concrete_terms_20260824/tables/01_concrete_GO_membership.csv",
    )
    legacy_rho_path = os.path.join(
        publication_root, "data/publication_input/go/13_tissue_pathway_Spearman_rho_matrix.csv"
    )
    for path in (metadata_path, term_stats_path, membership_path, legacy_rho_path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required input missing: {path}")

    metadata = pd.read_csv(metadata_path)
    term_stats = pd.read_csv(term_stats_path)
    membership = pd.read_csv(membership_path)
    legacy_rho = pd.read_csv(legacy_rho_path)

    term_meta = term_stats[["term_key", "term_name", "go_id", "term_order"]].drop_duplicates("term_key")
    term_meta = term_meta.sort_values("term_order").reset_index(drop=True)
    if term_meta["term_key"].nunique() != EXPECTED_N_TERMS or len(term_meta) != EXPECTED_N_TERMS:
        raise ValueError(f"Expected {EXPECTED_N_TERMS} GO terms, found {len(term_meta)}")
    term_keys = term_meta["term_key"].tolist()

    membership_ids = normalize_ensembl(membership["ensembl_gene_id"])
    if membership_ids.duplicated().any():
        raise ValueError("GO membership contains duplicate Ensembl IDs")
    member_masks = {}
    for term_key in term_keys:
        column = f"in_{term_key}"
        if column not in membership.columns:
            raise ValueError(f"Membership lacks column {column}")
        member_masks[term_key] = membership[column].astype(bool).to_numpy()
    gene_order = membership_ids.tolist()

    if metadata["Group"].nunique() != EXPECTED_N_TISSUES:
        raise ValueError(f"Expected {EXPECTED_N_TISSUES} tissues, found {metadata['Group'].nunique()}")
    flight = metadata[metadata["sample_status"] == "Flight"]
    if flight["sample_column"].nunique() != EXPECTED_N_FLIGHT_SAMPLES:
        raise ValueError(
            f"Expected {EXPECTED_N_FLIGHT_SAMPLES} flight samples, found {flight['sample_column'].nunique()}"
        )
    if metadata["sample_column"].duplicated().any():
        raise ValueError("Metadata contains duplicate sample columns")

    return {
        "root": root,
        "metadata": metadata,
        "term_meta": term_meta,
        "term_keys": term_keys,
        "gene_order": gene_order,
        "member_masks": member_masks,
        "membership_path": membership_path,
        "metadata_path": metadata_path,
        "term_stats_path": term_stats_path,
        "legacy_rho": legacy_rho,
        "legacy_rho_path": legacy_rho_path,
    }


def build_sample_term_long(inputs: dict) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Read all 48 raw matrices and build the flight sample x term log2FC table."""
    metadata = inputs["metadata"]
    term_keys = inputs["term_keys"]
    gene_order = inputs["gene_order"]
    member_masks = inputs["member_masks"]
    scope_rows = []
    term_frames = []
    for accession, chunk in metadata.groupby("accession", sort=True):
        path = os.path.join(inputs["root"], chunk["input_relative_path"].iloc[0])
        if not os.path.exists(path):
            raise FileNotFoundError(f"Raw matrix missing: {path}")
        # `sample_column` is the exact column name in the frozen raw matrix;
        # `sample_name_normalized` intentionally drops assay suffixes and is
        # therefore not safe for file-column lookup (e.g. the OSD-462 mRNA
        # columns).
        expected_samples = sorted(chunk["sample_column"].unique())
        use_cols = ["ENSEMBL"] + expected_samples
        raw = pd.read_csv(path, usecols=lambda c: c in set(use_cols))
        raw["ENSEMBL"] = normalize_ensembl(raw["ENSEMBL"])
        missing_samples = [s for s in expected_samples if s not in raw.columns]
        if missing_samples:
            raise ValueError(f"{accession}: missing sample columns {missing_samples[:5]}")
        raw = raw.drop_duplicates("ENSEMBL").set_index("ENSEMBL")
        raw_gene_set = set(raw.index)
        available_indices = np.array(
            [idx for idx, gene in enumerate(gene_order) if gene in raw_gene_set], dtype=int
        )
        if available_indices.size == 0:
            raise ValueError(f"{accession}: no GO-membership Ensembl IDs are present")
        available_genes = [gene_order[idx] for idx in available_indices]
        raw = raw.loc[available_genes, expected_samples]
        expression = raw.to_numpy(dtype=float)
        if not np.isfinite(expression).all():
            raise ValueError(f"{accession}: non-finite expression values")
        local_member_masks = {
            term_key: mask[available_indices]
            for term_key, mask in member_masks.items()
        }

        for tissue, tissue_chunk in chunk.groupby("Group", sort=True):
            ground_samples = tissue_chunk.loc[
                tissue_chunk["sample_status"].eq("Ground"), "sample_column"
            ].tolist()
            flight_samples = tissue_chunk.loc[
                tissue_chunk["sample_status"].eq("Flight"), "sample_column"
            ].tolist()
            sample_position = {sample: idx for idx, sample in enumerate(expected_samples)}
            ground_idx = [sample_position[s] for s in ground_samples]
            flight_idx = [sample_position[s] for s in flight_samples]
            if not ground_idx or not flight_idx:
                scope_rows.append(
                    {
                        "accession": accession,
                        "analysis_tissue": tissue,
                        "n_ground_samples": len(ground_idx),
                        "n_flight_samples": len(flight_idx),
                        "usable": bool(ground_idx and flight_idx),
                        "note": "" if (ground_idx and flight_idx) else "missing ground or flight samples",
                    }
                )
                continue
            frame = compute_sample_term_log2fc(
                expression, local_member_masks, np.array(ground_idx), np.array(flight_idx)
            )
            frame.index = [expected_samples[idx] for idx in flight_idx]
            frame = frame.reset_index()
            frame.columns = ["sample_column"] + term_keys
            frame.insert(0, "analysis_tissue", tissue)
            frame.insert(0, "accession", accession)
            term_frames.append(frame)
            scope_rows.append(
                {
                    "accession": accession,
                    "analysis_tissue": tissue,
                    "n_ground_samples": len(ground_idx),
                    "n_flight_samples": len(flight_idx),
                    "usable": True,
                    "note": "",
                }
            )
    scope = pd.DataFrame(scope_rows)
    sample_term = pd.concat(term_frames, ignore_index=True)
    sample_term = sample_term.melt(
        id_vars=["analysis_tissue", "accession", "sample_column"],
        value_vars=term_keys,
        var_name="term_key",
        value_name="term_log2fc",
    )
    if not np.isfinite(sample_term["term_log2fc"]).all():
        raise ValueError("Sample-level term log2FC contains non-finite values")
    return sample_term, scope, term_keys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Main project root")
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory (default: 03_analysis_results/06_cross_tissue_integration/10_log2fc_per_tissue_pathway_spearman_20260902)",
    )
    args = parser.parse_args(argv)
    root = os.path.abspath(args.root)
    out_dir = args.out_dir or os.path.join(
        root,
        "03_analysis_results/06_cross_tissue_integration/"
        "10_log2fc_per_tissue_pathway_spearman_20260902",
    )
    os.makedirs(out_dir, exist_ok=True)

    inputs = load_inputs(root)
    term_keys = inputs["term_keys"]
    term_meta = inputs["term_meta"]
    metadata = inputs["metadata"]

    sample_term, scope, _ = build_sample_term_long(inputs)
    flight = metadata[metadata["sample_status"] == "Flight"]
    ground = metadata[metadata["sample_status"] == "Ground"]
    if sample_term["sample_column"].nunique() != EXPECTED_N_FLIGHT_SAMPLES:
        raise ValueError(
            "Usable sample-level log2FC rows do not cover the expected "
            f"{EXPECTED_N_FLIGHT_SAMPLES} unique flight samples"
        )
    tissues = sorted(metadata["Group"].unique())
    per_tissue_rho: dict[str, np.ndarray] = {}
    pair_rows = []
    n_by_tissue: dict[str, int] = {}
    for tissue in tissues:
        chunk = sample_term[sample_term["analysis_tissue"] == tissue]
        wide = chunk.pivot_table(index="sample_column", columns="term_key", values="term_log2fc", aggfunc="mean")
        wide = wide[term_keys]
        n_samples = int(len(wide))
        n_by_tissue[tissue] = n_samples
        rho = tissue_spearman_matrix(wide) if n_samples >= 2 else np.full((len(term_keys), len(term_keys)), np.nan)
        per_tissue_rho[tissue] = rho
        for i, key_i in enumerate(term_keys):
            for j in range(i + 1, len(term_keys)):
                pair_rows.append(
                    {
                        "analysis_tissue": tissue,
                        "n_flight_samples": n_samples,
                        "pathway_1_key": key_i,
                        "pathway_2_key": term_keys[j],
                        "spearman_rho": float(rho[i, j]) if np.isfinite(rho[i, j]) else np.nan,
                    }
                )
    pairs = pd.DataFrame(pair_rows)

    # Count only usable flight rows that actually contributed all 15 terms.
    # This prevents metadata-only samples with no matched ground reference
    # from entering the Fisher-z inclusion decision.
    included = [t for t in tissues if n_by_tissue.get(t, 0) >= MIN_FLIGHT_SAMPLES_FOR_OVERALL]
    excluded = [t for t in tissues if t not in included]
    overall, z_long, audit = fisher_z_aggregate(per_tissue_rho, term_keys, included)
    for index in range(len(overall)):
        overall.iat[index, index] = 1.0  # pandas 3 exposes read-only NumPy views
    if overall.isna().any().any():
        raise ValueError("Overall Fisher-z matrix contains undefined pairs")
    overall = (overall + overall.T) / 2.0

    order = average_linkage_order(overall.to_numpy(), term_keys)
    order_frame = pd.DataFrame(
        {"cluster_position": range(1, len(order) + 1), "term_key": order}
    )

    inclusion_rows = []
    for tissue in tissues:
        is_in = tissue in included
        tissue_scope = scope[scope["analysis_tissue"] == tissue]
        n_ground = int(tissue_scope.loc[tissue_scope["usable"], "n_ground_samples"].sum())
        rho = per_tissue_rho[tissue]
        n_pairs_na = int(np.sum(~np.isfinite(rho[np.triu_indices_from(rho, 1)])))
        inclusion_rows.append(
            {
                "analysis_tissue": tissue,
                "n_flight_samples": int(n_by_tissue.get(tissue, 0)),
                "n_ground_samples": n_ground,
                "included_in_overall": bool(is_in),
                "exclusion_reason": "" if is_in else f"n_flight < {MIN_FLIGHT_SAMPLES_FOR_OVERALL}",
                "n_pairs_na": audit.get(tissue, {"n_pairs_na": n_pairs_na})["n_pairs_na"],
                "n_pairs_z_capped": audit.get(tissue, {"n_pairs_capped": 0})["n_pairs_capped"],
            }
        )
    inclusion = pd.DataFrame(inclusion_rows)

    legacy = inputs["legacy_rho"].set_index("term_key")
    legacy = legacy[term_keys].astype(float)
    legacy = (legacy + legacy.T) / 2.0
    comparison_rows = []
    for i, key_i in enumerate(term_keys):
        for j in range(i + 1, len(term_keys)):
            comparison_rows.append(
                {
                    "pathway_1_key": key_i,
                    "pathway_2_key": term_keys[j],
                    "legacy_nes_rho": float(legacy.loc[key_i, term_keys[j]]),
                    "new_log2fc_fisher_z_rho": float(overall.loc[key_i, term_keys[j]]),
                    "delta": float(overall.loc[key_i, term_keys[j]] - legacy.loc[key_i, term_keys[j]]),
                }
            )
    comparison = pd.DataFrame(comparison_rows)

    pair_audit = (
        z_long.groupby(["pathway_1_key", "pathway_2_key"], sort=False)
        .agg(
            n_tissues_used=("analysis_tissue", "nunique"),
            fisher_z_mean=("fisher_z", "mean"),
            n_tissues_z_capped=("z_capped", "sum"),
        )
        .reset_index()
    )
    pair_audit["overall_rho"] = [
        float(overall.loc[row.pathway_1_key, row.pathway_2_key])
        for row in pair_audit.itertuples(index=False)
    ]
    if len(pair_audit) != len(term_keys) * (len(term_keys) - 1) // 2:
        raise ValueError("Overall pair audit does not contain all upper-triangle pairs")

    overall_freeze = overall.reset_index().rename(columns={"index": "term_key"})

    provenance_rows = []
    provenance_sources = {
        "mouse_metadata": inputs["metadata_path"],
        "term_statistics": inputs["term_stats_path"],
        "go_membership": inputs["membership_path"],
        "legacy_nes_rho_matrix": inputs["legacy_rho_path"],
    }
    for name, path in provenance_sources.items():
        provenance_rows.append({"input_role": name, "relative_path": os.path.relpath(path, root), "sha256": sha256_file(path)})
    for accession, chunk in metadata.groupby("accession", sort=True):
        raw_path = os.path.join(root, chunk["input_relative_path"].iloc[0])
        provenance_rows.append(
            {
                "input_role": "raw_expression_matrix",
                "relative_path": chunk["input_relative_path"].iloc[0],
                "sha256": sha256_file(raw_path),
            }
        )
    provenance = pd.DataFrame(provenance_rows)

    scope.to_csv(os.path.join(out_dir, "01_sample_scope_audit.csv"), index=False)
    sample_term.to_csv(os.path.join(out_dir, "02_sample_term_log2fc_matrix.csv.gz"), index=False, compression="gzip")
    pairs.to_csv(os.path.join(out_dir, "03_tissue_pairwise_spearman_log2fc_long.csv"), index=False)
    overall_freeze.to_csv(os.path.join(out_dir, "04_overall_fisher_z_spearman_rho_matrix.csv"), index=False)
    z_long.to_csv(os.path.join(out_dir, "05_fisher_z_per_tissue_long.csv"), index=False)
    inclusion.to_csv(os.path.join(out_dir, "06_fisher_z_inclusion_audit.csv"), index=False)
    order_frame.to_csv(os.path.join(out_dir, "07_cluster_order_overall_log2fc.csv"), index=False)
    comparison.to_csv(os.path.join(out_dir, "08_comparison_vs_legacy_nes_rho.csv"), index=False)
    provenance.to_csv(os.path.join(out_dir, "09_input_provenance_sha256.csv"), index=False)
    term_meta.to_csv(os.path.join(out_dir, "10_term_metadata.csv"), index=False)
    pair_audit.to_csv(os.path.join(out_dir, "11_overall_pairwise_audit.csv"), index=False)

    n_capped_total = int(z_long["z_capped"].sum())
    readme = f"""# 样本级 log2FC 通路关系：组织内 Spearman + Fisher-z 跨组织汇总（2026-09-02）

老师指令（2026-09-02）：用 log2FC 作为输入；先计算每个组织内部 15 条 GO term
两两之间的 Spearman 关系（逐组织出图）；再根据每个组织的相关性汇总出整体
跨组织的 15×15 矩阵。本目录是主项目的日期化计算产物；审阅后可将其复制为
`publication_release/data/publication_input/go/` 下的带日期冻结副本，供
`final_result/Main/Fig_2` 与逐组织图渲染使用。

## 计算对象与方法

- 输入：48 个 OSD 数据集的原始表达矩阵（`02_processed_data/osdr_processed/`），
  样本-组织-状态映射用冻结的 761 样本元数据
  （`publication_release/data/publication_input/mouse_metadata/03_mouse_level_metadata.csv`，
  Flight {int(flight['sample_column'].nunique())} / Ground {int(ground['sample_column'].nunique())}），
  15 条 GO term 的 Ensembl 基因集用冻结的 concrete-term membership。
- 样本级 log2FC：每个 flight 样本对「同数据集、同组织」的 ground 样本 log2 均值
  计算 log2(x_sample) − mean[log2(x_ground)]。源矩阵已带正值偏移，
  本次不再额外添加 pseudocount（pseudocount = {PSEUDOCOUNT}）；ground 参照按数据集内取，
  避免跨批次参照。
- term 值：每个样本 × term = 该 term 成员基因 log2FC 的算术平均
  （分析键始终为 Ensembl Gene ID，Symbol 只用于显示）。
- 组织内 Spearman：每个组织内，15 条 term 的两两 Spearman rho，
  观测单位 = flight 样本（每组织 n = {min(n_by_tissue.values())}–{max(n_by_tissue.values())}）。
  数据集/mission 粒度不可行：26 个组织中 21 个不足 3 个 OSD 数据集。
- 跨组织汇总：对每组织 rho 做 Fisher z = atanh(rho)，{len(included)} 个
  n_flight ≥ {MIN_FLIGHT_SAMPLES_FOR_OVERALL} 的组织等权平均后 tanh 反变换。
  排除 {', '.join(f'{t} (n={int(n_by_tissue[t])})' for t in excluded)}；
  它们保留在逐组织表与图中。|z| > {Z_CAP:g}（含 rho = ±1）截断为 ±{Z_CAP:g}，
  本次共截断 {n_capped_total} 个组织×对。
- 展示顺序：对整体矩阵做 distance = (1 - rho)/2 的 average linkage（UPGMA）聚类。

## 输出

- `01_sample_scope_audit.csv`：每数据集×组织的样本数与可用性审计。
- `02_sample_term_log2fc_matrix.csv.gz`：flight 样本 × term 的 log2FC 长表
  （{len(sample_term)} 行 = {int(flight['sample_column'].nunique())} 样本 × {EXPECTED_N_TERMS} term）。
- `03_tissue_pairwise_spearman_log2fc_long.csv`：{len(tissues)} 组织 × 105 对的
  组织内 Spearman rho。
- `04_overall_fisher_z_spearman_rho_matrix.csv`：整体 15×15 Fisher-z 汇总矩阵。
- `05_fisher_z_per_tissue_long.csv`：参与汇总的组织×对的 z 值（含截断标记）。
- `06_fisher_z_inclusion_audit.csv`：组织纳入/排除审计。
- `07_cluster_order_overall_log2fc.csv`：整体矩阵的 average-linkage 展示顺序。
- `08_comparison_vs_legacy_nes_rho.csv`：与旧 mission-equal-NES 整体矩阵的逐对对比。
- `09_input_provenance_sha256.csv`：全部输入的 SHA-256。
- `10_term_metadata.csv`：15 条 GO term 的稳定 key、GO ID 和显示名称。
- `11_overall_pairwise_audit.csv`：每个整体 pair 的纳入组织数、Fisher-z 均值、
  截断数和最终 rho。

## 可解释边界

- 同一组织的样本可能来自不同 OSD 数据集（批次），批次结构会保留在样本值中，
  可能放大组织内相关；rho 描述共响应，不是因果或直接调控。
- GO term 之间存在基因重叠，相关也可能来自共享基因。
- 小样本组织（n < {MIN_FLIGHT_SAMPLES_FOR_OVERALL}）的 rho 不稳定，只展示不汇总。
"""
    with open(os.path.join(out_dir, "README.md"), "w", encoding="utf-8") as handle:
        handle.write(readme)

    print(
        f"LOG2FC_PER_TISSUE_SPEARMAN_PASS: tissues={len(tissues)}; "
        f"flight_samples={int(flight['sample_column'].nunique())}; "
        f"included_in_overall={len(included)}; z_capped={n_capped_total}; "
        f"overall_range=[{overall.min().min():.3f}, {overall.max().max():.3f}]"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
