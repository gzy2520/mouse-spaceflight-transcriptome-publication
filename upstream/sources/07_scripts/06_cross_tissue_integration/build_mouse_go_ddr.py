#!/usr/bin/env python3
"""Rebuild current DDR heatmaps with stable-ID-mapped mouse GO annotations.

The validated 2026-07-28 integration, direction, mission-clustering, final-FDR,
and plotting implementation is reused. The candidate family uses Ensembl Gene
IDs mapped from the TSV UniProtKB accessions; Symbol is display-only.
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE_PATH = (
    ROOT
    / "07_scripts/06_cross_tissue_integration"
    / "_ddr_heatmap_support.py"
)
DEFAULT_OUT = (
    ROOT
    / "03_analysis_results/06_cross_tissue_integration"
    / "04_mouse_go_ddr"
)
OUT = Path(os.environ.get("MOUSE_GO_DDR_OUTPUT_DIR", DEFAULT_OUT))
if not OUT.is_absolute():
    OUT = ROOT / OUT
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
ANNOTATION_DIR = ROOT / "08_GO_annotation"
DEFAULT_MAPPED_ANNOTATIONS = (
    ANNOTATION_DIR
    / "mgi_mod_gaf_20260804"
    / "final_manual_curated"
    / "03_mouse_GO_annotations_with_gene_ids.tsv.gz"
)
MAPPED_ANNOTATIONS = Path(
    os.environ.get(
        "MOUSE_GO_MAPPING_FILE",
        DEFAULT_MAPPED_ANNOTATIONS,
    )
)
if not MAPPED_ANNOTATIONS.is_absolute():
    MAPPED_ANNOTATIONS = ROOT / MAPPED_ANNOTATIONS
ANALYSIS_VARIANT = os.environ.get(
    "MOUSE_GO_ANALYSIS_VARIANT",
    "uniprot_tsv_stable_id",
)

GO_CONFIG = {
    "DNA Damage Response (GO:0006974)": (
        "DNA_damage_response_GO0006974",
        "in_DNA_damage_response_GO0006974",
    ),
    "DNA Repair (GO:0006281)": (
        "DNA_repair_GO0006281",
        "in_DNA_repair_GO0006281",
    ),
}


def load_base_module():
    spec = importlib.util.spec_from_file_location(
        "direct_ddr_final_fdr01_base",
        BASE_PATH,
    )
    if spec is None or spec.loader is None:
        raise ImportError(BASE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collapse_values(values: pd.Series) -> str:
    collected: set[str] = set()
    for value in values.dropna().astype(str):
        collected.update(part for part in value.split(";") if part)
    return ";".join(sorted(collected))

def count_collapsed_values(values: pd.Series) -> int:
    collapsed = collapse_values(values)
    return len(collapsed.split(";")) if collapsed else 0


def parse_target_go_mouse() -> tuple[dict[str, set[str]], pd.DataFrame]:
    table = pd.read_csv(
        MAPPED_ANNOTATIONS,
        sep="\t",
        encoding="utf-8-sig",
        low_memory=False,
    )
    required = {
        "requested_term_key",
        "gene_product_db",
        "gene_product_id",
        "SYMBOL",
        "ensembl_gene_id",
        "entrezgene_id",
        "mgi_id",
        "included_in_gene_set",
        "mapping_decision",
        "analysis_id_type",
        "symbol_used_as_analysis_key",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f"{MAPPED_ANNOTATIONS.name} missing columns: {sorted(missing)}"
        )
    included = (
        table["included_in_gene_set"]
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1"})
    )
    accepted = table.loc[included].copy()
    accepted["ensembl_gene_id"] = (
        accepted["ensembl_gene_id"].astype(str).str.strip().str.upper()
    )
    accepted = accepted.loc[
        accepted["ensembl_gene_id"].str.match(r"^ENSMUSG\d+$")
    ]
    if not accepted["analysis_id_type"].eq("Ensembl Gene ID").all():
        raise ValueError("Mapped GO annotations are not Ensembl-keyed")
    symbol_key_flag = (
        accepted["symbol_used_as_analysis_key"]
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1"})
    )
    if symbol_key_flag.any():
        raise ValueError("Symbol is incorrectly marked as an analysis key")

    go_sets: dict[str, set[str]] = {}
    audit_rows: list[dict[str, object]] = []
    for term_name, (term_key, _) in GO_CONFIG.items():
        annotations = accepted.loc[
            accepted["requested_term_key"].eq(term_key)
        ]
        genes = set(annotations["ensembl_gene_id"])
        if not genes:
            raise ValueError(f"No accepted Ensembl Gene IDs for {term_key}")
        go_sets[term_name] = genes
        audit_rows.append(
            {
                "term_name": term_name,
                "source_file": str(MAPPED_ANNOTATIONS.relative_to(ROOT)),
                "gene_set_rule": os.environ.get(
                    "MOUSE_GO_GENE_SET_RULE",
                    (
                        "accepted Ensembl Gene IDs mapped from TSV UniProtKB "
                        "accessions through MGI/NCBI Gene; Symbol not used as "
                        "analysis merge key"
                    ),
                ),
                "positive_annotation_rows": len(annotations),
                "unique_ensembl_gene_ids": len(genes),
                "unique_entrez_gene_ids": count_collapsed_values(
                    annotations["entrezgene_id"]
                ),
                "analysis_id_type": "Ensembl Gene ID",
                "symbol_used_as_analysis_key": False,
                "ensembl_release": 116,
                "ensembl_assembly": "GRCm39",
            }
        )

    damage = go_sets["DNA Damage Response (GO:0006974)"]
    repair = go_sets["DNA Repair (GO:0006281)"]
    if repair - damage:
        raise ValueError("Mapped DNA repair set is not a subset of damage response")
    lookup = (
        accepted.loc[
            accepted["ensembl_gene_id"].isin(damage | repair),
            [
                "ensembl_gene_id",
                "SYMBOL",
                "external_gene_name",
                "entrezgene_id",
                "mgi_id",
                "uniprot_base",
            ],
        ]
        .groupby("ensembl_gene_id", as_index=False)
        .agg(
            {
                "SYMBOL": collapse_values,
                "external_gene_name": collapse_values,
                "entrezgene_id": collapse_values,
                "mgi_id": collapse_values,
                "uniprot_base": collapse_values,
            }
        )
    )
    membership = pd.DataFrame({"ensembl_id": sorted(damage | repair)})
    membership = membership.merge(
        lookup,
        left_on="ensembl_id",
        right_on="ensembl_gene_id",
        how="left",
        validate="one_to_one",
    ).drop(columns="ensembl_gene_id")
    membership["in_DNA_damage_response_GO0006974"] = membership[
        "ensembl_id"
    ].isin(damage)
    membership["in_DNA_repair_GO0006281"] = membership["ensembl_id"].isin(
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
    TABLES.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audit_rows).to_csv(
        TABLES / "00_mouse_GO_source_audit.csv",
        index=False,
        encoding="utf-8-sig",
    )
    return go_sets, membership


def main() -> None:
    base = load_base_module()
    base.OUT = OUT
    base.TABLES = TABLES
    base.FIGURES = FIGURES
    base.GO_SET_ID_TYPE = "ensembl_gene_id"
    base.parse_target_go = parse_target_go_mouse
    base.main()

    explanation = (OUT / "00_DDR直接热图_最终FDR01说明.md").read_text(
        encoding="utf-8"
    )
    if ANALYSIS_VARIANT == "uniprot_tsv_stable_id":
        preface = """# 小鼠GO注释替换说明（2026-08-04）

本目录复用已验证的2026-07-28 DDR整合算法，只替换候选gene set来源和匹配主键：
使用用户提供的小鼠TSV，经
`UniProtKB accession → MGI/NCBI Gene → Ensembl Gene ID`
映射后，以Ensembl Gene ID直接匹配公共检测背景。Symbol仅作显示和歧义候选
验证，不作为分析merge键。ComplexPortal和RNAcentral条目不直接进入gene set。
单比较FDR、任务内处理、mission cluster等权、最终BH-FDR<0.10及Top 15
展示规则均未改变。旧版结果保留，不被覆盖。

"""
    elif ANALYSIS_VARIANT == "mgi_mod_gaf_final_manual_curated":
        preface = """# MGI MOD-GAF小鼠GO最终映射重分析说明（2026-08-04）

本目录复用已验证的DDR整合、方向、任务聚类和最终FDR算法，只替换候选gene set
来源。GO注释直接来自MGI MOD-GAF，以MGI ID为注释对象；先经完整Ensembl
release 116 BioMart、MGI MRK_ENSEMBL和MGI Entrez→Ensembl116路线自动核对，
再对5个严格映射例外逐项使用MGI/NCBI/Ensembl官方记录复核。Gcna、Csprs和
旧ID Rnf212b获得官方ID链确认；1700040F15Rik与Obox4因没有Ensembl116 Gene
ID而保留排除。Symbol仅作显示和检索入口，不作为候选裁决证据或分析merge键。
单比较FDR、任务内处理、mission cluster等权、最终BH-FDR<0.10及Top 15展示
规则均未改变。本变体是后续默认主版本；UniProt TSV稳定ID版保留为映射来源
敏感性对照。

"""
    else:
        preface = f"""# MGI MOD-GAF小鼠GO映射来源重分析说明（2026-08-04）

本目录复用已验证的DDR整合、方向、任务聚类和最终FDR算法，只替换候选gene set
来源。当前分析变体为`{ANALYSIS_VARIANT}`：GO注释直接来自MGI MOD-GAF，
以MGI ID为注释对象，经Ensembl release 116 BioMart、MGI MRK_ENSEMBL和
MGI Entrez→Ensembl116多路线交叉核对后，以Ensembl Gene ID匹配公共检测背景。
Symbol仅作显示，不用于一对多候选裁决。单比较FDR、任务内处理、mission
cluster等权、最终BH-FDR<0.10及Top 15展示规则均未改变。UniProt TSV稳定ID版
保留为映射来源敏感性对照，不被覆盖。

"""
    (OUT / "00_DDR直接热图_最终FDR01说明.md").write_text(
        preface + explanation,
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
