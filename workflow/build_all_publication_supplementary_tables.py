#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build all publication Supplementary Tables aligned with Teacher's Manuscript & Method citations.
Outputs both .xlsx workbooks and .csv files to /Users/gzy2520/Desktop/Supplementary_Tables/
and repo release directory.
"""

import os
import gzip
import pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DESKTOP_DIR = "/Users/gzy2520/Desktop/Supplementary_Tables"
CSV_DIR = os.path.join(DESKTOP_DIR, "csv")
os.makedirs(CSV_DIR, exist_ok=True)

FONT_FAMILY = "Calibri"
HEADER_FONT = Font(name=FONT_FAMILY, size=11, bold=True, color="FFFFFF")
HEADER_FILL = PatternFill(start_color="1F4E78", end_color="1F4E78", fill_type="solid")
HEADER_ALIGN = Alignment(horizontal="center", vertical="center", wrap_text=True)

DATA_FONT = Font(name=FONT_FAMILY, size=10, color="000000")
ZEBRA_FILL = PatternFill(start_color="F2F6F9", end_color="F2F6F9", fill_type="solid")
WHITE_FILL = PatternFill(start_color="FFFFFF", end_color="FFFFFF", fill_type="solid")

BORDER_THIN = Side(border_style="thin", color="D9D9D9")
CELL_BORDER = Border(left=BORDER_THIN, right=BORDER_THIN, top=BORDER_THIN, bottom=BORDER_THIN)

def format_worksheet(ws, df, is_first_col_bold=False):
    """Apply publication-quality formatting to openpyxl worksheet."""
    ws.views.sheetView[0].showGridLines = True
    ws.freeze_panes = "A2"
    
    ws.row_dimensions[1].height = 26
    for col_num in range(1, len(df.columns) + 1):
        cell = ws.cell(row=1, column=col_num)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = HEADER_ALIGN
        cell.border = CELL_BORDER

    num_rows = len(df)
    num_cols = len(df.columns)
    for row_idx in range(num_rows):
        row_num = row_idx + 2
        ws.row_dimensions[row_num].height = 19
        fill = ZEBRA_FILL if (row_idx % 2 == 1) else WHITE_FILL
        
        for col_idx in range(num_cols):
            cell = ws.cell(row=row_num, column=col_idx + 1)
            cell.font = DATA_FONT
            cell.fill = fill
            cell.border = CELL_BORDER
            
            val = cell.value
            if isinstance(val, (int,)):
                cell.alignment = Alignment(horizontal="right", vertical="center")
                cell.number_format = "#,##0"
            elif isinstance(val, (float,)):
                cell.alignment = Alignment(horizontal="right", vertical="center")
                if abs(val) < 0.001 and val != 0:
                    cell.number_format = "0.00E+00"
                else:
                    cell.number_format = "0.0000"
            else:
                if col_idx == 0 and is_first_col_bold:
                    cell.font = Font(name=FONT_FAMILY, size=10, bold=True, color="000000")
                cell.alignment = Alignment(horizontal="left", vertical="center")
                
    for col in ws.columns:
        col_letter = get_column_letter(col[0].column)
        max_len = 0
        for cell in col:
            v_str = str(cell.value or '')
            if len(v_str) > max_len:
                max_len = len(v_str)
        col_width = max(max_len + 3, 11)
        if col_width > 42:
            col_width = 42
        ws.column_dimensions[col_letter].width = col_width

def write_excel_sheets(output_xlsx_path, sheets_dict):
    """Write multiple dataframes to a styled Excel workbook."""
    wb = openpyxl.Workbook()
    default_sheet = wb.active
    wb.remove(default_sheet)
    
    for sheet_name, df in sheets_dict.items():
        ws = wb.create_sheet(title=sheet_name[:31])
        for col_idx, col_name in enumerate(df.columns, start=1):
            ws.cell(row=1, column=col_idx, value=col_name)
            
        for row_idx, row_data in enumerate(df.itertuples(index=False), start=2):
            for col_idx, value in enumerate(row_data, start=1):
                if pd.isna(value):
                    ws.cell(row=row_idx, column=col_idx, value="")
                else:
                    ws.cell(row=row_idx, column=col_idx, value=value)
                    
        format_worksheet(ws, df)
        
    wb.save(output_xlsx_path)
    print(f"Generated: {output_xlsx_path}")

def save_csvs(table_prefix, sheets_dict):
    """Save each sheet as a standalone CSV file."""
    for sheet_name, df in sheets_dict.items():
        clean_name = sheet_name.replace(" ", "_").replace("/", "_")
        csv_path = os.path.join(CSV_DIR, f"{table_prefix}_{clean_name}.csv")
        df.to_csv(csv_path, index=False)
        print(f"  CSV: {csv_path}")

def main():
    print("=== Step 1: Generating Table S1 (Animal Cohort & Sample Metadata) ===")
    df_s1_cohorts = pd.read_csv(os.path.join(REPO_ROOT, "reproduced_final_result_20260908_tree_axis_order_complete/provenance/Fig_1b_mouse_metadata_grouped.csv"))
    df_s1_samples = pd.read_csv(os.path.join(REPO_ROOT, "reproduced_final_result_20260908_tree_axis_order_complete/provenance/Fig_6_mouse_sample_metadata_complete_audit.csv"))
    
    s1_sheets = {
        "Cohort_Summary_Fig1b": df_s1_cohorts,
        "Sample_Metadata_Audit_761": df_s1_samples
    }
    s1_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S1_Cohort_and_Sample_Metadata.xlsx")
    write_excel_sheets(s1_path, s1_sheets)
    save_csvs("Table_S1", s1_sheets)

    print("\n=== Step 2: Generating Table S2 (DDR Overlap & Common DEGs) ===")
    df_s2_thymus_deg = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/common_direction/B_thymus_DDR_common_direction_top30_min_abs_log2fc_direction_ordered.csv"))
    df_s2_thymus_upset = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/venn/04_membership_matrix_thymus_four_comparisons.csv"))
    df_s2_kidney_deg = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/common_direction/A_kidney_DDR_common_direction_top30_min_abs_log2fc_direction_ordered.csv"))
    df_s2_kidney_upset = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/venn/04_membership_matrix_kidney_five_datasets.csv"))
    df_s2_up_deg = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/common_direction/D_NES_gt1_excl_Kidney_DDR_common_direction_top30_min_abs_log2fc_direction_ordered.csv"))
    df_s2_up_upset = pd.read_csv(os.path.join(REPO_ROOT, "results/tables/04_membership_matrix_up_tissues.csv"))
    df_s2_down_deg = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/common_direction/C_NES_lt_minus1_excl_Lung_Thymus_DDR_common_direction_top30_min_abs_log2fc_direction_ordered.csv"))
    df_s2_down_upset = pd.read_csv(os.path.join(REPO_ROOT, "results/tables/04_membership_matrix_down_tissues.csv"))

    s2_sheets = {
        "Thymus_Top30_DEGs": df_s2_thymus_deg,
        "Thymus_DDR_Overlap_893": df_s2_thymus_upset,
        "Kidney_Common12_DEGs": df_s2_kidney_deg,
        "Kidney_DDR_Overlap_893": df_s2_kidney_upset,
        "Up_Cohort_Common5_DEGs": df_s2_up_deg,
        "Up_Cohort_DDR_Overlap_893": df_s2_up_upset,
        "Down_Cohort_Common11_DEGs": df_s2_down_deg,
        "Down_Cohort_DDR_Overlap_893": df_s2_down_upset
    }
    s2_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S2_DDR_Overlap_and_Common_DEGs.xlsx")
    write_excel_sheets(s2_path, s2_sheets)
    save_csvs("Table_S2", s2_sheets)

    print("\n=== Step 3: Generating Table S3 (Seven DNA Repair Pathways Repertoires & Detection Audit) ===")
    df_s3_sizes = pd.read_csv(os.path.join(REPO_ROOT, "release/tables/08_Table_S3_pathway_repertoire_sizes.csv"))
    df_s3_tissue = pd.read_csv(os.path.join(REPO_ROOT, "release/tables/09_Table_S3_per_tissue_sample_summary.csv"))
    df_s3_samples = pd.read_csv(os.path.join(REPO_ROOT, "release/tables/01_seven_pathway_member_gene_counts_per_source_sample_wide.csv"))

    s3_sheets = {
        "Pathway_Repertoire_Sizes": df_s3_sizes,
        "Tissue_Detection_Summary": df_s3_tissue,
        "Sample_Pathway_Counts_761": df_s3_samples
    }
    s3_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S3_Seven_Pathways_Gene_Repertoires_and_Detection.xlsx")
    write_excel_sheets(s3_path, s3_sheets)
    save_csvs("Table_S3", s3_sheets)

    print("\n=== Step 4: Generating Table S4 (Core Repair Genes Expression & ANOVA) ===")
    df_s4_gene_anova = pd.read_csv(os.path.join(REPO_ROOT, "release/tables/log/02_gene_one_way_anova_log_summary.csv"))
    df_s4_pathway_anova = pd.read_csv(os.path.join(REPO_ROOT, "release/tables/log/01_pathway_one_way_anova_log_summary.csv"))
    df_s4_tissue_mean = pd.read_csv(os.path.join(REPO_ROOT, "release/tables/log/03_tissue_gene_mean_log_qsmooth_summary.csv"))
    with gzip.open(os.path.join(REPO_ROOT, "release/tables/log/04_sample_log_qsmooth_values_used.csv.gz"), "rt") as f:
        df_s4_samples = pd.read_csv(f)

    s4_sheets = {
        "ANOVA_Core_Genes": df_s4_gene_anova,
        "ANOVA_Pathways": df_s4_pathway_anova,
        "Tissue_Gene_Mean_Log2": df_s4_tissue_mean,
        "Sample_Log2_Values": df_s4_samples
    }
    s4_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S4_Repair_Genes_Expression_and_ANOVA.xlsx")
    write_excel_sheets(s4_path, s4_sheets)
    save_csvs("Table_S4", s4_sheets)

    print("\n=== Step 5: Generating Master Consolidated Table S1-S4 (Teacher Citations) ===")
    master_sheets = {
        "S1_Cohort_Summary": df_s1_cohorts,
        "S1_Sample_Metadata_761": df_s1_samples,
        "S2_Thymus_Top30_DEGs": df_s2_thymus_deg,
        "S2_Kidney_Common12_DEGs": df_s2_kidney_deg,
        "S2_Up_Cohort_5_DEGs": df_s2_up_deg,
        "S2_Down_Cohort_11_DEGs": df_s2_down_deg,
        "S3_Pathway_Sizes": df_s3_sizes,
        "S3_Tissue_Detection": df_s3_tissue,
        "S4_ANOVA_Core_Genes": df_s4_gene_anova,
        "S4_Tissue_Gene_Mean_Log2": df_s4_tissue_mean
    }
    master_path = os.path.join(DESKTOP_DIR, "Supplementary_Tables_S1_to_S4_Teacher_Citations.xlsx")
    write_excel_sheets(master_path, master_sheets)

    print("\n=== Step 6: Generating Additional Manuscript Supplementary Tables (S5 to S10) ===")
    # Table S5: GO Enrichment Matrix
    df_s5 = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/go/04_tissue_statistics_concrete_terms_and_context.csv"))
    s5_sheets = {"GO_Enrichment_Statistics": df_s5}
    s5_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S5_Cross_Tissue_GO_Enrichment_Matrix.xlsx")
    write_excel_sheets(s5_path, s5_sheets)
    save_csvs("Table_S5", s5_sheets)

    # Table S6: Hallmark GSEA
    df_s6 = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/hallmark/01_global_top_term_figure_selection.csv"))
    s6_sheets = {"Hallmark_Top12_Pathways": df_s6}
    s6_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S6_MSigDB_Hallmark_GSEA_Overview.xlsx")
    write_excel_sheets(s6_path, s6_sheets)
    save_csvs("Table_S6", s6_sheets)

    # Table S7: Cross-tissue Fisher-z
    df_s7 = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/go/log2fc_spearman_20260902/04_overall_fisher_z_spearman_rho_matrix.csv"))
    s7_sheets = {"Cross_Tissue_FisherZ_Matrix": df_s7}
    s7_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S7_Cross_Tissue_Pathway_Association_FisherZ.xlsx")
    write_excel_sheets(s7_path, s7_sheets)
    save_csvs("Table_S7", s7_sheets)

    # Table S8: Per-tissue Spearman
    df_s8 = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/go/log2fc_spearman_20260902/03_tissue_pairwise_spearman_log2fc_long.csv"))
    s8_sheets = {"Per_Tissue_Spearman_Long": df_s8}
    s8_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S8_Per_Tissue_Pathway_Spearman_Profiles.xlsx")
    write_excel_sheets(s8_path, s8_sheets)
    save_csvs("Table_S8", s8_sheets)

    # Table S9: DSB Meta-Analysis Matrix
    df_s9 = pd.read_csv(os.path.join(REPO_ROOT, "data/publication_input/meta/08_essential_components_counts_log2FC_tissue_matrix.csv"))
    s9_sheets = {"DSB_Meta_Analysis_26x29": df_s9}
    s9_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S9_DSB_Repair_Meta_Analysis_Matrix.xlsx")
    write_excel_sheets(s9_path, s9_sheets)
    save_csvs("Table_S9", s9_sheets)

    # Table S10: SSB Meta-Analysis Matrix
    df_s10 = pd.read_csv(os.path.join(REPO_ROOT, "release/tables/01_SSB_tissue_meta_log2FC_pvalue_matrix_without_Neil2.csv"))
    s10_sheets = {"SSB_Meta_Analysis_26x45": df_s10}
    s10_path = os.path.join(DESKTOP_DIR, "Supplementary_Table_S10_SSB_Repair_Meta_Analysis_Matrix.xlsx")
    write_excel_sheets(s10_path, s10_sheets)
    save_csvs("Table_S10", s10_sheets)

    print("\n=== Manifest generation ===")
    manifest_rows = [
        {"Table_ID": "Supplementary Table S1", "Title": "Animal Cohort Compilation and Multidimensional Biospecimen Metadata Atlas", "Teacher_Citation": "Fig. 1a-b, Table S1", "Sheets": "Cohort_Summary_Fig1b (52 cohorts), Sample_Metadata_Audit_761 (761 samples × 54 cols)", "File": "Supplementary_Table_S1_Cohort_and_Sample_Metadata.xlsx"},
        {"Table_ID": "Supplementary Table S2", "Title": "Spaceflight-Induced Co-Directionally Regulated DDR Genes and Repertoire Overlap", "Teacher_Citation": "Fig. 3a, Fig. S2a, Fig. S3ab, Fig. S4a, Table S2", "Sheets": "Thymus_Top30_DEGs (30), Thymus_Overlap (893), Kidney_DEGs (12), Kidney_Overlap (893), Up_DEGs (5), Up_Overlap (893), Down_DEGs (11), Down_Overlap (893)", "File": "Supplementary_Table_S2_DDR_Overlap_and_Common_DEGs.xlsx"},
        {"Table_ID": "Supplementary Table S3", "Title": "Seven DNA Damage Repair Pathways Gene Repertoire Sizes and Multi-Tissue Detection Audit", "Teacher_Citation": "Fig. S5, Table S3", "Sheets": "Pathway_Repertoire_Sizes (7 pathways), Tissue_Detection_Summary (182 rows), Sample_Pathway_Counts_761 (761 samples)", "File": "Supplementary_Table_S3_Seven_Pathways_Gene_Repertoires_and_Detection.xlsx"},
        {"Table_ID": "Supplementary Table S4", "Title": "Baseline Expression Profiles and ANOVA Heterogeneity of Core Repair Genes Across Tissues", "Teacher_Citation": "Fig. 4b,c, Fig. 5c, Fig. 4/5 Legends, Table S4", "Sheets": "ANOVA_Core_Genes (26 genes), ANOVA_Pathways (7 pathways), Tissue_Gene_Mean_Log2 (676 rows), Sample_Log2_Values (9,360 rows)", "File": "Supplementary_Table_S4_Repair_Genes_Expression_and_ANOVA.xlsx"},
        {"Table_ID": "Supplementary Table S5", "Title": "Cross-Tissue Gene Ontology Enrichment Profiling Matrix", "Teacher_Citation": "Fig. 1c", "Sheets": "GO_Enrichment_Statistics (390 rows: 26 tissues × 15 GO terms)", "File": "Supplementary_Table_S5_Cross_Tissue_GO_Enrichment_Matrix.xlsx"},
        {"Table_ID": "Supplementary Table S6", "Title": "MSigDB Hallmark Pathway GSEA Overview Across Missions", "Teacher_Citation": "Fig. S1", "Sheets": "Hallmark_Top12_Pathways (48 rows)", "File": "Supplementary_Table_S6_MSigDB_Hallmark_GSEA_Overview.xlsx"},
        {"Table_ID": "Supplementary Table S7", "Title": "Cross-Tissue Consensus Pathway Association Matrix (Fisher-z)", "Teacher_Citation": "Fig. 2", "Sheets": "Cross_Tissue_FisherZ_Matrix (15 × 15 correlation matrix)", "File": "Supplementary_Table_S7_Cross_Tissue_Pathway_Association_FisherZ.xlsx"},
        {"Table_ID": "Supplementary Table S8", "Title": "Per-Tissue Sample-Level Pathway Spearman Correlation Profiles", "Teacher_Citation": "Methods Section 3", "Sheets": "Per_Tissue_Spearman_Long (2,730 pairwise correlations)", "File": "Supplementary_Table_S8_Per_Tissue_Pathway_Spearman_Profiles.xlsx"},
        {"Table_ID": "Supplementary Table S9", "Title": "Double-Strand Break (DSB) Repair Multi-Mission Spaceflight Meta-Analysis Matrix", "Teacher_Citation": "Fig. 6a", "Sheets": "DSB_Meta_Analysis_26x29 (754 rows: 26 tissues × 29 components)", "File": "Supplementary_Table_S9_DSB_Repair_Meta_Analysis_Matrix.xlsx"},
        {"Table_ID": "Supplementary Table S10", "Title": "Single-Strand Break (SSB) Repair Multi-Mission Spaceflight Meta-Analysis Matrix", "Teacher_Citation": "Fig. 6b", "Sheets": "SSB_Meta_Analysis_26x45 (1,170 rows: 26 tissues × 45 components)", "File": "Supplementary_Table_S10_SSB_Repair_Meta_Analysis_Matrix.xlsx"},
    ]
    df_manifest = pd.DataFrame(manifest_rows)
    manifest_path = os.path.join(DESKTOP_DIR, "Supplementary_Tables_Manifest.csv")
    df_manifest.to_csv(manifest_path, index=False)
    print(f"Manifest written: {manifest_path}")

if __name__ == "__main__":
    main()
