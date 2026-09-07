#!/usr/bin/env python3
"""Generate the publication Methods Word document (.docx) organized strictly in Figure order."""

import os
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_ALIGN_VERTICAL
from docx.oxml import parse_xml
from docx.oxml.ns import nsdecls

def set_cell_background(cell, fill_hex):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{fill_hex}"/>')
    tcPr.append(shd)

def set_cell_margins(cell, top=100, bottom=100, left=140, right=140):
    tcPr = cell._tc.get_or_add_tcPr()
    tcMar = parse_xml(
        f'<w:tcMar {nsdecls("w")}>'
        f'  <w:top w:w="{top}" w:type="dxa"/>'
        f'  <w:bottom w:w="{bottom}" w:type="dxa"/>'
        f'  <w:left w:w="{left}" w:type="dxa"/>'
        f'  <w:right w:w="{right}" w:type="dxa"/>'
        f'</w:tcMar>'
    )
    tcPr.append(tcMar)

def set_three_line_borders(table):
    tblPr = table._tbl.tblPr
    borders = parse_xml(
        f'<w:tblBorders {nsdecls("w")}>'
        f'  <w:top w:val="single" w:sz="12" w:space="0" w:color="1F4E78"/>'
        f'  <w:bottom w:val="single" w:sz="12" w:space="0" w:color="1F4E78"/>'
        f'  <w:insideH w:val="none"/>'
        f'  <w:insideV w:val="none"/>'
        f'  <w:left w:val="none"/>'
        f'  <w:right w:val="none"/>'
        f'</w:tblBorders>'
    )
    tblPr.append(borders)
    for cell in table.rows[0].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tcBorders = parse_xml(
            f'<w:tcBorders {nsdecls("w")}>'
            f'  <w:bottom w:val="single" w:sz="6" w:space="0" w:color="2E75B6"/>'
            f'</w:tcBorders>'
        )
        tcPr.append(tcBorders)

def add_callout(doc, text_content, label="METHODOLOGICAL SPECIFICATION"):
    table = doc.add_table(rows=1, cols=1)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    table.columns[0].width = Inches(6.5)
    
    cell = table.cell(0, 0)
    set_cell_background(cell, "F2F5F9")
    set_cell_margins(cell, top=130, bottom=130, left=180, right=180)
    
    tcPr = cell._tc.get_or_add_tcPr()
    tcBorders = parse_xml(
        f'<w:tcBorders {nsdecls("w")}>'
        f'  <w:left w:val="single" w:sz="24" w:space="0" w:color="1F4E78"/>'
        f'  <w:top w:val="none"/>'
        f'  <w:right w:val="none"/>'
        f'  <w:bottom w:val="none"/>'
        f'</w:tcBorders>'
    )
    tcPr.append(tcBorders)
    
    p = cell.paragraphs[0]
    p.paragraph_format.line_spacing = 1.15
    p.paragraph_format.space_after = Pt(2)
    run_label = p.add_run(f"[{label}]\n")
    run_label.font.name = "Arial"
    run_label.font.size = Pt(9.5)
    run_label.font.bold = True
    run_label.font.color.rgb = RGBColor(0x1F, 0x4E, 0x78)
    
    run_text = p.add_run(text_content)
    run_text.font.name = "Calibri"
    run_text.font.size = Pt(9.5)
    run_text.font.italic = True
    run_text.font.color.rgb = RGBColor(0x33, 0x33, 0x33)
    
    p_after = doc.add_paragraph()
    p_after.paragraph_format.space_before = Pt(0)
    p_after.paragraph_format.space_after = Pt(3)

def build_docx(output_path):
    doc = Document()
    
    # Page setup: Standard 1 inch margins
    for section in doc.sections:
        section.top_margin = Inches(1.0)
        section.bottom_margin = Inches(1.0)
        section.left_margin = Inches(1.0)
        section.right_margin = Inches(1.0)
        
        header = section.header
        p_hdr = header.paragraphs[0]
        p_hdr.alignment = WD_ALIGN_PARAGRAPH.RIGHT
        r_hdr = p_hdr.add_run("Methods: Murine Spaceflight Transcriptomic Atlas & DNA Damage Repair Profiling")
        r_hdr.font.name = "Arial"
        r_hdr.font.size = Pt(8.5)
        r_hdr.font.color.rgb = RGBColor(0x88, 0x88, 0x88)
        
        footer = section.footer
        p_ftr = footer.paragraphs[0]
        p_ftr.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r_ftr = p_ftr.add_run("Methods Section — NASA OSDR Multi-Tissue Spaceflight Dataset (Figure-Aligned)")
        r_ftr.font.name = "Arial"
        r_ftr.font.size = Pt(8.5)
        r_ftr.font.color.rgb = RGBColor(0x88, 0x88, 0x88)
        
    COLOR_PRIMARY = RGBColor(0x1F, 0x4E, 0x78)
    COLOR_SECONDARY = RGBColor(0x2E, 0x75, 0xB6)
    COLOR_BODY = RGBColor(0x26, 0x26, 0x26)
    
    # Title
    p_title = doc.add_paragraph()
    p_title.paragraph_format.space_before = Pt(0)
    p_title.paragraph_format.space_after = Pt(4)
    r_title = p_title.add_run("Methods")
    r_title.font.name = "Arial"
    r_title.font.size = Pt(22)
    r_title.font.bold = True
    r_title.font.color.rgb = COLOR_PRIMARY
    
    # Subtitle
    p_sub = doc.add_paragraph()
    p_sub.paragraph_format.space_before = Pt(0)
    p_sub.paragraph_format.space_after = Pt(12)
    r_sub = p_sub.add_run("Multi-Tissue Murine Spaceflight Transcriptomic Atlas and DNA Damage Repair Pathway Profiling")
    r_sub.font.name = "Arial"
    r_sub.font.size = Pt(13)
    r_sub.font.bold = True
    r_sub.font.color.rgb = COLOR_SECONDARY
    
    # Governance Callout
    add_callout(
        doc,
        "Organization & Rigor Contract: This Methods section is structured strictly in the sequential order of publication figures "
        "(Figure 1b to Figure 6, Supplementary Figures S1 to S6, and statistical audit tables). The text faithfully mirrors the actual "
        "computational implementation in the verified repository code, including the formal renumbering of DNA repair linear boxplots "
        "into Fig. 4b–d (DSB pathways: NHEJ, HR, A-EJ) and Fig. 5b–e (SSB pathways: BER, NER, MMR, FA), followed by the respective meta-dotplots "
        "(Fig. 4e and Fig. 5f). All analyses adhere to immutable Ensembl Gene IDs (GRCm38/GRCm39) as computational primary keys and "
        "deterministic pseudo-random seed locking at seed = 25.",
        "FIGURE-ALIGNED METHODS SPECIFICATION"
    )
    
    def add_sec_heading(text):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(14)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.keep_with_next = True
        r = p.add_run(text)
        r.font.name = "Arial"
        r.font.size = Pt(12)
        r.font.bold = True
        r.font.color.rgb = COLOR_PRIMARY
        return p

    def add_subsec_heading(text):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(9)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.keep_with_next = True
        r = p.add_run(text)
        r.font.name = "Arial"
        r.font.size = Pt(11)
        r.font.bold = True
        r.font.color.rgb = COLOR_SECONDARY
        return p

    def add_body_p(text, bold_prefix=None):
        p = doc.add_paragraph()
        p.paragraph_format.line_spacing = 1.22
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(4.5)
        if bold_prefix:
            r_pre = p.add_run(bold_prefix)
            r_pre.font.name = "Calibri"
            r_pre.font.size = Pt(10)
            r_pre.font.bold = True
            r_pre.font.color.rgb = COLOR_BODY
        r = p.add_run(text)
        r.font.name = "Calibri"
        r.font.size = Pt(10)
        r.font.color.rgb = COLOR_BODY
        return p

    # ---------------------------------------------------------------------------
    # Section 1: Figure 1b & Figure 6
    # ---------------------------------------------------------------------------
    add_sec_heading("1. Experimental Design, Animal Cohort Compilation, and Multidimensional Metadata (Figure 1b, Figure 6)")
    
    add_body_p(
        "To establish an organism-wide, multi-tissue atlas of transcriptional adaptations and DNA damage repair alterations under spaceflight "
        "environmental stressors, transcriptomic profiles and experimental metadata were curated from the NASA Open Science Data Repository "
        "(OSDR, formerly GeneLab; https://osdr.nasa.gov/). The analytical cohort encompasses 48 independent transcriptomic accessions (OSD datasets) "
        "across 59 primary analysis units, comprising 761 murine RNA-seq biospecimens spanning 26 distinct anatomical tissues (Figure 1b, Figure 6).",
        bold_prefix="Cohort Scope and Biospecimen Selection: "
    )
    
    add_body_p(
        "Biological and study metadata were programmatically extracted directly from the official NASA OSDR RESTful API and ISA-Tab archive "
        "ZIP files (R/data/01_refresh_osdr_mouse_metadata.R). The unified cohort includes 360 spaceflight mice (Flight / μG) and 401 ground-control "
        "mice (Ground control, Vivarium control, or Ground centrifuge / 1G baseline control). Exact donor animal identities (accession::Source Name), "
        "flight mission assignments across 13 mission clusters (SpaceX-4 to SpaceX-21, NG-11), biological sex (female and male), and donor chronological "
        "age (15 discrete age brackets, ranging from 8 to 32 weeks) were audited to resolve legacy ISA syntax discrepancies (such as confirming the "
        "current API mission assignment SpaceX-8 for OSD-162) without artificial age conversion or continuous imputation. Technical replicates were "
        "collapsed to ensure biological independence across all downstream tests.",
        bold_prefix="Metadata Audit and Conflict Governance: "
    )
    
    add_body_p(
        "Experimental design and animal allocation across 26 tissues were visualized using a stratified bubble plot (Figure 1b, Figure 6). "
        "The horizontal axis follows the canonical cross-tissue DNA damage response (DDR) activation rank. Vertically, each mission cluster "
        "is bifurcated into paired female and male tracks, with point area scaling with unique donor mouse counts. Age brackets are mapped to a "
        "fixed chronological 15-color palette (yellow -> green -> blue -> purple) without continuous interpolation. The single-sample audit table "
        "Fig_6_mouse_sample_metadata_complete_audit.csv contains 761 rows and 54 audited attributes verifying complete sample-to-study traceability.",
        bold_prefix="Metadata Bubble Visualization and Complete Audit: "
    )

    # ---------------------------------------------------------------------------
    # Section 2: Figure 1c & Supplementary Figure S1
    # ---------------------------------------------------------------------------
    add_sec_heading("2. Functional Pathway Enrichment and Cross-Tissue Gene Ontology Profiling (Figure 1c, Supplementary Figure S1)")
    
    add_body_p(
        "Functional pathway enrichment was evaluated via pre-ranked Gene Set Enrichment Analysis using fgseaMultilevel (fgsea package; R 4.4.3), "
        "with pseudo-random seeds locked at seed = 25. Differential expression was modeled within each accession using negative binomial GLMs "
        "in DESeq2, with contrast orientation audited to ensure Spaceflight vs. Ground directionality. Genes within each comparison were ranked by "
        "signed differential Wald statistics: S = sign(log2FC) · (-log10 P_nominal).",
        bold_prefix="Pre-Ranked fgsea Multilevel Modeling: "
    )
    
    add_body_p(
        "Analyses queried the MSigDB Hallmark collection (2026.1.Mm; Supplementary Figure S1) and 15 curated Gene Ontology Biological Process (GO:BP) "
        "terms (Figure 1c). The GO panel incorporates the core root term for DNA damage response (GO:0006974), sub-pathway branches (DNA repair GO:0006281; "
        "DNA damage signal transduction GO:0042770; Intrinsic apoptotic signaling GO:0008630; DNA damage tolerance GO:0006301; Telomere maintenance in "
        "response to DNA damage GO:0043247; Chromosome telomeric region GO:0000781), physiological context terms (Cell cycle GO:0007049; DNA replication "
        "GO:0006260; Transcription GO:0006351; Translation GO:0006412; Mitochondrion GO:0005739; Cytoskeleton GO:0005856), and comparative mechanosensory "
        "processes (Sensory perception of mechanical stimulus GO:0050954; Cell death GO:0008219).",
        bold_prefix="Curated Pathway Repertoire: "
    )
    
    add_body_p(
        "To prevent large missions from exerting disproportionate leverage, Normalized Enrichment Scores (NES) were aggregated via unweighted "
        "arithmetic averaging across contributing missions within each tissue. Empirical significance was verified via exact sign-flip two-sided "
        "permutation testing across all 2^m sign configurations, coupled with Stouffer meta-Z transformation and Benjamini-Hochberg False Discovery "
        "Rate (FDR) multiple testing correction. Figure 1c presents the complete 26 × 15 matrix with a fine blue–white–red continuous gradient.",
        bold_prefix="Mission-Equal NES and Exact Permutation Statistics: "
    )

    # ---------------------------------------------------------------------------
    # Section 3: Figure 2 & Supplementary Figure S6
    # ---------------------------------------------------------------------------
    add_sec_heading("3. Sample-Level Pathway Co-Regulation and Cross-Tissue Fisher-z Integration (Figure 2, Supplementary Figure S6)")
    
    add_body_p(
        "To capture coordinated biological regulation across the 15 GO processes without suffering from ecological aggregation fallacies, "
        "a sample-level correlation framework was executed (python/plot_log2fc_per_tissue_pathway_spearman_20260902.py; Figure 2, Supplementary Figure S6).",
        bold_prefix="Analytical Rationale: "
    )
    
    add_body_p(
        "1. Ground-Referenced Sample-Level log2FC: For each flight animal k within a given accession and tissue, sample-level expression was referenced "
        "against matched ground controls:\n"
        "        Sample_log2FC_k = log2(Expr_flight_k + 1) - <log2(Expr_ground + 1)>_accession, tissue\n"
        "The sample-level pathway activation score was defined as the mean log2FC of constituent member genes.\n"
        "2. Intra-Tissue Spearman Correlation: Within each of the 26 tissues, taking individual flight mice (n = 360) as independent observation units, "
        "pairwise non-parametric Spearman rank correlation coefficients (ρ_ij) were computed across the 15 GO terms, generating 26 distinct 15 × 15 "
        "correlation matrices (Supplementary Figure S6).\n"
        "3. Fisher-z Meta-Aggregation: Intra-tissue correlation coefficients were transformed into normally distributed z-scores:\n"
        "        z = arctanh(ρ) = 0.5 · ln[(1 + ρ) / (1 - ρ)]\n"
        "Correlations were bounded within [-1 + 10^(-12), 1 - 10^(-12)], and z-scores were capped at |z| ≤ 5.0.\n"
        "4. Cross-Tissue Pooling Criteria: An unweighted arithmetic mean z̄ was calculated across the 23 eligible tissues satisfying the sample-size "
        "threshold (n_flight ≥ 5; Colon [n = 4], Heart [n = 3], and Femoral skin [n = 2] were retained in individual panels but excluded from cross-tissue "
        "pooling). The pooled z̄ was back-transformed to derive the overall meta-analytic Spearman rank correlation:\n"
        "        ρ̄ = tanh(z̄) = [exp(2z̄) - 1] / [exp(2z̄) + 1]\n"
        "The resulting 15 × 15 consensus matrix is presented in Figure 2.",
        bold_prefix="Mathematical Formulation and Pooling Criteria: "
    )

    # ---------------------------------------------------------------------------
    # Section 4: Figure 3a, 3b & Supplementary Figures S2–S4
    # ---------------------------------------------------------------------------
    add_sec_heading("4. Organ-Specific Repertoires, UpSet Intersections, and Common-Direction DEG Prioritization (Figure 3a, 3b, Supplementary Figures S2–S4)")
    
    add_body_p(
        "To examine shared and tissue-specific responsive elements within the core DDR machinery (GO:0006974), overlapping DEG repertoires "
        "were evaluated across multi-dataset organs, specifically Thymus (4 flight comparisons, Figure 3a) and Kidney (5 independent OSDR datasets, "
        "Supplementary Figure S2a) using ggVennDiagram. Organ-wide intersection topology across NES > 1 tissues (excluding Kidney, Supplementary Figure S3a) "
        "and NES < -1 tissues (excluding Lung and Thymus, Supplementary Figure S4a) was resolved using UpSet intersection geometry (R/figures/06_go_upsets.R).",
        bold_prefix="Repertoire Overlap and UpSet Topology: "
    )
    
    add_body_p(
        "To prioritize core genes demonstrating conserved responsiveness independent of mission idiosyncrasies, common-direction filtering "
        "was executed (python/plot_common_direction_heatmaps.py; Figure 3b, Supplementary Figure S2b, S3b, S4b). Candidate DEGs were required to demonstrate "
        "uniform directional concordance (all log2FC > 0 or all log2FC < 0) across all parallel dataset comparisons at nominal P < 0.05 (or P < 0.10 for global "
        "NES-stratified cohorts). Qualifying genes were ranked by their minimum absolute fold change across comparisons:\n"
        "        Score_g = min_c |log2FC_{g, c}|\n"
        "The top 30 ranked common-direction genes were displayed as heatmaps with standardized blue–white–red color gradients and exact log2FC text overlays.",
        bold_prefix="Common-Direction Screening Formulation: "
    )

    # ---------------------------------------------------------------------------
    # Section 5: Figure 4a–e (DSB Module)
    # ---------------------------------------------------------------------------
    add_sec_heading("5. Double-Strand Break (DSB) Repair Architecture: Baseline Clustering, Individual Linear Boxplots, and Spaceflight Meta-Analysis (Figure 4a–e)")
    
    add_subsec_heading("5.1 Baseline Quantile Smoothing (YARN / qsmooth) and Hierarchical Clustering (Figure 4a)")
    add_body_p(
        "To mitigate technical batch artifacts while preserving bona fide tissue-specific transcriptional baselines, cross-tissue expression normalization "
        "was performed using the YARN framework implementing tissue-aware quantile smoothing (qsmooth). For the DSB repair module, 29 curated components "
        "spanning NHEJ (7 components), HR (14 components), HR–A-EJ shared factors (3 components), and A-EJ (5 components) were assembled. "
        "Pairwise topological dissimilarity between tissues Ta and Tb was computed on baseline normalized expression (YARNNormalizedLog2) via Spearman "
        "rank distance: D(Ta, Tb) = 1 - ρ_Spearman(Ta, Tb). Hierarchical agglomerative clustering was executed using average linkage "
        "(hclust(as.dist(D), method = 'average')), defining the stable dendrogram order for Figure 4a.",
        bold_prefix="Normalization and Dendrogram Clustering: "
    )
    
    add_subsec_heading("5.2 Tissue-Stratified Individual Linear Abundance Modeling and One-Way ANOVA (Figure 4b–d)")
    add_body_p(
        "To assess physiological baseline expression and organ-level heterogeneity across the mammalian DSB repair machinery, key biological marker "
        "genes were selected and displayed on a linear expression scale (Expression_linear = 2^(YARNNormalizedLog2) - 1) across the 26 canonical tissues:\n"
        "• Figure 4b (NHEJ, 4 genes): Nhej1 / XLF (ENSMUSG00000026162), Paxx (ENSMUSG00000047617), Xrcc6 / Ku70 (ENSMUSG00000022471), Xrcc5 / Ku80 (ENSMUSG00000026187).\n"
        "• Figure 4c (HR, 4 genes): Brca1 (ENSMUSG00000017146), Bard1 (ENSMUSG00000026196), Blm (ENSMUSG00000030528), Rad51 (ENSMUSG00000027323).\n"
        "• Figure 4d (A-EJ, 4 genes): Parp1 (ENSMUSG00000026496), Polq (ENSMUSG00000034206), Lig1 (ENSMUSG00000056394), Lig3 (ENSMUSG00000020697).\n"
        "Boxplot Construction: The lower hinge represents Q1 (25%); the center solid line represents the sample mean (μ); the upper hinge represents Q3 (75%). "
        "Whiskers extend to the furthest observations within 1.5 × IQR. Individual flight mice biospecimens (n = 360) are overlaid as jittered points "
        "(dodge.width = 0.80, jitter.width = 0.11, seed = 25). Dashed lines connect tissue-level mean values for each gene.\n"
        "One-Way ANOVA: Tissue heterogeneity was evaluated via ordinary One-way ANOVA. Both gene-level tests (df_group = 25, df_residual = 334) and "
        "pathway-level tests (modeled on per-mouse mean pathway abundance) were conducted, with exact F-statistics and scientific P-values printed "
        "in figure subtitles: NHEJ (Pathway F = 47.99, P = 6.82 × 10^(-95)), HR (Pathway F = 39.53, P = 2.37 × 10^(-84)), and A-EJ (Pathway F = 42.60, P = 2.37 × 10^(-88)).",
        bold_prefix="Linear Scale Formulation, Boxplot Geometry, and ANOVA Testing: "
    )
    
    add_subsec_heading("5.3 Spaceflight-Induced Meta-Differential Expression Analysis (Figure 4e)")
    add_body_p(
        "To characterize spaceflight-induced perturbations across all 29 DSB components, multi-mission meta-analysis was performed (Figure 4e). "
        "Effect sizes are presented as the median of mission-level median log2FC values. Meta-significance across missions was integrated via signed "
        "Stouffer's Z-trend test: Z_meta = [∑ sign(log2FC_i) · Φ^(-1)(1 - P_i / 2)] / √m. Dot color reflects median log2FC and dot diameter scales with -log10(P_meta).",
        bold_prefix="Meta-Dotplot Construction: "
    )

    # ---------------------------------------------------------------------------
    # Section 6: Figure 5a–f (SSB Module)
    # ---------------------------------------------------------------------------
    add_sec_heading("6. Single-Strand Break (SSB) Repair Architecture: Baseline Clustering, Individual Linear Boxplots, and Spaceflight Meta-Analysis (Figure 5a–f)")
    
    add_subsec_heading("6.1 Baseline Quantile Smoothing, Neil2 Exclusion Audit, and Hierarchical Clustering (Figure 5a)")
    add_body_p(
        "Cross-tissue baseline normalization for the single-strand break and base excision machinery encompassed 45 curated components spanning "
        "BER (6 components), NER (18 components), MMR (5 components), and FA (16 components). Data completeness audits across all flight biospecimens "
        "identified recurrent sequencing dropout and high zero-count sparsity for Neil2 (ENSMUSG00000035121); Neil2 was therefore systematically excluded "
        "from the 45-component SSB matrix and audited in 02_SSB_component_stable_id_mapping_without_Neil2.csv. "
        "Pairwise topological distance was computed on baseline YARNNormalizedLog2 via Spearman rank distance (D = 1 - ρ), followed by average-linkage "
        "hierarchical clustering to generate Figure 5a.",
        bold_prefix="Quality Exclusion Audit and Dendrogram Topology: "
    )
    
    add_subsec_heading("6.2 Tissue-Stratified Individual Linear Abundance Modeling and One-Way ANOVA (Figure 5b–e)")
    add_body_p(
        "Key marker genes across the four SSB pathways were evaluated on the linear expression scale (2^(log2) - 1) across 360 flight mice in 26 tissues:\n"
        "• Figure 5b (BER, 3 genes): Ung (ENSMUSG00000029591), Ogg1 (ENSMUSG00000030271), Neil1 (ENSMUSG00000032298).\n"
        "• Figure 5c (NER, 7 genes): Xpc (ENSMUSG00000030094), Rad23b (ENSMUSG00000028426), Cetn2 (ENSMUSG00000031347), Ddb1 (ENSMUSG00000024740), "
        "Ddb2 (ENSMUSG00000002109), Ercc6 (ENSMUSG00000054051), Ercc8 (ENSMUSG00000021694).\n"
        "  NER Glyphic Disambiguation: Recognition complexes are encoded by paired color-and-shape glyphs: XPC/RAD23B/CETN2 (amber; cross ✕, triangle ▲, circle ●); "
        "DDB1/DDB2 (teal; cross ✕, triangle ▲); ERCC6/ERCC8 (blue; cross ✕, triangle ▲).\n"
        "• Figure 5d (MMR, 2 genes): Msh2 (ENSMUSG00000024151), Msh3 (ENSMUSG00000014850).\n"
        "• Figure 5e (FA, 2 genes): Fancd2 (ENSMUSG00000034023), Fanci (ENSMUSG00000039187). (df_group = 25, df_residual = 330 due to 4 missing observations).\n"
        "One-Way ANOVA: Cross-tissue expression heterogeneity was highly significant across all pathways: BER (Pathway F = 32.44, P = 3.49 × 10^(-74)), "
        "NER (Pathway F = 89.47, P = 7.98 × 10^(-132)), MMR (Pathway F = 56.20, P = 7.71 × 10^(-104)), and FA (Pathway F = 31.47, P = 2.86 × 10^(-72)).",
        bold_prefix="Marker Gene Selection, Glyphic Coding, and ANOVA Testing: "
    )
    
    add_subsec_heading("6.3 Spaceflight-Induced Meta-Differential Expression Analysis (Figure 5f)")
    add_body_p(
        "To assess spaceflight transcriptional shifts across all 45 SSB components, multi-mission meta-analysis was performed (Figure 5f). "
        "Effect sizes represent the median of mission-level median log2FC values, and statistical significance was pooled using signed Stouffer's "
        "Z-trend integration, scaled by nominal P-value significance.",
        bold_prefix="Meta-Dotplot Construction: "
    )

    # ---------------------------------------------------------------------------
    # Section 7: Supplementary Figure S5
    # ---------------------------------------------------------------------------
    add_sec_heading("7. Seven-Pathway Repertoire Library Representation and Gene Coverage Auditing (Supplementary Figure S5)")
    
    add_body_p(
        "To verify that pathway-level transcriptomic analyses were not compromised by library dropouts or incomplete capture across sequencing accessions, "
        "per-sample gene detection was audited against a master catalog of 335 DNA repair genes: BER (n = 54), NER (n = 53), MMR (n = 27), FA (n = 39), "
        "HR (n = 135), A-EJ (n = 9), and NHEJ (n = 18). Supplementary Figure S5 tracks the count of detected member genes across all 761 source samples "
        "(stratified by Ground and Flight), confirming comprehensive transcript coverage across all 26 anatomical tissues.",
        bold_prefix="Catalog Completeness and Library Representation: "
    )

    # ---------------------------------------------------------------------------
    # Section 8: Environment and Reproducibility
    # ---------------------------------------------------------------------------
    add_sec_heading("8. Computational Environment, Code Governance, and Reproducibility Standards")
    
    add_body_p(
        "All statistical modeling, data transformations, and figure generation were executed in R (version 4.4.3) and Python (version 3.13). "
        "Core R dependencies included data.table (v1.16.4), ggplot2 (v3.5.1), patchwork (v1.3.0), scales (v1.3.0), and ggVennDiagram (v1.5.2), "
        "managed via renv.lock. Python workflows utilized matplotlib (v3.10.0), pandas (v2.2.3), and numpy (v2.2.3).",
        bold_prefix="Software Architecture: "
    )
    
    add_body_p(
        "All stochastic operations, point jittering, and permutation tests were locked to random seed 25. Complete analytical workflows, scripts, "
        "and data matrices are version-controlled via GitHub (git). All released figure panels, data tables, and input manifests were independently "
        "audited and locked with cryptographic SHA-256 checksums (provenance/reference_sha256.csv), ensuring 100% deterministic reproducibility.",
        bold_prefix="Seed Locking and Version Control: "
    )

    # ---------------------------------------------------------------------------
    # Section 9: Reference Matrix
    # ---------------------------------------------------------------------------
    add_sec_heading("9. Production Figure and Analytical Workflow Cross-Reference Matrix")
    
    add_body_p(
        "The following audited reference matrix connects every publication figure and statistical table in sequential order to its underlying "
        "mathematical model, sample scope, and analytical script:",
        bold_prefix="Provenance Cross-Reference: "
    )
    
    table_data = [
        ["Figure / Table", "Analytical Scope", "Input Cohort / Sample Size", "Statistical Model / Algorithm", "Source Script"],
        ["Fig. 1b / Fig. 6", "Experimental Design Atlas", "761 samples, 48 accessions, 26 tissues", "Stratified categorical bubble matrix (15 discrete age colors)", "R/figures/08_mouse_metadata_bubble.R"],
        ["Fig. 1c", "Cross-Tissue GO Profiling", "26 tissues × 15 GO biological terms", "Mission-equal unweighted mean NES, exact sign-flip P", "R/figures/01_go_overview.R"],
        ["Fig. 2", "Cross-Tissue Co-Regulation", "23 eligible tissues (n_flight ≥ 5)", "Sample-level log2FC, intra-tissue Spearman ρ, Fisher-z pooling", "python/plot_log2fc_per_tissue_pathway_spearman_20260902.py"],
        ["Fig. 3a", "Thymus Repertoire Overlap", "Thymus (4 flight comparisons)", "Set intersection membership, ggVennDiagram", "R/figures/03_kidney_thymus_venn.R"],
        ["Fig. 3b", "Thymus Common-Direction DEGs", "Thymus (4 flight comparisons)", "Concordant directionality, ranked by min |log2FC|", "python/plot_common_direction_heatmaps.py"],
        ["Fig. 4a", "DSB Baseline Clustering", "26 tissues, 29 DSB components", "YARN qsmooth baseline, 1 - Spearman ρ, average linkage", "R/figures/04_qsmooth_tree_heatmaps.R"],
        ["Fig. 4b", "NHEJ Linear Boxplot", "360 flight mice × 4 genes (26 tissues)", "Linear scale 2^log2 - 1, Mean hinge, Gene & Pathway ANOVA", "R/figures/10_pathway_gene_tissue_boxplots_linear.R"],
        ["Fig. 4c", "HR Linear Boxplot", "360 flight mice × 4 genes (26 tissues)", "Linear scale 2^log2 - 1, Mean hinge, Gene & Pathway ANOVA", "R/figures/10_pathway_gene_tissue_boxplots_linear.R"],
        ["Fig. 4d", "A-EJ Linear Boxplot", "360 flight mice × 4 genes (26 tissues)", "Linear scale 2^log2 - 1, Mean hinge, Gene & Pathway ANOVA", "R/figures/10_pathway_gene_tissue_boxplots_linear.R"],
        ["Fig. 4e", "DSB Spaceflight Meta-Analysis", "26 tissues × 29 DSB components", "Median of mission medians log2FC, signed Stouffer meta-Z", "R/figures/05_meta_log2fc_heatmaps.R"],
        ["Fig. 5a", "SSB Baseline Clustering", "26 tissues, 45 SSB components", "YARN qsmooth baseline (Neil2 excluded), 1 - ρ, average linkage", "R/figures/04_qsmooth_tree_heatmaps.R"],
        ["Fig. 5b", "BER Linear Boxplot", "360 flight mice × 3 genes (26 tissues)", "Linear scale 2^log2 - 1, Mean hinge, Gene & Pathway ANOVA", "R/figures/10_pathway_gene_tissue_boxplots_linear.R"],
        ["Fig. 5c", "NER Linear Boxplot", "360 flight mice × 7 genes (26 tissues)", "Linear scale 2^log2 - 1, Glyphic coding, Gene & Pathway ANOVA", "R/figures/10_pathway_gene_tissue_boxplots_linear.R"],
        ["Fig. 5d", "MMR Linear Boxplot", "360 flight mice × 2 genes (26 tissues)", "Linear scale 2^log2 - 1, Mean hinge, Gene & Pathway ANOVA", "R/figures/10_pathway_gene_tissue_boxplots_linear.R"],
        ["Fig. 5e", "FA Linear Boxplot", "360 flight mice × 2 genes (26 tissues)", "Linear scale 2^log2 - 1, Mean hinge, Gene & Pathway ANOVA", "R/figures/10_pathway_gene_tissue_boxplots_linear.R"],
        ["Fig. 5f", "SSB Spaceflight Meta-Analysis", "26 tissues × 45 SSB components", "Median of mission medians log2FC, signed Stouffer meta-Z", "R/figures/05_meta_log2fc_heatmaps.R"],
        ["Fig. S1", "Hallmark GSEA Overview", "13 missions with complete coverage", "Top 12 global Hallmark terms (6 positive, 6 negative NES)", "R/figures/02_hallmark_gsea.R"],
        ["Fig. S2a / S2b", "Kidney Overlap & Common-Dir", "Kidney (5 OSDR accessions)", "Venn intersection & common-direction top 30 DEGs", "R/figures/03_kidney_thymus_venn.R"],
        ["Fig. S3a / S3b", "NES > 1 UpSet & Common-Dir", "NES > 1 tissues (excl. Kidney)", "UpSet binary membership & common-direction top 30", "R/figures/06_go_upsets.R"],
        ["Fig. S4a / S4b", "NES < -1 UpSet & Common-Dir", "NES < -1 tissues (excl. Lung/Thymus)", "UpSet binary membership & common-direction top 30", "R/figures/06_go_upsets.R"],
        ["Fig. S5", "Pathway Gene Coverage", "761 samples × 7 pathways (335 genes)", "Per-sample member gene library detection tracking", "R/figures/07_seven_pathway_counts.R"],
        ["Fig. S6", "Per-Tissue Co-Regulation", "26 individual tissue matrices", "Sample-level 15 × 15 Spearman rank correlation matrices", "python/plot_log2fc_per_tissue_pathway_spearman_20260902.py"],
        ["Tables 01–04", "Statistical Audit Tables", "9,360 sample observations, 26 tissues", "Pathway ANOVA, Gene ANOVA, Tissue Means, Sample Linear Values", "results/pathway_gene_tissue_boxplots_linear/tables/"]
    ]
    
    table = doc.add_table(rows=len(table_data), cols=5)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    
    col_widths = [Inches(1.2), Inches(1.3), Inches(1.3), Inches(1.5), Inches(1.2)]
    for i, col in enumerate(table.columns):
        col.width = col_widths[i]
        
    set_three_line_borders(table)
    
    for row_idx, row_content in enumerate(table_data):
        row = table.rows[row_idx]
        is_header = (row_idx == 0)
        
        trPr = row._tr.get_or_add_trPr()
        if is_header:
            trPr.append(parse_xml(f'<w:tblHeader {nsdecls("w")}/>'))
            
        for col_idx, text in enumerate(row_content):
            cell = row.cells[col_idx]
            cell.width = col_widths[col_idx]
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            
            if is_header:
                set_cell_background(cell, "EBF1F5")
                set_cell_margins(cell, top=100, bottom=100, left=70, right=70)
            else:
                if row_idx % 2 == 1:
                    set_cell_background(cell, "FAFCFD")
                set_cell_margins(cell, top=55, bottom=55, left=70, right=70)
                
            p = cell.paragraphs[0]
            p.paragraph_format.line_spacing = 1.15
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(0)
            
            run = p.add_run(text)
            run.font.name = "Calibri"
            if is_header:
                run.font.size = Pt(8.5)
                run.font.bold = True
                run.font.color.rgb = COLOR_PRIMARY
            else:
                run.font.size = Pt(8.0)
                run.font.color.rgb = COLOR_BODY

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    doc.save(output_path)
    print(f"Successfully generated Figure-Ordered Methods document: {output_path}")

if __name__ == "__main__":
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    out1 = os.path.join(repo_root, "manuscript_methods_mouse_spaceflight_transcriptome.docx")
    out2 = os.path.join(repo_root, "final_figures_acceptance_20260906/05_Statistical_Tables_and_Audits/manuscript_methods_mouse_spaceflight_transcriptome.docx")
    build_docx(out1)
    build_docx(out2)
