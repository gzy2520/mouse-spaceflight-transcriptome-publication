#!/usr/bin/env python3
"""Generate the publication Methods Word document (.docx) for the mouse spaceflight transcriptome study."""

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

def set_cell_margins(cell, top=100, bottom=100, left=150, right=150):
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
    set_cell_margins(cell, top=140, bottom=140, left=200, right=200)
    
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
    run_text.font.size = Pt(10)
    run_text.font.italic = True
    run_text.font.color.rgb = RGBColor(0x33, 0x33, 0x33)
    
    p_after = doc.add_paragraph()
    p_after.paragraph_format.space_before = Pt(0)
    p_after.paragraph_format.space_after = Pt(4)

def build_docx(output_path):
    doc = Document()
    
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
        r_ftr = p_ftr.add_run("Methods Section — NASA OSDR Multi-Tissue Spaceflight Dataset")
        r_ftr.font.name = "Arial"
        r_ftr.font.size = Pt(8.5)
        r_ftr.font.color.rgb = RGBColor(0x88, 0x88, 0x88)
        
    COLOR_PRIMARY = RGBColor(0x1F, 0x4E, 0x78)
    COLOR_SECONDARY = RGBColor(0x2E, 0x75, 0xB6)
    COLOR_BODY = RGBColor(0x26, 0x26, 0x26)
    
    p_title = doc.add_paragraph()
    p_title.paragraph_format.space_before = Pt(0)
    p_title.paragraph_format.space_after = Pt(4)
    r_title = p_title.add_run("Methods")
    r_title.font.name = "Arial"
    r_title.font.size = Pt(22)
    r_title.font.bold = True
    r_title.font.color.rgb = COLOR_PRIMARY
    
    p_sub = doc.add_paragraph()
    p_sub.paragraph_format.space_before = Pt(0)
    p_sub.paragraph_format.space_after = Pt(12)
    r_sub = p_sub.add_run("Multi-Tissue Murine Spaceflight Transcriptomic Atlas and DNA Damage Repair Pathway Profiling")
    r_sub.font.name = "Arial"
    r_sub.font.size = Pt(13)
    r_sub.font.bold = True
    r_sub.font.color.rgb = COLOR_SECONDARY
    
    add_callout(
        doc,
        "This Methods document is compiled with strict fidelity to the audited publication analysis pipeline, "
        "frozen multi-omics matrices, and production figures (Main Figures 1b–6, Supplementary Figures S1–S6, "
        "and the 7-pathway linear boxplots). Analytical governance enforces immutable Ensembl Gene IDs (GRCm38/GRCm39) "
        "as computational keys, R/Python reproducibility, and deterministic random seed locking at seed = 25.",
        "GOVERNANCE & REPRODUCIBILITY CONTRACT"
    )
    
    def add_sec_heading(text):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(14)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.keep_with_next = True
        r = p.add_run(text)
        r.font.name = "Arial"
        r.font.size = Pt(12.5)
        r.font.bold = True
        r.font.color.rgb = COLOR_PRIMARY
        return p

    def add_body_p(text, bold_prefix=None):
        p = doc.add_paragraph()
        p.paragraph_format.line_spacing = 1.25
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(5)
        if bold_prefix:
            r_pre = p.add_run(bold_prefix)
            r_pre.font.name = "Calibri"
            r_pre.font.size = Pt(10.5)
            r_pre.font.bold = True
            r_pre.font.color.rgb = COLOR_BODY
        r = p.add_run(text)
        r.font.name = "Calibri"
        r.font.size = Pt(10.5)
        r.font.color.rgb = COLOR_BODY
        return p

    # Section 1
    add_sec_heading("1. Spaceflight Transcriptome Cohort Compilation and Metadata Curation")
    add_body_p(
        "To establish an organism-wide, multi-tissue atlas of spaceflight-induced transcriptional adaptations and DNA damage repair "
        "alterations, transcriptomic profiles and experimental metadata were curated from the NASA Open Science Data Repository (OSDR, "
        "formerly GeneLab; https://osdr.nasa.gov/). The analytical cohort encompasses 48 independent transcriptomic accessions (OSD datasets) "
        "across 59 primary analysis units, comprising 761 murine RNA-seq biospecimens spanning 26 distinct anatomical tissues (Fig. 1b, Fig. 6).",
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
        "Experimental design and animal allocation across 26 tissues were visualized using a stratified bubble plot (Fig. 1b, Fig. 6). "
        "The horizontal axis follows the canonical cross-tissue DNA damage response (DDR) activation rank. Vertically, each mission cluster "
        "is bifurcated into paired female and male tracks, with point area scaling with unique donor mouse counts. Age brackets are mapped to a "
        "fixed chronological 15-color palette (yellow -> green -> blue -> purple) without continuous interpolation.",
        bold_prefix="Metadata Bubble Visualization: "
    )

    # Section 2
    add_sec_heading("2. Gene Annotation Governance and Stable Identifier Mapping")
    add_body_p(
        "To eliminate gene symbol instability, historical alias drift, and erroneous cross-version joins, all computational analyses strictly "
        "adopted immutable Ensembl Gene Identifiers (Mus musculus Ensembl release GRCm38/GRCm39, ENSMUSG*) as primary keys. Official Mouse Genome "
        "Informatics (MGI) symbols and Entrez IDs were maintained strictly as display annotations in final visualization layers.",
        bold_prefix="Immutable Identifier Governance: "
    )
    add_body_p(
        "A foundational catalog of 74 curated core components across seven DNA repair pathways was established (Sheet DSB: 29 components spanning "
        "NHEJ, HR, and A-EJ; Sheet SSB: 45 components spanning BER, NER, MMR, and FA). Data completeness and sequencing depth audits across all flight "
        "biospecimens were conducted. Neil2 (ENSMUSG00000035121) exhibited recurrent sequencing dropout and zero-count sparsity across multiple flight "
        "tissue libraries, and was therefore systematically excluded from the final 45-component SSB analytical matrix and documented in a dedicated "
        "quality exclusion audit table (02_SSB_component_stable_id_mapping_without_Neil2.csv).",
        bold_prefix="Core Component Curation and Quality Exclusion: "
    )
    add_body_p(
        "In parallel, full pathway repertoire coverage was characterized across 335 total cataloged DNA repair genes: Base Excision Repair (BER, n = 54), "
        "Nucleotide Excision Repair (NER, n = 53), Mismatch Repair (MMR, n = 27), Fanconi Anemia (FA, n = 39), Homologous Recombination (HR, n = 135), "
        "Alternative End Joining (A-EJ, n = 9), and Non-Homologous End Joining (NHEJ, n = 18). Per-sample gene coverage was tracked to confirm robust "
        "library representation across all 761 source samples (Fig. S5).",
        bold_prefix="Seven-Pathway Repertoire Benchmarking: "
    )

    # Section 3
    add_sec_heading("3. Cross-Tissue Quantile Smoothing Normalization (YARN / qsmooth)")
    add_body_p(
        "To mitigate technical batch artifacts, library depth disparities, and platform variances across the 48 independent OSDR accessions while "
        "rigorously preserving bona fide tissue-specific transcriptional architectures, cross-tissue expression normalization was conducted using "
        "the YARN framework implementing tissue-aware quantile smoothing (qsmooth; Bioconductor package yarn / qsmooth).",
        bold_prefix="Algorithmic Framework: "
    )
    add_body_p(
        "Conventional quantile normalization forcibly enforces identical empirical distribution shapes across all libraries, which erroneously eliminates "
        "biological differences when profiling highly heterogeneous organ systems (e.g., spleen vs brain). In contrast, qsmooth computes quantile-specific "
        "weights w_q based on the ratio of between-tissue variance to total (between- plus within-tissue) variance across all quantiles q ∈ [0, 1]:\n"
        "        w_q = Var_between(q) / [Var_between(q) + Var_within(q)]\n"
        "At quantiles where within-tissue technical noise dominates (w_q → 0), values are normalized toward the global cross-tissue empirical distribution. "
        "Conversely, at quantiles where authentic between-tissue biological divergence dominates (w_q → 1), tissue-specific quantile distributions are preserved.",
        bold_prefix="Quantile Weight Formulation: "
    )
    add_body_p(
        "Normalized counts were log2-transformed with a unit pseudocount: YARNNormalizedLog2 = log2(Normalized Count + 1). For linear-scale abundance "
        "evaluations and publication boxplots, values were back-transformed to the linear expression scale via: Expression_linear = 2^(YARNNormalizedLog2) - 1.",
        bold_prefix="Logarithmic and Linear Scale Transformation: "
    )

    # Section 4
    add_sec_heading("4. Differential Expression Analysis and Mission-Level Effect Summarization")
    add_body_p(
        "Differential expression between spaceflight mice and matched ground controls was estimated within each OSDR accession using negative binomial "
        "generalized linear models in DESeq2. Contrast orientations were audited via orientation_multiplier flags to guarantee uniform directional "
        "consistency (Spaceflight vs. Ground).",
        bold_prefix="Negative Binomial Modeling: "
    )
    add_body_p(
        "To capture individual flight animal variability relative to baseline, sample-level ground-referenced log2 fold changes were computed. "
        "For each flight mouse k within a given accession and tissue:\n"
        "        Sample_log2FC_k = log2(Expr_flight_k + 1) - <log2(Expr_ground + 1)>_accession, tissue\n"
        "where <·> represents the arithmetic mean log2 baseline of all matched ground-control mice from the identical study accession and tissue.",
        bold_prefix="Sample-Level Ground-Referenced log2FC: "
    )
    add_body_p(
        "To synthesize flight impact across heterogeneous missions without distortion by statistical outliers, mission-level median log2FC values were "
        "computed, and the cross-mission consensus was determined as the median of mission medians (Fig. 4b, Fig. 5b). Statistical meta-significance across "
        "m missions was evaluated using signed Stouffer's Z-trend integration:\n"
        "        Z_meta = [∑_{i=1}^m sign(log2FC_i) · Φ^(-1)(1 - P_i / 2)] / √m\n"
        "where Φ^(-1) is the standard normal inverse cumulative distribution function. The resulting Z_meta and two-sided nominal P-values capture both the "
        "magnitude and directional consistency of spaceflight-induced transcriptional shifts.",
        bold_prefix="Multi-Mission Meta-Aggregation and Signed Stouffer Test: "
    )

    # Section 5
    add_sec_heading("5. Gene Set Enrichment Analysis (GSEA) and Gene Ontology Pathway Profiling")
    add_body_p(
        "Functional pathway enrichment was evaluated via pre-ranked Gene Set Enrichment Analysis using fgseaMultilevel (fgsea package; R 4.4.3), "
        "with pseudo-random seeds locked at seed = 25. Genes within each comparison were ranked by signed differential Wald statistics: "
        "S = sign(log2FC) · (-log10 P_nominal).",
        bold_prefix="Pre-Ranked fgsea Multilevel Modeling: "
    )
    add_body_p(
        "Analyses queried the MSigDB Hallmark collection (2026.1.Mm, Fig. S1) and 15 curated Gene Ontology Biological Process (GO:BP) terms (Fig. 1c). "
        "The GO panel incorporates the core root term for DNA damage response (GO:0006974), sub-pathway branches (DNA repair GO:0006281; DNA damage signal "
        "transduction GO:0042770; Intrinsic apoptotic signaling GO:0008630; DNA damage tolerance GO:0006301; Telomere maintenance in response to DNA damage "
        "GO:0043247; Chromosome telomeric region GO:0000781), physiological context terms (Cell cycle GO:0007049; DNA replication GO:0006260; "
        "Transcription GO:0006351; Translation GO:0006412; Mitochondrion GO:0005739; Cytoskeleton GO:0005856), and comparative mechanosensory processes "
        "(Sensory perception of mechanical stimulus GO:0050954; Cell death GO:0008219).",
        bold_prefix="Curated Pathway Repertoire: "
    )
    add_body_p(
        "To prevent large missions from exerting disproportionate leverage, Normalized Enrichment Scores (NES) were aggregated via unweighted "
        "arithmetic averaging across contributing missions within each tissue. Empirical significance was verified via exact sign-flip two-sided "
        "permutation testing across all 2^m sign configurations, coupled with Stouffer meta-Z transformation and Benjamini-Hochberg False Discovery "
        "Rate (FDR) multiple testing correction.",
        bold_prefix="Mission-Equal NES and Exact Permutation Statistics: "
    )

    # Section 6
    add_sec_heading("6. Sample-Level Cross-Tissue Pathway Co-Regulation and Fisher-z Integration")
    add_body_p(
        "To quantify coordinated pathway co-regulation across the 15 GO biological processes, we implemented a sample-level correlation framework "
        "rather than relying on aggregated tissue means (build_log2fc_per_tissue_pathway_spearman_20260902.py, Fig. 2, Fig. S6).",
        bold_prefix="Analytical Rationale: "
    )
    add_body_p(
        "1. Intra-Tissue Correlation: Within each of the 26 tissues, taking individual flight mice (n = 360) as independent observation units, "
        "the sample-level pathway score was defined as the mean member-gene log2FC. Pairwise non-parametric Spearman rank correlation coefficients "
        "(ρ_ij) were computed across the 15 GO terms, generating 26 distinct 15 × 15 correlation matrices (Fig. S6).\n"
        "2. Fisher-z Meta-Aggregation: To synthesize an organism-wide pathway co-regulation matrix (Fig. 2) without sample-size distortion, "
        "intra-tissue correlation coefficients were transformed into normally distributed z-scores via Fisher's z-transformation:\n"
        "        z = arctanh(ρ) = 0.5 · ln[(1 + ρ) / (1 - ρ)]\n"
        "To guarantee numerical stability, correlations were bounded within [-1 + 10^(-12), 1 - 10^(-12)], and z-scores were capped at |z| ≤ 5.0.\n"
        "3. Cross-Tissue Pooling Criteria: An unweighted arithmetic mean z̄ was calculated across the 23 eligible tissues satisfying the sample-size "
        "threshold (n_flight ≥ 5; Colon [n = 4], Heart [n = 3], and Femoral skin [n = 2] were retained in individual panels but excluded from cross-tissue "
        "pooling). The pooled z̄ was back-transformed to derive the overall meta-analytic Spearman rank correlation:\n"
        "        ρ̄ = tanh(z̄) = [exp(2z̄) - 1] / [exp(2z̄) + 1]",
        bold_prefix="Mathematical Formulation: "
    )

    # Section 7
    add_sec_heading("7. Intersection Repertoires, UpSet Matrices, and Common-Direction DEG Prioritization")
    add_body_p(
        "To examine shared and tissue-specific responsive elements within the core DDR machinery (GO:0006974), overlapping DEG repertoires "
        "were evaluated across multi-dataset organs, specifically Thymus (4 flight comparisons, Fig. 3a) and Kidney (5 independent OSDR datasets, "
        "Fig. S2a) using ggVennDiagram. Organ-wide intersection topology across NES > 1 tissues (excluding Kidney, Fig. S3a) and NES < -1 tissues "
        "(excluding Lung and Thymus, Fig. S4a) was resolved using UpSet intersection geometry.",
        bold_prefix="Repertoire Overlap and UpSet Topology: "
    )
    add_body_p(
        "To prioritize core genes demonstrating conserved responsiveness independent of mission idiosyncrasies, common-direction filtering "
        "was executed (plot_common_direction_heatmaps.py; Fig. 3b, Fig. S2b, Fig. S3b, Fig. S4b). Candidate DEGs were required to demonstrate uniform "
        "directional concordance (all log2FC > 0 or all log2FC < 0) across all parallel dataset comparisons at nominal P < 0.05 (or P < 0.10 for global "
        "NES-stratified cohorts). Qualifying genes were ranked by their minimum absolute fold change across comparisons:\n"
        "        Score_g = min_c |log2FC_{g, c}|\n"
        "The top 30 ranked common-direction genes were displayed as heatmaps with standardized color gradients and exact log2FC text overlays.",
        bold_prefix="Common-Direction Screening Formulation: "
    )

    # Section 8
    add_sec_heading("8. Hierarchical Clustering and Topological Dendrogram Analysis")
    add_body_p(
        "To establish baseline tissue groupings based on DNA repair pathway expression, the 74 core components were partitioned into two functional "
        "modules: the Double-Strand Break (DSB) module (29 components across NHEJ, HR, and A-EJ; Fig. 4a) and the Single-Strand Break (SSB) module "
        "(45 components across BER, NER, MMR, and FA; Fig. 5a).",
        bold_prefix="Functional Module Partitioning: "
    )
    add_body_p(
        "Pairwise topological dissimilarity between tissues Ta and Tb was defined on normalized baseline YARN qsmooth expression matrices using "
        "Spearman rank distance: D(Ta, Tb) = 1 - ρ_Spearman(Ta, Tb). Hierarchical agglomerative clustering was executed using average linkage "
        "(hclust(as.dist(D), method = 'average')). The resulting dendrogram branch topologies were preserved to structure heatmap rows, while individual "
        "matrix cells were rendered as clean, high-contrast color blocks with English anatomical labels.",
        bold_prefix="Clustering Algorithm and Topology Governance: "
    )

    # Section 9
    add_sec_heading("9. Tissue-Stratified Linear-Scale Boxplots and One-Way ANOVA Modeling")
    add_body_p(
        "To evaluate physiological baseline expression and organ-level heterogeneity across the mammalian DNA repair machinery, 26 key biological "
        "marker genes across all seven canonical pathways were selected based on mechanistic centrality and manuscript narrative:\n"
        "• Nucleotide Excision Repair (NER, 7 genes): Xpc (ENSMUSG00000030094), Rad23b (ENSMUSG00000028426), Cetn2 (ENSMUSG00000031347), "
        "Ddb1 (ENSMUSG00000024740), Ddb2 (ENSMUSG00000002109), Ercc6 (ENSMUSG00000054051), Ercc8 (ENSMUSG00000021694).\n"
        "• Homologous Recombination (HR, 4 genes): Brca1 (ENSMUSG00000017146), Bard1 (ENSMUSG00000026196), Blm (ENSMUSG00000030528), Rad51 (ENSMUSG00000027323).\n"
        "• Non-Homologous End Joining (NHEJ, 4 genes): Nhej1 / XLF (ENSMUSG00000026162), Paxx (ENSMUSG00000047617), Xrcc6 / Ku70 (ENSMUSG00000022471), Xrcc5 / Ku80 (ENSMUSG00000026187).\n"
        "• Alternative End Joining (A-EJ, 4 genes): Parp1 (ENSMUSG00000026496), Polq (ENSMUSG00000034206), Lig1 (ENSMUSG00000056394), Lig3 (ENSMUSG00000020697).\n"
        "• Base Excision Repair (BER, 3 genes): Ung (ENSMUSG00000029591), Ogg1 (ENSMUSG00000030271), Neil1 (ENSMUSG00000032298).\n"
        "• Mismatch Repair (MMR, 2 genes): Msh2 (ENSMUSG00000024151), Msh3 (ENSMUSG00000014850).\n"
        "• Fanconi Anemia (FA, 2 genes): Fancd2 (ENSMUSG00000034023), Fanci (ENSMUSG00000039187).",
        bold_prefix="Marker Gene Selection and Immutable Identifier Registry: "
    )
    add_body_p(
        "Individual flight mice biospecimens (n = 360) were stratified across 26 canonical tissues (arranged along the horizontal axis in canonical "
        "DDR rank order). Linear expression values were derived as Expression_linear = 2^(YARNNormalizedLog2) - 1.\n"
        "• Boxplot Geometry: The lower hinge indicates the first quartile (Q1, 25%); the center solid line represents the sample mean (μ); the upper "
        "hinge indicates the third quartile (Q3, 75%). Whiskers extend to the furthest observations within 1.5 × IQR from the hinges.\n"
        "• Individual Biospecimen Scatter: All 360 individual flight mice were overlaid as jittered points using position_jitterdodge "
        "(dodge.width = 0.80, jitter.width = 0.11, seed = 25).\n"
        "• Cross-Tissue Mean Trajectories: Dashed lines connect tissue-level mean values for each gene to illustrate organ-specific expression gradients.\n"
        "• NER Subcomplex Color and Shape Coding: NER recognition complexes were disambiguated via paired color-and-shape glyphs: XPC/RAD23B/CETN2 "
        "(amber; cross ✕, triangle ▲, circle ●); DDB1/DDB2 (teal; cross ✕, triangle ▲); and ERCC6/ERCC8 (blue; cross ✕, triangle ▲).",
        bold_prefix="Linear Boxplot Structural Specifications: "
    )
    add_body_p(
        "Tissue-dependent expression heterogeneity was formally evaluated using ordinary One-way Analysis of Variance (ANOVA):\n"
        "1. Gene-Level One-Way ANOVA: Evaluated cross-tissue variation for each gene independently:\n"
        "        y_ijk = μ + α_j + ε_ijk,    ε_ijk ~ N(0, σ^2)\n"
        "yielding between-group degrees of freedom df_group = 25 and residual degrees of freedom df_residual = 334 (or 330 for FA), with exact "
        "F-statistics and scientific-notation P-values printed directly in figure subtitles.\n"
        "2. Pathway-Level One-Way ANOVA: Evaluated composite pathway heterogeneity by modeling the per-mouse mean pathway expression:\n"
        "        ȳ_ik = (1 / G) ∑_{g=1}^G y_igk\n"
        "testing overall pathway abundance divergence across the 26 tissues. Exact F-statistics and P-values are reported in subtitles and audited "
        "in comprehensive CSV summary tables (01_pathway_one_way_anova_linear_summary.csv, 02_gene_one_way_anova_linear_summary.csv).",
        bold_prefix="One-Way ANOVA Statistical Modeling: "
    )

    # Section 10
    add_sec_heading("10. Computational Environment, Code Governance, and Reproducibility")
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

    # Section 11: Summary Table
    add_sec_heading("11. Production Figure and Analytical Workflow Reference Matrix")
    add_body_p(
        "The following audited reference matrix connects every publication figure and statistical table to its underlying mathematical model, "
        "sample scope, and analytical script:",
        bold_prefix="Provenance Cross-Reference: "
    )
    
    table_data = [
        ["Figure / Table", "Analytical Scope", "Input Cohort / Sample Size", "Statistical Model / Algorithm", "Source Script"],
        ["Fig. 1b / Fig. 6", "Experimental Design Atlas", "761 samples, 48 accessions, 26 tissues", "Stratified categorical bubble matrix (15 discrete age colors)", "R/figures/08_mouse_metadata_bubble.R"],
        ["Fig. 1c", "Cross-Tissue GO Profiling", "26 tissues × 15 GO biological terms", "Mission-equal unweighted mean NES, exact sign-flip P", "R/figures/01_go_overview.R"],
        ["Fig. 2", "Cross-Tissue Co-Regulation", "23 eligible tissues (n_flight ≥ 5)", "Sample-level log2FC, intra-tissue Spearman ρ, Fisher-z pooling", "python/plot_log2fc_per_tissue_pathway_spearman_20260902.py"],
        ["Fig. 3a / S2a", "Tissue Repertoire Overlap", "Thymus (4 comp.) / Kidney (5 OSD)", "Set intersection membership, ggVennDiagram", "R/figures/03_kidney_thymus_venn.R"],
        ["Fig. 3b / S2b", "Common-Direction DEGs", "Thymus / Kidney flight comparisons", "Concordant directionality, ranked by min |log2FC|", "python/plot_common_direction_heatmaps.py"],
        ["Fig. 4a / 5a", "Repair Baseline Clustering", "26 tissues, 74 components (DSB 29, SSB 45)", "YARN qsmooth baseline, 1 - Spearman ρ, average linkage", "R/figures/04_qsmooth_tree_heatmaps.R"],
        ["Fig. 4b / 5b", "Spaceflight Meta-Analysis", "26 tissues × core components", "Median of mission medians log2FC, signed Stouffer meta-Z", "R/figures/05_meta_log2fc_heatmaps.R"],
        ["Fig. S1", "Hallmark GSEA Overview", "13 missions with complete coverage", "Top 12 global Hallmark terms (6 positive, 6 negative NES)", "R/figures/02_hallmark_gsea.R"],
        ["Fig. S3a / S4a", "Multi-Tissue UpSet Repertoire", "NES > 1 tissues / NES < -1 tissues", "Binary membership pattern matrix, UpSet intersections", "R/figures/06_go_upsets.R"],
        ["Fig. S3b / S4b", "Common-Direction DEGs", "NES > 1 / NES < -1 multi-tissue sets", "Directional concordance at P < 0.10, ranked by min |log2FC|", "python/plot_common_direction_heatmaps.py"],
        ["Fig. S5", "Pathway Gene Coverage", "761 samples × 7 pathways (335 genes)", "Per-sample member gene library detection tracking", "R/figures/07_seven_pathway_counts.R"],
        ["Fig. S6", "Per-Tissue Co-Regulation", "26 individual tissue matrices", "Sample-level 15 × 15 Spearman rank correlation matrices", "python/plot_log2fc_per_tissue_pathway_spearman_20260902.py"],
        ["Fig_A-EJ~NHEJ (Linear)", "Linear Pathway Boxplots", "360 flight mice × 26 marker genes", "Linear expression 2^log2 - 1, Mean hinge, Gene & Pathway ANOVA", "R/figures/10_pathway_gene_tissue_boxplots_linear.R"],
        ["Tables 01–04 (Linear)", "Statistical Audit Tables", "9,360 sample observations, 26 tissues", "Pathway ANOVA, Gene ANOVA, Tissue Means, Sample Linear Values", "results/pathway_gene_tissue_boxplots_linear/tables/"]
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
                set_cell_margins(cell, top=100, bottom=100, left=80, right=80)
            else:
                if row_idx % 2 == 1:
                    set_cell_background(cell, "FAFCFD")
                set_cell_margins(cell, top=60, bottom=60, left=80, right=80)
                
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
    print(f"Successfully generated Methods document: {output_path}")

if __name__ == "__main__":
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    out1 = os.path.join(repo_root, "manuscript_methods_mouse_spaceflight_transcriptome.docx")
    out2 = os.path.join(repo_root, "final_figures_acceptance_20260906/05_Statistical_Tables_and_Audits/manuscript_methods_mouse_spaceflight_transcriptome.docx")
    build_docx(out1)
    build_docx(out2)
