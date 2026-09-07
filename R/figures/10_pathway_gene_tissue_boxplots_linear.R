#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Script: R/figures/10_pathway_gene_tissue_boxplots_linear.R
# Purpose: Generate publication-grade tissue-stratified boxplots for 7 DNA repair
#          pathways on linear expression scale (2^log2 - 1) highlighting key
#          genes from manuscript text/images.
# Features:
#   - X-axis: 26 tissues (canonical DDR GO order)
#   - Y-axis: mouse-level YARN qsmooth linear expression (2^log2 - 1)
#   - Points: individual flight mice (n = 360)
#   - Boxplot: Q1, Mean (middle solid bar), Q3, 1.5 * IQR whiskers
#   - Dashed lines: connect tissue-level mean values for each gene
#   - One-way ANOVA: test tissue heterogeneity on linear expression
#   - Custom coloring: NER grouped (XPC/RAD23B/CETN2; DDB1/2; ERCC6/8);
#                     all other pathways have unique gene colors.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(25) # Personal symbol random seed

# Resolve repository root
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
} else {
  repo_root <- getwd()
}

out_dir <- file.path(repo_root, "results/pathway_gene_tissue_boxplots_linear")
fig_dir <- file.path(out_dir, "figures")
tbl_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

format_sci <- function(p) {
  if (is.na(p)) return("NA")
  if (p >= 0.01) return(sprintf("%.3f", p))
  sci_str <- sprintf("%.2e", p)
  parts <- strsplit(sci_str, "e")[[1]]
  base_part <- parts[1]
  exp_int <- as.integer(parts[2])
  digits <- c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3", "4" = "\u2074",
              "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078", "9" = "\u2079", "-" = "\u207b")
  exp_str <- as.character(exp_int)
  exp_unicode <- paste0(sapply(strsplit(exp_str, "")[[1]], function(ch) digits[[ch]]), collapse = "")
  sprintf("%s \u00d7 10%s", base_part, exp_unicode)
}

# Input data paths
qsmooth_path <- file.path(repo_root, "data/publication_input/qsmooth/06_component_flight_sample_yarn_qsmooth_values.csv.gz")
upstream_fallback <- "/Users/gzy2520/Desktop/graduate design/03_analysis_results/22_essential_components_yarn_qsmooth_no_Neil2_attached_row_dendrogram_20260828/tables/06_component_flight_sample_yarn_qsmooth_values.csv.gz"
if (!file.exists(qsmooth_path) && file.exists(upstream_fallback)) {
  qsmooth_path <- upstream_fallback
}
go_path <- file.path(repo_root, "data/publication_input/go/04_tissue_statistics_concrete_terms_and_context.csv")

if (!file.exists(qsmooth_path)) stop("Missing qsmooth sample table: ", qsmooth_path, call. = FALSE)
if (!file.exists(go_path)) stop("Missing GO tissue table: ", go_path, call. = FALSE)

dt <- fread(qsmooth_path, showProgress = FALSE)
go <- fread(go_path, showProgress = FALSE)

# Canonical 26-tissue order
tissue_order <- go[term_key == "DNA_damage_response_GO0006974" & term_order == 1]$analysis_tissue
stopifnot(length(tissue_order) == 26L, setequal(unique(dt$Group), tissue_order))

# Compute linear values: 2^log2 - 1
dt[, YARNLinear := 2^YARNNormalizedLog2 - 1]

# Define pathway gene definitions with stable Ensembl IDs, mapped to final publication figure panels
pathway_defs <- list(
  # Double-Strand Break (DSB) Repair -> Figure 4 (panels b to d)
  NHEJ = list(
    fig_panel = "Fig_4b",
    fig_letter = "b",
    fig_group = "DSB",
    title = "Non-Homologous End Joining (NHEJ)",
    caption_note = "Genes: NHEJ1 (XLF), PAXX, XRCC6 (Ku70), and XRCC5 (Ku80).",
    genes = c(
      Nhej1 = "ENSMUSG00000026162",
      Paxx = "ENSMUSG00000047617",
      Xrcc6 = "ENSMUSG00000022471",
      Xrcc5 = "ENSMUSG00000026187"
    ),
    colors = c(
      Nhej1 = "#E69F00",
      Paxx = "#56B4E9",
      Xrcc6 = "#009E73",
      Xrcc5 = "#D55E00"
    ),
    linetypes = c(
      Nhej1 = "dashed",
      Paxx = "dashed",
      Xrcc6 = "dashed",
      Xrcc5 = "dashed"
    )
  ),
  HR = list(
    fig_panel = "Fig_4c",
    fig_letter = "c",
    fig_group = "DSB",
    title = "Homologous Recombination (HR)",
    caption_note = "Genes: BRCA1, BARD1, BLM, and RAD51.",
    genes = c(
      Brca1 = "ENSMUSG00000017146",
      Bard1 = "ENSMUSG00000026196",
      Blm = "ENSMUSG00000030528",
      Rad51 = "ENSMUSG00000027323"
    ),
    colors = c(
      Brca1 = "#D55E00",
      Bard1 = "#E69F00",
      Blm = "#009E73",
      Rad51 = "#0072B2"
    ),
    linetypes = c(
      Brca1 = "dashed",
      Bard1 = "dashed",
      Blm = "dashed",
      Rad51 = "dashed"
    )
  ),
  "A-EJ" = list(
    fig_panel = "Fig_4d",
    fig_letter = "d",
    fig_group = "DSB",
    title = "Alternative End Joining (A-EJ / MMEJ)",
    caption_note = "Genes: PARP1, POLQ, LIG1, and LIG3.",
    genes = c(
      Parp1 = "ENSMUSG00000026496",
      Polq = "ENSMUSG00000034206",
      Lig1 = "ENSMUSG00000056394",
      Lig3 = "ENSMUSG00000020697"
    ),
    colors = c(
      Parp1 = "#E69F00",
      Polq = "#56B4E9",
      Lig1 = "#009E73",
      Lig3 = "#CC79A7"
    ),
    linetypes = c(
      Parp1 = "dashed",
      Polq = "dashed",
      Lig1 = "dashed",
      Lig3 = "dashed"
    )
  ),

  # Single-Strand Break (SSB) Repair -> Figure 5 (panels b to e)
  BER = list(
    fig_panel = "Fig_5b",
    fig_letter = "b",
    fig_group = "SSB",
    title = "Base Excision Repair (BER)",
    caption_note = "Genes: UNG (uracil sensor), OGG1 (8-oxoG glycosylase), and NEIL1 (oxidized base repair).",
    genes = c(
      Ung = "ENSMUSG00000029591",
      Ogg1 = "ENSMUSG00000030271",
      Neil1 = "ENSMUSG00000032298"
    ),
    colors = c(
      Ung = "#E69F00",
      Ogg1 = "#009E73",
      Neil1 = "#CC79A7"
    ),
    linetypes = c(
      Ung = "dashed",
      Ogg1 = "dashed",
      Neil1 = "dashed"
    )
  ),
  NER = list(
    fig_panel = "Fig_5c",
    fig_letter = "c",
    fig_group = "SSB",
    title = "Nucleotide Excision Repair (NER)",
    caption_note = "Shapes: Cross (\u2715), Triangle (\u25b2), Circle (\u25cf). Color coding: XPC/RAD23B/CETN2 (GG-NER surveillance; Amber), DDB1/2 (UV lesion detection; Teal), ERCC6/8 (TC-NER; Blue).",
    genes = c(
      Xpc = "ENSMUSG00000030094",
      Rad23b = "ENSMUSG00000028426",
      Cetn2 = "ENSMUSG00000031347",
      Ddb1 = "ENSMUSG00000024740",
      Ddb2 = "ENSMUSG00000002109",
      Ercc6 = "ENSMUSG00000054051",
      Ercc8 = "ENSMUSG00000021694"
    ),
    colors = c(
      Xpc = "#E69F00", Rad23b = "#E69F00", Cetn2 = "#E69F00",
      Ddb1 = "#009E73", Ddb2 = "#009E73",
      Ercc6 = "#0072B2", Ercc8 = "#0072B2"
    ),
    linetypes = c(
      Xpc = "dashed", Rad23b = "dotted", Cetn2 = "dotdash",
      Ddb1 = "dashed", Ddb2 = "dotted",
      Ercc6 = "dashed", Ercc8 = "dotted"
    ),
    shapes = c(
      Xpc = 4, Rad23b = 17, Cetn2 = 16,
      Ddb1 = 4, Ddb2 = 17,
      Ercc6 = 4, Ercc8 = 17
    ),
    labels = c(
      Xpc = "Xpc (\u2715)", Rad23b = "Rad23b (\u25b2)", Cetn2 = "Cetn2 (\u25cf)",
      Ddb1 = "Ddb1 (\u2715)", Ddb2 = "Ddb2 (\u25b2)",
      Ercc6 = "Ercc6 (\u2715)", Ercc8 = "Ercc8 (\u25b2)"
    )
  ),
  MMR = list(
    fig_panel = "Fig_5d",
    fig_letter = "d",
    fig_group = "SSB",
    title = "Mismatch Repair (MMR)",
    caption_note = "Genes: MSH2 and MSH3 (MutS\u03b2 complex).",
    genes = c(
      Msh2 = "ENSMUSG00000024151",
      Msh3 = "ENSMUSG00000014850"
    ),
    colors = c(
      Msh2 = "#D55E00",
      Msh3 = "#0072B2"
    ),
    linetypes = c(
      Msh2 = "dashed",
      Msh3 = "dashed"
    )
  ),
  FA = list(
    fig_panel = "Fig_5e",
    fig_letter = "e",
    fig_group = "SSB",
    title = "Fanconi Anemia (FA) Pathway",
    caption_note = "Genes: FANCD2 and FANCI (central ID complex essential for interstrand crosslink repair).",
    genes = c(
      Fancd2 = "ENSMUSG00000034023",
      Fanci = "ENSMUSG00000039187"
    ),
    colors = c(
      Fancd2 = "#E69F00",
      Fanci = "#56B4E9"
    ),
    linetypes = c(
      Fancd2 = "dashed",
      Fanci = "dashed"
    )
  )
)

# Background alternating bands data for 26 tissues
bands <- data.table(
  xmin = seq(0.5, 25.5, by = 2),
  xmax = seq(1.5, 26.5, by = 2),
  ymin = -Inf,
  ymax = Inf
)

anova_gene_list <- list()
anova_pw_list <- list()
summary_mean_list <- list()
used_data_list <- list()

for (pw_name in names(pathway_defs)) {
  pw_info <- pathway_defs[[pw_name]]
  gene_ids <- pw_info$genes
  gene_symbols <- names(gene_ids)
  
  sub_dt <- dt[EnsemblID %in% gene_ids]
  sub_dt[, Gene := factor(OfficialMouseSymbol, levels = gene_symbols)]
  sub_dt[, Group := factor(Group, levels = tissue_order)]
  sub_dt[, Pathway := pw_name]
  
  used_data_list[[pw_name]] <- sub_dt[, .(Pathway = pw_name, EnsemblID, OfficialMouseSymbol, SampleID, Group, accession, mission_cluster, YARNNormalizedLog2, YARNLinear)]
  
  # Compute boxplot statistics with mean as the center line on linear values
  box_stats <- sub_dt[, {
    valid_vals <- YARNLinear[!is.na(YARNLinear)]
    q <- quantile(valid_vals, c(0.25, 0.75))
    iqr <- q[2] - q[1]
    m <- mean(valid_vals)
    sd_val <- sd(valid_vals)
    n_mice <- length(valid_vals)
    ymin <- max(min(valid_vals), q[1] - 1.5 * iqr)
    ymax <- min(max(valid_vals), q[2] + 1.5 * iqr)
    .(ymin = ymin, lower = q[1], middle = m, upper = q[2], ymax = ymax, SD = sd_val, N = n_mice)
  }, by = .(Group, Gene)]
  
  summary_mean_list[[pw_name]] <- box_stats[, .(Pathway = pw_name, Group, Gene, Mean = middle, SD, Q1 = lower, Q3 = upper, Ymin = ymin, Ymax = ymax, N)]
  
  # One-way ANOVA per gene on linear values
  gene_p_strs <- character()
  for (sym in gene_symbols) {
    g_data <- sub_dt[Gene == sym & !is.na(YARNLinear)]
    fit_g <- aov(YARNLinear ~ Group, data = g_data)
    s_g <- summary(fit_g)[[1]]
    f_g <- s_g["Group", "F value"]
    p_g <- s_g["Group", "Pr(>F)"]
    df_g <- s_g["Group", "Df"]
    df_r <- s_g["Residuals", "Df"]
    sig_g <- fifelse(p_g < 0.001, "***", fifelse(p_g < 0.01, "**", fifelse(p_g < 0.05, "*", "ns")))
    gene_p_strs <- c(gene_p_strs, sprintf("%s (P = %s)", sym, format_sci(p_g)))
    
    anova_gene_list[[length(anova_gene_list) + 1]] <- data.table(
      Pathway = pw_name, EnsemblID = gene_ids[[sym]], Symbol = sym,
      Df_group = df_g, Df_residual = df_r, F_value = f_g, P_value = p_g, Significance = sig_g
    )
  }
  
  # Overall Pathway One-way ANOVA on linear values
  pw_sample <- sub_dt[!is.na(YARNLinear), .(MeanExpr = mean(YARNLinear)), by = .(Group, SampleID)]
  fit_pw <- aov(MeanExpr ~ Group, data = pw_sample)
  s_pw <- summary(fit_pw)[[1]]
  f_pw <- s_pw["Group", "F value"]
  p_pw <- s_pw["Group", "Pr(>F)"]
  df_pw <- s_pw["Group", "Df"]
  df_pwr <- s_pw["Residuals", "Df"]
  sig_pw <- fifelse(p_pw < 0.001, "***", fifelse(p_pw < 0.01, "**", fifelse(p_pw < 0.05, "*", "ns")))
  
  anova_pw_list[[length(anova_pw_list) + 1]] <- data.table(
    Pathway = pw_name, Df_group = df_pw, Df_residual = df_pwr,
    F_value = f_pw, P_value = p_pw, Significance = sig_pw
  )
  
  # Layout dodge width
  n_genes <- length(gene_symbols)
  dodge_w <- if (n_genes >= 6) 0.82 else if (n_genes >= 4) 0.78 else 0.70
  box_w <- if (n_genes >= 6) 0.74 else if (n_genes >= 4) 0.70 else 0.65
  dodge <- position_dodge(width = dodge_w)
  
  sub_dt_points <- sub_dt[!is.na(YARNLinear)]
  
  # Assign shapes and labels for unified point legend across all pathways
  pw_shapes <- if (!is.null(pw_info$shapes)) pw_info$shapes else setNames(rep(16, n_genes), gene_symbols)
  pw_labels <- if (!is.null(pw_info$labels)) pw_info$labels else setNames(gene_symbols, gene_symbols)
  
  # Plot assembly
  p <- ggplot() +
    geom_rect(
      data = bands,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#F6F8FA", inherit.aes = FALSE
    ) +
    geom_vline(xintercept = seq(1.5, 25.5, by = 1), colour = "#E5E9EE", linewidth = 0.35) +
    geom_boxplot(
      data = box_stats,
      aes(x = Group, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = Gene, color = Gene),
      stat = "identity", position = dodge, width = box_w, alpha = 0.38, linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_point(
      data = sub_dt_points,
      aes(x = Group, y = YARNLinear, color = Gene, shape = Gene),
      position = position_jitterdodge(jitter.width = 0.11, dodge.width = dodge_w, seed = 25),
      size = if (!is.null(pw_info$shapes)) 1.45 else 1.20,
      alpha = 0.65,
      show.legend = TRUE
    ) +
    geom_line(
      data = box_stats,
      aes(x = Group, y = middle, group = Gene, color = Gene, linetype = Gene),
      position = dodge, linewidth = 0.82,
      show.legend = FALSE
    )
  
  if (!is.null(pw_info$shapes)) {
    p <- p +
      geom_point(
        data = box_stats,
        aes(x = Group, y = middle, color = Gene, shape = Gene),
        position = dodge, size = 3.0, stroke = 1.0,
        show.legend = FALSE
      )
  }
  
  p <- p +
    scale_fill_manual(values = pw_info$colors, guide = "none") +
    scale_color_manual(values = pw_info$colors, labels = pw_labels) +
    scale_shape_manual(values = pw_shapes, labels = pw_labels) +
    scale_linetype_manual(values = pw_info$linetypes, guide = "none") +
    guides(
      color = guide_legend(
        title = "Gene",
        override.aes = list(
          shape = pw_shapes,
          size = 5.0,
          stroke = 1.2,
          alpha = 1,
          linetype = 0
        ),
        nrow = if (n_genes >= 6) 2 else 1
      ),
      shape = "none"
    ) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      title = paste0(pw_info$title, " Expression across 26 Tissues (Linear Scale)"),
      subtitle = sprintf(
        "One-way ANOVA (across 26 tissues): Pathway F(%d, %d) = %.2f, P = %s\nIndividual gene ANOVA: %s",
        df_pw, df_pwr, f_pw, format_sci(p_pw), paste(gene_p_strs, collapse = "  |  ")
      ),
      caption = paste0(
        "Boxes represent Q1, Mean (middle solid bar), and Q3; whiskers extend to 1.5 * IQR. Points represent individual flight mice (n = 360).\n",
        "Dashed lines connect tissue-level mean values for each gene. Values shown on linear scale (2^log2 - 1). ", pw_info$caption_note
      ),
      x = NULL,
      y = "YARN qsmooth normalized expression (linear scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#E8ECEF", linewidth = 0.5),
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, face = "bold", colour = "#222222", size = 11.5),
      axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
      plot.title = element_text(face = "bold", size = 17, colour = "#1A202C"),
      plot.subtitle = element_text(colour = "#334155", size = 11.5, lineheight = 1.3, margin = margin(b = 8)),
      plot.caption = element_text(colour = "#718096", size = 10, lineheight = 1.25, margin = margin(t = 10)),
      legend.position = "top",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11.5),
      plot.margin = margin(12, 16, 12, 16)
    )
  
  # File export (PNG and PDF)
  safe_name <- gsub("[^A-Za-z0-9_]", "-", pw_name)
  panel_id <- pw_info$fig_panel
  
  # 1. Output to results figures
  png_file <- file.path(fig_dir, sprintf("%s.png", panel_id))
  pdf_file <- file.path(fig_dir, sprintf("%s.pdf", panel_id))
  ggsave(png_file, p, width = 18, height = 8.5, dpi = 300, bg = "white")
  ggsave(pdf_file, p, width = 18, height = 8.5, device = grDevices::cairo_pdf, bg = "white")
  
  # Also write descriptive names in results
  file.copy(png_file, file.path(fig_dir, sprintf("%s_%s_tissue_boxplot_linear.png", panel_id, safe_name)), overwrite = TRUE)
  file.copy(pdf_file, file.path(fig_dir, sprintf("%s_%s_tissue_boxplot_linear.pdf", panel_id, safe_name)), overwrite = TRUE)
  file.copy(png_file, file.path(fig_dir, sprintf("Fig_%s_tissue_boxplot_linear.png", safe_name)), overwrite = TRUE)
  file.copy(pdf_file, file.path(fig_dir, sprintf("Fig_%s_tissue_boxplot_linear.pdf", safe_name)), overwrite = TRUE)
  
  # 2. Sync to final_figures_acceptance_20260906/01_Main_Figures
  accept_main_dir <- file.path(repo_root, "final_figures_acceptance_20260906/01_Main_Figures")
  if (dir.exists(accept_main_dir)) {
    file.copy(png_file, file.path(accept_main_dir, sprintf("%s.png", panel_id)), overwrite = TRUE)
    file.copy(pdf_file, file.path(accept_main_dir, sprintf("%s.pdf", panel_id)), overwrite = TRUE)
  }
  
  # 3. Sync to final_figures_acceptance_20260906/03_Pathway_Tissue_Boxplots_Linear
  accept_linear_dir <- file.path(repo_root, "final_figures_acceptance_20260906/03_Pathway_Tissue_Boxplots_Linear")
  if (dir.exists(accept_linear_dir)) {
    file.copy(png_file, file.path(accept_linear_dir, sprintf("%s.png", panel_id)), overwrite = TRUE)
    file.copy(pdf_file, file.path(accept_linear_dir, sprintf("%s.pdf", panel_id)), overwrite = TRUE)
    file.copy(png_file, file.path(accept_linear_dir, sprintf("%s_%s_tissue_boxplot_linear.png", panel_id, safe_name)), overwrite = TRUE)
    file.copy(pdf_file, file.path(accept_linear_dir, sprintf("%s_%s_tissue_boxplot_linear.pdf", panel_id, safe_name)), overwrite = TRUE)
    file.copy(png_file, file.path(accept_linear_dir, sprintf("Fig_%s_tissue_boxplot_linear.png", safe_name)), overwrite = TRUE)
    file.copy(pdf_file, file.path(accept_linear_dir, sprintf("Fig_%s_tissue_boxplot_linear.pdf", safe_name)), overwrite = TRUE)
  }
  
  # 4. Sync to final_result_integrated_20260905/Main
  integrated_main_dir <- file.path(repo_root, "final_result_integrated_20260905/Main")
  if (dir.exists(integrated_main_dir)) {
    file.copy(png_file, file.path(integrated_main_dir, sprintf("%s.png", panel_id)), overwrite = TRUE)
    file.copy(pdf_file, file.path(integrated_main_dir, sprintf("%s.pdf", panel_id)), overwrite = TRUE)
  }
  
  message("Generated ", panel_id, " (", pw_name, "): ", png_file)
}

# Export summary tables
anova_gene_dt <- rbindlist(anova_gene_list)
anova_pw_dt <- rbindlist(anova_pw_list)
summary_mean_dt <- rbindlist(summary_mean_list)
used_data_dt <- rbindlist(used_data_list)

fwrite(anova_pw_dt, file.path(tbl_dir, "01_pathway_one_way_anova_linear_summary.csv"))
fwrite(anova_gene_dt, file.path(tbl_dir, "02_gene_one_way_anova_linear_summary.csv"))
fwrite(summary_mean_dt, file.path(tbl_dir, "03_tissue_gene_mean_linear_qsmooth_summary.csv"))
fwrite(used_data_dt, file.path(tbl_dir, "04_sample_linear_qsmooth_values_used.csv.gz"))

# Clean up prototype if present
proto_file <- file.path(fig_dir, "prototype_HR_linear.png")
if (file.exists(proto_file)) file.remove(proto_file)

message("All 7 linear pathway figures and tables successfully generated under: ", out_dir)
