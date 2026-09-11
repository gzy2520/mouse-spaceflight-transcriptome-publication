#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Script: R/figures/09_pathway_gene_tissue_boxplots_log.R
# Purpose: Generate publication-grade tissue-stratified boxplots for 7 DNA repair
#          pathways on logarithmic expression scale (YARN qsmooth log2 values)
#          matching the exact design, single-row NER legend, tree tissue ordering,
#          and colorblind-safe palettes of the linear version (Fig. 4b-d, 5b-e).
# Features:
#   - X-axis: 26 tissues in matching Fig. 4a DSB or Fig. 5a SSB tree order
#   - Y-axis: mouse-level YARN qsmooth log2 expression
#   - Points: individual flight mice (n = 360), shape 16, seed 25
#   - Boxplot: Q1, Mean (middle solid bar), Q3, 1.5 * IQR whiskers
#   - Dashed lines: connect tissue-level mean values for each gene
#   - Significance is supplied separately as the approved external software table.
#   - Colorblind-safe palette: strict avoidance of red+green and blue+yellow
#   - Legend: single-row (nrow = 1), clean gene symbols without parentheses
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

args <- commandArgs(trailingOnly = TRUE)
pub_dir <- Sys.getenv("PUBLICATION_OUTPUT_DIR", unset = "")

if (nzchar(pub_dir)) {
  out_dir <- normalizePath(pub_dir, mustWork = TRUE)
  fig_dir <- file.path(out_dir, "Main")
  tbl_dir <- file.path(out_dir, "tables/log")
} else if (length(args)) {
  out_dir <- if (!grepl("^/", args[[1L]])) file.path(repo_root, args[[1L]]) else args[[1L]]
  fig_dir <- file.path(out_dir, "figures")
  tbl_dir <- file.path(out_dir, "tables")
} else if (nzchar(Sys.getenv("LOG_OUTPUT_DIR"))) {
  out_dir <- Sys.getenv("LOG_OUTPUT_DIR")
  if (!grepl("^/", out_dir)) out_dir <- file.path(repo_root, out_dir)
  fig_dir <- file.path(out_dir, "figures")
  tbl_dir <- file.path(out_dir, "tables")
} else {
  out_dir <- file.path(repo_root, "results/log_scale")
  fig_dir <- file.path(out_dir, "figures")
  tbl_dir <- file.path(out_dir, "tables")
}
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

# Input data paths
qsmooth_path <- file.path(repo_root, "data/publication_input/qsmooth/06_component_flight_sample_yarn_qsmooth_values.csv.gz")
tree_order_path <- file.path(repo_root, "data/publication_input/qsmooth/02_tissue_spearman_cluster_order.csv")

if (!file.exists(qsmooth_path)) stop("Missing qsmooth sample table: ", qsmooth_path, call. = FALSE)
if (!file.exists(tree_order_path)) stop("Missing Fig. 4a/5a tree-order contract: ", tree_order_path, call. = FALSE)

dt <- fread(qsmooth_path, showProgress = FALSE)
tree_order <- fread(tree_order_path)
required_tree_columns <- c("Figure", "OrderTopToBottom", "Tissue", "ComponentN")
stopifnot(all(required_tree_columns %in% names(tree_order)))

# Match Fig. 4a (DSB) or Fig. 5a (SSB) tree order from left to right on x-axis
tree_order_by_group <- lapply(c("DSB", "SSB"), function(figure_group) {
  ordered <- tree_order[Figure == figure_group][order(OrderTopToBottom)]
  stopifnot(
    nrow(ordered) == 26L,
    identical(ordered$OrderTopToBottom, seq_len(26L)),
    !anyDuplicated(ordered$Tissue),
    setequal(ordered$Tissue, unique(dt$Group))
  )
  ordered$Tissue
})
names(tree_order_by_group) <- c("DSB", "SSB")

display_tissue_name <- function(x) {
  x <- as.character(x)
  x[x == "Heart / Heart right ventricle"] <- "Heart right ventricle"
  x
}

# Define pathway gene definitions with stable Ensembl IDs, mapped to final publication figure panels
# Strict avoidance of red+green and blue+yellow combinations (only red, only blue/purple)
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
      Nhej1 = "#0072B2",
      Paxx = "#222222",
      Xrcc6 = "#785EF0",
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
      Bard1 = "#785EF0",
      Blm = "#222222",
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
      Parp1 = "#D55E00",
      Polq = "#785EF0",
      Lig1 = "#222222",
      Lig3 = "#0072B2"
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
      Ung = "#0072B2",
      Ogg1 = "#785EF0",
      Neil1 = "#D55E00"
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
    caption_note = "Shapes: Cross (\u2715), Triangle (\u25b2), Circle (\u25cf). Color coding: XPC/RAD23B/CETN2 (GG-NER surveillance; Blue), DDB1/2 (UV lesion detection; Purple), ERCC6/8 (TC-NER; Vermilion).",
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
      Xpc = "#0072B2",
      Rad23b = "#0072B2",
      Cetn2 = "#0072B2",
      Ddb1 = "#785EF0",
      Ddb2 = "#785EF0",
      Ercc6 = "#D55E00",
      Ercc8 = "#D55E00"
    ),
    linetypes = c(
      Xpc = "dashed",
      Rad23b = "dotted",
      Cetn2 = "dotdash",
      Ddb1 = "dashed",
      Ddb2 = "dotted",
      Ercc6 = "dashed",
      Ercc8 = "dotted"
    ),
    shapes = c(
      Xpc = 4,
      Rad23b = 17,
      Cetn2 = 16,
      Ddb1 = 4,
      Ddb2 = 17,
      Ercc6 = 4,
      Ercc8 = 17
    ),
    labels = c(
      Xpc = "Xpc",
      Rad23b = "Rad23b",
      Cetn2 = "Cetn2",
      Ddb1 = "Ddb1",
      Ddb2 = "Ddb2",
      Ercc6 = "Ercc6",
      Ercc8 = "Ercc8"
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
      Fancd2 = "#D55E00",
      Fanci = "#0072B2"
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

summary_mean_list <- list()
used_data_list <- list()

axis_order_audit <- rbindlist(lapply(names(pathway_defs), function(pw_name) {
  pw_info <- pathway_defs[[pw_name]]
  tissue_order <- tree_order_by_group[[pw_info$fig_group]]
  data.table(
    FigureGroup = pw_info$fig_group,
    Pathway = pw_name,
    FigurePanel = pw_info$fig_panel,
    OrderLeftToRight = seq_along(tissue_order),
    Tissue = tissue_order,
    DisplayTissue = display_tissue_name(tissue_order)
  )
}))

for (pw_name in names(pathway_defs)) {
  pw_info <- pathway_defs[[pw_name]]
  gene_ids <- pw_info$genes
  gene_symbols <- names(gene_ids)
  tissue_order <- tree_order_by_group[[pw_info$fig_group]]
  
  sub_dt <- dt[EnsemblID %in% gene_ids]
  stopifnot(setequal(unique(sub_dt$EnsemblID), unname(gene_ids)))
  sub_dt[, Gene := factor(names(gene_ids)[match(EnsemblID, gene_ids)], levels = gene_symbols)]
  sub_dt[, Group := factor(as.character(Group), levels = tissue_order)]
  sub_dt[, Pathway := pw_name]
  
  used_data_list[[pw_name]] <- sub_dt[, .(Pathway = pw_name, EnsemblID, OfficialMouseSymbol, SampleID, Group, accession, mission_cluster, YARNNormalizedLog2)]
  
  # Compute boxplot statistics with mean as the center line on log2 values
  box_stats <- sub_dt[, {
    valid_vals <- YARNNormalizedLog2[!is.na(YARNNormalizedLog2)]
    q <- quantile(valid_vals, c(0.25, 0.75))
    iqr <- q[2] - q[1]
    m <- mean(valid_vals)
    sd_val <- sd(valid_vals)
    n_mice <- length(valid_vals)
    ymin <- max(min(valid_vals), q[1] - 1.5 * iqr)
    ymax <- min(max(valid_vals), q[2] + 1.5 * iqr)
    .(ymin = ymin, lower = q[1], middle = m, upper = q[2], ymax = ymax, SD = sd_val, N = n_mice)
  }, by = .(Group, EnsemblID, Gene)]
  
  summary_mean_list[[pw_name]] <- box_stats[, .(Pathway = pw_name, Group, EnsemblID, Gene, Mean = middle, SD, Q1 = lower, Q3 = upper, Ymin = ymin, Ymax = ymax, N)]
  
  # Layout dodge width
  n_genes <- length(gene_symbols)
  dodge_w <- if (n_genes >= 6) 0.82 else if (n_genes >= 4) 0.78 else 0.70
  box_w <- if (n_genes >= 6) 0.74 else if (n_genes >= 4) 0.70 else 0.65
  dodge <- position_dodge(width = dodge_w)
  
  sub_dt_points <- sub_dt[!is.na(YARNNormalizedLog2)]
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
      aes(x = Group, y = YARNNormalizedLog2, color = Gene, shape = Gene),
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
    scale_x_discrete(labels = display_tissue_name) +
    guides(
      color = guide_legend(
        title = "Gene",
        override.aes = list(
          shape = pw_shapes,
          size = 8.5,
          stroke = 1.8,
          alpha = 1,
          linetype = 0
        ),
        nrow = 1
      ),
      shape = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.06))) +
    labs(
      title = paste0(pw_info$title, " Expression across 26 Tissues (Log2 Scale)"),
      caption = paste0(
        "Boxes represent Q1, Mean (middle solid bar), and Q3; whiskers extend to 1.5 * IQR. Points represent individual flight mice (n = 360).\n",
        "Dashed lines connect tissue-level mean values for each gene. Values shown on log2 scale (YARN qsmooth). ", pw_info$caption_note
      ),
      x = NULL,
      y = "YARN qsmooth normalized expression (log2 scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#E8ECEF", linewidth = 0.5),
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, face = "bold", colour = "#222222", size = 11.5),
      axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
      plot.title = element_text(face = "bold", size = 22, colour = "#1A202C", margin = margin(b = 14)),
      plot.caption = element_text(colour = "#718096", size = 10, lineheight = 1.25, margin = margin(t = 10)),
      legend.position = "top",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = 23),
      legend.text = element_text(size = 23),
      legend.key.size = unit(1.0, "cm"),
      plot.margin = margin(16, 20, 16, 20)
    )
  
  # File export (PNG and PDF)
  panel_id <- pw_info$fig_panel
  safe_name <- gsub("[^A-Za-z0-9_]", "-", pw_name)
  
  # Save primary panel name
  png_file <- file.path(fig_dir, sprintf("%s.png", panel_id))
  pdf_file <- file.path(fig_dir, sprintf("%s.pdf", panel_id))
  ggsave(png_file, p, width = 18, height = 8.5, dpi = 300, bg = "white")
  ggsave(pdf_file, p, width = 18, height = 8.5, device = grDevices::cairo_pdf, bg = "white")
  
  message("Generated ", panel_id, " (", pw_name, " log2): ", png_file)
}

# Export summary tables
summary_mean_dt <- rbindlist(summary_mean_list)
used_data_dt <- rbindlist(used_data_list)

fwrite(summary_mean_dt, file.path(tbl_dir, "03_tissue_gene_mean_log_qsmooth_summary.csv"))
fwrite(used_data_dt, file.path(tbl_dir, "04_sample_log_qsmooth_values_used.csv.gz"))

if (nzchar(pub_dir)) {
  provenance_dir <- file.path(out_dir, "provenance")
  dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(axis_order_audit, file.path(provenance_dir, "pathway_boxplot_x_axis_order.csv"))
  
}

message("Generated seven log2 expression panels and expression summaries: ", out_dir)
