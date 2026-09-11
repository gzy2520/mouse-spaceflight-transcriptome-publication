#!/usr/bin/env Rscript

# Teacher-defined essential DNA-repair components: an expression + response
# companion to the percentile heatmap. Stable Ensembl IDs are the sole merge
# key. Official mouse symbols are display-only.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]]) else normalizePath(".")

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

outdir <- file.path(
  root,
  "03_analysis_results",
  "11_essential_components_percentile_heatmap_20260825"
)
figure_dir <- file.path(outdir, "figures")
table_dir <- file.path(outdir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

component_mapping_path <- file.path(
  table_dir,
  "01_teacher_component_to_stable_id_mapping_audit.csv"
)
tissue_rna_path <- file.path(
  root,
  "03_analysis_results/10_full_stable_id_rerun_20260823",
  "04_seven_pathway_response_26_tissues_20260824/prepared_inputs",
  "02_selected_tissue_RNA_meta_Ensembl116.tsv.gz"
)
assert(file.exists(component_mapping_path), paste("Missing component stable-ID mapping:", component_mapping_path))
assert(file.exists(tissue_rna_path), paste("Missing 26-tissue RNA meta table:", tissue_rna_path))

components <- fread(component_mapping_path)
required_component_columns <- c(
  "RequestedComponent", "Pathway", "ComponentOrder", "EnsemblID",
  "MGI_ID", "EntrezID", "OfficialMouseSymbol"
)
assert(all(required_component_columns %in% names(components)), "Component mapping audit schema has changed.")
setorder(components, ComponentOrder)
assert(
  nrow(components) == 29L && uniqueN(components$EnsemblID) == 29L &&
    !anyNA(components$OfficialMouseSymbol),
  "The component mapping must contain exactly 27 uniquely mapped stable Ensembl IDs."
)

rna <- fread(tissue_rna_path)
required_rna_columns <- c(
  "Group", "EnsemblID", "log2FoldChange", "tissue_meta_p", "padj",
  "FlightExpressionLog2NormalizedCount", "n_missions", "n_analysis_units",
  "expression_n_missions", "expression_n_analysis_units",
  "expression_n_flight_samples"
)
assert(all(required_rna_columns %in% names(rna)), "26-tissue RNA meta table schema has changed.")

plot_table <- merge(
  components[, .(
    RequestedComponent, Pathway, ComponentOrder, EnsemblID, MGI_ID,
    EntrezID, OfficialMouseSymbol
  )],
  rna[, ..required_rna_columns],
  by = "EnsemblID",
  all.x = TRUE,
  sort = FALSE
)
assert(
  nrow(plot_table) == 26L * 29L && uniqueN(plot_table$Group) == 26L &&
    !anyDuplicated(plot_table[, .(Group, EnsemblID)]),
  "The companion matrix must contain one row for each of 26 tissues x 27 components."
)
assert(
  all(is.finite(plot_table$FlightExpressionLog2NormalizedCount)) &&
    all(is.finite(plot_table$log2FoldChange)) &&
    all(is.finite(plot_table$tissue_meta_p)) &&
    all(plot_table$tissue_meta_p >= 0 & plot_table$tissue_meta_p <= 1),
  "Expression, log2FC, or nominal tissue-meta P values are invalid."
)

# The sample-column values are already GeneLab's normalized-count scale with
# a +1 pseudocount. The source workflow therefore applies log2 directly
# (without adding another pseudocount), averages analysis units within mission,
# and averages missions equally within tissue.
nominal_p_cutoff <- 0.05
plot_table[, `:=`(
  Log2FCDisplayed = fifelse(tissue_meta_p < nominal_p_cutoff, log2FoldChange, NA_real_),
  NominalPDisplayed = tissue_meta_p < nominal_p_cutoff,
  ExpressionLabel = sprintf("%.1f", FlightExpressionLog2NormalizedCount),
  Log2FCLabel = fifelse(
    tissue_meta_p < nominal_p_cutoff,
    sprintf("%+.2f", log2FoldChange),
    ""
  ),
  PathwayDisplay = fifelse(Pathway == "HR_A-EJ", "HR–A-EJ", Pathway)
)]
plot_table[, `:=`(
  ExpressionWithinTissueScaled = {
    lower <- quantile(FlightExpressionLog2NormalizedCount, probs = 0.05, names = FALSE)
    upper <- quantile(FlightExpressionLog2NormalizedCount, probs = 0.95, names = FALSE)
    if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
      rep(0.5, .N)
    } else {
      pmin(1, pmax(0, (FlightExpressionLog2NormalizedCount - lower) / (upper - lower)))
    }
  }
), by = Group]

tissue_labels <- c(
  "Adrenal gland" = "肾上腺", "Bone marrow" = "骨髓", "Cecum" = "盲肠", "Cerebellum" = "小脑",
  "Colon" = "结肠", "Dorsal skin" = "背部皮肤", "Extensor digitorum longus" = "趾长伸肌", "Eye" = "眼",
  "Femoral lateral skin" = "股外侧皮肤", "Femoral skin" = "股部皮肤", "Gastrocnemius" = "腓肠肌",
  "Heart" = "心脏", "Heart / Heart right ventricle" = "心脏 / 右心室", "Kidney" = "肾脏",
  "Left lobe of the liver" = "肝左叶", "Liver" = "肝脏", "Lung" = "肺", "Mammary gland" = "乳腺",
  "Optic nerve" = "视神经", "Quadriceps femoris" = "股四头肌", "Retina" = "视网膜", "Soleus" = "比目鱼肌",
  "Spleen" = "脾脏", "Spleen-distal" = "脾脏（远端）", "Thymus" = "胸腺", "Tibialis anterior" = "胫骨前肌"
)
tissue_order <- sort(unique(plot_table$Group))
assert(identical(sort(names(tissue_labels)), tissue_order), "The bilingual tissue-label contract does not cover exactly the 26 tissues.")
pathway_order <- unique(components$Pathway)
tissue_display_labels <- setNames(paste0(tissue_order, "\n", tissue_labels[tissue_order]), tissue_order)

plot_table[, `:=`(
  TissueLabel = factor(
    paste0(Group, "\n", tissue_labels[Group]),
    levels = unname(tissue_display_labels[rev(tissue_order)])
  ),
  DisplayComponent = factor(OfficialMouseSymbol, levels = components$OfficialMouseSymbol),
  PathwayDisplay = factor(
    PathwayDisplay,
    levels = ifelse(pathway_order == "HR_A-EJ", "HR–A-EJ", pathway_order)
  )
)]

fc_limit <- quantile(abs(plot_table$Log2FCDisplayed), probs = 0.95, na.rm = TRUE, names = FALSE)
if (!is.finite(fc_limit) || fc_limit <= 0) fc_limit <- 1
assert(sum(plot_table$NominalPDisplayed) > 0L, "No component has nominal tissue-meta P < 0.05; the log2FC panel would be empty.")

shared_theme <- theme_minimal(base_size = 8.5, base_family = "Arial Unicode MS") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 6.7, lineheight = 0.85),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", colour = NA),
    strip.text.x.bottom = element_text(face = "bold", size = 9),
    strip.placement = "outside",
    panel.spacing.x = grid::unit(0.45, "lines"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5),
    plot.caption = element_text(size = 7.2, hjust = 0, colour = "#444444"),
    legend.position = "right",
    plot.margin = margin(6, 10, 6, 7)
  )

expression_plot <- ggplot(
  plot_table,
  aes(x = DisplayComponent, y = TissueLabel, fill = ExpressionWithinTissueScaled)
) +
  geom_tile(colour = "white", linewidth = 0.18) +
  geom_text(aes(label = ExpressionLabel), size = 2.25, colour = "#2F2F2F") +
  facet_grid(. ~ PathwayDisplay, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FFF3E0", "#FDBB84", "#FC8D59", "#B2182B"),
    values = rescale(c(0, 0.20, 0.50, 0.80, 1)),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("low", "mid", "high"),
    name = "Within-tissue\nexpression position"
  ) +
  labs(
    title = "Flight expression: log2 GeneLab sample-column value",
    subtitle = "Cell labels are mission-equal mean log2 expression; colour is independently scaled within each tissue (5th–95th percentile).",
    x = NULL, y = NULL
  ) +
  shared_theme

logfc_plot <- ggplot(
  plot_table,
  aes(x = DisplayComponent, y = TissueLabel, fill = Log2FCDisplayed)
) +
  geom_tile(colour = "white", linewidth = 0.18) +
  geom_text(aes(label = Log2FCLabel), size = 2.25, colour = "#2F2F2F") +
  facet_grid(. ~ PathwayDisplay, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0,
    limits = c(-fc_limit, fc_limit), oob = squish,
    na.value = "#EEF1F4", name = "Meta log2FC\n(P < 0.05)"
  ) +
  labs(
    title = "Flight response: tissue meta-log2FC",
    subtitle = "Only cells with unadjusted tissue-meta P < 0.05 are coloured and labelled; blue = lower in flight, red = higher in flight.",
    x = NULL, y = NULL,
    caption = "Meta-log2FC is the median of mission-level median log2FC values. Nominal tissue-meta P is from mission-level signed Stouffer aggregation."
  ) +
  shared_theme

# Unfiltered descriptive log2FC version: every tissue x component cell is
# rendered, irrespective of nominal or adjusted P value. P values remain in
# the audit table and are not encoded as a display filter in this figure.
fc_full_limit <- quantile(abs(plot_table$log2FoldChange), probs = 0.95, na.rm = TRUE, names = FALSE)
if (!is.finite(fc_full_limit) || fc_full_limit <= 0) fc_full_limit <- 1
plot_table[, Log2FCFullLabel := sprintf("%+.2f", log2FoldChange)]
logfc_full_plot <- ggplot(
  plot_table,
  aes(x = DisplayComponent, y = TissueLabel, fill = log2FoldChange)
) +
  geom_tile(colour = "white", linewidth = 0.18) +
  geom_text(aes(label = Log2FCFullLabel), size = 2.25, colour = "#2F2F2F") +
  facet_grid(. ~ PathwayDisplay, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0,
    limits = c(-fc_full_limit, fc_full_limit), oob = squish,
    name = "Meta log2FC\n(all cells)"
  ) +
  labs(
    title = "Flight response: tissue meta-log2FC (all components)",
    subtitle = "All 26 tissues × 27 components are shown; no P-value filtering is applied. Blue = lower in flight, red = higher in flight.",
    x = NULL, y = NULL,
    caption = "Meta-log2FC is the median of mission-level median log2FC values. Tissue-meta P values are retained in the audit table but are not used to hide cells."
  ) +
  shared_theme

expression_stem <- "02a_essential_components_flight_log2_expression_matrix"
ggsave(
  file.path(figure_dir, paste0(expression_stem, ".png")), expression_plot,
  width = 20, height = 10.5, dpi = 320, bg = "white"
)
ggsave(
  file.path(figure_dir, paste0(expression_stem, ".pdf")), expression_plot,
  width = 20, height = 10.5, device = cairo_pdf, bg = "white"
)
logfc_stem <- "02b_essential_components_tissue_meta_log2FC_P_lt_0.05_matrix"
ggsave(
  file.path(figure_dir, paste0(logfc_stem, ".png")), logfc_plot,
  width = 20, height = 10.5, dpi = 320, bg = "white"
)
ggsave(
  file.path(figure_dir, paste0(logfc_stem, ".pdf")), logfc_plot,
  width = 20, height = 10.5, device = cairo_pdf, bg = "white"
)
logfc_full_stem <- "02c_essential_components_tissue_meta_log2FC_all_cells_matrix"
ggsave(
  file.path(figure_dir, paste0(logfc_full_stem, ".png")), logfc_full_plot,
  width = 20, height = 10.5, dpi = 320, bg = "white"
)
ggsave(
  file.path(figure_dir, paste0(logfc_full_stem, ".pdf")), logfc_full_plot,
  width = 20, height = 10.5, device = cairo_pdf, bg = "white"
)

setcolorder(plot_table, c(
  "Group", "TissueLabel", "Pathway", "PathwayDisplay", "ComponentOrder",
  "RequestedComponent", "OfficialMouseSymbol", "EnsemblID", "MGI_ID", "EntrezID",
  "FlightExpressionLog2NormalizedCount", "ExpressionWithinTissueScaled",
  "log2FoldChange", "tissue_meta_p", "padj", "Log2FCDisplayed", "NominalPDisplayed",
  "n_missions", "n_analysis_units", "expression_n_missions",
  "expression_n_analysis_units", "expression_n_flight_samples"
))
fwrite(plot_table, file.path(table_dir, "08_essential_components_counts_log2FC_tissue_matrix.csv"))
fwrite(
  data.table(
    Field = c(
    "Component membership source", "Stable-ID analysis key", "RNA tissue-meta input",
      "Tissue scope", "Expression value", "Expression colour", "log2FC aggregation",
      "P-value display rule", "Random seed", "GO pathway status"
    ),
    Value = c(
      paste0(component_mapping_path, "; RPA2/RPA3 appended per teacher supplementary instruction"), "Ensembl Gene ID", tissue_rna_path,
      "26 original-Material-Type tissues; 48 OSDR accessions; 59 analysis units",
      "Mission-equal mean flight log2(GeneLab sample-column normalized-count value; +1 pseudocount already present in source)",
      "Independently scaled within tissue from 5th to 95th percentile across the 27 teacher-defined components; numeric cell labels retain log2 expression value",
      "Median of mission-level median log2FC values",
      "Show and colour log2FC only when unadjusted tissue-meta P < 0.05; P from mission-level signed Stouffer aggregation",
      "25", "Previous GO seven-pathway membership is not used; the teacher DSB worksheet defines component membership and grouping."
    )
  ),
  file.path(table_dir, "09_essential_components_counts_log2FC_metadata.csv")
)

readme <- c(
  "# Essential components: expression and log2FC companion matrix",
  "",
  "- `figures/02a_essential_components_flight_log2_expression_matrix.png` is the standalone expression figure: mission-equal mean flight `log2(GeneLab sample-column value)`. Input values already contain the GeneLab +1 pseudocount, so no additional +1 is added. Colour is scaled independently within tissue and labels show the log2 expression values.",
  "- `figures/02b_essential_components_tissue_meta_log2FC_P_lt_0.05_matrix.png` is the standalone response figure: tissue meta-log2FC. A cell is shown only for unadjusted tissue-meta `P < 0.05`; its numeric label is log2FC. The P value is the signed Stouffer meta P from mission-level statistics, and log2FC is the median of mission-level median log2FC values.",
  "- `figures/02c_essential_components_tissue_meta_log2FC_all_cells_matrix.png` is the unfiltered response figure: all 702 tissue–component cells are shown, regardless of nominal or adjusted P value. It is descriptive and should not be interpreted as evidence of significance from colour alone.",
  "- All joins use Ensembl Gene ID. Official mouse symbols are display-only. The `PRA1` label in the teacher source is interpreted as RPA1/Rpa1 (`ENSMUSG00000000751`); see `01_teacher_component_to_stable_id_mapping_audit.csv`."
)
writeLines(readme, file.path(outdir, "02_counts_logFC_README.md"), useBytes = TRUE)

message("Created expression + log2FC companion matrix in: ", outdir)
