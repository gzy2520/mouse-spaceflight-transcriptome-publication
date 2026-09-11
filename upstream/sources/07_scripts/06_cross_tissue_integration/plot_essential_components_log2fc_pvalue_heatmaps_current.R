#!/usr/bin/env Rscript

# Separate DSB and SSB tissue meta-log2FC heatmaps.
#
# This is the teacher-requested log2FC + P-value-display version.  The
# historical display rule is retained: only cells with unadjusted
# tissue_meta_p < 0.05 are coloured and labelled; all cells and their exact
# log2FC, nominal P, and padj values are retained in the audit table.
# Stable Ensembl Gene IDs are the only merge key.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]]) else normalizePath(".")

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

normalize_ensembl <- function(x) {
  sub("\\.[0-9]+$", "", trimws(as.character(x)))
}

outdir <- file.path(
  root,
  "03_analysis_results",
  "13_essential_components_log2fc_pvalue_heatmaps_20260825"
)
figure_dir <- file.path(outdir, "figures")
table_dir <- file.path(outdir, "tables")
input_dir <- file.path(outdir, "input")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

mapping_path <- file.path(
  root,
  "03_analysis_results/12_essential_components_expression_log2_heatmap_20260825",
  "tables/01_included_component_stable_id_mapping_audit.csv"
)
scope_path <- file.path(
  root,
  "03_analysis_results/12_essential_components_expression_log2_heatmap_20260825",
  "tables/00_teacher_component_scope_audit.csv"
)
unit_manifest_path <- file.path(
  root,
  "03_analysis_results/10_full_stable_id_rerun_20260823",
  "04_seven_pathway_response_26_tissues_20260824/prepared_inputs",
  "01_selected_tissue_analysis_units.csv"
)
expression_rna_path <- file.path(
  root,
  "03_analysis_results/10_full_stable_id_rerun_20260823",
  "04_seven_pathway_response_26_tissues_20260824/prepared_inputs",
  "02_selected_tissue_RNA_meta_Ensembl116.tsv.gz"
)

required_files <- c(mapping_path, scope_path, unit_manifest_path, expression_rna_path)
assert(
  all(file.exists(required_files)),
  paste("Missing required input(s):", paste(required_files[!file.exists(required_files)], collapse = "; "))
)

scope <- fread(scope_path, showProgress = FALSE)
components <- fread(mapping_path, showProgress = FALSE)
unit_manifest <- fread(unit_manifest_path, showProgress = FALSE)
assert(nrow(scope) == 85L && sum(scope$Included) == 75L && sum(scope$Included == FALSE) == 10L,
       "The teacher component scope audit is not the expected 85-row / 75-included contract.")
assert(nrow(components) == 75L && uniqueN(components$EnsemblID) == 75L,
       "The included component mapping must contain 75 unique stable Ensembl IDs.")
assert(nrow(unit_manifest) == 59L && uniqueN(unit_manifest$Group) == 26L,
       "The selected manifest must contain 59 analysis units across 26 tissues.")
setorder(components, ComponentOrder)

rna <- fread(expression_rna_path, showProgress = FALSE)
required_rna_columns <- c(
  "Group", "EnsemblID", "log2FoldChange", "tissue_meta_p", "padj",
  "n_missions", "n_analysis_units", "expression_n_missions",
  "expression_n_analysis_units", "expression_n_flight_samples"
)
assert(all(required_rna_columns %in% names(rna)), "The RNA meta table schema has changed.")
rna[, EnsemblID := normalize_ensembl(EnsemblID)]
rna_component <- rna[EnsemblID %chin% components$EnsemblID, ..required_rna_columns]
assert(!anyDuplicated(rna_component[, .(Group, EnsemblID)]),
       "The RNA meta table contains duplicated tissue + Ensembl rows.")
assert(nrow(rna_component) == 26L * 75L && uniqueN(rna_component$Group) == 26L,
       "The selected component RNA table must contain exactly 26 x 75 rows.")
assert(
  all(is.finite(rna_component$log2FoldChange)) &&
    all(is.finite(rna_component$tissue_meta_p)) &&
    all(rna_component$tissue_meta_p >= 0 & rna_component$tissue_meta_p <= 1) &&
    all(is.finite(rna_component$padj)) &&
    all(rna_component$padj >= 0 & rna_component$padj <= 1),
  "The selected component log2FC, P, or padj values are invalid."
)

plot_table <- merge(
  components[, .(
    Pathway, PathwayOrder, ComponentOrderWithinPathway, ComponentOrder,
    RequestedComponent, MappingStatus, LookupMouseSymbol, DisplayComponent,
    EnsemblID, MGI_ID, EntrezID, OfficialMouseSymbol, PathwayDisplay
  )],
  rna_component,
  by = "EnsemblID",
  all.x = TRUE,
  sort = FALSE,
  allow.cartesian = TRUE
)
assert(
  nrow(plot_table) == 26L * 75L && uniqueN(plot_table$Group) == 26L &&
    !anyDuplicated(plot_table[, .(Group, EnsemblID)]),
  "The final log2FC/P matrix must contain one row for each of 26 tissues x 75 components."
)
plot_table[, `:=`(
  P_lt_0_05 = tissue_meta_p < 0.05,
  Log2FCDisplayed = fifelse(tissue_meta_p < 0.05, log2FoldChange, NA_real_),
  Log2FCLabel = fifelse(tissue_meta_p < 0.05, sprintf("%+.2f", log2FoldChange), ""),
  PathwayDisplay = fifelse(Pathway == "HR_A-EJ", "HR–A-EJ", Pathway)
)]

tissue_labels <- c(
  "Adrenal gland" = "肾上腺", "Bone marrow" = "骨髓", "Cecum" = "盲肠", "Cerebellum" = "小脑",
  "Colon" = "结肠", "Dorsal skin" = "背部皮肤", "Extensor digitorum longus" = "趾长伸肌", "Eye" = "眼",
  "Femoral lateral skin" = "股外侧皮肤", "Femoral skin" = "股部皮肤", "Gastrocnemius" = "腓肠肌",
  "Heart" = "心脏", "Heart / Heart right ventricle" = "心脏 / 右心室", "Kidney" = "肾脏",
  "Left lobe of the liver" = "肝左叶", "Liver" = "肝脏", "Lung" = "肺", "Mammary gland" = "乳腺",
  "Optic nerve" = "视神经", "Quadriceps femoris" = "股四头肌", "Retina" = "视网膜", "Soleus" = "比目鱼肌",
  "Spleen" = "脾脏", "Spleen-distal" = "脾脏（远端）", "Thymus" = "胸腺", "Tibialis anterior" = "胫骨前肌"
)
tissue_order <- sort(unique(unit_manifest$Group))
assert(identical(sort(names(tissue_labels)), tissue_order), "The bilingual tissue-label contract does not cover exactly the 26 tissues.")
tissue_display_labels <- setNames(paste0(tissue_order, "\n", tissue_labels[tissue_order]), tissue_order)
plot_table[, `:=`(
  TissueLabel = paste0(Group, "\n", unname(tissue_labels[Group])),
  TissueOrder = match(Group, tissue_order)
)]
setorder(plot_table, TissueOrder, ComponentOrder)
plot_table[, TissueLabel := factor(TissueLabel, levels = unname(tissue_display_labels[rev(tissue_order)]))]

displayed_fc <- plot_table[P_lt_0_05 == TRUE, abs(log2FoldChange)]
fc_limit <- as.numeric(quantile(displayed_fc, probs = 0.95, names = FALSE, type = 7))
if (!is.finite(fc_limit) || fc_limit <= 0) fc_limit <- max(displayed_fc, na.rm = TRUE)
if (!is.finite(fc_limit) || fc_limit <= 0) fc_limit <- 1
plot_table[, `:=`(
  DisplayComponent = as.character(DisplayComponent),
  PathwayDisplay = as.character(PathwayDisplay)
)]

shared_theme <- theme_minimal(base_size = 8.8, base_family = "Arial Unicode MS") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.5, colour = "#3D3D3D"),
    axis.text.y = element_text(size = 7.1, lineheight = 0.84, colour = "#3D3D3D"),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", colour = NA),
    strip.text.x.bottom = element_text(face = "bold", size = 9.5),
    strip.placement = "outside",
    panel.spacing.x = grid::unit(0.42, "lines"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 7.5, hjust = 0, colour = "#444444"),
    legend.position = "right",
    plot.margin = margin(8, 12, 10, 8)
  )

make_fc_plot <- function(plot_data, title, subtitle, text_size = 2.0) {
  d <- copy(plot_data)
  component_levels <- unique(d[order(ComponentOrder), DisplayComponent])
  pathway_levels <- unique(d[order(PathwayOrder), PathwayDisplay])
  d[, `:=`(
    DisplayComponent = factor(DisplayComponent, levels = component_levels),
    PathwayDisplay = factor(PathwayDisplay, levels = pathway_levels)
  )]
  ggplot(d, aes(x = DisplayComponent, y = TissueLabel, fill = Log2FCDisplayed)) +
    geom_tile(colour = "white", linewidth = 0.16) +
    geom_text(aes(label = Log2FCLabel), size = text_size, colour = "#2F2F2F", family = "Arial Unicode MS") +
    facet_grid(. ~ PathwayDisplay, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0,
      limits = c(-fc_limit, fc_limit), oob = squish,
      na.value = "#EEF1F4",
      name = "Meta log2FC\n(P < 0.05)"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      caption = paste(
        "Only unadjusted tissue_meta_p < 0.05 cells are coloured and labelled; grey cells are not displayed as significant.",
        "log2FC is the median of mission-level median log2FC values; P is from mission-level signed Stouffer aggregation."
      )
    ) +
    shared_theme
}

dsb_plot <- make_fc_plot(
  plot_table[Pathway %chin% c("NHEJ", "HR", "HR_A-EJ", "A-EJ")],
  "DSB essential components: tissue meta-log2FC (P < 0.05)",
  "26 mouse tissues; 29 DSB components; blue = lower in flight, red = higher in flight",
  text_size = 2.05
)
ssb_plot <- make_fc_plot(
  plot_table[Pathway %chin% c("BER", "NER", "MMR", "FA")],
  "SSB essential components: tissue meta-log2FC (P < 0.05)",
  "26 mouse tissues; 46 SSB components; shared rows excluded; blue = lower in flight, red = higher in flight",
  text_size = 1.9
)

dsb_stem <- "01_DSB_essential_components_tissue_meta_log2FC_P_lt_0.05_heatmap"
ssb_stem <- "02_SSB_essential_components_tissue_meta_log2FC_P_lt_0.05_heatmap"
ggsave(file.path(figure_dir, paste0(dsb_stem, ".png")), dsb_plot, width = 20, height = 12.5, dpi = 320, bg = "white", limitsize = FALSE)
ggsave(file.path(figure_dir, paste0(dsb_stem, ".pdf")), dsb_plot, width = 20, height = 12.5, device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(file.path(figure_dir, paste0(ssb_stem, ".png")), ssb_plot, width = 28, height = 12.5, dpi = 320, bg = "white", limitsize = FALSE)
ggsave(file.path(figure_dir, paste0(ssb_stem, ".pdf")), ssb_plot, width = 28, height = 12.5, device = cairo_pdf, bg = "white", limitsize = FALSE)

setcolorder(plot_table, c(
  "Group", "TissueLabel", "TissueOrder", "Pathway", "PathwayDisplay", "PathwayOrder",
  "ComponentOrder", "ComponentOrderWithinPathway", "RequestedComponent", "MappingStatus",
  "LookupMouseSymbol", "OfficialMouseSymbol", "DisplayComponent", "EnsemblID", "MGI_ID", "EntrezID",
  "log2FoldChange", "tissue_meta_p", "padj", "P_lt_0_05", "Log2FCDisplayed", "Log2FCLabel",
  "n_missions", "n_analysis_units", "expression_n_missions", "expression_n_analysis_units", "expression_n_flight_samples"
))
fwrite(plot_table, file.path(table_dir, "01_DSB_SSB_tissue_meta_log2FC_pvalue_matrix.csv"))
fwrite(components, file.path(table_dir, "02_included_component_stable_id_mapping_audit.csv"))
fwrite(unit_manifest, file.path(table_dir, "03_26_tissue_analysis_unit_manifest.csv"))

metadata <- data.table(
  Field = c(
    "Component scope", "Tissue scope", "Stable-ID analysis key", "RNA input",
    "Effect value", "P-value value", "Display rule", "P-value adjustment retained",
    "Common colour limit", "Figures", "Random seed"
  ),
  Value = c(
    "75 included components: DSB 29 + SSB 46; shared 10 excluded",
    "26 tissues; 48 OSDR accessions; 59 analysis units",
    "Ensembl Gene ID",
    expression_rna_path,
    "tissue-level log2FoldChange = median of mission-level median log2FC values",
    "tissue_meta_p = two-sided P from mission-level signed Stouffer aggregation",
    "Only unadjusted tissue_meta_p < 0.05 cells are coloured and labelled; no padj filtering",
    "padj is retained in the matrix table but is not used to hide cells",
    sprintf("Shared symmetric scale: +/- %.4f (95th percentile of absolute displayed log2FC across all 75 components)", fc_limit),
    "Separate DSB and SSB PNG/PDF figures",
    "25"
  )
)
fwrite(metadata, file.path(table_dir, "04_run_metadata.csv"))

file.copy(mapping_path, file.path(input_dir, "included_component_stable_id_mapping_audit.csv"), overwrite = TRUE)
file.copy(scope_path, file.path(input_dir, "teacher_component_scope_audit.csv"), overwrite = TRUE)

readme_lines <- c(
  "# DSB / SSB tissue meta-log2FC + P-value display heatmaps",
  "",
  "- 按老师要求，将 DSB 和 SSB 分成两张图。DSB 为 29 个成分；SSB 为 46 个成分；`shared` 10 个成分不纳入。",
  "- 两张图均使用当前 26 个组织、48 个 OSDR accession、59 个 analysis unit 的固定 RNA meta 输入。",
  "",
  "## 图中显示规则",
  "",
  "- 只有 `tissue_meta_p < 0.05` 的格子才着色并显示 log2FC 数值；其余格子为灰色，不解释为显著变化。这里使用的是未校正 tissue-meta P，沿用之前确认的 0.05 显示规则。",
  "- 颜色使用 DSB+SSB 全部 75 个成分共同计算的对称范围，因此两张图可以直接比较颜色深浅。蓝色表示飞行中较低，红色表示飞行中较高。",
  "- 完整的 log2FC、`tissue_meta_p`、`padj` 和显示标记保存在 `tables/01_DSB_SSB_tissue_meta_log2FC_pvalue_matrix.csv`；padj 不用于本图隐藏格子。",
  "- 所有匹配使用 Ensembl Gene ID；官方 mouse Symbol 只用于横轴显示。",
  "",
  "## 输出",
  "",
  "- `figures/01_DSB_essential_components_tissue_meta_log2FC_P_lt_0.05_heatmap.png/pdf`：DSB 独立图。",
  "- `figures/02_SSB_essential_components_tissue_meta_log2FC_P_lt_0.05_heatmap.png/pdf`：SSB 独立图。",
  "- `tables/01_DSB_SSB_tissue_meta_log2FC_pvalue_matrix.csv`：26 × 75 = 1950 个组织--成分格子的完整数值审计。"
)
writeLines(readme_lines, file.path(outdir, "README.md"), useBytes = TRUE)

message("Created separate DSB and SSB tissue meta-log2FC/P<0.05 heatmaps in: ", outdir)
