#!/usr/bin/env Rscript

# Re-render the SSB tissue meta-log2FC heatmap after the teacher-requested
# removal of Neil2.  The source values and the P < 0.05 display rule are
# unchanged from the frozen 26-tissue meta-log2FC table; only the component
# row is excluded and the SSB panel is re-laid out.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

root <- normalizePath(Sys.getenv("GRADUATE_DESIGN_ROOT", unset = getwd()), mustWork = TRUE)
source_dir <- file.path(root, "03_analysis_results/13_essential_components_log2fc_pvalue_heatmaps_20260825")
out_dir <- Sys.getenv(
  "SSB_META_OUTPUT_DIR",
  unset = file.path(root, "03_analysis_results/24_teacher_selected_figure_release_20260829/ssb_meta_no_neil2")
)
if (!grepl("^/", out_dir)) out_dir <- file.path(root, out_dir)
figure_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
input_dir <- file.path(out_dir, "input")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
normalize_ensembl <- function(x) sub("\\.[0-9]+$", "", toupper(trimws(as.character(x))))

matrix_path <- file.path(source_dir, "tables/01_DSB_SSB_tissue_meta_log2FC_pvalue_matrix.csv")
mapping_path <- file.path(source_dir, "tables/02_included_component_stable_id_mapping_audit.csv")
scope_path <- file.path(source_dir, "input/teacher_component_scope_audit.csv")
assert(all(file.exists(c(matrix_path, mapping_path, scope_path))), "Missing frozen meta-log2FC inputs")

full <- fread(matrix_path, showProgress = FALSE)
mapping <- fread(mapping_path, showProgress = FALSE)
scope <- fread(scope_path, showProgress = FALSE)
full[, EnsemblID := normalize_ensembl(EnsemblID)]
mapping[, EnsemblID := normalize_ensembl(EnsemblID)]
assert(nrow(full) == 1950L && uniqueN(full[, .(Group, EnsemblID)]) == 1950L, "Frozen table is not 26 x 75")
assert(nrow(mapping) == 75L && uniqueN(mapping$EnsemblID) == 75L, "Frozen mapping is not 75 unique Ensembl IDs")

neil2_id <- "ENSMUSG00000035121"
excluded <- mapping[
  EnsemblID == neil2_id |
    toupper(trimws(RequestedComponent)) == "NEIL2" |
    toupper(trimws(DisplayComponent)) == "NEIL2" |
    toupper(trimws(OfficialMouseSymbol)) == "NEIL2"
]
assert(nrow(excluded) == 1L && excluded$EnsemblID[[1L]] == neil2_id, "Neil2 mapping was not uniquely resolved")
keep_ids <- setdiff(mapping$EnsemblID, neil2_id)
ssb <- full[
  Pathway %chin% c("BER", "NER", "MMR", "FA") & EnsemblID %chin% keep_ids
]
ssb_map <- mapping[Pathway %chin% c("BER", "NER", "MMR", "FA") & EnsemblID %chin% keep_ids]
assert(nrow(ssb) == 26L * 45L && uniqueN(ssb$EnsemblID) == 45L, "Filtered SSB table is not 26 x 45")
assert(!neil2_id %chin% ssb$EnsemblID, "Neil2 remains in filtered SSB table")
assert(!anyDuplicated(ssb[, .(Group, EnsemblID)]), "Filtered SSB table has duplicate tissue/component rows")

ssb[, `:=`(
  P_lt_0_05 = is.finite(tissue_meta_p) & tissue_meta_p < 0.05,
  Log2FCDisplayed = fifelse(is.finite(tissue_meta_p) & tissue_meta_p < 0.05, log2FoldChange, NA_real_),
  Log2FCLabel = fifelse(is.finite(tissue_meta_p) & tissue_meta_p < 0.05, sprintf("%+.2f", log2FoldChange), ""),
  PathwayDisplay = fifelse(Pathway == "HR_A-EJ", "HR–A-EJ", Pathway)
)]
ssb[, `:=`(TissueOrder = as.integer(TissueOrder), ComponentOrder = as.integer(ComponentOrder), PathwayOrder = as.integer(PathwayOrder))]
setorder(ssb, TissueOrder, PathwayOrder, ComponentOrder)
tissue_label_levels <- full[order(TissueOrder), unique(TissueLabel)]
assert(length(tissue_label_levels) == 26L, "Frozen tissue labels do not contain the 26-tissue order")
# ggplot draws the last discrete level at the top; reversing the established
# first-heatmap order therefore keeps Adrenal gland at the top and Tibialis
# anterior at the bottom, rather than silently reverting to A--Z order.
ssb[, TissueLabel := factor(as.character(TissueLabel), levels = rev(tissue_label_levels))]

displayed_fc <- full[is.finite(tissue_meta_p) & tissue_meta_p < 0.05, abs(log2FoldChange)]
fc_limit <- as.numeric(quantile(displayed_fc, probs = 0.95, names = FALSE, type = 7))
if (!is.finite(fc_limit) || fc_limit <= 0) fc_limit <- max(displayed_fc, na.rm = TRUE)
if (!is.finite(fc_limit) || fc_limit <= 0) fc_limit <- 1

component_levels <- unique(ssb[order(PathwayOrder, ComponentOrder), DisplayComponent])
pathway_levels <- c("BER", "NER", "MMR", "FA")
ssb[, DisplayComponent := factor(as.character(DisplayComponent), levels = component_levels)]
ssb[, PathwayDisplay := factor(as.character(PathwayDisplay), levels = pathway_levels)]

shared_theme <- theme_minimal(base_size = 8.8, base_family = "Arial Unicode MS") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.5, colour = "#3D3D3D"),
    axis.text.y = element_text(size = 7.1, lineheight = 0.84, colour = "#3D3D3D"),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", colour = NA),
    strip.text.x.bottom = element_text(face = "bold", size = 9.5),
    strip.placement = "outside", panel.spacing.x = grid::unit(0.42, "lines"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9), plot.caption = element_text(size = 7.5, hjust = 0, colour = "#444444"),
    legend.position = "right", plot.margin = margin(8, 12, 10, 8)
  )

plot <- ggplot(ssb, aes(x = DisplayComponent, y = TissueLabel, fill = Log2FCDisplayed)) +
  geom_tile(colour = "white", linewidth = 0.16) +
  geom_text(aes(label = Log2FCLabel), size = 1.9, colour = "#2F2F2F", family = "Arial Unicode MS") +
  facet_grid(. ~ PathwayDisplay, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0,
    limits = c(-fc_limit, fc_limit), oob = squish, na.value = "#EEF1F4", name = "Meta log2FC\n(P < 0.05)"
  ) +
  labs(
    title = "SSB essential components: tissue meta-log2FC (P < 0.05)",
    subtitle = "26 mouse tissues; 45 SSB components; Neil2 removed to match Fig 5a; blue = lower in flight, red = higher in flight",
    x = NULL, y = NULL,
    caption = paste(
      "Only unadjusted tissue_meta_p < 0.05 cells are coloured and labelled; grey cells are not displayed as significant.",
      "log2FC is the median of mission-level median log2FC values; P is from mission-level signed Stouffer aggregation."
    )
  ) + shared_theme

stem <- "Fig_5b_SSB_essential_components_tissue_meta_log2FC_P_lt_0.05_without_Neil2"
png_path <- file.path(figure_dir, paste0(stem, ".png"))
pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
ggsave(png_path, plot, width = 28, height = 12.5, units = "in", dpi = 320, bg = "white", limitsize = FALSE)
ggsave(pdf_path, plot, width = 28, height = 12.5, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)

fwrite(ssb, file.path(table_dir, "01_SSB_tissue_meta_log2FC_pvalue_matrix_without_Neil2.csv"), bom = TRUE)
fwrite(ssb_map, file.path(table_dir, "02_SSB_component_stable_id_mapping_without_Neil2.csv"), bom = TRUE)
fwrite(excluded[, .(RequestedComponent, DisplayComponent, OfficialMouseSymbol, EnsemblID, MGI_ID, EntrezID, ExclusionReason = "Teacher-requested removal; match Fig 5a")],
       file.path(table_dir, "03_excluded_Neil2_audit.csv"), bom = TRUE)
fwrite(data.table(
  component_count = uniqueN(ssb$EnsemblID), tissue_count = uniqueN(ssb$Group), cell_count = nrow(ssb),
  excluded_component = "Neil2", excluded_ensembl_id = neil2_id, display_rule = "unadjusted tissue_meta_p < 0.05",
  common_symmetric_fc_limit = fc_limit, png = file.path("figures", basename(png_path)), pdf = file.path("figures", basename(pdf_path))
), file.path(table_dir, "04_run_audit.csv"), bom = TRUE)

file.copy(mapping_path, file.path(input_dir, basename(mapping_path)), overwrite = TRUE)
file.copy(scope_path, file.path(input_dir, basename(scope_path)), overwrite = TRUE)
writeLines(c(
  "# Fig 5b：去除 Neil2 的 SSB tissue meta-log2FC 图", "",
  "- 使用冻结的 26 组织 × 75 成分 meta-log2FC/P 表；仅从 SSB 的 BER、NER、MMR、FA 中去除 Neil2（Ensembl `ENSMUSG00000035121`）。",
  "- 数值、组织顺序、颜色范围和未校正 `tissue_meta_p < 0.05` 的显示规则与原 Fig 5b 一致；未重新计算上游 meta 统计。",
  "- 过滤后的表和排除审计在 `tables/`；原始 75 成分版本仍保留在 13 号目录。"
), file.path(out_dir, "README.md"), useBytes = TRUE)
cat("Wrote Neil2-excluded SSB meta-log2FC figure to", figure_dir, "\n")
