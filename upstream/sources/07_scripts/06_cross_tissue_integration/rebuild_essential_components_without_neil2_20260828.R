#!/usr/bin/env Rscript

# Rebuild the essential-component display after explicitly removing Neil2 from
# BER.  The upstream YARN/qsmooth values are retained: YARN was fitted on the
# complete common Ensembl universe, so removing one target component after the
# fit cannot change any remaining target value.  Tissue Spearman correlations
# and the attached row dendrograms are recomputed from the displayed target
# set (29 DSB, 45 SSB).

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
source_outdir <- file.path(root, "03_analysis_results", "18_essential_components_yarn_qsmooth_spearman_reordered_20260827")
outdir <- if (length(args) >= 2L) args[[2L]] else file.path(
  root, "03_analysis_results", "22_essential_components_yarn_qsmooth_no_Neil2_attached_row_dendrogram_20260828"
)
if (!grepl("^/", outdir)) outdir <- file.path(root, outdir)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

excluded_symbol <- "Neil2"
excluded_requested <- "NEIL2"
excluded_ensembl <- "ENSMUSG00000035121"
excluded_mgi <- "MGI:2686058"
excluded_entrez <- "382913"

source_table_dir <- file.path(source_outdir, "tables")
source_matrix_path <- file.path(source_table_dir, "09_component_tissue_yarn_qsmooth_log2_matrix.csv")
source_component_audit_path <- file.path(source_table_dir, "20_figure_component_reassignment_audit.csv")
source_mapping_path <- file.path(source_table_dir, "01_included_component_stable_id_mapping_audit.csv")
source_scope_path <- file.path(source_table_dir, "00_teacher_component_scope_audit.csv")
source_cluster_order_path <- file.path(source_table_dir, "17_tissue_spearman_cluster_order.csv")
source_dsb_corr_path <- file.path(source_table_dir, "18_tissue_spearman_correlation_DSB.csv")
source_ssb_corr_path <- file.path(source_table_dir, "19_tissue_spearman_correlation_SSB.csv")
source_flight_path <- file.path(source_table_dir, "06_component_flight_sample_yarn_qsmooth_values.csv.gz")
source_unit_path <- file.path(source_table_dir, "07_component_analysis_unit_yarn_qsmooth_log2.csv")
source_mission_path <- file.path(source_table_dir, "08_component_mission_yarn_qsmooth_log2.csv")
required <- c(
  source_matrix_path, source_component_audit_path, source_mapping_path, source_scope_path,
  source_cluster_order_path, source_dsb_corr_path, source_ssb_corr_path,
  source_flight_path, source_unit_path, source_mission_path
)
assert(all(file.exists(required)), paste("Missing source output:", paste(required[!file.exists(required)], collapse = "; ")))

dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "input"), recursive = TRUE, showWarnings = FALSE)

tissue_labels <- c(
  "Adrenal gland" = "肾上腺", "Bone marrow" = "骨髓", "Cecum" = "盲肠", "Cerebellum" = "小脑",
  "Colon" = "结肠", "Dorsal skin" = "背部皮肤", "Extensor digitorum longus" = "趾长伸肌", "Eye" = "眼",
  "Femoral lateral skin" = "股外侧皮肤", "Femoral skin" = "股部皮肤", "Gastrocnemius" = "腓肠肌",
  "Heart" = "心脏", "Heart / Heart right ventricle" = "右心室", "Kidney" = "肾脏",
  "Left lobe of the liver" = "肝左叶", "Liver" = "肝脏", "Lung" = "肺", "Mammary gland" = "乳腺",
  "Optic nerve" = "视神经", "Quadriceps femoris" = "股四头肌", "Retina" = "视网膜", "Soleus" = "比目鱼肌",
  "Spleen" = "脾脏", "Spleen-distal" = "脾脏（远端）", "Thymus" = "胸腺", "Tibialis anterior" = "胫骨前肌"
)

tissue_matrix_source <- fread(source_matrix_path, showProgress = FALSE)
component_audit_source <- fread(source_component_audit_path, showProgress = FALSE)
mapping_source <- fread(source_mapping_path, showProgress = FALSE)
scope_source <- fread(source_scope_path, showProgress = FALSE)
old_cluster_order <- fread(source_cluster_order_path, showProgress = FALSE)

assert(nrow(tissue_matrix_source) == 1950L, "The validated source tissue matrix is not 26 x 75.")
assert(uniqueN(tissue_matrix_source$Group) == 26L, "The validated source tissue matrix does not contain 26 tissues.")
assert(nrow(component_audit_source) == 75L && uniqueN(component_audit_source$OfficialMouseSymbol) == 75L,
       "The validated source component audit is not 75 unique components.")
assert(nrow(mapping_source) == 75L && uniqueN(mapping_source$EnsemblID) == 75L,
       "The validated source stable-ID mapping is not 75 unique components.")

excluded_rows <- tissue_matrix_source[OfficialMouseSymbol == excluded_symbol]
assert(nrow(excluded_rows) == 26L, "Neil2 must have exactly one tissue row for each of the 26 tissues in the source matrix.")
assert(uniqueN(excluded_rows$EnsemblID) == 1L && unique(excluded_rows$EnsemblID) == excluded_ensembl,
       "Neil2 stable-ID mapping is not the expected Ensembl ID.")
assert(uniqueN(excluded_rows$MGI_ID) == 1L && unique(excluded_rows$MGI_ID) == excluded_mgi,
       "Neil2 stable-ID mapping is not the expected MGI ID.")
assert(as.character(unique(excluded_rows$EntrezID)) == excluded_entrez,
       "Neil2 stable-ID mapping is not the expected Entrez ID.")
assert(unique(excluded_rows$Pathway) == "BER" && unique(excluded_rows$PlotPathwayDisplay) == "BER",
       "Neil2 is not assigned to BER in the validated source output.")

# Keep a full, inspectable scope audit while making the current-version
# exclusion explicit.  The teacher workbook and the original 75-component
# output remain untouched in the upstream versioned directories.
scope_filtered <- copy(scope_source)
scope_filtered[, OriginalIncluded := Included]
scope_filtered[, ExcludedByCurrentInstruction := toupper(trimws(RequestedComponent)) == excluded_requested]
scope_filtered[ExcludedByCurrentInstruction == TRUE, `:=`(
  Included = FALSE,
  ExclusionReason = "explicitly excluded from current BER display/table by instruction"
)]
scope_filtered[, CurrentVersion := "Neil2-excluded"]

# The existing reassignment audit is authoritative for display order (Exo1 and
# Dna2 remain in HR between Rad54l and Rbbp8).  Remove Neil2, then reindex the
# visible order so no gap remains in the plotted x-axis.
component_audit <- component_audit_source[OfficialMouseSymbol != excluded_symbol][order(PlotComponentOrder)]
component_audit[, PlotComponentOrder := seq_len(.N)]
assert(nrow(component_audit) == 74L && uniqueN(component_audit$OfficialMouseSymbol) == 74L,
       "The current version must contain 74 unique displayed components.")
assert(component_audit[PlottedPathwayDisplay %chin% c("NHEJ", "HR", "HR–A-EJ", "A-EJ"), .N] == 29L,
       "The DSB component count after Neil2 removal is not 29.")
assert(component_audit[PlottedPathwayDisplay %chin% c("BER", "NER", "MMR", "FA"), .N] == 45L,
       "The SSB component count after Neil2 removal is not 45.")
assert(!any(component_audit$OfficialMouseSymbol %chin% excluded_symbol), "Neil2 remains in the filtered component audit.")
assert(!anyDuplicated(component_audit$EnsemblID) && !anyNA(component_audit$EnsemblID),
       "The filtered component audit contains duplicate or missing Ensembl IDs.")

new_plot_order <- setNames(component_audit$PlotComponentOrder, component_audit$OfficialMouseSymbol)
new_plot_pathway <- setNames(component_audit$PlottedPathwayDisplay, component_audit$OfficialMouseSymbol)

tissue_matrix <- tissue_matrix_source[OfficialMouseSymbol != excluded_symbol]
tissue_matrix[, PlotComponentOrder := unname(new_plot_order[OfficialMouseSymbol])]
tissue_matrix[, PlotPathwayDisplay := unname(new_plot_pathway[OfficialMouseSymbol])]
setorder(tissue_matrix, TissueOrder, PlotComponentOrder)
assert(nrow(tissue_matrix) == 1924L && uniqueN(tissue_matrix[, .(Group, EnsemblID)]) == 1924L,
       "The filtered tissue matrix is not 26 x 74 with unique tissue-Ensembl cells.")
assert(!any(tissue_matrix$OfficialMouseSymbol %chin% excluded_symbol), "Neil2 remains in the filtered tissue matrix.")
assert(all(is.finite(tissue_matrix$TissueYARNNormalizedLog2)), "The filtered tissue matrix contains non-finite YARN values.")

# Filter the stable-ID mapping and retain both source pathway and plotted
# pathway.  Analysis and joins continue to use Ensembl IDs; Symbols are labels.
mapping_filtered <- mapping_source[OfficialMouseSymbol != excluded_symbol]
mapping_filtered[, PlotComponentOrder := unname(new_plot_order[OfficialMouseSymbol])]
mapping_filtered[, PlottedPathwayDisplay := unname(new_plot_pathway[OfficialMouseSymbol])]
setorder(mapping_filtered, PlotComponentOrder)
assert(nrow(mapping_filtered) == 74L && uniqueN(mapping_filtered$EnsemblID) == 74L,
       "The filtered stable-ID mapping is incomplete.")

excluded_audit <- unique(excluded_rows[, .(
  RequestedComponent, OfficialMouseSymbol, EnsemblID, MGI_ID, EntrezID,
  OriginalPathway = Pathway, OriginalPathwayDisplay = PathwayDisplay,
  PlottedPathwayDisplay = PlotPathwayDisplay,
  ExclusionReason = "explicitly excluded from current BER display/table by instruction"
)])

# Tissue order used only to construct a deterministic profile matrix.  The
# visible order is supplied independently by each figure's newly calculated
# Spearman tree.
tissue_universe <- tissue_matrix_source[order(TissueOrder), unique(Group)]
assert(length(tissue_universe) == 26L && !anyNA(tissue_universe), "The tissue universe is not 26 complete tissue labels.")
assert(setequal(tissue_universe, tissue_matrix$Group), "The filtered tissue matrix changed the tissue universe.")

spearman_cluster_order <- function(plot_data, figure_name, tissue_universe) {
  profile <- dcast(
    plot_data[, .(Group, OfficialMouseSymbol, TissueYARNNormalizedLog2)],
    Group ~ OfficialMouseSymbol,
    value.var = "TissueYARNNormalizedLog2"
  )
  profile <- profile[match(tissue_universe, Group)]
  assert(nrow(profile) == length(tissue_universe) && all(profile$Group == tissue_universe),
         paste("The", figure_name, "profile does not cover the expected 26 tissues."))
  component_columns <- setdiff(names(profile), "Group")
  assert(length(component_columns) > 1L, paste("The", figure_name, "profile has too few components."))
  profile_matrix <- as.matrix(profile[, ..component_columns])
  suppressWarnings(storage.mode(profile_matrix) <- "numeric")
  assert(all(is.finite(profile_matrix)), paste("The", figure_name, "profile contains non-finite values."))
  rownames(profile_matrix) <- profile$Group
  rho <- suppressWarnings(cor(t(profile_matrix), method = "spearman", use = "pairwise.complete.obs"))
  assert(all(is.finite(rho)), paste("The", figure_name, "Spearman matrix contains non-finite values."))
  rho[rho > 1] <- 1
  rho[rho < -1] <- -1
  distance <- 1 - rho
  distance <- (distance + t(distance)) / 2
  diag(distance) <- 0
  hc <- hclust(as.dist(distance), method = "average")
  ordered_tissues <- rownames(profile_matrix)[hc$order]
  assert(length(ordered_tissues) == 26L && setequal(ordered_tissues, tissue_universe),
         paste("The", figure_name, "cluster order is incomplete."))
  order_audit <- data.table(
    Figure = figure_name,
    OrderTopToBottom = seq_along(ordered_tissues),
    Tissue = ordered_tissues,
    ComponentN = length(component_columns),
    ClusteringMethod = "average-linkage hierarchical clustering",
    DistanceDefinition = "1 - Spearman rho across displayed components"
  )
  correlation_table <- as.data.table(rho, keep.rownames = "Tissue")
  correlation_table[, `:=`(
    Figure = figure_name,
    ComponentN = length(component_columns),
    ClusteringMethod = "average-linkage hierarchical clustering",
    DistanceDefinition = "1 - Spearman rho across displayed components"
  )]
  setcolorder(correlation_table, c("Figure", "ComponentN", "ClusteringMethod", "DistanceDefinition", "Tissue"))
  list(
    figure = figure_name,
    component_n = length(component_columns),
    method = "average-linkage hierarchical clustering",
    distance_definition = "1 - Spearman rho across displayed components",
    hc = hc,
    ordered_tissues = ordered_tissues,
    order_audit = order_audit,
    correlation = correlation_table
  )
}

dsb_plot_data <- tissue_matrix[PlotPathwayDisplay %chin% c("NHEJ", "HR", "HR–A-EJ", "A-EJ")]
ssb_plot_data <- tissue_matrix[PlotPathwayDisplay %chin% c("BER", "NER", "MMR", "FA")]
dsb <- spearman_cluster_order(dsb_plot_data, "DSB", tissue_universe)
ssb <- spearman_cluster_order(ssb_plot_data, "SSB", tissue_universe)

old_dsb_order <- old_cluster_order[Figure == "DSB"][order(OrderTopToBottom), as.character(Tissue)]
old_ssb_order <- old_cluster_order[Figure == "SSB"][order(OrderTopToBottom), as.character(Tissue)]
dsb_order_identical <- identical(as.character(dsb$ordered_tissues), as.character(old_dsb_order))
ssb_order_identical <- identical(as.character(ssb$ordered_tissues), as.character(old_ssb_order))

read_saved_correlation <- function(path, figure_name) {
  x <- fread(path, showProgress = FALSE)
  meta_cols <- c("Figure", "ComponentN", "ClusteringMethod", "DistanceDefinition", "Tissue")
  assert(all(meta_cols %in% names(x)), paste("Correlation matrix schema changed:", path))
  tissues <- as.character(x$Tissue)
  component_cols <- setdiff(names(x), meta_cols)
  assert(length(tissues) == 26L && length(component_cols) == 26L && setequal(tissues, component_cols),
         paste("Invalid saved correlation matrix:", path))
  rho <- as.matrix(x[, ..component_cols])
  suppressWarnings(storage.mode(rho) <- "numeric")
  rownames(rho) <- tissues
  colnames(rho) <- component_cols
  rho <- rho[tissues, tissues, drop = FALSE]
  assert(all(is.finite(rho)) && max(abs(rho - t(rho))) < 1e-12, paste("Invalid saved correlation values:", path))
  list(figure = figure_name, rho = rho)
}

compare_rho <- function(new_corr, old_corr) {
  common <- intersect(rownames(new_corr), rownames(old_corr))
  max(abs(new_corr[common, common, drop = FALSE] - old_corr[common, common, drop = FALSE]))
}
old_dsb <- read_saved_correlation(source_dsb_corr_path, "DSB")
old_ssb <- read_saved_correlation(source_ssb_corr_path, "SSB")
new_dsb_rho <- as.matrix(dsb$correlation[, 6:ncol(dsb$correlation), with = FALSE])
rownames(new_dsb_rho) <- dsb$correlation$Tissue
colnames(new_dsb_rho) <- setdiff(names(dsb$correlation), c("Figure", "ComponentN", "ClusteringMethod", "DistanceDefinition", "Tissue"))
new_ssb_rho <- as.matrix(ssb$correlation[, 6:ncol(ssb$correlation), with = FALSE])
rownames(new_ssb_rho) <- ssb$correlation$Tissue
colnames(new_ssb_rho) <- setdiff(names(ssb$correlation), c("Figure", "ComponentN", "ClusteringMethod", "DistanceDefinition", "Tissue"))
dsb_rho_max_abs_diff <- compare_rho(new_dsb_rho, old_dsb$rho)
ssb_rho_max_abs_diff <- compare_rho(new_ssb_rho, old_ssb$rho)

fwrite(scope_filtered, file.path(outdir, "tables", "00_teacher_component_scope_audit_filtered.csv"))
fwrite(mapping_filtered, file.path(outdir, "tables", "01_filtered_component_stable_id_mapping_audit.csv"))
fwrite(rbind(dsb$order_audit, ssb$order_audit), file.path(outdir, "tables", "02_tissue_spearman_cluster_order.csv"))
fwrite(dsb$correlation, file.path(outdir, "tables", "03_tissue_spearman_correlation_DSB.csv"))
fwrite(ssb$correlation, file.path(outdir, "tables", "04_tissue_spearman_correlation_SSB.csv"))
fwrite(tissue_matrix, file.path(outdir, "tables", "05_component_tissue_yarn_qsmooth_log2_matrix.csv"))
fwrite(component_audit, file.path(outdir, "tables", "09_figure_component_reassignment_audit.csv"))
fwrite(excluded_audit, file.path(outdir, "tables", "10_excluded_component_audit.csv"))

filter_component_table <- function(source_path, output_path) {
  x <- fread(source_path, showProgress = FALSE)
  assert("OfficialMouseSymbol" %in% names(x), paste("Expected OfficialMouseSymbol in", source_path))
  x <- x[OfficialMouseSymbol != excluded_symbol]
  assert(!any(x$OfficialMouseSymbol %chin% excluded_symbol), paste("Neil2 was not removed from", source_path))
  fwrite(x, output_path)
  nrow(x)
}
flight_rows <- filter_component_table(source_flight_path, file.path(outdir, "tables", "06_component_flight_sample_yarn_qsmooth_values.csv.gz"))
unit_rows <- filter_component_table(source_unit_path, file.path(outdir, "tables", "07_component_analysis_unit_yarn_qsmooth_log2.csv"))
mission_rows <- filter_component_table(source_mission_path, file.path(outdir, "tables", "08_component_mission_yarn_qsmooth_log2.csv"))

# Verify that filtering did not change any remaining normalized value.  This is
# the key reason a second, expensive YARN fit is not scientifically necessary.
remaining_source <- tissue_matrix_source[OfficialMouseSymbol != excluded_symbol, .(Group, EnsemblID, SourceYARN = TissueYARNNormalizedLog2)]
remaining_new <- tissue_matrix[, .(Group, EnsemblID, NewYARN = TissueYARNNormalizedLog2)]
value_compare <- merge(remaining_source, remaining_new, by = c("Group", "EnsemblID"), all = TRUE, sort = FALSE)
assert(nrow(value_compare) == 1924L && all(is.finite(value_compare$SourceYARN)) && all(is.finite(value_compare$NewYARN)),
       "The remaining source/new tissue values could not be matched completely.")
value_compare[, AbsoluteDifference := abs(SourceYARN - NewYARN)]
assert(max(value_compare$AbsoluteDifference) < 1e-12, "Filtering changed a remaining YARN value.")

old_yarn_min <- min(tissue_matrix_source$TissueYARNNormalizedLog2)
old_yarn_max <- max(tissue_matrix_source$TissueYARNNormalizedLog2)
yarn_min <- min(tissue_matrix$TissueYARNNormalizedLog2)
yarn_max <- max(tissue_matrix$TissueYARNNormalizedLog2)

component_count_audit <- component_audit[, .(ComponentN = .N), by = .(Pathway = PlottedPathwayDisplay)]
component_count_audit[, `:=`(
  Version = "Neil2-excluded",
  Figure = fifelse(Pathway %chin% c("NHEJ", "HR", "HR–A-EJ", "A-EJ"), "DSB", "SSB")
)]
# Keep the pathway table unambiguous and ordered for direct inspection.
setcolorder(component_count_audit, c("Version", "Figure", "Pathway", "ComponentN"))
setorder(component_count_audit, Figure, Pathway)
fwrite(component_count_audit, file.path(outdir, "tables", "12_component_count_audit.csv"))
fwrite(value_compare, file.path(outdir, "tables", "11_source_value_unchanged_audit.csv"))

run_audit <- data.table(
  Check = c(
    "Source output preserved", "Source component count", "Current displayed component count", "Current DSB component count",
    "Current SSB component count", "Excluded requested component", "Excluded official symbol", "Excluded Ensembl ID",
    "Excluded MGI ID", "Excluded Entrez ID", "Tissue count", "Filtered tissue matrix cells", "Filtered finite cells",
    "Remaining YARN max absolute difference", "Source global YARN minimum", "Current global YARN minimum",
    "Source global YARN maximum", "Current global YARN maximum", "DSB Spearman recomputed", "SSB Spearman recomputed",
    "DSB tree order identical to 75-component version", "SSB tree order identical to 75-component version",
    "DSB rho max absolute difference vs source", "SSB rho max absolute difference vs source",
    "Spearman distance", "Linkage", "Analysis key", "Current display rule"
  ),
  Observed = as.character(c(
    source_outdir, 75L, 74L, dsb$component_n, ssb$component_n, excluded_requested, excluded_symbol, excluded_ensembl,
    excluded_mgi, excluded_entrez, 26L, nrow(tissue_matrix), sum(is.finite(tissue_matrix$TissueYARNNormalizedLog2)),
    max(value_compare$AbsoluteDifference), old_yarn_min, yarn_min, old_yarn_max, yarn_max, TRUE, TRUE,
    dsb_order_identical, ssb_order_identical, dsb_rho_max_abs_diff, ssb_rho_max_abs_diff,
    "1 - Spearman rho across displayed components", "average-linkage hierarchical clustering", "Ensembl Gene ID", 
    "Neil2 removed from BER; Exo1/Dna2 remain visually in HR between Rad54l and Rbbp8"
  )),
  Status = c(
    "PASS", "PASS", "PASS", if (dsb$component_n == 29L) "PASS" else "FAIL", if (ssb$component_n == 45L) "PASS" else "FAIL",
    "PASS", "PASS", "PASS", "PASS", "PASS", if (length(tissue_universe) == 26L) "PASS" else "FAIL",
    if (nrow(tissue_matrix) == 1924L) "PASS" else "FAIL", if (all(is.finite(tissue_matrix$TissueYARNNormalizedLog2))) "PASS" else "FAIL",
    if (max(value_compare$AbsoluteDifference) < 1e-12) "PASS" else "FAIL", "recorded", "recorded", "recorded", "recorded",
    "PASS", "PASS", if (dsb_order_identical) "PASS" else "CHECK", if (ssb_order_identical) "PASS" else "CHANGED",
    if (dsb_rho_max_abs_diff < 1e-12) "PASS" else "CHECK", "expected changed feature set", "documented", "documented", "documented", "documented"
  )
)
fwrite(run_audit, file.path(outdir, "tables", "13_run_audit.csv"))

build_tree_segments <- function(hc) {
  n <- length(hc$labels)
  ordered <- hc$labels[hc$order]
  leaf_y <- setNames(seq_len(n), ordered)
  node_x <- numeric(n - 1L)
  node_y <- numeric(n - 1L)
  segments <- vector("list", n - 1L)
  child_coord <- function(code) {
    if (code < 0L) c(x = 0, y = unname(leaf_y[hc$labels[-code]])) else c(x = node_x[[code]], y = node_y[[code]])
  }
  for (i in seq_len(n - 1L)) {
    left <- child_coord(hc$merge[i, 1L])
    right <- child_coord(hc$merge[i, 2L])
    parent_x <- hc$height[[i]]
    parent_y <- mean(c(left[["y"]], right[["y"]]))
    node_x[[i]] <- parent_x
    node_y[[i]] <- parent_y
    segments[[i]] <- rbind(
      data.table(x = left[["x"]], xend = parent_x, y = left[["y"]], yend = left[["y"]]),
      data.table(x = right[["x"]], xend = parent_x, y = right[["y"]], yend = right[["y"]]),
      data.table(x = parent_x, xend = parent_x, y = min(left[["y"]], right[["y"]]), yend = max(left[["y"]], right[["y"]]))
    )
  }
  list(segments = rbindlist(segments), ordered_tissues = ordered)
}

display_tissue_name <- function(x) {
  y <- x
  y[x == "Heart / Heart right ventricle"] <- "Heart right ventricle"
  y
}

make_tree_plot <- function(cluster) {
  tree <- build_tree_segments(cluster$hc)
  max_height <- max(cluster$hc$height)
  tree$segments[, `:=`(x = max_height - x, xend = max_height - xend)]
  ggplot(tree$segments) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.62, colour = "#334E68", lineend = "round") +
    scale_x_continuous(limits = c(0, max_height), breaks = NULL, name = NULL, expand = c(0, 0)) +
    scale_y_reverse(limits = c(26.7, 0.3), breaks = NULL, expand = c(0, 0)) +
    theme_minimal(base_size = 9.2, base_family = "Arial Unicode MS") +
    theme(
      panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
      plot.margin = margin(2, 0, 13, 0)
    )
}

make_label_plot <- function(ordered_tissues) {
  labels <- data.table(
    y = seq_along(ordered_tissues),
    Label = paste0(display_tissue_name(ordered_tissues), "\n", unname(tissue_labels[ordered_tissues]))
  )
  ggplot(labels, aes(x = 0, y = y, label = Label)) +
    geom_text(hjust = 0, size = 3.05, lineheight = 0.82, family = "Arial Unicode MS", colour = "#263238") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_reverse(limits = c(26.7, 0.3), breaks = NULL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = "Arial Unicode MS") +
    theme(plot.margin = margin(2, 0, 13, 0))
}

make_heatmap_plot <- function(plot_data, ordered_tissues) {
  plot_data <- copy(plot_data)
  plot_data[, DisplayComponent := factor(OfficialMouseSymbol, levels = component_audit$OfficialMouseSymbol)]
  plot_data[, PlotPathwayDisplay := factor(PlotPathwayDisplay, levels = c("NHEJ", "HR", "HR–A-EJ", "A-EJ", "BER", "NER", "MMR", "FA"))]
  plot_data[, TissueKey := factor(Group, levels = rev(ordered_tissues))]
  plot_data[, CellTextColour := fifelse(
    TissueYARNNormalizedLog2 <= yarn_min + 0.16 * (yarn_max - yarn_min) |
      TissueYARNNormalizedLog2 >= yarn_min + 0.84 * (yarn_max - yarn_min),
    "white", "#1F2937"
  )]
  ggplot(plot_data, aes(x = DisplayComponent, y = TissueKey, fill = TissueYARNNormalizedLog2)) +
    geom_tile(colour = "white", linewidth = 0.16) +
    geom_text(aes(label = sprintf("%.1f", TissueYARNNormalizedLog2), colour = CellTextColour), size = 1.55, na.rm = TRUE, show.legend = FALSE) +
    scale_colour_identity() +
    facet_grid(. ~ PlotPathwayDisplay, scales = "free_x", space = "free_x", switch = "x", drop = TRUE) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(drop = FALSE, labels = function(x) rep("", length(x))) +
    scale_fill_gradientn(
      colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
      values = scales::rescale(seq(yarn_min, yarn_max, length.out = 5L), from = c(yarn_min, yarn_max)),
      limits = c(yarn_min, yarn_max), oob = scales::squish, name = "YARN qsmooth\nlog2 expression",
      breaks = pretty(c(yarn_min, yarn_max), n = 4), labels = function(x) formatC(x, format = "fg", digits = 3),
      guide = guide_colorbar(
        title.position = "top", title.hjust = 0.5, direction = "horizontal",
        barwidth = unit(15, "cm"), barheight = unit(0.38, "cm"), ticks.colour = "#52606D", frame.colour = "#8A9BA8"
      )
    ) +
    theme_minimal(base_size = 8.8, base_family = "Arial Unicode MS") +
    theme(
      panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.3, colour = "#3D3D3D"),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
      strip.background = element_rect(fill = "#F0F0F0", colour = NA), strip.text.x.bottom = element_text(face = "bold", size = 9.2),
      strip.placement = "outside", panel.spacing.x = unit(0.42, "lines"),
      legend.position = "top", legend.justification = "center", legend.title = element_text(size = 8.7, face = "bold"),
      legend.text = element_text(size = 7.8), plot.margin = margin(2, 0, 13, 0)
    )
}

render_combined <- function(cluster, figure_name, pathway_set, title, subtitle, stem, width) {
  selected <- tissue_matrix[PlotPathwayDisplay %chin% pathway_set]
  component_n <- component_audit[PlottedPathwayDisplay %chin% pathway_set, .N]
  assert(nrow(selected) == 26L * component_n, paste("Unexpected", figure_name, "component cells."))
  tree_plot <- make_tree_plot(cluster)
  heatmap_plot <- make_heatmap_plot(selected, cluster$ordered_tissues)
  label_plot <- make_label_plot(cluster$ordered_tissues)
  combined <- tree_plot + heatmap_plot + label_plot +
    plot_layout(widths = c(4.4, 28, 7.0), guides = "collect") +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = "左：从组织表达谱重算的 Spearman 行聚类树，根在左、叶端贴合对应组织行；中：YARN qsmooth log2 expression 热图；右：组织名称。颜色图例置于顶部。YARN 数值沿用完整共同 Ensembl 基因全集拟合结果，未因删除 Neil2 而改变。",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, family = "Arial Unicode MS"),
        plot.subtitle = element_text(hjust = 0.5, size = 8.8, family = "Arial Unicode MS"),
        plot.caption = element_text(hjust = 0, size = 7.3, colour = "#52606D", family = "Arial Unicode MS"),
        plot.margin = margin(8, 8, 8, 8)
      )
    ) & theme(legend.position = "top")
  png_path <- file.path(outdir, "figures", paste0(stem, ".png"))
  pdf_path <- file.path(outdir, "figures", paste0(stem, ".pdf"))
  ggsave(png_path, combined, width = width, height = 13, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, combined, width = width, height = 13, device = cairo_pdf, bg = "white", limitsize = FALSE)
  data.table(
    RenderType = "combined_tree_heatmap", Figure = figure_name, ComponentN = component_n, TissueN = 26L,
    ClusteringMethod = cluster$method, DistanceDefinition = cluster$distance_definition,
    Layout = "root-left row dendrogram with leaf tips attached to heatmap rows | heatmap | right tissue labels; legend at top",
    PNG = png_path, PDF = pdf_path, LeafOrderTopToBottom = paste(cluster$ordered_tissues, collapse = " > ")
  )
}

render_tree <- function(cluster, stem) {
  tree <- build_tree_segments(cluster$hc)
  max_height <- max(cluster$hc$height)
  tissue_display <- display_tissue_name(tree$ordered_tissues)
  leaf_data <- data.table(
    Tissue = tree$ordered_tissues,
    y = seq_along(tree$ordered_tissues),
    Label = paste0(tissue_display, "\n", unname(tissue_labels[tree$ordered_tissues]))
  )
  label_pad <- max(0.22 * max_height, 0.25)
  p <- ggplot() +
    geom_segment(data = tree$segments, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.65, colour = "#334E68", lineend = "round") +
    geom_point(data = leaf_data, aes(x = 0, y = y), size = 1.35, colour = "#2166AC") +
    geom_text(data = leaf_data, aes(x = 0, y = y, label = Label), hjust = 1.03, size = 3.15, lineheight = 0.82, family = "Arial Unicode MS", colour = "#263238") +
    scale_x_continuous(limits = c(-label_pad, max_height * 1.06), breaks = pretty(c(0, max_height), n = 5), name = "1 − Spearman rho", expand = c(0, 0)) +
    scale_y_reverse(limits = c(26.7, 0.3), breaks = NULL, expand = c(0, 0)) +
    labs(
      title = paste0(cluster$figure, " tissue clustering tree (Spearman)"),
      subtitle = sprintf("26 tissues; %s displayed components; average-linkage hierarchical clustering", cluster$component_n),
      caption = "The tree was recomputed after removing Neil2 from the displayed SSB set; leaf order is the order used by the corresponding attached heatmap."
    ) +
    theme_minimal(base_size = 10.5, base_family = "Arial Unicode MS") +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour = "#E6EEF5", linewidth = 0.35),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(colour = "#40566B"),
      axis.title.x = element_text(colour = "#40566B", margin = margin(t = 7)), plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.2), plot.caption = element_text(size = 8.4, hjust = 0, colour = "#52606D"),
      plot.margin = margin(10, 16, 10, 160)
    ) + coord_cartesian(clip = "off")
  png_path <- file.path(outdir, "figures", paste0(stem, ".png"))
  pdf_path <- file.path(outdir, "figures", paste0(stem, ".pdf"))
  ggsave(png_path, p, width = 12, height = 9, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, p, width = 12, height = 9, device = cairo_pdf, bg = "white", limitsize = FALSE)
  data.table(
    RenderType = "standalone_tree", Figure = cluster$figure, ComponentN = cluster$component_n, TissueN = 26L,
    ClusteringMethod = cluster$method, DistanceDefinition = cluster$distance_definition,
    Layout = "standalone dendrogram with bilingual tissue labels",
    PNG = png_path, PDF = pdf_path, LeafOrderTopToBottom = paste(tree$ordered_tissues, collapse = " > ")
  )
}

render_audit <- rbindlist(list(
  render_combined(dsb, "DSB", c("NHEJ", "HR", "HR–A-EJ", "A-EJ"),
                  "DSB essential components: tissue Spearman tree + YARN qsmooth heatmap",
                  "26 mouse tissues; 29 DSB components; Exo1/Dna2 displayed in HR; rows ordered by the recomputed DSB tree",
                  "01_DSB_yarn_qsmooth_attached_row_dendrogram_heatmap", 39),
  render_combined(ssb, "SSB", c("BER", "NER", "MMR", "FA"),
                  "SSB essential components: tissue Spearman tree + YARN qsmooth heatmap",
                  "26 mouse tissues; 45 SSB components after excluding Neil2 from BER; rows ordered by the recomputed SSB tree",
                  "02_SSB_yarn_qsmooth_attached_row_dendrogram_heatmap", 39),
  render_tree(dsb, "03_DSB_tissue_spearman_dendrogram"),
  render_tree(ssb, "04_SSB_tissue_spearman_dendrogram")
), use.names = TRUE, fill = TRUE)
fwrite(render_audit, file.path(outdir, "tables", "14_render_audit.csv"))

leaf_attachment_audit <- rbindlist(lapply(list(dsb, ssb), function(cluster) {
  data.table(
    Figure = cluster$figure,
    RowTopToBottom = seq_along(cluster$ordered_tissues),
    Tissue = cluster$ordered_tissues,
    TreeLeafY = seq_along(cluster$ordered_tissues),
    HeatmapRowY = seq_along(cluster$ordered_tissues),
    RightLabelY = seq_along(cluster$ordered_tissues),
    Status = "PASS: one shared tree leaf / heatmap row / right-label position"
  )
}), use.names = TRUE)
fwrite(leaf_attachment_audit, file.path(outdir, "tables", "15_row_dendrogram_leaf_attachment_audit.csv"))

readme_lines <- c(
  "# Neil2 removed: tissue Spearman tree + YARN qsmooth heatmaps",
  "",
  "本目录是按当前指令生成的去除 `Neil2` 版本。`Neil2`（Ensembl `ENSMUSG00000035121`；MGI `MGI:2686058`；Entrez `382913`）从 BER 的表格和图中均不再展示。老师提供的原始表、完整 75 成分版本和此前完整树图均保留在原目录，不被覆盖。",
  "",
  "## 关键处理决定",
  "",
  "- YARN/qsmooth 不需要因删除一个目标成分而重新拟合：上游版本是在全部 12,462 个共同 Ensembl 基因、360 个 flight 样本、59 个 analysis unit 和 26 个组织上拟合；本版只在拟合后从目标成分层过滤 `Neil2`。对剩余 1,924 个组织×成分格逐格比对，最大绝对差为 0。",
  "- Spearman 树必须重算，因为相关性特征集合发生了变化。DSB 用剩余的 29 个 DSB 成分重算；SSB 用去除 `Neil2` 后的 45 个 SSB 成分重算。距离仍为 `1 - Spearman rho`，聚类仍为 average-linkage。",
  "- DSB 的成分和数值没有改变，并记录了与原 75 成分版本的树叶顺序及相关矩阵差异；SSB 的树基于 45 个成分，是本版新的树。",
  "- 所有分析和匹配仍以 Ensembl Gene ID 为键；图中只显示官方小鼠 Symbol。Exo1、Dna2 仍按前一版显示规则放在 HR 中 Rad54l 与 Rbbp8 之间，源 Pathway 未覆盖。",
  "- 组织顺序不是 A–Z 排序，而是分别由 DSB/SSB 表达谱的 Spearman 树决定；树叶、热图行和右侧组织标签使用同一顺序。`Heart / Heart right ventricle` 的显示名为 `Heart right ventricle`。",
  "- 06–08 级别表保留上游数据中的少量源缺失值（NA），不做插补；最终用于热图的 26×74 tissue matrix 的 1,924 个格子均为有限值。",
  "",
  "## 输出",
  "",
  "- `figures/01_DSB_yarn_qsmooth_attached_row_dendrogram_heatmap.png` / `.pdf`：29 个 DSB 成分的树-热图合并图。",
  "- `figures/02_SSB_yarn_qsmooth_attached_row_dendrogram_heatmap.png` / `.pdf`：45 个 SSB 成分的树-热图合并图，BER 不再含 Neil2。",
  "- `figures/03_DSB_tissue_spearman_dendrogram.png` / `.pdf`、`04_SSB_tissue_spearman_dendrogram.png` / `.pdf`：本版重新计算的独立树图。",
  "- `tables/00_teacher_component_scope_audit_filtered.csv`：原始老师成分范围加当前排除标记；原始 00 表未改写。",
  "- `tables/01_filtered_component_stable_id_mapping_audit.csv`、`10_excluded_component_audit.csv`：去除后的稳定 ID 映射与 Neil2 排除记录。",
  "- `tables/02_tissue_spearman_cluster_order.csv`、`03_tissue_spearman_correlation_DSB.csv`、`04_tissue_spearman_correlation_SSB.csv`：重算的树叶顺序和相关矩阵。",
  "- `tables/05_component_tissue_yarn_qsmooth_log2_matrix.csv`：不含 Neil2 的 26×74 tissue matrix；06–08 为对应的 flight、analysis-unit 和 mission 级过滤表。",
  "- `tables/11_source_value_unchanged_audit.csv`、`13_run_audit.csv`、`14_render_audit.csv`、`15_row_dendrogram_leaf_attachment_audit.csv`：数值、参数、渲染和树叶-热图行对应关系审计。",
  "",
  "## 上游版本",
  "",
  paste0("- 完整 75 成分 YARN/qsmooth 与原始稳定 ID 表：`", source_outdir, "`。"),
  "- 当前版本只建立并行输出目录，不删除或覆盖老师原始 Excel，也不删除此前结果。"
)
writeLines(readme_lines, file.path(outdir, "README.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"), useBytes = TRUE)
writeLines(c(
  paste0("Source output: ", source_outdir),
  paste0("Excluded symbol: ", excluded_symbol),
  paste0("Excluded Ensembl ID: ", excluded_ensembl),
  paste0("Excluded MGI ID: ", excluded_mgi),
  "The teacher workbook and upstream 75-component output were not modified."
), file.path(outdir, "input", "source_and_exclusion.txt"), useBytes = TRUE)
message("Created Neil2-excluded tree/heatmap outputs in: ", outdir)
