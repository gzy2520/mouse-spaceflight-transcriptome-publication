#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(png)
})

set.seed(25)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate tests/validate_final_release.R", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1L]] else file.path(root, "final_result")
if (length(args) < 2L || !args[[2L]] %in% c("full", "compact")) {
  stop("Usage: Rscript tests/validate_final_release.R <output_dir> <full|compact>", call. = FALSE)
}
validation_mode <- args[[2L]]
if (!grepl("^/", target)) target <- file.path(root, target)
target <- normalizePath(target, mustWork = TRUE)

assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
sha256_file <- function(path) {
  answer <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  assert(length(answer) > 0L, paste("Unable to hash", path))
  sub("[[:space:]].*$", "", answer[[1L]])
}

input_manifest <- fread(file.path(root, "provenance/publication_input_sha256.csv"))
input_paths <- file.path(root, "data/publication_input", input_manifest$relative_path)
assert(all(file.exists(input_paths)), "A frozen publication input is missing")
input_hashes <- vapply(input_paths, sha256_file, character(1L))
assert(identical(unname(input_hashes), unname(input_manifest$sha256)), paste(
  "Frozen publication input SHA-256 mismatch:",
  paste(input_manifest$relative_path[input_hashes != input_manifest$sha256], collapse = "; ")
))

main_stems <- paste0("Fig_", c("1b", "1c", "2", "3a", "3b", "4a", "4b", "5a", "5b"))
suppl_stems <- paste0("Fig_S", c("1", "2a", "2b", "3a", "3b", "4a", "4b", "5"))
expected_png <- c(file.path("Main", paste0(main_stems, ".png")),
                  file.path("Suppl", paste0(suppl_stems, ".png")))
expected_pdf <- sub("[.]png$", ".pdf", expected_png)

observed_figures <- sort(c(
  file.path("Main", list.files(file.path(target, "Main"), full.names = FALSE, all.files = FALSE)),
  file.path("Suppl", list.files(file.path(target, "Suppl"), full.names = FALSE, all.files = FALSE))
))
expected_figures <- sort(c(expected_png, if (validation_mode == "full") expected_pdf else character()))
assert(identical(observed_figures, expected_figures), paste0(
  "Final figure file set differs. Missing: ", paste(setdiff(expected_figures, observed_figures), collapse = "; "),
  " | Extra: ", paste(setdiff(observed_figures, expected_figures), collapse = "; ")
))
assert(all(file.info(file.path(target, expected_figures))$size > 1000), "A figure file is empty or truncated")

png_info <- lapply(file.path(target, expected_png), function(path) {
  attr(readPNG(path, native = TRUE, info = TRUE), "info")$dim
})
names(png_info) <- expected_png
assert(all(vapply(png_info, function(x) all(x >= 1000L), logical(1L))), "Unexpectedly small PNG dimensions")
fig3a_dim <- png_info[["Main/Fig_3a.png"]]
assert(fig3a_dim[[2L]] <= 2200L && fig3a_dim[[1L]] / fig3a_dim[[2L]] >= 3.2,
       "Fig. 3a is not the compressed half-height layout")

excluded_tables <- c("Fig_6_mouse_metadata_grouped.csv", "Fig_6_mouse_sample_metadata_complete_audit.csv")
source_tables <- list.files(file.path(root, "results/tables"), pattern = "[.]csv$", full.names = TRUE)
source_tables <- sort(source_tables[!basename(source_tables) %in% excluded_tables])
observed_tables <- sort(list.files(file.path(target, "tables"), pattern = "[.]csv$", full.names = FALSE))
assert(length(source_tables) == 27L && identical(observed_tables, sort(basename(source_tables))),
       "Final table file set is not the frozen 27-table set")
source_hash <- vapply(source_tables, sha256_file, character(1L))
target_hash <- vapply(file.path(target, "tables", basename(source_tables)), sha256_file, character(1L))
reference_manifest <- fread(file.path(root, "provenance/reference_sha256.csv"))
table_relative_paths <- file.path("tables", basename(source_tables))
frozen_table_hash <- reference_manifest$sha256[match(table_relative_paths, reference_manifest$relative_path)]
assert(!anyNA(frozen_table_hash), "The frozen reference manifest does not cover all 27 final tables")
assert(identical(unname(source_hash), unname(frozen_table_hash)), "A source manuscript table drifted from its frozen hash")
assert(identical(unname(target_hash), unname(frozen_table_hash)), "A copied final manuscript table changed")

provenance_expected <- sort(c(
  "Fig_1b_mouse_metadata_grouped.csv", "Fig_3b_column_label_audit.csv",
  "Fig_S5_scope_audit.csv", "GO_display_label_audit.csv",
  "GO_figure_scope_audit.csv", "expression_tree_structure_audit.csv",
  "figure_palette_audit.csv", "final_figure_sources.csv", "final_palette_contract.csv",
  "meta_dotplot_scope_audit.csv",
  "publication_table_sha256.csv", "upset_render_audit.csv"
))
observed_provenance <- sort(list.files(file.path(target, "provenance"), pattern = "[.]csv$"))
assert(identical(observed_provenance, provenance_expected), "Final provenance file set differs")

palette_contract <- fread(file.path(target, "provenance/final_palette_contract.csv"), na.strings = NULL)
palette_actual <- fread(file.path(target, "provenance/figure_palette_audit.csv"), na.strings = NULL)
required_palette_columns <- c("figure", "role", "hex")
assert(identical(names(palette_contract), required_palette_columns) &&
         identical(names(palette_actual), required_palette_columns),
       "Invalid final palette audit columns")
assert(!anyDuplicated(palette_contract[, .(figure, role)]) &&
         !anyDuplicated(palette_actual[, .(figure, role)]),
       "A figure palette role is duplicated")
assert(identical(sort(unique(palette_contract$figure)), sort(c(main_stems, suppl_stems))),
       "The retained palette contract does not cover all 17 figures")
palette_contract[, hex := toupper(hex)]
palette_actual[, hex := toupper(hex)]
setorder(palette_contract, figure, role)
setorder(palette_actual, figure, role)
assert(all(grepl("^#[0-9A-F]{6}$", palette_contract$hex)), "Invalid retained palette colour")
assert(identical(palette_actual, palette_contract),
       "The colours actually supplied to the final plotting scales differ from the retained palettes")

table_manifest <- fread(file.path(target, "provenance/publication_table_sha256.csv"))
assert(nrow(table_manifest) == 27L && !anyDuplicated(table_manifest$relative_path), "Invalid table hash manifest")
assert(all(table_manifest$sha256 == frozen_table_hash[match(table_manifest$relative_path, table_relative_paths)]),
       "Table hash manifest differs from copied tables")

go_source <- fread(file.path(root, "data/publication_input/go/04_tissue_statistics_concrete_terms_and_context.csv"))
go_labels <- fread(file.path(target, "provenance/GO_display_label_audit.csv"))
go_scope <- fread(file.path(target, "provenance/GO_figure_scope_audit.csv"))
assert(nrow(go_source) == 390L && uniqueN(go_source$analysis_tissue) == 26L && uniqueN(go_source$term_key) == 15L,
       "Frozen GO matrix scope changed")
assert(go_labels[term_order == 4L, display_name] == "Intrinsic apoptotic signaling",
       "Fourth GO display label is incorrect")
assert(go_labels[term_order == 7L, display_name] == "Telomeric region",
       "Seventh GO display label is incorrect")
assert(!grepl("pathway", go_labels[term_order == 4L, display_name], ignore.case = TRUE) &&
         !grepl("chromosome", go_labels[term_order == 7L, display_name], ignore.case = TRUE),
       "Removed GO display words remain")
assert(identical(go_scope$n_cells, c(390L, 120L)), "GO figure cell counts changed")

mouse <- fread(file.path(root, "data/publication_input/mouse_metadata/03_mouse_level_metadata.csv"))
mouse_grouped <- fread(file.path(target, "provenance/Fig_1b_mouse_metadata_grouped.csv"))
assert(nrow(mouse) == 761L && uniqueN(mouse$accession) == 48L && uniqueN(mouse$Group) == 26L,
       "Frozen mouse metadata scope changed")
assert(sum(mouse_grouped$n_samples) == 761L && uniqueN(mouse_grouped$Group) == 26L,
       "Fig. 1b grouped metadata scope changed")

fig3b_labels <- fread(file.path(target, "provenance/Fig_3b_column_label_audit.csv"))
expected_fig3b <- c(
  "OSD-289 MHU-1 / Flight-Ground",
  "OSD-289 MHU-2 / Flight-Ground (vivarium)",
  "OSD-421 / Flight-Ground",
  "OSD-515 / Flight-Ground"
)
assert(identical(fig3b_labels$display_label, expected_fig3b), "Fig. 3b display labels changed")

group_order <- fread(file.path(root, "data/publication_input/common_direction/06_group_column_order.csv"))
common_direction_specs <- data.table(
  group_id = c(
    "A_kidney", "B_thymus", "C_NES_lt_minus1_excl_Lung_Thymus",
    "D_NES_gt1_excl_Kidney"
  ),
  n_genes = c(12L, 30L, 11L, 5L)
)
for (row in seq_len(nrow(common_direction_specs))) {
  selected_group_id <- common_direction_specs$group_id[[row]]
  expected_n <- common_direction_specs$n_genes[[row]]
  matrix_path <- file.path(
    root, "data/publication_input/common_direction",
    paste0(selected_group_id, "_DDR_common_direction_top30_min_abs_log2fc_direction_ordered.csv")
  )
  matrix <- fread(matrix_path)
  expected_columns <- paste0(
    "log2fc__",
    group_order[group_order[["group_id"]] == selected_group_id][order(column_position), column_id]
  )
  assert(nrow(matrix) == expected_n && uniqueN(matrix$ensembl_id) == expected_n &&
           !anyNA(matrix$ensembl_id) && all(grepl("^ENSMUSG[0-9]+$", matrix$ensembl_id)),
         paste(selected_group_id, "stable-ID gene membership changed"))
  assert(length(expected_columns) > 0L && all(expected_columns %in% names(matrix)),
         paste(selected_group_id, "display matrix columns changed"))
  numeric_values <- as.matrix(matrix[, ..expected_columns])
  storage.mode(numeric_values) <- "numeric"
  assert(all(is.finite(numeric_values)) && all(is.finite(matrix$final_fdr)) &&
           all(is.finite(matrix$ranking_score)),
         paste(selected_group_id, "contains a non-finite plotted value"))
  assert(all(matrix$common_direction %chin% c("common_up", "common_down")),
         paste(selected_group_id, "common-direction classification changed"))
}

upset <- fread(file.path(target, "provenance/upset_render_audit.csv"))
assert(identical(upset$figure, c("Fig_3a", "Fig_S2a", "Fig_S3a", "Fig_S4a")), "UpSet audit order changed")
assert(identical(upset$union_n, c(884L, 886L, 886L, 888L)) &&
         identical(upset$core_n, c(855L, 870L, 834L, 850L)),
       "UpSet union/core membership counts changed")

component_audit <- fread(file.path(root, "data/publication_input/qsmooth/01_filtered_component_stable_id_mapping_audit.csv"))
tissue_matrix <- fread(file.path(root, "data/publication_input/qsmooth/05_component_tissue_yarn_qsmooth_log2_matrix.csv"))
setorder(component_audit, PlotComponentOrder)
tissue_universe <- tissue_matrix[order(TissueOrder), unique(Group)]
rebuild_tree <- function(pathways) {
  selected <- tissue_matrix[PlotPathwayDisplay %in% pathways]
  profile <- dcast(selected[, .(Group, OfficialMouseSymbol, TissueYARNNormalizedLog2)],
                   Group ~ OfficialMouseSymbol, value.var = "TissueYARNNormalizedLog2")
  profile <- profile[match(tissue_universe, Group)]
  matrix <- as.matrix(profile[, setdiff(names(profile), "Group"), with = FALSE])
  storage.mode(matrix) <- "numeric"
  rownames(matrix) <- profile$Group
  rho <- suppressWarnings(cor(t(matrix), method = "spearman", use = "pairwise.complete.obs"))
  rho[rho > 1] <- 1
  rho[rho < -1] <- -1
  distance <- (1 - rho + t(1 - rho)) / 2
  diag(distance) <- 0
  hclust(as.dist(distance), method = "average")
}
tree_audit <- fread(file.path(target, "provenance/expression_tree_structure_audit.csv"))
formal_tree_order <- fread(file.path(
  root, "data/publication_input/qsmooth/02_tissue_spearman_cluster_order.csv"
))
tree_specs <- list(DSB = c("NHEJ", "HR", "HR–A-EJ", "A-EJ"),
                   SSB = c("BER", "NER", "MMR", "FA"))
for (figure in names(tree_specs)) {
  hc <- rebuild_tree(tree_specs[[figure]])
  figure_name <- figure
  observed <- tree_audit[figure == figure_name]
  assert(nrow(observed) == 25L && identical(observed$left_child, hc$merge[, 1L]) &&
           identical(observed$right_child, hc$merge[, 2L]) &&
           max(abs(observed$height - hc$height)) < 1e-12 &&
           unique(observed$ordered_tissues) == paste(hc$labels[hc$order], collapse = "|"),
         paste(figure, "tree distances or topology changed"))
  frozen_order <- formal_tree_order[Figure == figure][order(OrderTopToBottom)]
  expected_component_n <- if (figure == "DSB") 29L else 45L
  assert(nrow(frozen_order) == 26L && unique(frozen_order$ComponentN) == expected_component_n &&
           identical(hc$labels[hc$order], frozen_order$Tissue) &&
           unique(observed$ordered_tissues) == paste(frozen_order$Tissue, collapse = "|"),
         paste(figure, "tree order differs from the frozen formal order"))
}

meta_audit <- fread(file.path(target, "provenance/meta_dotplot_scope_audit.csv"))
dsb <- fread(file.path(root, "data/publication_input/meta/08_essential_components_counts_log2FC_tissue_matrix.csv"))
ssb <- fread(file.path(root, "results/tables/01_SSB_tissue_meta_log2FC_pvalue_matrix_without_Neil2.csv"))
expected_sig <- c(sum(is.finite(dsb$tissue_meta_p) & dsb$tissue_meta_p < 0.05 & is.finite(dsb$Log2FCDisplayed)),
                  sum(is.finite(ssb$tissue_meta_p) & ssb$tissue_meta_p < 0.05 & is.finite(ssb$Log2FCDisplayed)))
assert(identical(meta_audit$n_rows, c(754L, 1170L)) &&
         identical(meta_audit$n_components, c(29L, 45L)) &&
         identical(meta_audit$n_significant_dots, as.integer(expected_sig)),
       "Meta-dotplot data scope changed")
assert(!"ENSMUSG00000035121" %in% ssb$EnsemblID, "Neil2 remains in the SSB display matrix")

s5 <- fread(file.path(target, "provenance/Fig_S5_scope_audit.csv"))
assert(s5$n_rows == 5327L && s5$n_tissues == 26L && s5$n_pathways == 7L &&
         s5$n_samples == 761L && s5$seed == 25L,
       "Fig. S5 scope or random seed changed")

profile <- if (validation_mode == "full") "17 PNG + 17 PDF" else "17 PNG compact snapshot"
message("FINAL_RELEASE_CONTRACT_PASS: ", profile,
        "; 27 tables unchanged; palette, numeric and tree contracts passed")
