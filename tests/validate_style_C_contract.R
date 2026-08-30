#!/usr/bin/env Rscript

# Contract for the C-story visual variant.  It validates the release boundary
# and frozen inputs; it deliberately does not compare image bytes because the
# PDFs contain renderer metadata and the C variant is a new visual design.

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate validation script", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1L]] else file.path(root, "results_style_C")
if (!grepl("^/", target)) target <- file.path(root, target)
target <- normalizePath(target, mustWork = TRUE)

assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
sha256_file <- function(path) {
  answer <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", answer[[1L]])
}

manifest <- fread(file.path(root, "provenance/reference_sha256.csv"))
table_paths <- sort(manifest$relative_path[grepl("^tables/", manifest$relative_path)])
assert(length(table_paths) == 27L, "Expected 27 publication tables in the reference manifest")

main_names <- paste0("Main/Fig_", c("1", "2", "3a", "3b", "4a", "4b", "5a", "5b"))
suppl_names <- paste0("Suppl/Fig_S", c("1", "2a", "2b", "3a", "3b", "4a", "4b", "5"))
figure_paths <- sort(c(
  paste0(main_names, ".png"), paste0(main_names, ".pdf"),
  paste0(suppl_names, ".png"), paste0(suppl_names, ".pdf")
))
expected <- sort(c(figure_paths, table_paths))
observed <- list.files(target, recursive = TRUE, full.names = FALSE, all.files = FALSE)
observed <- sort(observed[file.info(file.path(target, observed))$isdir %in% FALSE])
assert(identical(observed, expected), paste0(
  "Style C output file set differs. Missing: ",
  paste(setdiff(expected, observed), collapse = "; "),
  " | Extra: ", paste(setdiff(observed, expected), collapse = "; ")
))

input_manifest <- fread(file.path(root, "provenance/publication_input_sha256.csv"))
input_paths <- file.path(root, "data/publication_input", input_manifest$relative_path)
assert(all(file.exists(input_paths)), "Frozen publication input is missing")
input_hashes <- vapply(input_paths, sha256_file, character(1L))
assert(identical(unname(input_hashes), unname(input_manifest$sha256)), "Frozen publication input SHA-256 mismatch")

target_table_hashes <- vapply(file.path(target, table_paths), sha256_file, character(1L))
reference_table_hashes <- manifest$sha256[match(table_paths, manifest$relative_path)]
assert(identical(unname(target_table_hashes), unname(reference_table_hashes)), "Publication table bytes changed")

figure_files <- file.path(target, figure_paths)
assert(all(file.info(figure_files)$size > 0), "One or more C-story figure files are empty")
assert(sum(grepl("^Main/", figure_paths)) == 16L, "Expected 8 main PNGs and 8 main PDFs")
assert(sum(grepl("^Suppl/", figure_paths)) == 16L, "Expected 8 supplementary PNGs and 8 supplementary PDFs")

# Frozen scope checks used by every C renderer.
go <- fread(file.path(root, "data/publication_input/go/04_tissue_statistics_concrete_terms_and_context.csv"))
rho <- fread(file.path(root, "data/publication_input/go/13_tissue_pathway_Spearman_rho_matrix.csv"))
meta <- fread(file.path(root, "data/publication_input/meta/08_essential_components_counts_log2FC_tissue_matrix.csv"))
qsmooth <- fread(file.path(root, "data/publication_input/qsmooth/05_component_tissue_yarn_qsmooth_log2_matrix.csv"))
assert(nrow(go) == 26L * 15L && uniqueN(go$analysis_tissue) == 26L, "GO scope changed")
assert(nrow(rho) == 15L && ncol(rho) == 16L, "Correlation scope changed")
assert(nrow(meta) == 26L * 29L && uniqueN(meta$EnsemblID) == 29L, "DSB scope changed")
assert(nrow(qsmooth) == 26L * 74L && uniqueN(qsmooth$EnsemblID) == 74L, "qsmooth scope changed")

message("STYLE_C_CONTRACT_PASS: 59 files; frozen inputs, tables and analysis scope preserved")
