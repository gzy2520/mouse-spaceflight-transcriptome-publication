#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate validation script", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1L]] else file.path(root, "results")
if (!grepl("^/", target)) target <- file.path(root, target)
target <- normalizePath(target, mustWork = TRUE)

assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
sha256_file <- function(path) {
  answer <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", answer[[1L]])
}

manifest <- fread(file.path(root, "provenance/reference_sha256.csv"))
expected <- sort(manifest$relative_path)

input_manifest <- fread(file.path(root, "provenance/publication_input_sha256.csv"))
input_paths <- file.path(root, "data/publication_input", input_manifest$relative_path)
assert(all(file.exists(input_paths)), paste(
  "Frozen publication input is missing:",
  paste(input_manifest$relative_path[!file.exists(input_paths)], collapse = "; ")
))
input_hashes <- vapply(input_paths, sha256_file, character(1L))
assert(identical(unname(input_hashes), unname(input_manifest$sha256)), paste(
  "Frozen publication input SHA-256 mismatch:",
  paste(input_manifest$relative_path[input_hashes != input_manifest$sha256], collapse = "; ")
))

observed <- list.files(target, recursive = TRUE, full.names = FALSE, all.files = FALSE)
observed <- sort(observed[file.info(file.path(target, observed))$isdir %in% FALSE])
assert(identical(observed, expected), paste0(
  "Output file set differs from the 24-directory release. Missing: ",
  paste(setdiff(expected, observed), collapse = "; "),
  " | Extra: ", paste(setdiff(observed, expected), collapse = "; ")
))

hashes <- vapply(file.path(target, expected), sha256_file, character(1L))
reference_hashes <- manifest$sha256[match(expected, manifest$relative_path)]
assert(identical(unname(hashes), unname(reference_hashes)), paste(
  "SHA-256 mismatch:", paste(expected[hashes != reference_hashes], collapse = "; ")
))

assert(sum(grepl("^Main/Fig_.*[.]png$", expected)) == 8L, "Expected 8 main figures")
assert(sum(grepl("^Suppl/Fig_S.*[.]png$", expected)) == 8L, "Expected 8 supplementary PNG figures")
assert(sum(grepl("^Suppl/Fig_S5[.]pdf$", expected)) == 1L, "Expected the Fig S5 PDF")
assert(sum(grepl("^tables/.*[.]csv$", expected)) == 27L, "Expected 27 tables")

tables <- file.path(target, "tables")
ssb_map <- fread(file.path(tables, "02_SSB_component_stable_id_mapping_without_Neil2.csv"))
assert(nrow(ssb_map) == 45L && uniqueN(ssb_map$EnsemblID) == 45L, "SSB map is not 45 unique Ensembl IDs")
assert(!"ENSMUSG00000035121" %chin% ssb_map$EnsemblID, "Neil2 remains in the SSB publication map")
ssb <- fread(file.path(tables, "01_SSB_tissue_meta_log2FC_pvalue_matrix_without_Neil2.csv"))
assert(nrow(ssb) == 26L * 45L && !anyDuplicated(ssb[, .(Group, EnsemblID)]), "Invalid SSB matrix")

seven <- fread(file.path(tables, "02_seven_pathway_member_gene_counts_per_source_sample_long.csv"))
assert(uniqueN(seven$Group) == 26L && uniqueN(seven$Pathway) == 7L, "Seven-pathway table scope changed")
assert(uniqueN(seven[, .(analysis_unit_id, sample_column, sample_status)]) == 761L, "Source-sample count changed")
assert(!anyDuplicated(seven[, .(analysis_unit_id, sample_column, sample_status, Pathway)]), "Duplicate seven-pathway sample rows")

for (family in c("up_tissues", "down_tissues")) {
  membership <- fread(file.path(tables, paste0("04_membership_matrix_", family, ".csv")))
  assert(!anyDuplicated(membership$ensembl_id), paste("Duplicate Ensembl IDs in", family))
}

message("PUBLICATION_CONTRACT_PASS: 44 files; SHA-256, stable-ID and scope checks passed")
