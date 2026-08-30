#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate validation script", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1L]] else file.path(root, "results_style_A")
if (!grepl("^/", target)) target <- file.path(root, target)
target <- normalizePath(target, mustWork = TRUE)

assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
sha256_file <- function(path) {
  answer <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", answer[[1L]])
}

manifest <- fread(file.path(root, "provenance/reference_sha256.csv"))
expected <- sort(manifest$relative_path)
observed <- list.files(target, recursive = TRUE, full.names = FALSE, all.files = FALSE)
observed <- sort(observed[file.info(file.path(target, observed))$isdir %in% FALSE])
assert(identical(observed, expected), paste0(
  "Style A output file set differs from the 24-directory release. Missing: ",
  paste(setdiff(expected, observed), collapse = "; "),
  " | Extra: ", paste(setdiff(observed, expected), collapse = "; ")
))

input_manifest <- fread(file.path(root, "provenance/publication_input_sha256.csv"))
input_paths <- file.path(root, "data/publication_input", input_manifest$relative_path)
assert(all(file.exists(input_paths)), "Frozen publication input is missing")
input_hashes <- vapply(input_paths, sha256_file, character(1L))
assert(identical(unname(input_hashes), unname(input_manifest$sha256)), "Frozen publication input SHA-256 mismatch")

# Tables are copied byte-for-byte from the approved release. This is the
# strongest available guard that the visual-only renderer did not alter data.
table_paths <- expected[grepl("^tables/", expected)]
target_table_hashes <- vapply(file.path(target, table_paths), sha256_file, character(1L))
reference_table_hashes <- manifest$sha256[match(table_paths, manifest$relative_path)]
assert(identical(unname(target_table_hashes), unname(reference_table_hashes)), "Publication table bytes changed")

assert(sum(grepl("^Main/Fig_.*[.]png$", expected)) == 8L, "Expected 8 main figures")
assert(sum(grepl("^Suppl/Fig_S.*[.]png$", expected)) == 8L, "Expected 8 supplementary PNG figures")
assert(sum(grepl("^Suppl/Fig_S5[.]pdf$", expected)) == 1L, "Expected the Fig S5 PDF")
assert(sum(grepl("^tables/.*[.]csv$", expected)) == 27L, "Expected 27 tables")

png_paths <- file.path(target, expected[grepl("[.]png$", expected)])
assert(all(file.info(png_paths)$size > 0), "One or more Style A PNGs are empty")
message("STYLE_A_CONTRACT_PASS: 44 files; frozen inputs, tables and scope preserved")
