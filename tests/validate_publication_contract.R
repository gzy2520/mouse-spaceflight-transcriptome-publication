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

assert(sum(grepl("^Main/Fig_.*[.]png$", expected)) == 9L, "Expected 9 main figures")
assert(sum(grepl("^Main/Fig_.*[.]pdf$", expected)) == 1L, "Expected the Fig. 6 vector export")
assert(sum(grepl("^Suppl/Fig_S.*[.]png$", expected)) == 8L, "Expected 8 supplementary PNG figures")
assert(sum(grepl("^Suppl/Fig_S5[.]pdf$", expected)) == 1L, "Expected the Fig S5 PDF")
assert(sum(grepl("^tables/.*[.]csv$", expected)) == 29L, "Expected 29 tables")

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

mouse_dir <- file.path(root, "data/publication_input/mouse_metadata")
mouse <- fread(file.path(mouse_dir, "03_mouse_level_metadata.csv"))
assert(nrow(mouse) == 761L && !anyDuplicated(mouse$sample_column), "Mouse metadata input scope changed")
assert(uniqueN(mouse$accession) == 48L && uniqueN(mouse$analysis_unit_id) == 59L, "Mouse metadata accession/unit scope changed")
assert(uniqueN(mouse$Group) == 26L && uniqueN(mouse$mission_cluster) == 12L, "Mouse metadata mission/tissue scope changed")
assert(!"SpaceX-8+SpaceX-9" %chin% mouse$mission_cluster, "The erroneous composite OSD-162 mission remains")
assert(unique(mouse[accession == "OSD-162", mission_cluster]) == "SpaceX-8", "OSD-162 is not assigned to the current OSDR mission")
assert(all(mouse$sex %chin% c("Female", "Male")) && all(!is.na(mouse$age_order_weeks)), "Mouse metadata has missing sex or age")
assert(all(mouse$isa_status_matches_analysis), "Mouse metadata conflicts with the analysis Flight/Ground scope")
assert(!anyNA(mouse[, .(mouse_id, metadata_zip_sha256, age_source_field)]), "Mouse metadata provenance is incomplete")

source_core <- fread(file.path(mouse_dir, "01_osdr_isa_sample_source_core_all.csv"))
download_manifest <- fread(file.path(mouse_dir, "02_osdr_download_manifest.csv"))
join_audit <- fread(file.path(mouse_dir, "04_sample_join_audit.csv"))
conflict_audit <- fread(file.path(mouse_dir, "06_metadata_conflict_audit.csv"))
assert(nrow(source_core) == 1588L && uniqueN(source_core$accession) == 48L, "Complete OSDR ISA source table scope changed")
assert(nrow(download_manifest) == 48L && !anyDuplicated(download_manifest$accession), "OSDR download manifest scope changed")
assert(all(grepl("^[0-9a-f]{64}$", download_manifest$metadata_zip_sha256)), "Invalid OSDR ZIP SHA-256")
assert(all(mouse$metadata_zip_sha256 == download_manifest$metadata_zip_sha256[match(mouse$accession, download_manifest$accession)]), "Selected rows do not match the OSDR ZIP manifest")
assert(nrow(join_audit) == 48L && sum(join_audit$n_missing_age) == 0L && sum(join_audit$n_missing_sex) == 0L && sum(join_audit$n_status_conflicts) == 0L, "OSDR sample join audit failed")
assert(nrow(conflict_audit) == 5L && all(conflict_audit$accession == "OSD-162") && all(conflict_audit$canonical_value == "SpaceX-8"), "Expected OSD-162 legacy mission conflict audit changed")
mouse_grouped <- fread(file.path(tables, "Fig_6_mouse_metadata_grouped.csv"))
assert(!anyDuplicated(mouse_grouped[, .(mission_cluster, Group, sex, age_label)]), "Duplicate Fig. 6 strata")
assert(all(mouse_grouped$n_mice > 0L), "Fig. 6 contains a non-positive mouse count")
assert(all(c("source_age_values", "age_source_fields", "age_is_range", "accessions") %chin% names(mouse_grouped)), "Fig. 6 source-age provenance columns are missing")
age_map <- unique(mouse_grouped[, .(age_label, age_colour)])
assert(!anyNA(mouse_grouped$age_colour) && nrow(age_map) == uniqueN(mouse_grouped$age_label) && uniqueN(age_map$age_colour) == nrow(age_map), "Fig. 6 age colours are not one-to-one")

mouse_audit <- fread(file.path(tables, "Fig_6_mouse_sample_metadata_complete_audit.csv"))
assert(nrow(mouse_audit) == 761L && ncol(mouse_audit) == 54L, "Complete mouse sample audit dimensions changed")
assert(!anyDuplicated(mouse_audit$audit_row_id) && all(mouse_audit$audit_status == "PASS"), "Complete mouse sample audit has duplicate or failed rows")
assert(uniqueN(mouse_audit$accession) == 48L && uniqueN(mouse_audit$analysis_unit_id) == 59L && uniqueN(mouse_audit$tissue) == 26L, "Complete mouse sample audit scope changed")
assert(sum(mouse_audit$age_is_range) == 312L && sum(mouse_audit$mission_conflict_detected) == 5L, "Complete mouse sample audit exception counts changed")
assert(all(mouse_audit$status_matches_analysis) && all(grepl("^[0-9a-f]{64}$", mouse_audit$metadata_zip_sha256)), "Complete mouse sample audit status or ZIP provenance failed")

for (family in c("up_tissues", "down_tissues")) {
  membership <- fread(file.path(tables, paste0("04_membership_matrix_", family, ".csv")))
  assert(!anyDuplicated(membership$ensembl_id), paste("Duplicate Ensembl IDs in", family))
}

message("PUBLICATION_CONTRACT_PASS: 48 files; SHA-256, stable-ID and scope checks passed")
