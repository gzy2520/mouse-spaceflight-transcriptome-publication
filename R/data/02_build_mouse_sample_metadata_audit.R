#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate audit-table script", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

input_dir <- file.path(root, "data/publication_input/mouse_metadata")
output_root <- Sys.getenv("PUBLICATION_OUTPUT_DIR", unset = file.path(root, "results"))
output_dir <- file.path(output_root, "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "Fig_6_mouse_sample_metadata_complete_audit.csv")

scope <- fread(file.path(input_dir, "00_analysis_sample_scope.csv"))
source_core <- fread(file.path(input_dir, "01_osdr_isa_sample_source_core_all.csv"))
download_manifest <- fread(file.path(input_dir, "02_osdr_download_manifest.csv"))
selected <- fread(file.path(input_dir, "03_mouse_level_metadata.csv"))
conflicts <- fread(file.path(input_dir, "06_metadata_conflict_audit.csv"))

stopifnot(
  nrow(scope) == 761L,
  nrow(selected) == 761L,
  !anyDuplicated(scope$sample_column),
  !anyDuplicated(selected$sample_column),
  nrow(download_manifest) == 48L,
  !anyDuplicated(download_manifest$accession)
)

scope_order <- scope[, .(
  sample_column,
  scope_order = .I
)]

source_index <- source_core[, .(
  accession,
  sample_name_normalized = isa_sample_name,
  isa_sample_name,
  source_name,
  isa_status,
  mission_project_token
)]
stopifnot(!anyDuplicated(source_index[, .(accession, sample_name_normalized)]))

manifest_index <- download_manifest[, .(
  accession,
  osdr_study_url,
  dataset_api_url,
  files_api_url,
  metadata_download_url,
  manifest_metadata_zip_filename = metadata_zip_filename,
  manifest_metadata_zip_sha256 = metadata_zip_sha256,
  metadata_zip_size_bytes,
  metadata_downloaded_at_utc = downloaded_at_utc,
  api_project_identifiers,
  manifest_api_mission_names = api_mission_names,
  api_mission_start_dates,
  api_mission_end_dates
)]

conflict_index <- conflicts[, .(
  sample_column,
  mission_conflict_detected = TRUE,
  mission_conflict_issue_type = issue_type,
  mission_conflict_source_value = source_value,
  mission_conflict_canonical_value = canonical_value,
  mission_conflict_resolution = resolution
)]
stopifnot(!anyDuplicated(conflict_index$sample_column))

audit <- merge(selected, scope_order, by = "sample_column", all.x = TRUE, sort = FALSE)
audit <- merge(
  audit,
  source_index,
  by = c("accession", "sample_name_normalized"),
  all.x = TRUE,
  sort = FALSE
)
audit <- merge(audit, manifest_index, by = "accession", all.x = TRUE, sort = FALSE)
audit <- merge(audit, conflict_index, by = "sample_column", all.x = TRUE, sort = FALSE)
setorder(audit, scope_order)

audit[is.na(mission_conflict_detected), mission_conflict_detected := FALSE]
for (column in c(
  "mission_conflict_issue_type",
  "mission_conflict_source_value",
  "mission_conflict_canonical_value",
  "mission_conflict_resolution"
)) {
  set(audit, which(is.na(audit[[column]])), column, "")
}

audit[, accession_scoped_mouse_id := paste(accession, mouse_id, sep = "::")]
audit[, audit_row_id := paste(accession, analysis_unit_id, sample_column, sep = "::")]
audit[, age_value_status := fifelse(
  age_is_range,
  "RANGE_REPORTED_BY_SOURCE",
  "EXACT_VALUE_REPORTED_BY_SOURCE"
)]
audit[, sample_match_status := fifelse(
  metadata_match_method == "exact ISA Sample Name",
  "EXACT",
  "NORMALIZED_TRANSFORM"
)]

audit[, audit_flags := vapply(seq_len(.N), function(row_index) {
  flags <- character()
  if (isTRUE(age_is_range[[row_index]])) flags <- c(flags, "SOURCE_AGE_RANGE")
  if (isTRUE(sample_is_technical_replicate[[row_index]])) flags <- c(flags, "TECHNICAL_REPLICATE")
  if (sample_match_status[[row_index]] == "NORMALIZED_TRANSFORM") flags <- c(flags, "SAMPLE_NAME_NORMALIZED")
  if (isTRUE(mission_conflict_detected[[row_index]])) flags <- c(flags, "MISSION_CONFLICT_RESOLVED")
  if (!length(flags)) "NONE" else paste(flags, collapse = ";")
}, character(1L))]

audit[, audit_status := fifelse(
  !is.na(scope_order) &
    !is.na(isa_sample_name) &
    !is.na(source_name) &
    source_name == mouse_id &
    isa_status_matches_analysis &
    !is.na(sex) &
    !is.na(age_label) &
    !is.na(mission_cluster) &
    metadata_zip_filename == manifest_metadata_zip_filename &
    metadata_zip_sha256 == manifest_metadata_zip_sha256 &
    api_mission_names == manifest_api_mission_names,
  "PASS",
  "FAIL"
)]

result <- audit[, .(
  audit_row_id,
  audit_status,
  audit_flags,
  accession,
  analysis_unit_id,
  tissue = Group,
  analysis_input_relative_path = input_relative_path,
  analysis_sample_column = sample_column,
  normalized_sample_name = sample_name_normalized,
  isa_sample_name,
  source_name,
  accession_scoped_mouse_id,
  sample_match_method = metadata_match_method,
  sample_match_status,
  sample_is_technical_replicate,
  analysis_sample_status = sample_status,
  isa_spaceflight_raw,
  isa_status,
  status_matches_analysis = isa_status_matches_analysis,
  sex,
  sex_raw,
  age_source_field,
  age_value_raw,
  age_unit_raw,
  age_raw,
  age_label,
  age_order_weeks,
  age_is_range,
  age_value_status,
  organism_raw,
  source_material_type,
  mission_cluster,
  mission_source,
  api_mission_names,
  mission_project_token,
  source_space_mission_raw,
  launch_mission_raw,
  mission_conflict_detected,
  mission_conflict_issue_type,
  mission_conflict_source_value,
  mission_conflict_canonical_value,
  mission_conflict_resolution,
  isa_source_table,
  metadata_zip_filename,
  metadata_zip_sha256,
  metadata_zip_size_bytes,
  metadata_downloaded_at_utc,
  osdr_study_url,
  dataset_api_url,
  files_api_url,
  metadata_download_url,
  api_project_identifiers,
  api_mission_start_dates,
  api_mission_end_dates
)]

stopifnot(
  nrow(result) == 761L,
  !anyDuplicated(result$audit_row_id),
  all(result$audit_status == "PASS"),
  uniqueN(result$accession) == 48L,
  uniqueN(result$analysis_unit_id) == 59L,
  uniqueN(result$tissue) == 26L,
  sum(result$age_is_range) == 312L,
  sum(result$mission_conflict_detected) == 5L,
  all(grepl("^[0-9a-f]{64}$", result$metadata_zip_sha256))
)

fwrite(result, output_path, bom = FALSE, na = "")
message(
  "MOUSE_SAMPLE_AUDIT_PASS: ", output_path, "; ",
  nrow(result), " rows; ", ncol(result), " columns; all rows PASS"
)
