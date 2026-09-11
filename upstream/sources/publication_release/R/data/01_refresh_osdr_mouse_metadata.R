#!/usr/bin/env Rscript

# Refresh the mouse metadata used by Fig. 6 directly from NASA OSDR.
#
# The analysis sample scope is taken from the checked-in seven-pathway source
# table.  Age, sex, organism/source identity, material type, spaceflight status
# and sample-level mission fields are read from each accession's ISA metadata
# ZIP.  Canonical mission names come from the current OSDR dataset API; this is
# important because some legacy ISA sample rows contain stale launch labels.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate 01_refresh_osdr_mouse_metadata.R", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
reuse_cache <- "--reuse-cache" %chin% args

value_after_prefix <- function(prefix, default) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

scope_path <- value_after_prefix(
  "--scope=",
  file.path(repo_root, "results/tables/02_seven_pathway_member_gene_counts_per_source_sample_long.csv")
)
output_dir <- value_after_prefix(
  "--output-dir=",
  file.path(repo_root, "data/publication_input/mouse_metadata")
)
cache_dir <- value_after_prefix(
  "--cache-dir=",
  file.path(repo_root, "data/.cache/osdr_isa")
)
for (path_name in c("scope_path", "output_dir", "cache_dir")) {
  path_value <- get(path_name)
  if (!grepl("^/", path_value)) path_value <- file.path(repo_root, path_value)
  assign(path_name, normalizePath(path_value, mustWork = FALSE))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

sha256_file <- function(path) {
  result <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  assert(length(result) == 1L, paste("Could not hash", path))
  sub("[[:space:]].*$", "", result[[1L]])
}

read_json_retry <- function(url, attempts = 3L) {
  last_error <- NULL
  for (attempt in seq_len(attempts)) {
    answer <- tryCatch(
      fromJSON(url, simplifyVector = FALSE),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(answer)) return(answer)
    if (attempt < attempts) Sys.sleep(attempt)
  }
  stop("Failed to read OSDR JSON after ", attempts, " attempts: ", url,
       " | ", last_error, call. = FALSE)
}

download_retry <- function(url, destination, attempts = 3L) {
  last_error <- NULL
  partial <- paste0(destination, ".part")
  if (file.exists(partial)) unlink(partial)
  for (attempt in seq_len(attempts)) {
    ok <- tryCatch({
      download.file(url, partial, mode = "wb", method = "libcurl", quiet = TRUE)
      file.exists(partial) && file.info(partial)$size > 0
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      FALSE
    })
    if (isTRUE(ok)) {
      if (file.exists(destination)) unlink(destination)
      assert(file.rename(partial, destination), paste("Could not finalize", destination))
      return(invisible(destination))
    }
    if (file.exists(partial)) unlink(partial)
    if (attempt < attempts) Sys.sleep(attempt)
  }
  stop("Failed to download OSDR metadata after ", attempts, " attempts: ", url,
       " | ", last_error, call. = FALSE)
}

flatten_chr <- function(value) {
  answer <- as.character(unlist(value, recursive = TRUE, use.names = FALSE))
  answer <- trimws(answer)
  unique(answer[!is.na(answer) & nzchar(answer)])
}

split_csv_values <- function(value) {
  answer <- flatten_chr(value)
  if (!length(answer)) return(character())
  answer <- trimws(unlist(strsplit(answer, ",", fixed = TRUE), use.names = FALSE))
  unique(answer[nzchar(answer)])
}

first_col_index <- function(column_names, patterns) {
  for (pattern in patterns) {
    hit <- which(grepl(pattern, column_names, ignore.case = TRUE, perl = TRUE))
    if (length(hit)) return(hit[[1L]])
  }
  NA_integer_
}

column_or_na <- function(table, index) {
  if (is.na(index)) return(rep(NA_character_, nrow(table)))
  as.character(table[[index]])
}

unit_after <- function(table, index) {
  if (is.na(index) || index >= ncol(table)) return(rep(NA_character_, nrow(table)))
  if (!grepl("^Unit$", names(table)[index + 1L], ignore.case = TRUE)) {
    return(rep(NA_character_, nrow(table)))
  }
  as.character(table[[index + 1L]])
}

clean_text <- function(value) {
  value <- trimws(as.character(value))
  value[value %chin% c("", "NA", "N/A")] <- NA_character_
  value
}

canonical_project_token <- function(value) {
  value <- toupper(clean_text(value))
  value <- gsub("MOUSE HABITAT UNIT", "MHU", value, fixed = TRUE)
  answer <- rep(NA_character_, length(value))
  mhu <- regexpr("MHU[ -]?[0-9]+", value, perl = TRUE)
  has_mhu <- !is.na(value) & mhu > 0L
  if (any(has_mhu)) {
    hit <- regmatches(value, mhu)
    answer[has_mhu] <- paste0("MHU-", gsub("[^0-9]", "", hit[has_mhu]))
  }
  rr <- regexpr("RR[ -]?[0-9]+", value, perl = TRUE)
  has_rr <- is.na(answer) & !is.na(value) & rr > 0L
  if (any(has_rr)) {
    hit <- regmatches(value, rr)
    answer[has_rr] <- paste0("RR-", gsub("[^0-9]", "", hit[has_rr]))
  }
  answer
}

normalize_status <- function(value) {
  value <- tolower(clean_text(value))
  answer <- rep(NA_character_, length(value))
  answer[!is.na(value) & grepl("space[ -]?flight", value)] <- "Flight"
  answer[!is.na(value) & grepl("ground control|vivarium control", value)] <- "Ground"
  answer[!is.na(value) & grepl("basal control", value)] <- "Basal"
  answer
}

normalize_sex <- function(value) {
  value <- tolower(clean_text(value))
  answer <- rep(NA_character_, length(value))
  answer[value %chin% c("female", "f")] <- "Female"
  answer[value %chin% c("male", "m")] <- "Male"
  answer
}

age_label_one <- function(value, unit) {
  if (is.na(value) || !nzchar(trimws(value))) return(NA_character_)
  label <- trimws(value)
  label <- gsub("(?i)\\s+to\\s+", "–", label, perl = TRUE)
  label <- gsub("\\s*[-–—]\\s*", "–", label, perl = TRUE)
  label <- gsub("([0-9]+)[.]0(?=($|–))", "\\1", label, perl = TRUE)
  unit_clean <- tolower(trimws(ifelse(is.na(unit), "", unit)))
  suffix <- if (unit_clean %chin% c("week", "weeks")) {
    "wk"
  } else if (unit_clean %chin% c("month", "months")) {
    "mo"
  } else if (unit_clean %chin% c("day", "days")) {
    "d"
  } else if (unit_clean %chin% c("year", "years")) {
    "yr"
  } else if (!nzchar(unit_clean) || unit_clean %chin% c("not applicable", "not available")) {
    ""
  } else {
    unit_clean
  }
  paste0(label, if (nzchar(suffix)) paste0(" ", suffix) else "")
}

age_order_one <- function(value, unit) {
  if (is.na(value)) return(NA_real_)
  hits <- regmatches(value, gregexpr("[0-9]+(?:[.][0-9]+)?", value, perl = TRUE))[[1L]]
  if (!length(hits)) return(NA_real_)
  midpoint <- mean(as.numeric(hits))
  unit_clean <- tolower(trimws(ifelse(is.na(unit), "", unit)))
  if (unit_clean %chin% c("week", "weeks")) return(midpoint)
  if (unit_clean %chin% c("month", "months")) return(midpoint * 365.2425 / 12 / 7)
  if (unit_clean %chin% c("day", "days")) return(midpoint / 7)
  if (unit_clean %chin% c("year", "years")) return(midpoint * 365.2425 / 7)
  NA_real_
}

age_is_range_one <- function(value) {
  if (is.na(value)) return(FALSE)
  grepl("(?i)\\bto\\b|[0-9]\\s*[-–—]\\s*[0-9]", value, perl = TRUE)
}

meaningful_mission <- function(value) {
  value <- clean_text(value)
  !is.na(value) & !tolower(value) %chin% c("not applicable", "not available", "unknown")
}

scope_long <- fread(scope_path, na.strings = c("", "NA"), encoding = "UTF-8")
required_scope <- c(
  "Group", "analysis_unit_id", "accession", "input_relative_path",
  "sample_status", "sample_column", "sample_is_technical_replicate"
)
assert(all(required_scope %chin% names(scope_long)), "The seven-pathway scope table changed schema")
scope <- unique(scope_long[, ..required_scope])
setorder(scope, accession, analysis_unit_id, sample_status, sample_column)
assert(nrow(scope) == 761L, "Expected 761 unique analysis sample columns")
assert(uniqueN(scope$accession) == 48L, "Expected 48 OSDR accessions")
assert(!anyDuplicated(scope$sample_column), "Analysis sample columns are not globally unique")
fwrite(scope, file.path(output_dir, "00_analysis_sample_scope.csv"))

download_time <- format(Sys.time(), tz = "UTC", usetz = TRUE)
manifest_rows <- list()
source_rows <- list()

accessions <- sort(unique(scope$accession))
for (index in seq_along(accessions)) {
  accession <- accessions[[index]]
  message(sprintf("[%02d/%02d] NASA OSDR %s", index, length(accessions), accession))

  dataset_api <- sprintf(
    "https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/%s/?format=json",
    accession
  )
  files_api <- sprintf(
    "https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/%s/files/",
    accession
  )
  dataset_json <- read_json_retry(dataset_api)
  files_json <- read_json_retry(files_api)
  dataset <- dataset_json[[accession]]
  files <- files_json[[accession]][["files"]]
  assert(!is.null(dataset), paste("Dataset API did not return", accession))
  assert(length(files) > 0L, paste("Files API did not return files for", accession))

  metadata_names <- names(files)[grepl("metadata.*ISA[.]zip$", names(files), ignore.case = TRUE)]
  assert(length(metadata_names) == 1L, paste(
    accession, "metadata ISA ZIP count was", length(metadata_names),
    paste(metadata_names, collapse = ";")
  ))
  metadata_name <- metadata_names[[1L]]
  metadata_url <- files[[metadata_name]][["URL"]]
  assert(length(metadata_url) == 1L && nzchar(metadata_url), paste("Missing metadata URL for", accession))
  zip_path <- file.path(cache_dir, metadata_name)
  if (!reuse_cache || !file.exists(zip_path)) download_retry(metadata_url, zip_path)
  assert(file.exists(zip_path) && file.info(zip_path)$size > 0L, paste("Missing ZIP for", accession))

  extract_dir <- tempfile(paste0(accession, "_isa_"))
  dir.create(extract_dir)
  on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
  unzip(zip_path, exdir = extract_dir)
  source_files <- list.files(
    extract_dir, pattern = "^s_.*[.]txt$", full.names = TRUE, recursive = TRUE
  )
  assert(length(source_files) >= 1L, paste("No ISA source/sample table in", metadata_name))
  source_table <- rbindlist(lapply(source_files, function(path) {
    answer <- fread(
      path, sep = "\t", quote = "", check.names = FALSE, fill = TRUE,
      na.strings = c("", "NA"), encoding = "UTF-8"
    )
    answer[, isa_source_table := basename(path)]
    answer
  }), fill = TRUE, use.names = TRUE)

  column_names <- names(source_table)
  source_index <- first_col_index(column_names, "^Source Name$")
  sample_index <- first_col_index(column_names, "^Sample Name$")
  sex_index <- first_col_index(column_names, c(
    "^Characteristics\\[Sex\\]$", "^Characteristics\\[Gender\\]$"
  ))
  age_index <- first_col_index(column_names, c(
    "^Characteristics\\[Age at Launch\\]$",
    "^Factor Value\\[Age\\]$",
    "^Characteristics\\[Age\\]$",
    "^Parameter Value\\[Age at Euthanasia\\]$"
  ))
  organism_index <- first_col_index(column_names, "^Characteristics\\[Organism\\]$")
  material_index <- first_col_index(column_names, "^Characteristics\\[Material Type\\]$")
  flight_index <- first_col_index(column_names, c(
    "^Factor Value\\[Spaceflight\\]$", "^Characteristics\\[Spaceflight\\]$"
  ))
  launch_mission_index <- first_col_index(column_names, c(
    "^Characteristics\\[Launch Mission\\]$", "^Characteristics\\[Mission.*\\]$"
  ))
  space_mission_index <- first_col_index(column_names, c(
    "^Factor Value\\[Space Mission\\]$", "^Factor Value\\[Mission.*\\]$"
  ))
  assert(!is.na(source_index) && !is.na(sample_index), paste("Missing source/sample names in", accession))
  assert(!is.na(sex_index), paste("Missing sex field in", accession))
  assert(!is.na(age_index), paste("Missing age field in", accession))
  assert(!is.na(flight_index), paste("Missing spaceflight field in", accession))

  metadata <- dataset[["metadata"]]
  mission_names <- flatten_chr(metadata[["mission"]][["name"]])
  mission_starts <- flatten_chr(metadata[["mission"]][["start date"]])
  mission_ends <- flatten_chr(metadata[["mission"]][["end date"]])
  project_ids <- split_csv_values(metadata[["project identifier"]])
  assert(length(mission_names) >= 1L, paste("Dataset API mission missing for", accession))

  core <- data.table(
    accession = accession,
    source_name = clean_text(column_or_na(source_table, source_index)),
    isa_sample_name = clean_text(column_or_na(source_table, sample_index)),
    organism_raw = clean_text(column_or_na(source_table, organism_index)),
    source_material_type = clean_text(column_or_na(source_table, material_index)),
    sex_raw = clean_text(column_or_na(source_table, sex_index)),
    age_value_raw = clean_text(column_or_na(source_table, age_index)),
    age_unit_raw = clean_text(unit_after(source_table, age_index)),
    isa_spaceflight_raw = clean_text(column_or_na(source_table, flight_index)),
    launch_mission_raw = clean_text(column_or_na(source_table, launch_mission_index)),
    source_space_mission_raw = clean_text(column_or_na(source_table, space_mission_index)),
    isa_source_table = source_table$isa_source_table,
    metadata_zip_filename = metadata_name
  )
  core <- core[!is.na(source_name) | !is.na(isa_sample_name)]
  assert(!anyDuplicated(core$isa_sample_name[!is.na(core$isa_sample_name)]),
         paste("Duplicate ISA Sample Name in", accession))
  core[, sex := normalize_sex(sex_raw)]
  core[, isa_status := normalize_status(isa_spaceflight_raw)]
  core[, age_label := mapply(age_label_one, age_value_raw, age_unit_raw, USE.NAMES = FALSE)]
  core[, age_order_weeks := mapply(age_order_one, age_value_raw, age_unit_raw, USE.NAMES = FALSE)]
  core[, age_is_range := mapply(age_is_range_one, age_value_raw, USE.NAMES = FALSE)]
  core[, age_source_field := column_names[[age_index]]]
  core[, api_mission_names := paste(mission_names, collapse = ";")]
  core[, mission_source := NA_character_]

  if (length(mission_names) == 1L) {
    core[, mission_cluster := mission_names[[1L]]]
    core[, mission_source := "OSDR dataset API: single mission"]
  } else {
    assert(length(project_ids) == length(mission_names), paste(
      "Cannot align project identifiers and missions for", accession,
      "projects=", paste(project_ids, collapse = ";"),
      "missions=", paste(mission_names, collapse = ";")
    ))
    mission_map <- data.table(
      project_token = canonical_project_token(project_ids),
      mission_cluster = mission_names
    )
    assert(!anyNA(mission_map$project_token) && !anyDuplicated(mission_map$project_token),
           paste("Invalid multi-mission project map for", accession))
    core[, mission_project_token := canonical_project_token(source_space_mission_raw)]
    core[, mission_cluster := mission_map$mission_cluster[
      match(mission_project_token, mission_map$project_token)
    ]]
    core[!is.na(mission_cluster), mission_source :=
           "OSDR ISA Space Mission mapped to OSDR dataset API mission"]
  }

  zip_hash <- sha256_file(zip_path)
  core[, metadata_zip_sha256 := zip_hash]
  source_rows[[accession]] <- core

  manifest_rows[[accession]] <- data.table(
    accession = accession,
    osdr_study_url = sprintf("https://osdr.nasa.gov/bio/repo/data/studies/%s/", accession),
    dataset_api_url = dataset_api,
    files_api_url = files_api,
    metadata_download_url = metadata_url,
    metadata_zip_filename = metadata_name,
    metadata_zip_sha256 = zip_hash,
    metadata_zip_size_bytes = file.info(zip_path)$size,
    downloaded_at_utc = download_time,
    api_project_identifiers = paste(project_ids, collapse = ";"),
    api_mission_names = paste(mission_names, collapse = ";"),
    api_mission_start_dates = paste(mission_starts, collapse = ";"),
    api_mission_end_dates = paste(mission_ends, collapse = ";"),
    isa_source_tables = paste(basename(source_files), collapse = ";"),
    n_isa_source_rows = nrow(core),
    age_source_field = column_names[[age_index]],
    sex_source_field = column_names[[sex_index]],
    spaceflight_source_field = column_names[[flight_index]],
    launch_mission_source_field = if (is.na(launch_mission_index)) "" else column_names[[launch_mission_index]],
    space_mission_source_field = if (is.na(space_mission_index)) "" else column_names[[space_mission_index]]
  )
  unlink(extract_dir, recursive = TRUE)
}

source_core <- rbindlist(source_rows, fill = TRUE, use.names = TRUE)
download_manifest <- rbindlist(manifest_rows, fill = TRUE, use.names = TRUE)
setorder(source_core, accession, source_name, isa_sample_name)
setorder(download_manifest, accession)
fwrite(source_core, file.path(output_dir, "01_osdr_isa_sample_source_core_all.csv"))
fwrite(download_manifest, file.path(output_dir, "02_osdr_download_manifest.csv"))

matched_rows <- list()
for (accession in accessions) {
  accession_value <- accession
  accession_scope <- copy(scope[accession == accession_value])
  accession_source <- source_core[accession == accession_value]
  source_index <- match(accession_scope$sample_column, accession_source$isa_sample_name)
  match_method <- rep(NA_character_, nrow(accession_scope))
  match_method[!is.na(source_index)] <- "exact ISA Sample Name"

  unresolved <- which(is.na(source_index))
  if (length(unresolved)) {
    transformed <- gsub(
      "(?i)_mRNA(?=(_techrep[0-9]+)?$)", "",
      accession_scope$sample_column[unresolved], perl = TRUE
    )
    transformed_index <- match(transformed, accession_source$isa_sample_name)
    resolved <- !is.na(transformed_index)
    source_index[unresolved[resolved]] <- transformed_index[resolved]
    match_method[unresolved[resolved]] <- "remove expression-only _mRNA token"
  }

  appended <- accession_source[source_index]
  appended[, accession := NULL]
  matched <- cbind(accession_scope, appended)
  matched[, metadata_match_method := match_method]
  matched_rows[[accession]] <- matched
}

mouse <- rbindlist(matched_rows, fill = TRUE, use.names = TRUE)
setorder(mouse, accession, analysis_unit_id, sample_status, sample_column)
mouse[, mouse_id := source_name]
mouse[, sample_name_normalized := isa_sample_name]
mouse[, age_raw := paste0(age_value_raw, " {", age_unit_raw, "}")]
mouse[, isa_status_matches_analysis := isa_status == sample_status]
mouse[, launch_mission_conflict := meaningful_mission(launch_mission_raw) &
        !is.na(mission_cluster) & launch_mission_raw != mission_cluster]

hard_fail <- mouse[
  is.na(isa_sample_name) | is.na(mouse_id) | is.na(sex) | is.na(age_label) |
    is.na(age_order_weeks) | is.na(mission_cluster) | !isa_status_matches_analysis
]
assert(nrow(mouse) == 761L && !anyDuplicated(mouse$sample_column),
       "Refreshed selected metadata no longer has one row per analysis sample")
assert(nrow(hard_fail) == 0L, paste(
  "Refreshed metadata has unmatched or invalid selected samples:",
  paste(unique(hard_fail$accession), collapse = ";")
))
assert(all(mouse$sex %chin% c("Female", "Male")), "Unexpected selected-sample sex value")

mouse_output <- mouse[, .(
  sample_column,
  sample_name_normalized,
  mouse_id,
  sample_status,
  isa_spaceflight_raw,
  isa_status_matches_analysis,
  sample_is_technical_replicate,
  Group,
  analysis_unit_id,
  accession,
  input_relative_path,
  mission_cluster,
  mission_source,
  api_mission_names,
  source_space_mission_raw,
  launch_mission_raw,
  launch_mission_conflict,
  sex,
  sex_raw,
  age_value_raw,
  age_unit_raw,
  age_raw,
  age_label,
  age_order_weeks,
  age_is_range,
  age_source_field,
  organism_raw,
  source_material_type,
  metadata_match_method,
  isa_source_table,
  metadata_zip_filename,
  metadata_zip_sha256
)]
fwrite(mouse_output, file.path(output_dir, "03_mouse_level_metadata.csv"))

join_audit <- mouse_output[, .(
  n_selected_sample_columns = .N,
  n_unique_mice = uniqueN(mouse_id),
  n_exact_matches = sum(metadata_match_method == "exact ISA Sample Name"),
  n_transformed_matches = sum(metadata_match_method != "exact ISA Sample Name"),
  n_missing_age = sum(is.na(age_label)),
  n_age_range_samples = sum(age_is_range),
  n_missing_sex = sum(is.na(sex)),
  n_status_conflicts = sum(!isa_status_matches_analysis),
  n_launch_mission_conflicts = sum(launch_mission_conflict),
  selected_age_labels = paste(sort(unique(age_label)), collapse = ";"),
  selected_sexes = paste(sort(unique(sex)), collapse = ";"),
  selected_missions = paste(sort(unique(mission_cluster)), collapse = ";")
), by = accession]
setorder(join_audit, accession)
fwrite(join_audit, file.path(output_dir, "04_sample_join_audit.csv"))

source_value_audit <- mouse_output[, .(
  n_selected_sample_columns = .N,
  n_unique_mice = uniqueN(mouse_id),
  age_source_fields = paste(sort(unique(age_source_field)), collapse = ";"),
  age_values_raw = paste(sort(unique(paste0(age_value_raw, " {", age_unit_raw, "}"))), collapse = ";"),
  age_labels = paste(sort(unique(age_label)), collapse = ";"),
  n_age_range_samples = sum(age_is_range),
  sex_values_raw = paste(sort(unique(sex_raw)), collapse = ";"),
  api_missions = paste(sort(unique(api_mission_names)), collapse = ";"),
  canonical_missions = paste(sort(unique(mission_cluster)), collapse = ";"),
  source_launch_missions = paste(sort(unique(launch_mission_raw[meaningful_mission(launch_mission_raw)])), collapse = ";")
), by = accession]
setorder(source_value_audit, accession)
fwrite(source_value_audit, file.path(output_dir, "05_selected_source_value_audit.csv"))

conflict_audit <- rbindlist(list(
  mouse_output[launch_mission_conflict == TRUE, .(
    accession, analysis_unit_id, sample_column, mouse_id,
    issue_type = "ISA launch mission conflicts with current dataset API mission",
    source_value = launch_mission_raw,
    canonical_value = mission_cluster,
    resolution = "Use the current OSDR dataset API mission for plotting; retain the ISA value for audit"
  )],
  mouse_output[isa_status_matches_analysis == FALSE, .(
    accession, analysis_unit_id, sample_column, mouse_id,
    issue_type = "ISA spaceflight status conflicts with analysis sample status",
    source_value = isa_spaceflight_raw,
    canonical_value = sample_status,
    resolution = "Hard failure: resolve before plotting"
  )]
), fill = TRUE)
if (!nrow(conflict_audit)) {
  conflict_audit <- data.table(
    accession = character(), analysis_unit_id = character(), sample_column = character(),
    mouse_id = character(), issue_type = character(), source_value = character(),
    canonical_value = character(), resolution = character()
  )
}
setorder(conflict_audit, accession, analysis_unit_id, sample_column)
fwrite(conflict_audit, file.path(output_dir, "06_metadata_conflict_audit.csv"))

message(
  "OSDR_MOUSE_METADATA_REFRESH_PASS: ", nrow(download_manifest), " accessions; ",
  nrow(source_core), " complete ISA source/sample rows; ", nrow(mouse_output),
  " selected analysis sample columns; ",
  uniqueN(paste(mouse_output$accession, mouse_output$mouse_id)),
  " accession-scoped mouse identifiers"
)
