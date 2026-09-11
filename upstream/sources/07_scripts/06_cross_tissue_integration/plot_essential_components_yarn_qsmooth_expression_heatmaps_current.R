#!/usr/bin/env Rscript

# Tissue-aware expression heatmaps for the teacher's essential DNA-repair
# component set.
#
# The source OSDR/GeneLab tables contain sample-level normalized expression
# values with a documented +1 pseudocount.  They are not a common raw-count
# matrix, so TMM/CPM is not re-applied here.  Instead, the values are restored
# to the non-pseudocount scale (value - 1), assembled on a complete common
# Ensembl-gene universe, and passed to yarn::normalizeTissueAware() with its
# qsmooth method.  YARN's qsmooth(log=TRUE) then evaluates log2((value - 1)+1)
# and performs tissue-aware quantile-shrinkage normalization.  All joins use
# Ensembl Gene IDs; official symbols are display labels only.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)

# Prefer the project-scoped library created by setup_project_r_environment.R
# before loading YARN and its Bioconductor dependencies.
project_library <- file.path(root, "R", "library")
if (dir.exists(project_library)) .libPaths(unique(c(project_library, .libPaths())))

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(scales)
  library(Biobase)
  library(yarn)
  library(digest)
})

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

normalize_ensembl <- function(x) {
  sub("\\.[0-9]+$", "", trimws(as.character(x)))
}

strip_bom_names <- function(x) {
  sub("^\uFEFF", "", x)
}

safe_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

mean_or_na <- function(x) {
  keep <- is.finite(x)
  if (!any(keep)) NA_real_ else mean(x[keep])
}

summarise_vector <- function(x) {
  x <- as.numeric(x)
  q <- quantile(x, probs = c(0.01, 0.25, 0.50, 0.75, 0.99), names = FALSE, type = 7)
  data.table(
    Min = min(x), Q01 = q[[1L]], Q25 = q[[2L]], Median = q[[3L]],
    Q75 = q[[4L]], Q99 = q[[5L]], Max = max(x), Mean = mean(x), SD = sd(x)
  )
}

default_outdir <- file.path(
  root,
  "03_analysis_results",
  "18_essential_components_yarn_qsmooth_spearman_reordered_20260827"
)
outdir <- Sys.getenv("MOUSE_GO_QSMOOTH_OUTPUT_DIR", unset = default_outdir)
if (!grepl("^/", outdir)) outdir <- file.path(root, outdir)
figure_dir <- file.path(outdir, "figures")
table_dir <- file.path(outdir, "tables")
input_dir <- file.path(outdir, "input")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

teacher_xlsx_original <- file.path(root, "design/Essential_components_1_teacher_source.xlsx")
teacher_xlsx_frozen <- file.path(
  root,
  "03_analysis_results/16_essential_components_yarn_qsmooth_expression_heatmaps_20260826",
  "input/Essential_components_1_teacher_source.xlsx"
)
teacher_xlsx <- if (file.exists(teacher_xlsx_original)) teacher_xlsx_original else teacher_xlsx_frozen
unit_manifest_path <- file.path(
  root,
  "03_analysis_results/10_full_stable_id_rerun_20260823",
  "04_seven_pathway_response_26_tissues_20260824/prepared_inputs",
  "01_selected_tissue_analysis_units.csv"
)
expression_audit_path <- file.path(
  root,
  "03_analysis_results/06_cross_tissue_integration/05_proliferation_state/tables",
  "02_analysis_unit_expression_and_group_audit.csv"
)
component_mapping_path <- file.path(
  root,
  "03_analysis_results/14_essential_components_actb_normalized_expression_heatmaps_20260826",
  "tables/01_included_component_stable_id_mapping_audit.csv"
)
ensembl_mapping_path <- file.path(
  root,
  "08_GO_annotation/mgi_mod_gaf_20260804",
  "ensembl116_full_mgi_entrez_mapping.tsv.gz"
)
cross_tissue_order_path <- file.path(
  root,
  "03_analysis_results/06_cross_tissue_integration/03_mouse_go_tissue_extended",
  "tables/05_mouse_GO_tissue_twelve_sets_sorted_high_to_low.csv"
)

required_files <- c(
  teacher_xlsx, unit_manifest_path, expression_audit_path,
  component_mapping_path, ensembl_mapping_path, cross_tissue_order_path
)
assert(
  all(file.exists(required_files)),
  paste("Missing required input(s):", paste(required_files[!file.exists(required_files)], collapse = "; "))
)

# -------------------------------------------------------------------------
# 1. Freeze the teacher scope and the reviewed stable-ID map.
# -------------------------------------------------------------------------

source_sheets <- excel_sheets(teacher_xlsx)
assert(
  all(c("DSB", "SSB_others") %in% source_sheets),
  paste("The teacher workbook must contain DSB and SSB_others; found", paste(source_sheets, collapse = ", "))
)

read_component_sheet <- function(sheet_name, sheet_index) {
  x <- as.data.table(read_excel(teacher_xlsx, sheet = sheet_name, col_names = TRUE))
  setnames(x, strip_bom_names(trimws(names(x))))
  assert(
    identical(names(x), c("Pathway", "Essential components")),
    paste0("Worksheet ", sheet_name, " must contain exactly 'Pathway' and 'Essential components'.")
  )
  x[, `:=`(
    Sheet = sheet_name,
    SheetOrder = sheet_index,
    RowInSheet = seq_len(.N),
    Pathway = safe_text(Pathway),
    RequestedComponent = safe_text(`Essential components`)
  )]
  x[!is.na(Pathway) & nzchar(Pathway) & !is.na(RequestedComponent) & nzchar(RequestedComponent)]
}

teacher_all <- rbindlist(
  lapply(seq_along(source_sheets), function(i) read_component_sheet(source_sheets[[i]], i)),
  use.names = TRUE,
  fill = TRUE
)
assert(nrow(teacher_all) == 85L, "The teacher workbook must contain 85 nonempty rows.")
assert(teacher_all[Sheet == "DSB", .N] == 29L, "The DSB worksheet must contain 29 rows.")
assert(teacher_all[Sheet == "SSB_others", .N] == 56L, "The SSB_others worksheet must contain 56 rows.")
teacher_all[, `:=`(
  Included = !(Sheet == "SSB_others" & Pathway == "shared"),
  ExclusionReason = fifelse(
    Sheet == "SSB_others" & Pathway == "shared",
    "shared rows excluded according to the current instruction",
    ""
  )
)]
assert(teacher_all[Pathway == "shared", .N] == 10L, "Exactly 10 shared rows must be excluded.")
assert(teacher_all[Included == TRUE, .N] == 75L, "The included scope must contain 75 components.")

mapping_audit <- fread(component_mapping_path, showProgress = FALSE)
setnames(mapping_audit, strip_bom_names(names(mapping_audit)))
required_mapping_columns <- c(
  "ComponentOrder", "PathwayOrder", "ComponentOrderWithinPathway", "Pathway",
  "PathwayDisplay", "RequestedComponent", "MappingStatus", "LookupMouseSymbol",
  "DisplayComponent", "EnsemblID", "MGI_ID", "EntrezID", "OfficialMouseSymbol",
  "GeneBiotype", "Sheet", "RowInSheet"
)
assert(
  all(required_mapping_columns %in% names(mapping_audit)),
  "The reviewed component stable-ID audit schema has changed."
)
components <- mapping_audit
components[, EnsemblID := normalize_ensembl(EnsemblID)]
components[, `:=`(
  RequestedComponent = safe_text(RequestedComponent),
  Pathway = safe_text(Pathway),
  OfficialMouseSymbol = safe_text(OfficialMouseSymbol),
  MGI_ID = safe_text(MGI_ID),
  EntrezID = safe_text(EntrezID),
  DisplayComponent = safe_text(OfficialMouseSymbol),
  PathwayDisplay = fifelse(Pathway == "HR_A-EJ", "HR–A-EJ", Pathway)
)]
setorder(components, ComponentOrder)
assert(nrow(components) == 75L, "The reviewed stable-ID map must contain 75 included components.")
assert(!anyDuplicated(components$EnsemblID), "The reviewed stable-ID map contains duplicated Ensembl IDs.")
assert(
  all(grepl("^ENSMUSG[0-9]+$", components$EnsemblID)) &&
    all(nzchar(components$OfficialMouseSymbol)) && all(nzchar(components$MGI_ID)),
  "Every included component must have a valid Ensembl ID, MGI ID, and official mouse Symbol."
)

teacher_included <- teacher_all[Included == TRUE, .(Sheet, Pathway, RequestedComponent)]
map_scope <- components[, .(Sheet, Pathway, RequestedComponent)]
assert(
  nrow(merge(teacher_included, map_scope, by = c("Sheet", "Pathway", "RequestedComponent"))) == 75L,
  "The reviewed stable-ID map does not match all included teacher rows."
)

# -------------------------------------------------------------------------
# 2. Validate the fixed 26-tissue manifest and flight samples.
# -------------------------------------------------------------------------

unit_manifest <- fread(unit_manifest_path, showProgress = FALSE)
setnames(unit_manifest, strip_bom_names(names(unit_manifest)))
expression_audit <- fread(expression_audit_path, showProgress = FALSE)
setnames(expression_audit, strip_bom_names(names(expression_audit)))

assert(nrow(unit_manifest) == 59L, "The selected manifest must contain 59 analysis units.")
assert(uniqueN(unit_manifest$Group) == 26L, "The selected manifest must contain exactly 26 tissue groups.")
assert(!anyDuplicated(unit_manifest$analysis_unit_id), "Manifest analysis-unit IDs must be unique.")
assert(!anyDuplicated(expression_audit$analysis_unit_id), "Expression-audit analysis-unit IDs must be unique.")
assert(
  setequal(unit_manifest$analysis_unit_id, expression_audit$analysis_unit_id),
  "The manifest and flight-sample audit do not contain the same analysis units."
)
assert(
  all(c("input_relative_path", "flight_samples", "n_flight_samples") %in% names(expression_audit)),
  "Flight-sample audit schema has changed."
)
audit_by_unit <- setNames(seq_len(nrow(expression_audit)), expression_audit$analysis_unit_id)

unit_paths <- file.path(root, unit_manifest$input_relative_path)
assert(all(file.exists(unit_paths)), "At least one selected OSDR input file is missing.")

unit_profile_audit <- vector("list", nrow(unit_manifest))
sample_meta_list <- vector("list", nrow(unit_manifest))
common_universe <- NULL
total_sample_n <- 0L

for (index in seq_len(nrow(unit_manifest))) {
  unit <- unit_manifest[index]
  analysis_unit_id <- as.character(unit$analysis_unit_id)
  audit_index <- match(analysis_unit_id, expression_audit$analysis_unit_id)
  assert(length(audit_index) == 1L && !is.na(audit_index), paste("No unique expression audit row for", analysis_unit_id))
  audit_row <- expression_audit[audit_index]
  assert(
    identical(as.character(unit$input_relative_path), as.character(audit_row$input_relative_path)),
    paste("Manifest and expression-audit input paths disagree for", analysis_unit_id)
  )
  input_path <- file.path(root, as.character(unit$input_relative_path))
  flight_samples <- trimws(strsplit(as.character(audit_row$flight_samples), ";", fixed = TRUE)[[1L]])
  flight_samples <- flight_samples[nzchar(flight_samples)]
  assert(length(flight_samples) > 0L && !anyDuplicated(flight_samples), paste("Invalid flight sample list for", analysis_unit_id))
  assert(length(flight_samples) == as.integer(audit_row$n_flight_samples), paste("Flight sample count disagrees for", analysis_unit_id))

  header <- names(fread(input_path, nrows = 0L, check.names = FALSE, showProgress = FALSE))
  assert("ENSEMBL" %in% header, paste("Raw input has no ENSEMBL column for", analysis_unit_id))
  assert(all(flight_samples %chin% header), paste("Flight sample columns are missing for", analysis_unit_id))
  raw <- fread(input_path, select = c("ENSEMBL", flight_samples), check.names = FALSE, showProgress = FALSE)
  raw[, EnsemblID := normalize_ensembl(ENSEMBL)]
  valid_id <- grepl("^ENSMUSG[0-9]+$", raw$EnsemblID)
  assert(sum(valid_id) > 0L, paste("No valid mouse Ensembl IDs for", analysis_unit_id))
  assert(!anyDuplicated(raw$EnsemblID[valid_id]), paste("Duplicate valid Ensembl IDs for", analysis_unit_id))
  value_matrix <- as.matrix(raw[valid_id, ..flight_samples])
  suppressWarnings(storage.mode(value_matrix) <- "numeric")
  finite_row <- rowSums(!is.finite(value_matrix)) == 0L
  usable_ids <- raw$EnsemblID[valid_id][finite_row]
  assert(length(usable_ids) > 0L, paste("No finite Ensembl rows for", analysis_unit_id))
  assert(all(value_matrix[finite_row, , drop = FALSE] >= 1), paste("A source value below the +1 pseudocount floor was found for", analysis_unit_id))
  if (is.null(common_universe)) common_universe <- usable_ids else common_universe <- intersect(common_universe, usable_ids)

  actual_sha256 <- unname(digest(input_path, algo = "sha256", file = TRUE))
  manifest_sha256 <- as.character(unit$input_sha256)
  assert(identical(actual_sha256, manifest_sha256), paste("Input SHA256 disagrees with manifest for", analysis_unit_id))
  unit_profile_audit[[index]] <- data.table(
    Group = as.character(unit$Group), accession = as.character(unit$accession),
    analysis_unit_id = analysis_unit_id, mission_cluster = as.character(unit$mission_cluster),
    input_relative_path = as.character(unit$input_relative_path), input_sha256 = actual_sha256,
    raw_rows = nrow(raw), valid_unique_Ensembl_rows = sum(valid_id),
    finite_unique_Ensembl_rows = length(usable_ids), nonfinite_value_rows = sum(valid_id) - length(usable_ids),
    flight_sample_n = length(flight_samples), source_min = min(value_matrix), source_max = max(value_matrix),
    source_exactly_one_n = sum(value_matrix == 1), common_universe_rows_after_unit = length(common_universe)
  )
  sample_meta_list[[index]] <- data.table(
    Group = as.character(unit$Group), accession = as.character(unit$accession),
    analysis_unit_id = analysis_unit_id, mission_cluster = as.character(unit$mission_cluster),
    FlightSample = flight_samples
  )
  total_sample_n <- total_sample_n + length(flight_samples)
  if (index %% 10L == 0L || index == nrow(unit_manifest)) message(sprintf("Profiled %d/%d: %s; current common universe %d genes", index, nrow(unit_manifest), analysis_unit_id, length(common_universe)))
}

common_universe <- sort(unique(common_universe))
assert(length(common_universe) > 1000L, paste("The common Ensembl universe is unexpectedly small:", length(common_universe)))
unit_profile_audit <- rbindlist(unit_profile_audit, use.names = TRUE)
sample_meta <- rbindlist(sample_meta_list, use.names = TRUE)
sample_meta[, SampleID := make.unique(paste0(analysis_unit_id, "::", FlightSample))]
sample_meta[, SampleIndex := seq_len(.N)]
assert(nrow(sample_meta) == total_sample_n, "Sample metadata row count is inconsistent.")
assert(!anyDuplicated(sample_meta$SampleID), "Sample IDs must be unique after analysis-unit prefixing.")
assert(total_sample_n == sum(as.integer(expression_audit$n_flight_samples)), "Total flight sample count disagrees with expression audit.")

unit_profile_audit[, common_universe_rows_final := length(common_universe)]
unit_profile_audit[, retained_fraction_of_finite_rows := common_universe_rows_final / finite_unique_Ensembl_rows]

# -------------------------------------------------------------------------
# 3. Build a complete common-universe matrix and run YARN qsmooth.
# -------------------------------------------------------------------------

expression_after_pseudocount_removal <- matrix(
  NA_real_, nrow = length(common_universe), ncol = nrow(sample_meta),
  dimnames = list(common_universe, sample_meta$SampleID)
)
component_source_plus1 <- matrix(
  NA_real_, nrow = nrow(components), ncol = nrow(sample_meta),
  dimnames = list(components$EnsemblID, sample_meta$SampleID)
)

column_cursor <- 0L
for (index in seq_len(nrow(unit_manifest))) {
  unit <- unit_manifest[index]
  analysis_unit_id <- as.character(unit$analysis_unit_id)
  audit_index <- match(analysis_unit_id, expression_audit$analysis_unit_id)
  assert(length(audit_index) == 1L && !is.na(audit_index), paste("No unique expression audit row for", analysis_unit_id))
  audit_row <- expression_audit[audit_index]
  input_path <- file.path(root, as.character(unit$input_relative_path))
  flight_samples <- trimws(strsplit(as.character(audit_row$flight_samples), ";", fixed = TRUE)[[1L]])
  flight_samples <- flight_samples[nzchar(flight_samples)]
  raw <- fread(input_path, select = c("ENSEMBL", flight_samples), check.names = FALSE, showProgress = FALSE)
  raw[, EnsemblID := normalize_ensembl(ENSEMBL)]
  valid_id <- grepl("^ENSMUSG[0-9]+$", raw$EnsemblID)
  ids <- raw$EnsemblID[valid_id]
  assert(!anyDuplicated(ids), paste("Duplicate IDs while building matrix for", analysis_unit_id))
  value_matrix <- as.matrix(raw[valid_id, ..flight_samples])
  suppressWarnings(storage.mode(value_matrix) <- "numeric")
  row_index <- match(common_universe, ids)
  assert(all(!is.na(row_index)), paste("The common universe is not complete for", analysis_unit_id))
  source_values <- value_matrix[row_index, , drop = FALSE]
  assert(all(is.finite(source_values)) && all(source_values >= 1), paste("Invalid common-universe source value for", analysis_unit_id))
  col_range <- (column_cursor + 1L):(column_cursor + ncol(source_values))
  expression_after_pseudocount_removal[, col_range] <- source_values - 1

  # Target components that are not present in every table are retained as
  # explicit NA observations.  They do not enter the qsmooth training matrix,
  # but their source values are kept here for transparent coverage audits and
  # are used whenever that component exists in the current input table.
  target_index <- match(components$EnsemblID, ids)
  target_present <- !is.na(target_index)
  target_values <- matrix(NA_real_, nrow = nrow(components), ncol = ncol(value_matrix))
  if (any(target_present)) {
    target_values[target_present, ] <- value_matrix[target_index[target_present], , drop = FALSE]
    observed_target_values <- target_values[target_present, , drop = FALSE]
    assert(all(is.finite(observed_target_values)) && all(observed_target_values >= 1), paste("Invalid target component source value for", analysis_unit_id))
  }
  component_source_plus1[, col_range] <- target_values
  column_cursor <- column_cursor + ncol(source_values)
  if (index %% 10L == 0L || index == nrow(unit_manifest)) message(sprintf("Loaded %d/%d: %s", index, nrow(unit_manifest), analysis_unit_id))
}
assert(column_cursor == nrow(sample_meta), "Expression matrix columns do not match sample metadata.")
assert(all(is.finite(expression_after_pseudocount_removal)) && all(expression_after_pseudocount_removal >= 0), "The common expression matrix must be finite and non-negative after removing +1.")

input_log2_matrix <- log2(expression_after_pseudocount_removal + 1)
sample_groups <- sample_meta$Group
assert(length(sample_groups) == ncol(input_log2_matrix), "Tissue group labels do not match matrix columns.")
assert(uniqueN(sample_groups) == 26L, "YARN groups must contain exactly 26 tissues.")

set.seed(25)
eset <- ExpressionSet(
  assayData = expression_after_pseudocount_removal,
  phenoData = AnnotatedDataFrame(
    data = as.data.frame(sample_meta[, .(Group, accession, analysis_unit_id, mission_cluster, FlightSample, SampleIndex)], row.names = sample_meta$SampleID)
  )
)
assert(identical(colnames(exprs(eset)), sample_meta$SampleID), "ExpressionSet column names do not match sample metadata.")

normalization_method <- "qsmooth"
qsmooth_window <- 0.05
qsmooth_log <- TRUE
normalized_eset <- normalizeTissueAware(
  eset,
  groups = sample_groups,
  normalizationMethod = normalization_method,
  window = qsmooth_window,
  log = qsmooth_log
)
yarn_log2_matrix <- assayData(normalized_eset)[["normalizedMatrix"]]
assert(is.matrix(yarn_log2_matrix), "YARN did not return a normalized matrix.")
assert(identical(dim(yarn_log2_matrix), dim(input_log2_matrix)), "YARN matrix dimensions changed unexpectedly.")
assert(identical(rownames(yarn_log2_matrix), rownames(input_log2_matrix)) && identical(colnames(yarn_log2_matrix), colnames(input_log2_matrix)), "YARN matrix dimnames changed unexpectedly.")
assert(all(is.finite(yarn_log2_matrix)), "YARN normalized matrix contains non-finite values.")

# Recompute the qstats object for an inspectable weight audit.  This is the
# same deterministic qstats calculation used by yarn::qsmooth; no expression
# values are changed by this audit call.
qstats_audit_obj <- yarn:::qstats(input_log2_matrix, groups = sample_groups, window = qsmooth_window)
qsmooth_quantile_audit <- data.table(
  QuantileRank = seq_len(nrow(input_log2_matrix)),
  InputQuantileReference = qstats_audit_obj$Qref,
  SST = qstats_audit_obj$SST,
  SSB = qstats_audit_obj$SSB,
  RoughWeight = qstats_audit_obj$roughWeights,
  SmoothWeight = qstats_audit_obj$smoothWeights
)
assert(all(is.finite(qsmooth_quantile_audit$SmoothWeight)), "YARN qsmooth weights are non-finite.")

# Build the per-sample qsmooth reference quantiles once.  For a component that
# is absent from one or more source tables (and therefore cannot be part of the
# complete normalization universe), the same rank-based mapping is applied to
# its observed log2 value.  This is a transparent extension of the fitted
# qsmooth transform, not an imputation: missing source observations remain NA.
group_factor <- factor(as.character(sample_groups))
design_matrix <- model.matrix(~0 + group_factor)
qhat_matrix <- qstats_audit_obj$QBETAS %*% t(design_matrix)
qref_matrix <- matrix(qstats_audit_obj$Qref, nrow = nrow(input_log2_matrix), ncol = ncol(input_log2_matrix))
weight_matrix <- matrix(qstats_audit_obj$smoothWeights, nrow = nrow(input_log2_matrix), ncol = ncol(input_log2_matrix))
qsmooth_reference_quantiles <- weight_matrix * qref_matrix + (1 - weight_matrix) * qhat_matrix
assert(all(is.finite(qsmooth_reference_quantiles)), "The qsmooth reference quantiles are non-finite.")

# Sample-level distribution audit before and after qsmooth.
distribution_audit <- rbindlist(lapply(seq_len(ncol(input_log2_matrix)), function(j) {
  before <- summarise_vector(input_log2_matrix[, j])
  after <- summarise_vector(yarn_log2_matrix[, j])
  rbind(
    cbind(data.table(Stage = "input_log2_after_pseudocount_removal"), before),
    cbind(data.table(Stage = "yarn_qsmooth_log2"), after)
  )[, `:=`(
    SampleID = sample_meta$SampleID[[j]], Group = sample_meta$Group[[j]],
    accession = sample_meta$accession[[j]], analysis_unit_id = sample_meta$analysis_unit_id[[j]],
    mission_cluster = sample_meta$mission_cluster[[j]], FlightSample = sample_meta$FlightSample[[j]],
    SampleIndex = sample_meta$SampleIndex[[j]]
  )]
}), use.names = TRUE, fill = TRUE)
setcolorder(distribution_audit, c("SampleIndex", "SampleID", "Group", "accession", "analysis_unit_id", "mission_cluster", "FlightSample", "Stage", "Min", "Q01", "Q25", "Median", "Q75", "Q99", "Max", "Mean", "SD"))

# -------------------------------------------------------------------------
# 4. Extract components and aggregate sample -> unit -> mission -> tissue.
# -------------------------------------------------------------------------

component_row_index <- match(components$EnsemblID, common_universe)
component_in_normalization_universe <- !is.na(component_row_index)
component_input_log2 <- log2(component_source_plus1)
component_yarn_log2 <- matrix(
  NA_real_, nrow = nrow(components), ncol = ncol(yarn_log2_matrix),
  dimnames = list(components$EnsemblID, colnames(yarn_log2_matrix))
)
if (any(component_in_normalization_universe)) {
  component_yarn_log2[component_in_normalization_universe, ] <- yarn_log2_matrix[component_row_index[component_in_normalization_universe], , drop = FALSE]
}
assert(
  all(is.finite(component_source_plus1[is.finite(component_source_plus1)])) &&
    all(component_source_plus1[is.finite(component_source_plus1)] >= 1),
  "Observed component source values must be finite and >=1."
)

quantile_transfer_one <- function(x, sorted_input, normalized_reference) {
  if (!is.finite(x)) return(NA_real_)
  n <- length(sorted_input)
  left_n <- findInterval(x, sorted_input, left.open = TRUE)
  right_n <- findInterval(x, sorted_input)
  if (right_n > left_n) {
    mean(normalized_reference[(left_n + 1L):right_n])
  } else {
    normalized_reference[min(n, max(1L, left_n + 1L))]
  }
}

if (any(!component_in_normalization_universe)) {
  for (k in which(!component_in_normalization_universe)) {
    for (j in seq_len(ncol(component_yarn_log2))) {
      if (is.finite(component_input_log2[k, j])) {
        component_yarn_log2[k, j] <- quantile_transfer_one(
          component_input_log2[k, j], qstats_audit_obj$Q[, j], qsmooth_reference_quantiles[, j]
        )
      }
    }
  }
}

component_flight <- data.table(
  EnsemblID = rep(components$EnsemblID, each = ncol(component_yarn_log2)),
  SampleID = rep(colnames(component_yarn_log2), times = nrow(component_yarn_log2)),
  SourceExpressionPlus1 = as.vector(t(component_source_plus1)),
  ExpressionAfterPseudocountRemoval = as.vector(t(component_source_plus1 - 1)),
  InputExpressionLog2 = as.vector(t(component_input_log2)),
  YARNNormalizedLog2 = as.vector(t(component_yarn_log2))
)
component_flight[, `:=`(
  ComponentOrder = match(EnsemblID, components$EnsemblID),
  ComponentPresentInRaw = is.finite(SourceExpressionPlus1),
  ComponentPresentInNormalizationUniverse = rep(component_in_normalization_universe, each = ncol(component_yarn_log2)),
  YARNValueDerivation = ifelse(
    rep(component_in_normalization_universe, each = ncol(component_yarn_log2)),
    "direct_yarn_normalizeTissueAware_qsmooth",
    "rank_transfer_from_fitted_qsmooth_reference"
  )
)]
component_flight <- merge(
  component_flight,
  sample_meta[, .(SampleID, SampleIndex, Group, accession, analysis_unit_id, mission_cluster, FlightSample)],
  by = "SampleID", all.x = TRUE, sort = FALSE
)
component_flight <- merge(
  component_flight,
  components[, .(ComponentOrder, Pathway, PathwayOrder, ComponentOrderWithinPathway, PathwayDisplay, RequestedComponent, MappingStatus, DisplayComponent, OfficialMouseSymbol, EnsemblID, MGI_ID, EntrezID)],
  by = c("ComponentOrder", "EnsemblID"), all.x = TRUE, sort = FALSE
)
setorder(component_flight, SampleIndex, ComponentOrder)
assert(nrow(component_flight) == nrow(sample_meta) * nrow(components), "Component flight table has unexpected dimensions.")
assert(
  all(is.finite(component_flight$YARNNormalizedLog2[component_flight$ComponentPresentInNormalizationUniverse])),
  "Component-level YARN values are non-finite for components in the normalization universe."
)

unit_yarn <- component_flight[
  , .(
    UnitMeanYARNNormalizedLog2 = mean_or_na(YARNNormalizedLog2),
    UnitMeanInputExpressionLog2 = mean_or_na(InputExpressionLog2),
    UnitMeanSourceExpressionPlus1 = mean_or_na(SourceExpressionPlus1),
    FlightSampleN = .N,
    ObservedFlightSampleN = sum(is.finite(YARNNormalizedLog2))
  ),
  by = .(
    Group, accession, analysis_unit_id, mission_cluster, Pathway, PathwayOrder,
    PathwayDisplay, ComponentOrder, ComponentOrderWithinPathway, RequestedComponent,
    MappingStatus, DisplayComponent, OfficialMouseSymbol, EnsemblID, MGI_ID, EntrezID
  )
]
mission_yarn <- unit_yarn[
  , .(
    MissionMeanUnitYARNNormalizedLog2 = mean_or_na(UnitMeanYARNNormalizedLog2),
    MissionMeanUnitInputExpressionLog2 = mean_or_na(UnitMeanInputExpressionLog2),
    MissionMeanUnitSourceExpressionPlus1 = mean_or_na(UnitMeanSourceExpressionPlus1),
    MissionAnalysisUnitN = .N,
    MissionObservedAnalysisUnitN = sum(is.finite(UnitMeanYARNNormalizedLog2)),
    MissionFlightSampleN = sum(FlightSampleN)
  ),
  by = .(
    Group, mission_cluster, Pathway, PathwayOrder, PathwayDisplay, ComponentOrder,
    ComponentOrderWithinPathway, RequestedComponent, MappingStatus, DisplayComponent,
    OfficialMouseSymbol, EnsemblID, MGI_ID, EntrezID
  )
]
tissue_yarn <- mission_yarn[
  , .(
    TissueYARNNormalizedLog2 = mean_or_na(MissionMeanUnitYARNNormalizedLog2),
    TissueInputExpressionLog2 = mean_or_na(MissionMeanUnitInputExpressionLog2),
    TissueSourceExpressionPlus1 = mean_or_na(MissionMeanUnitSourceExpressionPlus1),
    MissionN = .N,
    ObservedMissionN = sum(is.finite(MissionMeanUnitYARNNormalizedLog2)),
    AnalysisUnitN = sum(MissionAnalysisUnitN),
    FlightSampleN = sum(MissionFlightSampleN)
  ),
  by = .(
    Group, Pathway, PathwayOrder, PathwayDisplay, ComponentOrder,
    ComponentOrderWithinPathway, RequestedComponent, MappingStatus, DisplayComponent,
    OfficialMouseSymbol, EnsemblID, MGI_ID, EntrezID
  )
]
assert(
  nrow(unit_yarn) == nrow(unit_manifest) * nrow(components) &&
    !anyDuplicated(unit_yarn[, .(analysis_unit_id, EnsemblID)]) &&
    nrow(mission_yarn) > 0L,
  "Unit-level YARN aggregation is incomplete."
)
assert(
  nrow(tissue_yarn) == 26L * 75L && uniqueN(tissue_yarn$Group) == 26L &&
    !anyDuplicated(tissue_yarn[, .(Group, EnsemblID)]) &&
    all(is.finite(tissue_yarn$TissueYARNNormalizedLog2)) &&
    all(tissue_yarn$ObservedMissionN > 0L),
  "The final YARN tissue matrix must contain exactly 26 x 75 finite cells."
)

# -------------------------------------------------------------------------
# 5. Tissue labels, audits, figures, and provenance.
# -------------------------------------------------------------------------

tissue_labels <- c(
  "Adrenal gland" = "肾上腺", "Bone marrow" = "骨髓", "Cecum" = "盲肠", "Cerebellum" = "小脑",
  "Colon" = "结肠", "Dorsal skin" = "背部皮肤", "Extensor digitorum longus" = "趾长伸肌", "Eye" = "眼",
  "Femoral lateral skin" = "股外侧皮肤", "Femoral skin" = "股部皮肤", "Gastrocnemius" = "腓肠肌",
  "Heart" = "心脏", "Heart / Heart right ventricle" = "右心室", "Kidney" = "肾脏",
  "Left lobe of the liver" = "肝左叶", "Liver" = "肝脏", "Lung" = "肺", "Mammary gland" = "乳腺",
  "Optic nerve" = "视神经", "Quadriceps femoris" = "股四头肌", "Retina" = "视网膜", "Soleus" = "比目鱼肌",
  "Spleen" = "脾脏", "Spleen-distal" = "脾脏（远端）", "Thymus" = "胸腺", "Tibialis anterior" = "胫骨前肌"
)
cross_tissue_order <- fread(cross_tissue_order_path, showProgress = FALSE)
setnames(cross_tissue_order, strip_bom_names(names(cross_tissue_order)))
assert(
  all(c("NES_sort_order", "analysis_tissue") %in% names(cross_tissue_order)) &&
    nrow(cross_tissue_order) == 26L &&
    !anyDuplicated(cross_tissue_order$analysis_tissue) &&
    !anyDuplicated(cross_tissue_order$NES_sort_order),
  "The cross-tissue heatmap order table must contain 26 unique ordered tissues."
)
setorder(cross_tissue_order, NES_sort_order)
tissue_order <- trimws(as.character(cross_tissue_order$analysis_tissue))
tissue_display_names <- setNames(tissue_order, tissue_order)
tissue_display_names["Heart / Heart right ventricle"] <- "Heart right ventricle"
assert(
  setequal(tissue_order, unique(unit_manifest$Group)),
  "The cross-tissue heatmap order does not cover exactly the 26 qsmooth tissues."
)
assert(setequal(names(tissue_labels), tissue_order), "The bilingual tissue-label map must cover exactly 26 tissues.")
tissue_yarn[, `:=`(
  TissueChinese = unname(tissue_labels[Group]),
  TissueLabel = paste0(unname(tissue_display_names[Group]), "\n", unname(tissue_labels[Group])),
  TissueOrder = match(Group, tissue_order)
)]
setorder(tissue_yarn, TissueOrder, ComponentOrder)

yarn_min <- min(tissue_yarn$TissueYARNNormalizedLog2)
yarn_max <- max(tissue_yarn$TissueYARNNormalizedLog2)
assert(yarn_max > yarn_min, "The YARN tissue matrix has no colour-scale span.")
global_colour_audit <- data.table(
  Matrix = "YARN qsmooth normalized log2 expression",
  CellN = nrow(tissue_yarn), TissueN = uniqueN(tissue_yarn$Group), ComponentN = uniqueN(tissue_yarn$EnsemblID),
  GlobalMinimum = yarn_min, GlobalMaximum = yarn_max,
  MinimumColour = "blue (#2166AC)", MidpointColour = "white (#F7F7F7)", MaximumColour = "red (#B2182B)",
  CellLabelDigits = 1L,
  ColourScope = "all 26 tissues x all 75 included components", DSBAndSSBShareScale = TRUE
)

component_context <- components[, .(
  ComponentOrder, PathwayOrder, ComponentOrderWithinPathway, Pathway, PathwayDisplay,
  RequestedComponent, MappingStatus, LookupMouseSymbol, DisplayComponent, EnsemblID,
  MGI_ID, EntrezID, OfficialMouseSymbol, GeneBiotype, Sheet, RowInSheet
)]
# The source pathway assignment remains authoritative.  For the DSB figure
# only, Exo1 and Dna2 are displayed with HR (between Rad54l and Rbbp8).  This
# is a presentation-level reassignment; no source annotation, ID, or
# expression value is changed.
reassigned_symbols <- component_context[
  toupper(RequestedComponent) %chin% c("EXO1", "DNA2"), OfficialMouseSymbol
]
assert(
  setequal(reassigned_symbols, c("Exo1", "Dna2")) && length(reassigned_symbols) == 2L,
  "The reviewed component map must contain exactly Exo1 and Dna2 for the requested display reassignment."
)
plot_component_order <- components$OfficialMouseSymbol
plot_component_order <- plot_component_order[!plot_component_order %chin% reassigned_symbols]
rad54l_position <- match("Rad54l", plot_component_order)
rbbp8_position <- match("Rbbp8", plot_component_order)
assert(
  !is.na(rad54l_position) && !is.na(rbbp8_position) && rbbp8_position == rad54l_position + 1L,
  "Rad54l and Rbbp8 must be adjacent before inserting Exo1 and Dna2."
)
plot_component_order <- append(plot_component_order, reassigned_symbols, after = rad54l_position)
assert(
  length(plot_component_order) == nrow(components) &&
    setequal(plot_component_order, components$OfficialMouseSymbol) &&
    all(plot_component_order[(rad54l_position + 1L):(rad54l_position + 2L)] == reassigned_symbols) &&
    match("Rbbp8", plot_component_order) == rad54l_position + 3L,
  "The plotted component order did not place Exo1/Dna2 between Rad54l and Rbbp8."
)
plot_component_context <- copy(component_context)
plot_component_context[, PlotPathwayDisplay := PathwayDisplay]
plot_component_context[OfficialMouseSymbol %chin% reassigned_symbols, PlotPathwayDisplay := "HR"]
plot_component_context[, PlotComponentOrder := match(OfficialMouseSymbol, plot_component_order)]
assert(!anyNA(plot_component_context$PlotComponentOrder), "Every component must have a plotted component order.")

# `tissue_yarn` already carries the fields used to aggregate the component
# values.  Add only the remaining mapping context here so that merge() does
# not create `.x`/`.y` duplicates for DisplayComponent or MappingStatus.
tissue_matrix <- merge(
  tissue_yarn,
  plot_component_context[, .(ComponentOrder, LookupMouseSymbol, GeneBiotype, Sheet, RowInSheet, PlotPathwayDisplay, PlotComponentOrder)],
  by = "ComponentOrder", all.x = TRUE, sort = FALSE
)
setorder(tissue_matrix, TissueOrder, ComponentOrder)
assert(
  nrow(tissue_matrix) == 26L * 75L && !anyNA(tissue_matrix$DisplayComponent) &&
    !anyNA(tissue_matrix$PlotPathwayDisplay) && !anyNA(tissue_matrix$PlotComponentOrder),
  "Tissue matrix stable-ID context is incomplete."
)

plot_component_reassignment_audit <- plot_component_context[, .(
  ComponentOrder, PlotComponentOrder, RequestedComponent, OfficialMouseSymbol, EnsemblID,
  OriginalPathway = Pathway, OriginalPathwayDisplay = PathwayDisplay,
  PlottedPathwayDisplay = PlotPathwayDisplay,
  DisplayOrderRule = fifelse(
    OfficialMouseSymbol %chin% reassigned_symbols,
    "visually moved into HR between Rad54l and Rbbp8; source mapping unchanged",
    "source component order"
  )
)]
setorder(plot_component_reassignment_audit, PlotComponentOrder)

missing_component_audit <- component_flight[
  ComponentPresentInRaw == FALSE,
  .(
    MissingFlightSampleN = .N,
    MissingAnalysisUnitN = uniqueN(analysis_unit_id),
    MissingTissueN = uniqueN(Group),
    MissingSampleNames = paste(unique(SampleID), collapse = ";")
  ),
  by = .(Group, analysis_unit_id, RequestedComponent, OfficialMouseSymbol, EnsemblID)
]
if (!nrow(missing_component_audit)) {
  missing_component_audit <- data.table(
    Group = character(), analysis_unit_id = character(), RequestedComponent = character(),
    OfficialMouseSymbol = character(), EnsemblID = character(), MissingFlightSampleN = integer(),
    MissingAnalysisUnitN = integer(), MissingTissueN = integer(), MissingSampleNames = character()
  )
}
missing_component_sample_cells <- sum(!component_flight$ComponentPresentInRaw)
partial_tissue_cells <- sum(tissue_matrix$ObservedMissionN < tissue_matrix$MissionN)

# Direct comparison to the Actb version is sensitivity-only; it does not
# affect the YARN values or colour scales.
actb_tissue_path <- file.path(root, "03_analysis_results/14_essential_components_actb_normalized_expression_heatmaps_20260826/tables/06_component_tissue_actb_ratio_matrix.csv")
actb_compare <- data.table()
actb_compare_summary <- data.table()
if (file.exists(actb_tissue_path)) {
  actb <- fread(actb_tissue_path, showProgress = FALSE)
  actb[, EnsemblID := normalize_ensembl(EnsemblID)]
  actb <- actb[, .(Group, EnsemblID, ActbNormalizedLog2 = TissueLog2RatioPlus1)]
  actb_compare <- merge(tissue_matrix[, .(Group, EnsemblID, YARNNormalizedLog2 = TissueYARNNormalizedLog2)], actb, by = c("Group", "EnsemblID"), all = FALSE)
  actb_compare_summary <- actb_compare[, .(
    CellN = .N,
    SpearmanRho = suppressWarnings(cor(YARNNormalizedLog2, ActbNormalizedLog2, method = "spearman")),
    PearsonR = suppressWarnings(cor(YARNNormalizedLog2, ActbNormalizedLog2, method = "pearson")),
    MeanAbsoluteDifference = mean(abs(YARNNormalizedLog2 - ActbNormalizedLog2))
  ), by = Group]
  actb_compare_summary <- rbind(
    actb_compare_summary,
    actb_compare[, .(
      Group = "ALL_TISSUES", CellN = .N,
      SpearmanRho = suppressWarnings(cor(YARNNormalizedLog2, ActbNormalizedLog2, method = "spearman")),
      PearsonR = suppressWarnings(cor(YARNNormalizedLog2, ActbNormalizedLog2, method = "pearson")),
      MeanAbsoluteDifference = mean(abs(YARNNormalizedLog2 - ActbNormalizedLog2))
    )],
    fill = TRUE
  )
}

# Reorder tissue rows independently for the DSB and SSB figures using the
# observed component profiles.  The distance is 1 - Spearman rho between two
# tissue profiles, and average-linkage hierarchical clustering supplies the
# leaf order.  This changes display order only; it does not alter any value or
# aggregation step above.
spearman_cluster_order <- function(plot_data, figure_name, tissue_universe) {
  profile <- dcast(
    plot_data[, .(Group, OfficialMouseSymbol, TissueYARNNormalizedLog2)],
    Group ~ OfficialMouseSymbol,
    value.var = "TissueYARNNormalizedLog2"
  )
  profile <- profile[match(tissue_universe, Group)]
  assert(
    nrow(profile) == length(tissue_universe) && all(profile$Group == tissue_universe),
    paste("The", figure_name, "profile does not cover the expected 26 tissues.")
  )
  component_columns <- setdiff(names(profile), "Group")
  assert(length(component_columns) > 1L, paste("The", figure_name, "profile has too few components for correlation clustering."))
  profile_matrix <- as.matrix(profile[, ..component_columns])
  suppressWarnings(storage.mode(profile_matrix) <- "numeric")
  assert(all(is.finite(profile_matrix)), paste("The", figure_name, "profile contains non-finite values."))
  rownames(profile_matrix) <- profile$Group
  rho <- suppressWarnings(cor(t(profile_matrix), method = "spearman", use = "pairwise.complete.obs"))
  assert(all(is.finite(rho)), paste("The", figure_name, "Spearman correlation matrix contains non-finite values."))
  # Keep the matrix dimensions and dimnames while guarding against tiny
  # floating-point excursions outside the correlation interval.
  rho[rho > 1] <- 1
  rho[rho < -1] <- -1
  distance <- 1 - rho
  distance <- (distance + t(distance)) / 2
  diag(distance) <- 0
  hc <- hclust(as.dist(distance), method = "average")
  ordered_tissues <- rownames(profile_matrix)[hc$order]
  assert(length(ordered_tissues) == 26L && setequal(ordered_tissues, tissue_universe), paste("The", figure_name, "cluster order is incomplete."))
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
  list(order = ordered_tissues, order_audit = order_audit, correlation = correlation_table)
}

dsb_plot_data <- tissue_matrix[PlotPathwayDisplay %chin% c("NHEJ", "HR", "HR–A-EJ", "A-EJ")]
ssb_plot_data <- tissue_matrix[PlotPathwayDisplay %chin% c("BER", "NER", "MMR", "FA")]
dsb_spearman <- spearman_cluster_order(dsb_plot_data, "DSB", tissue_order)
ssb_spearman <- spearman_cluster_order(ssb_plot_data, "SSB", tissue_order)
tissue_spearman_cluster_order_audit <- rbindlist(
  list(dsb_spearman$order_audit, ssb_spearman$order_audit), use.names = TRUE
)
fwrite(tissue_spearman_cluster_order_audit, file.path(table_dir, "17_tissue_spearman_cluster_order.csv"))
fwrite(dsb_spearman$correlation, file.path(table_dir, "18_tissue_spearman_correlation_DSB.csv"))
fwrite(ssb_spearman$correlation, file.path(table_dir, "19_tissue_spearman_correlation_SSB.csv"))
fwrite(plot_component_reassignment_audit, file.path(table_dir, "20_figure_component_reassignment_audit.csv"))

setcolorder(teacher_all, c("Sheet", "SheetOrder", "RowInSheet", "Pathway", "RequestedComponent", "Included", "ExclusionReason"))
fwrite(teacher_all, file.path(table_dir, "00_teacher_component_scope_audit.csv"))
fwrite(component_context, file.path(table_dir, "01_included_component_stable_id_mapping_audit.csv"))
fwrite(unit_profile_audit, file.path(table_dir, "03_common_ensembl_universe_unit_audit.csv"))
fwrite(data.table(EnsemblID = common_universe, CommonUniverseOrder = seq_along(common_universe)), file.path(table_dir, "04_common_ensembl_gene_universe.csv"))
fwrite(sample_meta, file.path(table_dir, "05_flight_sample_metadata.csv"))
fwrite(component_flight, file.path(table_dir, "06_component_flight_sample_yarn_qsmooth_values.csv.gz"))
fwrite(unit_yarn, file.path(table_dir, "07_component_analysis_unit_yarn_qsmooth_log2.csv"))
fwrite(mission_yarn, file.path(table_dir, "08_component_mission_yarn_qsmooth_log2.csv"))
fwrite(tissue_matrix, file.path(table_dir, "09_component_tissue_yarn_qsmooth_log2_matrix.csv"))
fwrite(global_colour_audit, file.path(table_dir, "10_global_colour_scale_audit.csv"))
fwrite(qsmooth_quantile_audit, file.path(table_dir, "11_qsmooth_quantile_weight_audit.csv"))
fwrite(distribution_audit, file.path(table_dir, "12_sample_distribution_before_after_qsmooth.csv"))
fwrite(missing_component_audit, file.path(table_dir, "13_missing_component_presence_audit.csv"))
if (nrow(actb_compare)) fwrite(actb_compare, file.path(table_dir, "14_yarn_vs_actb_tissue_cell_comparison.csv"))
if (nrow(actb_compare_summary)) fwrite(actb_compare_summary, file.path(table_dir, "15_yarn_vs_actb_tissue_correlation_summary.csv"))

integrity_checks <- data.table(
  Check = c(
    "teacher nonempty rows", "included component rows", "excluded shared rows", "analysis units", "tissue groups",
    "flight samples", "common Ensembl universe genes", "common matrix cells", "component flight cells",
    "tissue matrix cells", "unique tissue-Ensembl cells", "source values below 1", "source values exactly 1",
    "common matrix finite cells", "YARN matrix finite cells", "YARN smooth weight finite rows",
    "missing component sample cells", "tissue cells with zero mission coverage", "tissue cells with partial mission coverage",
    "YARN global minimum", "YARN global maximum"
  ),
  Observed = c(
    nrow(teacher_all), nrow(components), teacher_all[Included == FALSE, .N], nrow(unit_manifest), uniqueN(unit_manifest$Group),
    nrow(sample_meta), length(common_universe), length(expression_after_pseudocount_removal), nrow(component_flight),
    nrow(tissue_matrix), uniqueN(tissue_matrix[, .(Group, EnsemblID)]), sum(expression_after_pseudocount_removal < 0),
    sum(expression_after_pseudocount_removal == 0), sum(is.finite(expression_after_pseudocount_removal)), sum(is.finite(yarn_log2_matrix)),
    sum(is.finite(qsmooth_quantile_audit$SmoothWeight)), missing_component_sample_cells,
    sum(tissue_matrix$ObservedMissionN == 0L), partial_tissue_cells, yarn_min, yarn_max
  ),
  Expected = c(
    85L, 75L, 10L, 59L, 26L, nrow(sample_meta), length(common_universe), length(expression_after_pseudocount_removal), nrow(sample_meta) * 75L,
    1950L, 1950L, 0L, NA, length(expression_after_pseudocount_removal), length(yarn_log2_matrix), length(common_universe), NA,
    0L, NA, yarn_min, yarn_max
  ),
  Status = c(
    if (nrow(teacher_all) == 85L) "PASS" else "FAIL", if (nrow(components) == 75L) "PASS" else "FAIL",
    if (teacher_all[Included == FALSE, .N] == 10L) "PASS" else "FAIL", if (nrow(unit_manifest) == 59L) "PASS" else "FAIL",
    if (uniqueN(unit_manifest$Group) == 26L) "PASS" else "FAIL", if (nrow(sample_meta) == total_sample_n) "PASS" else "FAIL",
    "recorded", if (length(expression_after_pseudocount_removal) == length(common_universe) * nrow(sample_meta)) "PASS" else "FAIL",
    if (nrow(component_flight) == nrow(sample_meta) * 75L) "PASS" else "FAIL", if (nrow(tissue_matrix) == 1950L) "PASS" else "FAIL",
    if (uniqueN(tissue_matrix[, .(Group, EnsemblID)]) == 1950L) "PASS" else "FAIL", if (!any(expression_after_pseudocount_removal < 0)) "PASS" else "FAIL",
    "recorded", if (all(is.finite(expression_after_pseudocount_removal))) "PASS" else "FAIL", if (all(is.finite(yarn_log2_matrix))) "PASS" else "FAIL",
    if (sum(is.finite(qsmooth_quantile_audit$SmoothWeight)) == length(common_universe)) "PASS" else "FAIL",
    "recorded", if (!any(tissue_matrix$ObservedMissionN == 0L)) "PASS" else "FAIL", "recorded", "recorded", "recorded"
  )
)
fwrite(integrity_checks, file.path(table_dir, "15_integrity_checks.csv"))

metadata <- data.table(
  Field = c(
    "Teacher workbook", "Teacher worksheets used", "DSB rows", "SSB_others rows", "shared rows excluded", "Included components",
    "Tissue groups", "OSDR accessions", "Analysis units", "Flight samples", "Common Ensembl gene universe",
    "Stable-ID analysis key", "Stable-ID mapping source", "Raw expression source", "Source value contract",
    "Normalization method", "YARN version", "Bioconductor version", "R version", "YARN groups", "Tissue row order source",
    "Tissue row order rule", "Figure component display reassignment", "Spearman clustering method", "Spearman distance",
    "YARN log argument", "YARN qsmooth window", "Pseudocount handling", "Aggregation hierarchy",
    "Tissue aggregation", "Figure colour rule", "Figure labels", "Random seed", "TMM/CPM status"
  ),
  Value = c(
    teacher_xlsx, paste(source_sheets, collapse = "; "), "29", "56", "10", "75 (29 DSB + 46 SSB)",
    "26", uniqueN(unit_manifest$accession), nrow(unit_manifest), nrow(sample_meta), length(common_universe),
    "Ensembl Gene ID", component_mapping_path, "OSDR/GeneLab differential-expression tables; flight sample columns",
    "All selected source values are finite and >=1; source value 1 is treated as the documented +1 pseudocount floor",
    "yarn::normalizeTissueAware(..., normalizationMethod='qsmooth')", as.character(packageVersion("yarn")), as.character(BiocManager::version()), R.version.string,
    "26 tissue labels, one group per tissue", cross_tissue_order_path,
    "DSB and SSB independently ordered by the corresponding figure's tissue profile",
    "Exo1/Dna2 moved from HR-A-EJ display to HR between Rad54l/Rbbp8; source Pathway retained",
    "Average-linkage hierarchical clustering", "1 - Spearman rho across displayed components",
    "TRUE (YARN computes log2((value - 1) + 1))", as.character(qsmooth_window),
    "Subtract 1 before YARN; no second external +1; YARN's internal +1 restores the original log2 scale",
    "Flight samples mean within analysis unit -> analysis-unit mean within mission -> equal mission mean within tissue",
    "Equal-weight mean of mission-level normalized log2 values within tissue",
    sprintf("One shared blue-white-red scale for DSB and SSB: global minimum %.12g = blue; mid-range = white; global maximum %.12g = red across 26 x 75 tissue cells; each cell labelled to 1 decimal", yarn_min, yarn_max),
    "Official mouse Symbols only on x-axis; each finite cell has a one-decimal label; stable IDs remain in audit tables", "25",
    "Not run: complete raw-count matrices are not available for all 59 units; TMM/CPM is reserved for a future raw-count version"
  )
)
fwrite(metadata, file.path(table_dir, "16_run_metadata.csv"))

saveRDS(
  list(
    EnsemblID = common_universe,
    SampleMetadata = sample_meta,
    YARNNormalizedLog2 = yarn_log2_matrix,
    Normalization = list(method = normalization_method, groups = "26 tissue labels", log = qsmooth_log, window = qsmooth_window, seed = 25L)
  ),
  file.path(outdir, "yarn_qsmooth_common_universe_normalized_matrix.rds"), compress = "xz"
)

make_plot <- function(plot_data, title, subtitle, caption, plot_tissue_order) {
  plot_data <- copy(plot_data)
  assert(length(plot_tissue_order) == 26L && setequal(plot_tissue_order, tissue_order), "A figure tissue order must cover exactly the 26 tissues.")
  display_tissue_levels_for_plot <- setNames(
    paste0(unname(tissue_display_names[plot_tissue_order]), "\n", unname(tissue_labels[plot_tissue_order])), plot_tissue_order
  )
  plot_data[, DisplayComponent := factor(OfficialMouseSymbol, levels = plot_component_order)]
  plot_data[, PlotPathwayDisplay := factor(PlotPathwayDisplay, levels = c("NHEJ", "HR", "HR–A-EJ", "A-EJ", "BER", "NER", "MMR", "FA"))]
  plot_data[, TissueLabel := factor(
    paste0(unname(tissue_display_names[Group]), "\n", unname(tissue_labels[Group])),
    levels = unname(display_tissue_levels_for_plot[rev(plot_tissue_order)])
  )]
  # Use white text only on the darkest tails of the shared blue-white-red
  # scale; dark text remains readable in the pale/mid-range cells.
  plot_data[, CellTextColour := fifelse(
    TissueYARNNormalizedLog2 <= yarn_min + 0.16 * (yarn_max - yarn_min) |
      TissueYARNNormalizedLog2 >= yarn_min + 0.84 * (yarn_max - yarn_min),
    "white", "#1F2937"
  )]
  ggplot(plot_data, aes(x = DisplayComponent, y = TissueLabel, fill = TissueYARNNormalizedLog2)) +
    geom_tile(colour = "white", linewidth = 0.16) +
    geom_text(
      aes(label = sprintf("%.1f", TissueYARNNormalizedLog2), colour = CellTextColour),
      size = 1.55, fontface = "plain", na.rm = TRUE, show.legend = FALSE
    ) +
    scale_colour_identity() +
    facet_grid(. ~ PlotPathwayDisplay, scales = "free_x", space = "free_x", switch = "x", drop = TRUE) +
    scale_fill_gradientn(
      colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
      values = rescale(seq(yarn_min, yarn_max, length.out = 5L), from = c(yarn_min, yarn_max)),
      limits = c(yarn_min, yarn_max), oob = squish, name = "YARN qsmooth\nlog2 expression",
      breaks = pretty(c(yarn_min, yarn_max), n = 4), labels = function(x) formatC(x, format = "fg", digits = 3)
    ) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, caption = caption) +
    theme_minimal(base_size = 8.8, base_family = "Arial Unicode MS") +
    theme(
      panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.3, colour = "#3D3D3D"),
      axis.text.y = element_text(size = 7.1, lineheight = 0.84, colour = "#3D3D3D"), axis.ticks = element_blank(),
      strip.background = element_rect(fill = "#F0F0F0", colour = NA), strip.text.x.bottom = element_text(face = "bold", size = 9.2),
      strip.placement = "outside", panel.spacing.x = grid::unit(0.42, "lines"), plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 8.7), plot.caption = element_text(size = 7.2, hjust = 0, colour = "#444444"),
      legend.position = "right", plot.margin = margin(8, 12, 10, 8)
    )
}

common_caption <- paste(
  "Source flight values were reduced by 1 and normalized with YARN qsmooth using all finite genes in the common Ensembl universe.",
  "Values were averaged within analysis unit, mission, and then equally across missions within tissue.",
  "Blue = global minimum, white = mid-range, and red = global maximum across all 26 tissues x 75 components; each finite cell is labelled with one decimal place.",
  "Colours are not ranks or percentiles. Rows are ordered by average-linkage clustering of 1 - Spearman rho across the displayed components, separately for DSB and SSB.",
  "All joins use Ensembl Gene ID; x-axis labels are official mouse Symbols; Actb is not used as a denominator."
)

dsb_plot <- make_plot(
  dsb_plot_data,
  "DSB essential components: YARN qsmooth tissue-aware flight RNA expression",
  sprintf("26 mouse tissues; 29 DSB components; tissues ordered by Spearman tree; Exo1/Dna2 displayed in HR; global log2 range %.3g–%.3g", yarn_min, yarn_max), common_caption,
  dsb_spearman$order
)
ssb_plot <- make_plot(
  ssb_plot_data,
  "SSB essential components: YARN qsmooth tissue-aware flight RNA expression",
  sprintf("26 mouse tissues; 46 SSB components; shared rows excluded; tissues ordered by Spearman tree; global log2 range %.3g–%.3g", yarn_min, yarn_max), common_caption,
  ssb_spearman$order
)

figure_specs <- list(
  list(stem = "01_DSB_yarn_qsmooth_tissue_aware_log2_heatmap", plot = dsb_plot, width = 22),
  list(stem = "02_SSB_yarn_qsmooth_tissue_aware_log2_heatmap", plot = ssb_plot, width = 31)
)
for (spec in figure_specs) {
  ggsave(file.path(figure_dir, paste0(spec$stem, ".png")), spec$plot, width = spec$width, height = 12.5, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(file.path(figure_dir, paste0(spec$stem, ".pdf")), spec$plot, width = spec$width, height = 12.5, device = cairo_pdf, bg = "white", limitsize = FALSE)
}

file.copy(teacher_xlsx, file.path(input_dir, "Essential_components_1_teacher_source.xlsx"), overwrite = TRUE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"), useBytes = TRUE)

readme_lines <- c(
  "# YARN qsmooth tissue-aware essential-component expression heatmaps",
  "",
  "## 选择这个归一化方案的原因",
  "",
  "- 老师建议使用 `YARN normalizeTissueAware()` 或 TMM/CPM。当前 59 个 analysis unit 的输入是 GeneLab/OSDR differential-expression 表中的 sample-level normalized expression 列，而不是所有数据集都具备、且可合并的统一 raw-count 矩阵；因此本版选择 YARN 的 `qsmooth`，不把 TMM/CPM 误用于已经归一化的值。",
  "- Actb 除法版本保留在 `03_analysis_results/14_essential_components_actb_normalized_expression_heatmaps_20260826`，作为敏感性分析；本版不覆盖它，也不使用 Actb 作为分母。",
  "",
  "## 数据范围和稳定 ID",
  "",
  "- 使用老师 `Essential_components(1).xlsx` 的 DSB（29 行）和 SSB_others（56 行）。SSB_others 中 `shared` 的 10 行按当前指令排除，最终纳入 75 个成分（DSB 29 + SSB 46）。",
  "- 使用 26 个组织、48 个 OSDR accession、59 个 analysis unit 和全部 360 个 flight 样本。图中横轴只显示官方小鼠 Symbol；所有匹配、合并和统计都以 Ensembl Gene ID 为键。DSB 与 SSB 分别根据本图所含成分的组织表达谱，按 `1 - Spearman rho` 距离进行平均连接层次聚类并使用树的叶顺序；因此不按 A–Z 排列。跨组织热图的 NES 顺序仅作为完整组织清单和审计基准。",
  "- 归一化不是只在 75 个目标基因上估计，而是在 59 个输入表都能提供的完整共同 Ensembl gene universe 上估计；本版共同集合大小和每个单位的保留数见 `tables/03_common_ensembl_universe_unit_audit.csv` 与 `tables/04_common_ensembl_gene_universe.csv`。少数目标成分在旧表中没有对应行（共 38 个 sample-component 单元，集中于 OSD-419 和 OSD-686），这些观测保持 NA，不参与该单位/mission 的均值。",
  "",
  "## YARN qsmooth 计算合同",
  "",
  "1. 每个 OSDR 表的 flight 样本值均有限且不小于 1。按 GeneLab 数据约定，把源值 1 视为已包含的 +1 pseudocount，并对共同集合内每个值计算 `source_value - 1`。没有把缺失或低于 1 的值静默替换为 0。",
  "2. 把 59 个 analysis unit 的 360 个 flight 样本拼成 `Ensembl gene × sample` 矩阵，并以 26 个组织作为 YARN 分组。共同集合保证矩阵不含 NA。",
  "3. 调用 `yarn::normalizeTissueAware(eset, groups = tissue, normalizationMethod = \"qsmooth\", log = TRUE, window = 0.05)`。YARN 的 `qsmooth(log=TRUE)` 对第 1 步的值执行 `log2((source_value - 1) + 1)`，再按组织分组进行 quantile-shrinkage；所以没有发生“双加一”。随机种子固定为 25，以固定 qsmooth 对并列秩的处理。",
  "4. 对不在完整共同全集、但在部分表中存在的目标成分，使用同一 fitted qsmooth quantile reference 按其观测 log2 值进行 rank transfer；缺失的源行仍为 NA，绝不补 0、插值或推断。该扩展在样本表的 `YARNValueDerivation` 字段中标记。",
  "5. 归一化后的样本值先在 analysis unit 内取 flight 样本均值，再在同一 mission 内取 analysis-unit 均值，最后在组织内对 mission 等权平均。该层级与此前跨组织结果一致，避免样本较多的 mission 单独支配组织值。",
  "6. 图使用 26 × 75 组织格子的 YARN normalized log2 expression。DSB/SSB 共享蓝—白—红全局色标：最低值为蓝色，中间范围为白色，最高值为红色；每个有限格叠加一位小数。颜色不是 percentile、组织内排名或 log2FC。",
  "7. 为了反映老师提出的成分归类，DSB 图中 Exo1 和 Dna2 从 `HR–A-EJ` 的显示分组移入 `HR`，并插在 Rad54l 与 Rbbp8 之间；`tables/20_figure_component_reassignment_audit.csv` 同时保留原始 Pathway 和显示 Pathway，源映射没有被覆盖。",
  "8. 该版本用于跨组织表达可视化和敏感性比较；qsmooth 不能恢复原始文库大小因子，也不能保证消除所有实验室/批次差异。它不替代原始 DE 表中的 log2FC、P 值或 FDR，也不据此重新进行差异检验。",
  "",
  "## 交付图",
  "",
  "- `figures/01_DSB_yarn_qsmooth_tissue_aware_log2_heatmap.png` / `.pdf`：DSB 29 个成分。",
  "- `figures/02_SSB_yarn_qsmooth_tissue_aware_log2_heatmap.png` / `.pdf`：SSB 46 个成分。",
  "- 图中使用蓝—白—红色标并在每个有限格显示一位小数；完整数值保存在 `tables/09_component_tissue_yarn_qsmooth_log2_matrix.csv`。",
  "",
  "## 审计文件",
  "",
  "- `tables/03_common_ensembl_universe_unit_audit.csv`：每个输入表的行数、有效 Ensembl 行数、源值范围、SHA256 和最终共同集合覆盖。",
  "- `tables/05_flight_sample_metadata.csv`：360 个 flight 样本及其组织、analysis unit 和 mission。",
  "- `tables/06`、`07`、`08`、`09`：component 的样本、analysis-unit、mission 和组织层级结果。",
  "- `tables/11_qsmooth_quantile_weight_audit.csv`：YARN qsmooth 的 quantile 权重审计。",
  "- `tables/12_sample_distribution_before_after_qsmooth.csv`：每个样本归一化前后分布。",
  "- `tables/13_missing_component_presence_audit.csv`：缺失的目标成分源行（共 38 个 sample-component 单元），缺失值未被当作 0。",
  "- `tables/14_yarn_vs_actb_tissue_cell_comparison.csv`、`15_yarn_vs_actb_tissue_correlation_summary.csv`：如 Actb 版本存在，提供 YARN 与 Actb 组织格子的敏感性对照，不改变主结果。",
  "- `tables/17_tissue_spearman_cluster_order.csv`：DSB/SSB 的 26 组织树叶顺序和聚类参数；`18_tissue_spearman_correlation_DSB.csv` 与 `19_tissue_spearman_correlation_SSB.csv`：对应 Spearman 相关矩阵。",
  "- `tables/20_figure_component_reassignment_audit.csv`：Exo1/Dna2 的显示重排与其原始稳定 ID、原始 Pathway 对照。",
  "- `tables/15_integrity_checks.csv`、`16_run_metadata.csv` 和 `sessionInfo.txt`：完整性、参数和软件版本。旧版 `17_essential_components_yarn_qsmooth_expression_heatmaps_20260827` 保留不覆盖。"
)
writeLines(readme_lines, file.path(outdir, "README.md"), useBytes = TRUE)

message("Created YARN qsmooth tissue-aware DSB and SSB heatmaps in: ", outdir)
