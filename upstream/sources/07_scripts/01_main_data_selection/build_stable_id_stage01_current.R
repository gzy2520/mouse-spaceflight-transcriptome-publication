#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)
set.seed(25)

root <- normalizePath(
  Sys.getenv("GRADUATE_DESIGN_ROOT", unset = "."),
  mustWork = TRUE
)
run_root <- Sys.getenv(
  "STABLE_ID_RERUN_ROOT",
  unset = file.path(
    root,
    "03_analysis_results/10_full_stable_id_rerun_20260823"
  )
)
out_dir <- file.path(run_root, "01_stable_id_registry")
source_dir <- file.path(run_root, "00_frozen_inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(root, "07_scripts/config/dataset_manifest.csv")
manifest <- fread(manifest_path, fill = TRUE)
manifest[, contrast_id := fifelse(
  is.na(contrast_id),
  "",
  trimws(as.character(contrast_id))
)]
manifest[, analysis_unit_id := fifelse(
  nzchar(contrast_id),
  paste(accession, contrast_id, sep = "__"),
  accession
)]
stopifnot(
  nrow(manifest[analysis_role == "primary" & analysis_window == "main_30_40d"]) == 59L,
  uniqueN(manifest[
    analysis_role == "primary" & analysis_window == "main_30_40d",
    analysis_unit_id
  ]) == 59L
)

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    c("-a", "256", shQuote(path)),
    stdout = TRUE
  )
  sub(paste0("  ", path, "$"), "", output)
}

signed_z_from_p <- function(p_value, effect) {
  p_value <- as.numeric(p_value)
  effect <- as.numeric(effect)
  result <- rep(NA_real_, length(p_value))
  keep <- is.finite(p_value) & p_value >= 0 & p_value <= 1 &
    is.finite(effect)
  bounded <- pmin(pmax(p_value[keep], .Machine$double.xmin), 1)
  result[keep] <- qnorm(bounded / 2, lower.tail = FALSE) *
    sign(effect[keep])
  result
}

collapse_values <- function(x) {
  values <- sort(unique(trimws(as.character(
    x[!is.na(x) & nzchar(trimws(as.character(x)))]
  ))))
  paste(values, collapse = ";")
}

companion_column <- function(log2fc_column, prefix) {
  sub("^Log2fc_", paste0(prefix, "_"), log2fc_column)
}

read_unit <- function(index) {
  m <- manifest[index]
  input_path <- file.path(root, m$input_relative_path)
  if (!file.exists(input_path)) {
    stop("Missing raw DE input: ", input_path)
  }
  header <- names(fread(input_path, nrows = 0L, showProgress = FALSE))
  log_col <- m$raw_log2fc_column
  stat_col <- companion_column(log_col, "Stat")
  p_col <- companion_column(log_col, "P.value")
  adj_col <- companion_column(log_col, "Adj.p.value")
  required <- c("ENSEMBL", log_col, p_col, adj_col)
  if (!all(required %chin% header)) {
    stop(
      m$analysis_unit_id, " lacks: ",
      paste(setdiff(required, header), collapse = ";")
    )
  }
  select_columns <- unique(c(
    "ENSEMBL",
    if ("SYMBOL" %chin% header) "SYMBOL" else character(),
    log_col,
    if (stat_col %chin% header) stat_col else character(),
    p_col,
    adj_col
  ))
  raw <- fread(input_path, select = select_columns, showProgress = FALSE)
  if (!"SYMBOL" %chin% names(raw)) raw[, SYMBOL := ""]
  multiplier <- as.integer(m$orientation_multiplier)
  if (!multiplier %in% c(-1L, 1L)) {
    stop("Invalid orientation multiplier: ", m$analysis_unit_id)
  }

  prepared <- data.table(
    ensembl_id = sub(
      "\\.[0-9]+$",
      "",
      toupper(trimws(as.character(raw$ENSEMBL)))
    ),
    gene_symbol = trimws(as.character(raw$SYMBOL)),
    log2fc_spaceflight_vs_ground =
      as.numeric(raw[[log_col]]) * multiplier,
    raw_p_value = as.numeric(raw[[p_col]]),
    adjusted_p_value = as.numeric(raw[[adj_col]])
  )
  statistic_source <- if (stat_col %chin% names(raw)) {
    prepared[, stat_spaceflight_vs_ground :=
      as.numeric(raw[[stat_col]]) * multiplier]
    "OSDR_Stat"
  } else {
    prepared[, stat_spaceflight_vs_ground := signed_z_from_p(
      raw_p_value,
      log2fc_spaceflight_vs_ground
    )]
    "signed_Z_from_raw_P"
  }

  valid_id <- grepl("^ENSMUSG[0-9]+$", prepared$ensembl_id)
  finite_stat <- is.finite(prepared$stat_spaceflight_vs_ground)
  finite_effect <- is.finite(prepared$log2fc_spaceflight_vs_ground)
  kept <- prepared[valid_id & finite_stat & finite_effect]
  duplicate_rows <- kept[, .N, by = ensembl_id][N > 1L]

  stable <- kept[
    ,
    .(
      gene_symbol = collapse_values(gene_symbol),
      log2fc_spaceflight_vs_ground =
        median(log2fc_spaceflight_vs_ground),
      stat_spaceflight_vs_ground =
        median(stat_spaceflight_vs_ground),
      raw_p_value = median(raw_p_value, na.rm = TRUE),
      adjusted_p_value = median(adjusted_p_value, na.rm = TRUE),
      n_source_rows = .N
    ),
    by = ensembl_id
  ]
  stable[is.nan(raw_p_value), raw_p_value := NA_real_]
  stable[is.nan(adjusted_p_value), adjusted_p_value := NA_real_]
  stable[, `:=`(
    accession = m$accession,
    analysis_unit_id = m$analysis_unit_id,
    contrast_id = m$contrast_id,
    analysis_window = m$analysis_window,
    analysis_role = m$analysis_role,
    tissue = m$tissue,
    mission_raw = m$mission_raw,
    mission_cluster = m$mission_cluster,
    duration_days = m$duration_days,
    biological_experiment_id = m$biological_experiment_id,
    animal_cohort_id = m$animal_cohort_id,
    raw_log2fc_column = m$raw_log2fc_column,
    raw_left_label = m$raw_left_label,
    raw_right_label = m$raw_right_label,
    orientation_multiplier = multiplier,
    unified_contrast_label = m$unified_contrast_label,
    curation_note = m$curation_note,
    analysis_id_type = "Ensembl Gene ID",
    symbol_used_as_analysis_key = FALSE
  )]

  audit <- data.table(
    accession = m$accession,
    analysis_unit_id = m$analysis_unit_id,
    analysis_window = m$analysis_window,
    analysis_role = m$analysis_role,
    mission_cluster = m$mission_cluster,
    tissue = m$tissue,
    input_relative_path = m$input_relative_path,
    input_sha256 = sha256_file(input_path),
    orientation_multiplier = multiplier,
    unified_contrast_label = m$unified_contrast_label,
    statistic_source = statistic_source,
    raw_rows = nrow(raw),
    valid_ensembl_effect_stat_rows = nrow(kept),
    unique_stable_ensembl_rows = nrow(stable),
    invalid_ensembl_rows = sum(!valid_id),
    nonfinite_stat_rows = sum(valid_id & !finite_stat),
    nonfinite_effect_rows = sum(valid_id & !finite_effect),
    duplicate_stable_id_rows_collapsed = sum(duplicate_rows$N - 1L),
    semicolon_ensembl_cells = sum(grepl(";", stable$ensembl_id, fixed = TRUE)),
    analysis_id_type = "Ensembl Gene ID",
    symbol_used_as_analysis_key = FALSE
  )
  list(data = stable, audit = audit)
}

loaded <- lapply(seq_len(nrow(manifest)), read_unit)
stable_long <- rbindlist(lapply(loaded, `[[`, "data"), use.names = TRUE)
unit_audit <- rbindlist(lapply(loaded, `[[`, "audit"), use.names = TRUE)

stopifnot(
  !anyDuplicated(stable_long[, .(analysis_unit_id, ensembl_id)]),
  all(grepl("^ENSMUSG[0-9]+$", stable_long$ensembl_id)),
  all(stable_long$symbol_used_as_analysis_key == FALSE),
  unit_audit[
    analysis_role == "primary" & analysis_window == "main_30_40d",
    .N
  ] == 59L,
  unit_audit[, sum(semicolon_ensembl_cells)] == 0L
)

fwrite(
  manifest,
  file.path(out_dir, "01_dataset_manifest_validated.csv"),
  bom = TRUE
)
fwrite(
  unit_audit,
  file.path(out_dir, "02_stable_id_input_and_direction_audit.csv"),
  bom = TRUE
)
fwrite(
  stable_long,
  file.path(out_dir, "03_full_gene_statistics_ensembl_long.csv.gz")
)

display_audit <- stable_long[
  ,
  .(
    display_symbols = collapse_values(gene_symbol),
    n_display_symbols = uniqueN(gene_symbol[nzchar(gene_symbol)]),
    n_analysis_units = uniqueN(analysis_unit_id)
  ),
  by = ensembl_id
]
fwrite(
  display_audit,
  file.path(out_dir, "04_ensembl_to_display_symbol_audit.csv"),
  bom = TRUE
)

collections <- list(
  list(
    library = "MSigDB_2026.1.Mm_Hallmark",
    collection = "MH",
    subcollection = NA_character_
  ),
  list(
    library = "MSigDB_2026.1.Mm_GO_BP",
    collection = "M5",
    subcollection = "GO:BP"
  ),
  list(
    library = "MSigDB_2026.1.Mm_Reactome",
    collection = "M2",
    subcollection = "CP:REACTOME"
  ),
  list(
    library = "MSigDB_2026.1.Mm_WikiPathways",
    collection = "M2",
    subcollection = "CP:WIKIPATHWAYS"
  )
)

frozen_msigdb <- Sys.getenv("MSIGDB_MEMBERSHIP_FILE", unset = "")
if (nzchar(frozen_msigdb)) {
  msigdb <- fread(frozen_msigdb)
} else {
msig_rows <- lapply(collections, function(spec) {
  args <- list(
    db_species = "MM",
    species = "Mus musculus",
    collection = spec$collection
  )
  if (!is.na(spec$subcollection)) {
    args$subcollection <- spec$subcollection
  }
  x <- as.data.table(do.call(msigdbr, args))
  if (!all(x$db_target_species == "MM")) {
    stop("A non-mouse MSigDB target entered ", spec$library)
  }
  if (!all(x$db_ensembl_gene == x$ensembl_gene)) {
    stop("Mouse MSigDB native and requested Ensembl IDs differ in ", spec$library)
  }
  x <- x[
    grepl("^ENSMUSG[0-9]+$", db_ensembl_gene),
    .(
      library = spec$library,
      gene_set_id = gs_id,
      gene_set_name = gs_name,
      ensembl_gene_id = db_ensembl_gene,
      ncbi_gene_id = db_ncbi_gene,
      source_gene = source_gene,
      collection = gs_collection,
      subcollection = gs_subcollection,
      source_species = gs_source_species,
      db_version,
      db_target_species
    )
  ]
  x <- x[
    ,
    .(
      source_gene = collapse_values(source_gene),
      source_membership_rows = .N
    ),
    by = .(
      library,
      gene_set_id,
      gene_set_name,
      ensembl_gene_id,
      ncbi_gene_id,
      collection,
      subcollection,
      source_species,
      db_version,
      db_target_species
    )
  ]
  x
})
msigdb <- rbindlist(msig_rows, use.names = TRUE)
}
stopifnot(
  all(grepl("^ENSMUSG[0-9]+$", msigdb$ensembl_gene_id)),
  identical(unique(msigdb$db_version), "2026.1.Mm"),
  !anyDuplicated(msigdb[, .(library, gene_set_name, ensembl_gene_id)])
)
msigdb_path <- file.path(
  source_dir,
  "MSigDB_2026.1.Mm_mouse_Ensembl_membership.tsv.gz"
)
fwrite(msigdb, msigdb_path, sep = "\t")

msig_sizes <- msigdb[
  ,
  .(gene_set_size = uniqueN(ensembl_gene_id)),
  by = .(library, gene_set_name)
]
msig_audit <- msigdb[
  ,
  .(
    gene_sets = uniqueN(gene_set_name),
    unique_ensembl_gene_ids = uniqueN(ensembl_gene_id),
    gene_set_gene_pairs = .N,
    db_version = unique(db_version),
    target_species = unique(db_target_species),
    analysis_id_type = "Ensembl Gene ID",
    symbol_used_as_analysis_key = FALSE
  ),
  by = library
]
msig_audit <- merge(
  msig_audit,
  msig_sizes[
    ,
    .(
      minimum_gene_set_size = min(gene_set_size),
      median_gene_set_size = median(gene_set_size),
      maximum_gene_set_size = max(gene_set_size)
    ),
    by = library
  ],
  by = "library",
  sort = FALSE
)
fwrite(
  msig_audit,
  file.path(source_dir, "MSigDB_2026.1.Mm_collection_audit.csv"),
  bom = TRUE
)

provenance <- data.table(
  item = c(
    "dataset_manifest",
    "stable_id_long_table",
    "mouse_msigdb_membership"
  ),
  path = c(
    manifest_path,
    file.path(out_dir, "03_full_gene_statistics_ensembl_long.csv.gz"),
    msigdb_path
  )
)
provenance[, sha256 := vapply(path, sha256_file, character(1))]
provenance[, note := c(
  "59 primary units plus declared sensitivity rows; direction contract",
  "raw OSDR rows aggregated only within identical stable Ensembl Gene ID",
  paste0(
    "msigdbr ", as.character(packageVersion("msigdbr")),
    "; native mouse MSigDB 2026.1.Mm; no Symbol matching"
  )
)]
fwrite(
  provenance,
  file.path(out_dir, "05_input_output_provenance_sha256.csv"),
  bom = TRUE
)

summary <- data.table(
  metric = c(
    "primary_analysis_units",
    "stable_id_rows_primary",
    "median_unique_ensembl_per_primary_unit",
    "total_semicolon_ensembl_cells",
    "primary_duplicate_stable_id_rows_collapsed",
    "mouse_msigdb_version",
    "random_seed"
  ),
  value = c(
    "59",
    stable_long[
      analysis_role == "primary" & analysis_window == "main_30_40d",
      as.character(.N)
    ],
    unit_audit[
      analysis_role == "primary" & analysis_window == "main_30_40d",
      as.character(median(unique_stable_ensembl_rows))
    ],
    as.character(unit_audit[, sum(semicolon_ensembl_cells)]),
    unit_audit[
      analysis_role == "primary" & analysis_window == "main_30_40d",
      as.character(sum(duplicate_stable_id_rows_collapsed))
    ],
    unique(msigdb$db_version),
    "25"
  )
)
fwrite(
  summary,
  file.path(out_dir, "06_stage01_stable_id_validation_summary.csv"),
  bom = TRUE
)

message("Stable-ID Stage01 and mouse MSigDB snapshot completed: ", run_root)
