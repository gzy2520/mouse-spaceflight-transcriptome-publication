#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(25)
fgsea_nproc <- suppressWarnings(as.integer(Sys.getenv(
  "FGSEA_NPROC", unset = "4"
)))
if (!is.finite(fgsea_nproc) || fgsea_nproc < 1L) fgsea_nproc <- 1L

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
stage01_dir <- file.path(run_root, "01_stable_id_registry")
source_dir <- file.path(run_root, "00_frozen_inputs")
global_dir <- Sys.getenv(
  "STABLE_GLOBAL_OUTPUT_DIR",
  unset = file.path(run_root, "02_global_gsea")
)
go_dir <- Sys.getenv(
  "STABLE_GO_OUTPUT_DIR",
  unset = file.path(run_root, "03_mouse_go_concrete_terms_20260824")
)
global_cache <- file.path(global_dir, "dataset_library_cache")
dir.create(global_cache, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(go_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(go_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

stable_path <- file.path(
  stage01_dir,
  "03_full_gene_statistics_ensembl_long.csv.gz"
)
msigdb_path <- file.path(
  source_dir,
  "MSigDB_2026.1.Mm_mouse_Ensembl_membership.tsv.gz"
)
annotation_root <- Sys.getenv(
  "MOUSE_GO_CURRENT_ANNOTATION_ROOT",
  unset = file.path(
    root,
    "08_GO_annotation/mgi_mod_gaf_five_relation_concrete_terms_20260824"
  )
)
annotation_path <- file.path(
  annotation_root,
  "final_manual_curated/03_mouse_GO_annotations_with_gene_ids.tsv.gz"
)
tissue_mapping_path <- file.path(
  root,
  "03_analysis_results/04_dna_damage_repair_focus/current",
  "00_original_material_tissue_mapping_audit.csv"
)
required <- c(stable_path, msigdb_path, annotation_path, tissue_mapping_path)
if (!all(file.exists(required))) {
  stop("Missing stable-ID GSEA input(s): ", paste(required[!file.exists(required)], collapse = ";"))
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

exact_signflip_p <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  n <- length(values)
  if (!n) return(NA_real_)
  observed <- abs(mean(values))
  patterns <- 0:(2^n - 1)
  permuted <- vapply(patterns, function(pattern) {
    bits <- as.integer(intToBits(pattern))[seq_len(n)]
    signs <- ifelse(bits == 0L, -1, 1)
    abs(mean(signs * values))
  }, numeric(1))
  mean(permuted >= observed - 1e-12)
}

collapse_values <- function(x) {
  values <- sort(unique(as.character(x[
    !is.na(x) & nzchar(as.character(x))
  ])))
  paste(values, collapse = ";")
}

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    c("-a", "256", shQuote(path)),
    stdout = TRUE
  )
  sub(paste0("  ", path, "$"), "", output)
}

run_fgsea_safe <- function(pathways, stats, min_size, max_size) {
  set.seed(25)
  result <- suppressWarnings(fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = as.integer(min_size),
    maxSize = as.integer(max_size),
    eps = 1e-10,
    scoreType = "std",
    nproc = fgsea_nproc
  ))
  result <- as.data.table(result)
  result[, `:=`(
    fallback_used = FALSE,
    fgsea_score_type = "std"
  )]
  missing <- result[!is.finite(NES), .(pathway, ES)]
  run_directional <- function(keys, score_type) {
    if (!length(keys)) return(NULL)
    set.seed(25)
    x <- suppressWarnings(fgseaMultilevel(
      pathways = pathways[keys],
      stats = stats,
      minSize = as.integer(min_size),
      maxSize = as.integer(max_size),
      eps = 1e-10,
      scoreType = score_type,
      nproc = fgsea_nproc
    ))
    x <- as.data.table(x)
    x[, `:=`(
      fallback_used = TRUE,
      fgsea_score_type = score_type
    )]
    x
  }
  if (nrow(missing)) {
    fallback <- rbindlist(list(
      run_directional(missing[ES >= 0, pathway], "pos"),
      run_directional(missing[ES < 0, pathway], "neg")
    ), fill = TRUE)
    result <- rbindlist(list(
      result[!pathway %chin% missing$pathway],
      fallback
    ), fill = TRUE)
  }
  if (result[, any(!is.finite(NES))]) {
    stop(
      "fgsea returned non-finite NES: ",
      paste(result[!is.finite(NES), pathway], collapse = ";")
    )
  }
  result[, leadingEdge := vapply(
    leadingEdge,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
  )]
  result
}

stable <- fread(stable_path)
stable <- stable[
  analysis_role == "primary" &
    analysis_window == "main_30_40d"
]
stopifnot(
  uniqueN(stable$analysis_unit_id) == 59L,
  !anyDuplicated(stable[, .(analysis_unit_id, ensembl_id)]),
  all(stable$symbol_used_as_analysis_key == FALSE)
)

tissue_map <- fread(tissue_mapping_path)
tissue_map <- unique(tissue_map[
  included_in_exact_go_analysis == TRUE,
  .(accession, analysis_tissue, tissue_definition)
])
metadata <- unique(stable[, .(
  accession,
  analysis_unit_id,
  mission_cluster,
  biological_experiment_id,
  animal_cohort_id,
  tissue
)])
metadata <- merge(
  metadata,
  tissue_map,
  by = "accession",
  all.x = TRUE,
  sort = FALSE
)
if (metadata[, any(is.na(analysis_tissue) | !nzchar(analysis_tissue))]) {
  stop(
    "Missing normalized original-material tissue mapping: ",
    paste(
      metadata[is.na(analysis_tissue) | !nzchar(analysis_tissue), accession],
      collapse = ";"
    )
  )
}
stopifnot(!anyDuplicated(metadata$analysis_unit_id), nrow(metadata) == 59L)

annotation <- fread(annotation_path)
evidence_exclude <- trimws(strsplit(
  Sys.getenv("GO_EXCLUDE_EVIDENCE", unset = ""),
  ",",
  fixed = TRUE
)[[1L]])
evidence_exclude <- evidence_exclude[nzchar(evidence_exclude)]
accepted <- annotation[
  included_in_gene_set == TRUE &
    grepl("^ENSMUSG[0-9]+$", ensembl_gene_id)
]
if (length(evidence_exclude)) {
  accepted <- accepted[!EVIDENCE_CODE %chin% evidence_exclude]
}
if (accepted[, any(symbol_used_as_analysis_key == TRUE)]) {
  stop("GO annotation incorrectly marks Symbol as analysis key")
}
term_config <- fread(file.path(annotation_root, "00_extended_GO_term_config.csv"))
go_pathways <- setNames(
  lapply(term_config$term_key, function(key) {
    sort(unique(accepted[
      requested_term_key == key,
      ensembl_gene_id
    ]))
  }),
  term_config$term_key
)
if (any(lengths(go_pathways) == 0L)) {
  stop(
    "Empty GO sets: ",
    paste(names(go_pathways)[lengths(go_pathways) == 0L], collapse = ";")
  )
}
ddr_key <- "DNA_damage_response_GO0006974"
repair_key <- "DNA_repair_GO0006281"
signal_key <- "DNA_damage_signal_transduction_GO0042770"
tolerance_key <- "DNA_damage_tolerance_GO0006301"
apoptosis_key <- "intrinsic_apoptotic_signaling_GO0008630"
telomere_key <- "telomere_maintenance_GO0043247"
telomere_region_key <- "telomere_region_GO0000781"
cell_death_key <- "cell_death_GO0008219"
required_concrete_keys <- c(
  ddr_key, signal_key, repair_key, tolerance_key, apoptosis_key,
  telomere_key, telomere_region_key, cell_death_key
)
if (!all(required_concrete_keys %chin% names(go_pathways))) {
  stop(
    "Required explicit GO root(s) absent: ",
    paste(setdiff(required_concrete_keys, names(go_pathways)), collapse = ";")
  )
}
if (length(setdiff(go_pathways[[repair_key]], go_pathways[[ddr_key]]))) {
  stop("DNA repair is not a subset of DNA damage response")
}

# The source/configuration order is the original 1--15 order.  The display
# order is an explicit, auditable input because it is a teacher-facing choice
# and may be revised without changing any GSEA values.  Prefer the order file
# in the current output; when running a sensitivity output, fall back to the
# main GO output's order file, then the version-controlled project contract.
original_keys <- term_config$term_key
stopifnot(length(original_keys) == 15L)
local_order_file <- file.path(
  go_dir, "tables/00_teacher_requested_term_order.csv"
)
main_order_file <- file.path(
  run_root,
  "03_mouse_go_concrete_terms_20260824/tables/00_teacher_requested_term_order.csv"
)
project_order_file <- file.path(
  root, "07_scripts/config/mouse_full_go_term_display_order.csv"
)
order_file <- Sys.getenv(
  "GO_TERM_ORDER_FILE",
  unset = if (file.exists(local_order_file)) {
    local_order_file
  } else if (file.exists(main_order_file)) {
    main_order_file
  } else {
    project_order_file
  }
)
if (file.exists(order_file)) {
  order_audit_input <- fread(order_file)
  required_order_columns <- c(
    "original_term_order", "display_term_order", "term_key"
  )
  if (!all(required_order_columns %chin% names(order_audit_input))) {
    stop(
      "GO term order audit is missing required columns: ",
      paste(setdiff(required_order_columns, names(order_audit_input)), collapse = ",")
    )
  }
  if (nrow(order_audit_input) != length(original_keys) ||
      anyDuplicated(order_audit_input$term_key) ||
      !setequal(order_audit_input$term_key, original_keys) ||
      anyDuplicated(order_audit_input$display_term_order) ||
      !identical(sort(as.integer(order_audit_input$display_term_order)), seq_len(15L)) ||
      !identical(sort(as.integer(order_audit_input$original_term_order)), seq_len(15L))) {
    stop("GO term order audit must contain a one-to-one permutation of the 15 configured terms.")
  }
  config_original_order <- term_config$official_plot_order[
    match(order_audit_input$term_key, term_config$term_key)
  ]
  if (anyNA(config_original_order) ||
      !identical(as.integer(order_audit_input$original_term_order),
                 as.integer(config_original_order))) {
    stop("GO term order audit original_term_order does not match the configured term order.")
  }
  ordered_keys <- order_audit_input[
    order(display_term_order), term_key
  ]
} else {
  stop("No validated GO term display-order contract found: ", order_file)
}
go_pathways <- go_pathways[ordered_keys]
term_labels <- term_config[
  match(ordered_keys, term_key),
  .(term_key, term_name, go_id)
]
term_labels[, `:=`(
  original_term_order = match(term_key, original_keys),
  term_order = seq_len(.N),
  term_role = fifelse(
    term_key %chin% c(
      ddr_key, repair_key, signal_key, tolerance_key, apoptosis_key,
      telomere_key
    ),
    "official_root",
    fifelse(term_key == cell_death_key, "comparison_context", "context")
  )
)]

gene_set_audit <- term_labels[, {
  genes <- go_pathways[[term_key]]
  list(
    unique_ensembl_gene_ids = length(genes),
    analysis_id_type = "Ensembl Gene ID",
    symbol_used_as_analysis_key = FALSE,
    relation_scope = Sys.getenv(
      "GO_RELATION_LABEL",
      unset = paste(
        "root plus is_a, part_of, regulates,",
        "positively_regulates, negatively_regulates descendants"
      )
    ),
    evidence_code_filter = ifelse(
      length(evidence_exclude),
      paste0("excluded: ", paste(evidence_exclude, collapse = ",")),
      "none"
    ),
    taxon = "NCBITaxon:10090",
    qualifier_filter = "NOT excluded",
    source = paste0(
      "MGI MOUSE-mod.gaf.gz date-generated 2026-08-04; annotation root ",
      basename(annotation_root)
    )
  )
}, by = .(
  term_key, term_name, go_id, original_term_order, term_order, term_role
)]
fwrite(
  gene_set_audit,
  file.path(go_dir, "tables/00_concrete_GO_gene_set_audit.csv"),
  bom = TRUE
)
fwrite(
  term_labels[, .(
    original_term_order, display_term_order = term_order,
    term_key, term_name, go_id, term_role
  )],
  file.path(go_dir, "tables/00_teacher_requested_term_order.csv"),
  bom = TRUE
)

membership_ids <- sort(unique(unlist(go_pathways, use.names = FALSE)))
membership <- data.table(ensembl_gene_id = membership_ids)
for (key in names(go_pathways)) {
  membership[, (paste0("in_", key)) := ensembl_gene_id %chin% go_pathways[[key]]]
}
fwrite(
  membership,
  file.path(go_dir, "tables/01_concrete_GO_membership.csv"),
  bom = TRUE
)

go_cache_dir <- file.path(go_dir, "dataset_cache")
dir.create(go_cache_dir, recursive = TRUE, showWarnings = FALSE)
go_results <- vector("list", nrow(metadata))
for (i in seq_len(nrow(metadata))) {
  m <- metadata[i]
  cache_path <- file.path(
    go_cache_dir,
    paste0(gsub("[^A-Za-z0-9_-]", "_", m$analysis_unit_id), ".csv.gz")
  )
  if (file.exists(cache_path)) {
    result <- fread(cache_path)
    cache_valid <- "pathway" %chin% names(result) &&
      identical(sort(unique(result$pathway)), sort(names(go_pathways)))
  } else {
    cache_valid <- FALSE
  }
  if (!cache_valid) {
    rows <- stable[analysis_unit_id == m$analysis_unit_id]
    rank <- rows[
      is.finite(stat_spaceflight_vs_ground),
      .(rank_score = median(stat_spaceflight_vs_ground)),
      by = ensembl_id
    ]
    setorder(rank, -rank_score, ensembl_id)
    rank[, rank_score := rank_score - seq_len(.N) * 1e-12]
    stats <- rank$rank_score
    names(stats) <- rank$ensembl_id
    result <- run_fgsea_safe(
      go_pathways,
      stats,
      min_size = 5L,
      max_size = 5000L
    )
    result[, `:=`(
      analysis_unit_id = m$analysis_unit_id,
      accession = m$accession,
      mission_cluster = m$mission_cluster,
      biological_experiment_id = m$biological_experiment_id,
      analysis_tissue = m$analysis_tissue,
      tissue_definition = m$tissue_definition,
      n_ranked_genes = length(stats),
      signed_nominal_z = signed_z_from_p(pval, NES),
      seed = 25L
    )]
    fwrite(result, cache_path)
  }
  go_results[[i]] <- result
  message(sprintf(
    "[GO %02d/59] %s",
    i,
    m$analysis_unit_id
  ))
}
go_dataset <- rbindlist(go_results, fill = TRUE)
setnames(go_dataset, "pathway", "term_key")
go_dataset <- merge(
  go_dataset,
  term_labels,
  by = "term_key",
  all.x = TRUE,
  sort = FALSE
)
fwrite(
  go_dataset,
  file.path(go_dir, "tables/02_dataset_GSEA_concrete_terms_and_context.csv.gz")
)

go_mission <- go_dataset[
  ,
  .(
    mission_nes = mean(NES),
    mission_signed_z = median(signed_nominal_z),
    n_biological_experiments = uniqueN(biological_experiment_id),
    n_analysis_units = uniqueN(analysis_unit_id),
    analysis_units = collapse_values(analysis_unit_id)
  ),
  by = .(
    analysis_tissue,
    tissue_definition,
    term_key,
    term_name,
    go_id,
    term_order,
    term_role,
    mission_cluster
  )
]
fwrite(
  go_mission,
  file.path(go_dir, "tables/03_tissue_mission_GSEA_concrete_terms_and_context.csv"),
  bom = TRUE
)

go_tissue <- go_mission[
  ,
  {
    z <- sum(mission_signed_z) / sqrt(.N)
    list(
      n_missions = uniqueN(mission_cluster),
      n_biological_experiments = sum(n_biological_experiments),
      n_analysis_units = sum(n_analysis_units),
      mean_mission_nes = mean(mission_nes),
      median_mission_nes = median(mission_nes),
      n_positive_missions = sum(mission_nes > 0),
      n_negative_missions = sum(mission_nes < 0),
      exact_signflip_p_two_sided = exact_signflip_p(mission_nes),
      stouffer_meta_z = z,
      stouffer_p_two_sided = 2 * pnorm(abs(z), lower.tail = FALSE),
      missions = collapse_values(mission_cluster)
    )
  },
  by = .(
    analysis_tissue,
    tissue_definition,
    term_key,
    term_name,
    go_id,
    term_order,
    term_role
  )
]
go_tissue[, exact_signflip_FDR_all_cells := p.adjust(
  exact_signflip_p_two_sided,
  method = "BH"
)]
go_tissue[, stouffer_FDR_all_cells := p.adjust(
  stouffer_p_two_sided,
  method = "BH"
)]
setorder(go_tissue, term_order, -mean_mission_nes, analysis_tissue)
fwrite(
  go_tissue,
  file.path(go_dir, "tables/04_tissue_statistics_concrete_terms_and_context.csv"),
  bom = TRUE
)

focus_terms <- c(ddr_key, repair_key)
focus <- go_tissue[term_key %chin% focus_terms]
focus_wide <- dcast(
  focus,
  analysis_tissue + tissue_definition + n_missions +
    n_biological_experiments ~ term_key,
  value.var = "mean_mission_nes"
)
ddr_col <- ddr_key
repair_col <- repair_key
focus_wide[, `:=`(
  joint_effect_strength = pmin(abs(get(ddr_col)), abs(get(repair_col))),
  both_terms_direction_concordant =
    sign(get(ddr_col)) == sign(get(repair_col))
)]
loo <- go_mission[term_key %chin% focus_terms, {
  full_sign <- sign(mean(mission_nes))
  missions <- unique(mission_cluster)
  stable_direction <- if (length(missions) < 2L) {
    FALSE
  } else {
    all(vapply(missions, function(omitted) {
      sign(mean(mission_nes[mission_cluster != omitted])) == full_sign
    }, logical(1)))
  }
  list(loo_stable = stable_direction)
}, by = .(analysis_tissue, term_key)]
loo_wide <- dcast(
  loo,
  analysis_tissue ~ term_key,
  value.var = "loo_stable",
  fill = FALSE
)
loo_wide[, leave_one_mission_direction_stable :=
  get(ddr_key) & get(repair_key)]
focus_wide <- merge(
  focus_wide,
  loo_wide[, .(
    analysis_tissue,
    leave_one_mission_direction_stable
  )],
  by = "analysis_tissue",
  all.x = TRUE
)
focus_wide[, eligible_for_focus :=
  n_missions >= 2L &
  n_biological_experiments >= 2L &
  both_terms_direction_concordant &
  leave_one_mission_direction_stable]
setorderv(
  focus_wide,
  c(
    "eligible_for_focus",
    "joint_effect_strength",
    "n_missions",
    "n_biological_experiments",
    "analysis_tissue"
  ),
  c(-1L, -1L, -1L, -1L, 1L)
)
focus_wide[, selection_rank := seq_len(.N)]
focus_wide[, selected_focus := eligible_for_focus & selection_rank == 1L]
fwrite(
  focus_wide,
  file.path(go_dir, "tables/05_focus_tissue_ranking_stable_id.csv"),
  bom = TRUE
)

term_order_for_plot <- term_labels$term_key
tissue_order <- go_tissue[
  term_key == ddr_key,
  analysis_tissue[order(-mean_mission_nes, analysis_tissue)]
]
plot_data <- copy(go_tissue)
plot_data[, analysis_tissue := factor(
  analysis_tissue,
  levels = rev(tissue_order)
)]
plot_data[, term_display := paste0(term_name, "\n(", go_id, ")")]
plot_data[, term_display := factor(
  term_display,
  levels = paste0(
    term_labels$term_name[
    match(term_order_for_plot, term_labels$term_key)
    ],
    "\n(",
    term_labels$go_id[match(term_order_for_plot, term_labels$term_key)],
    ")"
  )
)]
limit <- max(abs(plot_data$mean_mission_nes), na.rm = TRUE)
plot_data[, `:=`(
  mean_NES_label = sprintf("%.2f", mean_mission_nes),
  tile_text_colour = ifelse(abs(mean_mission_nes) >= limit * 0.48,
                            "white", "#1A1A1A")
)]
plot_all <- ggplot(
  plot_data,
  aes(x = term_display, y = analysis_tissue, fill = mean_mission_nes)
) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(
    aes(label = mean_NES_label, colour = tile_text_colour),
    size = 2.25,
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-limit, limit),
    name = "Mission-equal\nmean NES"
  ) +
  labs(
    title = "Mouse tissue GO enrichment using stable Ensembl Gene IDs",
    subtitle = "Explicit GO root terms; no derived DDR-minus-repair gene set",
    x = NULL,
    y = NULL,
    caption = paste(
      "NES is enrichment direction, not pathway activation/inhibition.",
      "Exact mission sign-flip inference is reported in the tables."
    )
  ) +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 42,
      hjust = 1,
      vjust = 1,
      size = 7.2
    ),
    axis.text.y = element_text(size = 7.5),
    plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(go_dir, "figures/01_tissue_15_explicit_GO_terms.png"),
  plot_all,
  width = 19,
  height = 12,
  dpi = 300
)
ggsave(
  file.path(go_dir, "figures/01_tissue_15_explicit_GO_terms.pdf"),
  plot_all,
  width = 19,
  height = 12,
  device = cairo_pdf
)

core_keys <- c(
  ddr_key,
  repair_key,
  signal_key,
  tolerance_key,
  apoptosis_key,
  telomere_key,
  telomere_region_key,
  cell_death_key
)
plot_core <- plot_all %+% plot_data[term_key %chin% core_keys] +
  labs(
    title = "Mouse DDR tissue comparison with explicit GO root terms"
  )
ggsave(
  file.path(go_dir, "figures/02_tissue_DDR_related_explicit_GO_terms.png"),
  plot_core,
  width = 12,
  height = 12,
  dpi = 300
)
ggsave(
  file.path(go_dir, "figures/02_tissue_DDR_related_explicit_GO_terms.pdf"),
  plot_core,
  width = 12,
  height = 12,
  device = cairo_pdf
)

run_global <- tolower(Sys.getenv(
  "RUN_GLOBAL_MSIGDB",
  unset = "true"
)) %chin% c("1", "true", "yes")
global_summary <- NULL
if (run_global) {
  membership_msig <- fread(msigdb_path)
  library_names <- unique(membership_msig$library)
  global_parts <- list()
  part_index <- 0L
  for (library_name in library_names) {
    library_members <- membership_msig[library == library_name]
    pathways <- split(
      library_members$ensembl_gene_id,
      library_members$gene_set_name
    )
    pathways <- lapply(pathways, unique)
    for (i in seq_len(nrow(metadata))) {
      m <- metadata[i]
      cache_path <- file.path(
        global_cache,
        paste0(
          gsub("[^A-Za-z0-9_-]", "_", m$analysis_unit_id),
          "__",
          gsub("[^A-Za-z0-9_-]", "_", library_name),
          ".csv.gz"
        )
      )
      if (file.exists(cache_path)) {
        result <- fread(cache_path)
      } else {
        rows <- stable[analysis_unit_id == m$analysis_unit_id]
        rank <- rows[
          is.finite(stat_spaceflight_vs_ground),
          .(rank_score = median(stat_spaceflight_vs_ground)),
          by = ensembl_id
        ]
        setorder(rank, -rank_score, ensembl_id)
        rank[, rank_score := rank_score - seq_len(.N) * 1e-12]
        stats <- rank$rank_score
        names(stats) <- rank$ensembl_id
        result <- run_fgsea_safe(
          pathways,
          stats,
          min_size = 10L,
          max_size = 5000L
        )
        result[, `:=`(
          library = library_name,
          analysis_unit_id = m$analysis_unit_id,
          accession = m$accession,
          mission_cluster = m$mission_cluster,
          biological_experiment_id = m$biological_experiment_id,
          analysis_tissue = m$analysis_tissue,
          n_ranked_genes = length(stats),
          seed = 25L
        )]
        fwrite(result, cache_path)
      }
      part_index <- part_index + 1L
      global_parts[[part_index]] <- result
      message(sprintf(
        "[GLOBAL %s %02d/59] %s",
        library_name,
        i,
        m$analysis_unit_id
      ))
    }
  }
  global_dataset <- rbindlist(global_parts, fill = TRUE)
  setnames(global_dataset, "pathway", "gene_set_name")
  fwrite(
    global_dataset,
    file.path(global_dir, "01_dataset_GSEA_mouse_MSigDB_Ensembl.csv.gz")
  )
  global_mission <- global_dataset[
    ,
    .(
      mission_mean_NES_across_analysis_units = mean(NES),
      n_analysis_units = uniqueN(analysis_unit_id),
      n_tissues = uniqueN(analysis_tissue)
    ),
    by = .(library, gene_set_name, mission_cluster)
  ]
  fwrite(
    global_mission,
    file.path(global_dir, "02_mission_descriptive_GSEA_summary.csv.gz")
  )
  global_summary <- global_mission[
    ,
    .(
      n_missions = uniqueN(mission_cluster),
      mean_of_mission_mean_NES = mean(
        mission_mean_NES_across_analysis_units
      ),
      median_of_mission_mean_NES = median(
        mission_mean_NES_across_analysis_units
      ),
      positive_missions = sum(
        mission_mean_NES_across_analysis_units > 0
      ),
      negative_missions = sum(
        mission_mean_NES_across_analysis_units < 0
      ),
      inference_status =
        "descriptive_only_heterogeneous_tissues_not_common_effect_meta"
    ),
    by = .(library, gene_set_name)
  ]
  fwrite(
    global_summary,
    file.path(global_dir, "03_global_descriptive_term_summary.csv.gz")
  )
}

provenance <- data.table(
  item = c(
    "stable_id_long_table",
    "concrete_GO_annotations",
    "mouse_msigdb_membership",
    "tissue_mapping"
  ),
  path = c(stable_path, annotation_path, msigdb_path, tissue_mapping_path)
)
provenance[, sha256 := vapply(path, sha256_file, character(1))]
provenance[, note := c(
  "raw OSDR data; Ensembl-only aggregation",
  paste(
    "MGI MOUSE-mod 2026-08-04; GO 2026-07-26;",
    "five ontology relations; explicit GO roots only"
  ),
  "native mouse MSigDB 2026.1.Mm stable Ensembl IDs",
  "normalized original-material tissue labels; no gene values"
)]
fwrite(
  provenance,
  file.path(go_dir, "tables/06_input_provenance_sha256.csv"),
  bom = TRUE
)

selected <- focus_wide[selected_focus == TRUE, analysis_tissue]
readme <- c(
  "# Stable-ID mouse GO rerun with explicit GO root terms",
  "",
  "- Analysis key: Ensembl Gene ID; Symbol is display-only.",
  "- Mouse annotation: MGI MOUSE-mod GAF, date-generated 2026-08-04.",
  "- Ontology: go-basic.obo releases/2026-07-26.",
  paste0(
    "- Relations: is_a, part_of, regulates, positively_regulates, ",
    "negatively_regulates."
  ),
  "- Taxon: mouse 10090; NOT excluded; no evidence-code filter.",
  "- Rank: complete signed OSDR statistic; no gene-level significance cutoff.",
  "- Mission aggregation: same-mission units first, then missions equally weighted.",
  "- Primary cross-mission inference: exact two-sided mission sign-flip.",
  "- Stouffer results are secondary continuous evidence.",
  paste0(
    "- Teacher-requested display order is read from ",
    "tables/00_teacher_requested_term_order.csv; the original 1--15 order is retained."
  ),
  "- DNA damage signaling: GO:0042770 signal transduction in response to DNA damage.",
  "- Additional explicit roots: GO:0006301 DNA damage tolerance; GO:0008630 intrinsic apoptotic signaling; GO:0043247 telomere maintenance.",
  "- GO:0000781 chromosome, telomeric region is retained separately as a Cellular Component context set.",
  "- No derived DDR-minus-repair complement is used or interpreted as a GO pathway.",
  paste0("- Selected exploratory focus tissue: ", selected, "."),
  "- NES must not be described as pathway activation or inhibition.",
  "",
  "The global MSigDB summaries are descriptive across heterogeneous tissues;",
  "they are not common-effect meta-analysis estimates."
)
writeLines(readme, file.path(go_dir, "README.md"), useBytes = TRUE)

message("Stable-ID GO and global GSEA completed: ", run_root)
