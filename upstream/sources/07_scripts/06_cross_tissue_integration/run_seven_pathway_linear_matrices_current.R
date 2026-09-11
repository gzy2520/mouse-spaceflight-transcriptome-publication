#!/usr/bin/env Rscript

# Build one ontology-complete seven-pathway linear gene matrix per strict
# tissue in the primary analysis registry, using Ensembl Gene ID as the only
# analysis/merge key.
#
# RNA evidence:
#   - unified Space Flight vs Ground gene statistics
#   - mission-equal tissue summary
#   - median log2FC for the displayed tissue effect
#   - signed Stouffer evidence, followed by within-tissue BH adjustment
#
# GO evidence:
#   - frozen MGI MOD-GAF mouse annotations
#   - a mouse-supported seven-pathway GO catalog: BP terms are constrained to
#     DNA repair (GO:0006281), while pathway-specific MF/CC terms remain
#     eligible when directly annotated in mouse
#   - MGI -> Ensembl release 116 official ID mapping
#
# Gene symbols are retained only as optional display labels. They are never
# used to join RNA and GO tables or to choose an Ensembl candidate.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, datatable.allow.cartesian = TRUE)
set.seed(25)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

parse_args <- function(arguments) {
  defaults <- list(
    root = ".",
    `ignore-pvalue` = "false",
    `skill-root` = Sys.getenv(
      "PLOT_GO_PATHWAY_LINEAR_MATRIX_SKILL",
      unset = file.path("vendor", "plot-go-pathway-linear-matrix")
    ),
    outdir = paste0(
      "03_analysis_results/10_full_stable_id_rerun_20260823/",
      "04_seven_pathway_response_26_tissues_20260824"
    ),
    `tissue-mapping` = paste0(
      "03_analysis_results/04_dna_damage_repair_focus/current/",
      "00_original_material_tissue_mapping_audit.csv"
    ),
    rules = paste0(
      "07_scripts/config/",
      "mouse_seven_pathway_go_seed_rules_five_relation_20260822.csv"
    ),
    `dna-repair-relations` = paste(
      c(
        "is_a", "part_of", "regulates",
        "positively_regulates", "negatively_regulates"
      ),
      collapse = ","
    )
  )
  index <- 1L
  while (index <= length(arguments)) {
    key <- arguments[[index]]
    assert(startsWith(key, "--"), paste("Unexpected argument:", key))
    assert(index < length(arguments), paste("Missing value for", key))
    defaults[[sub("^--", "", key)]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  defaults
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
ignore_pvalue <- tolower(as.character(args[["ignore-pvalue"]])) %in%
  c("1", "true", "yes", "y")
root <- normalizePath(args$root, mustWork = TRUE)
skill_root <- normalizePath(args[["skill-root"]], mustWork = TRUE)
resolve_input_path <- function(path) {
  if (grepl("^/", path)) path else file.path(root, path)
}
outdir <- if (grepl("^/", args$outdir)) {
  args$outdir
} else {
  file.path(root, args$outdir)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, mustWork = TRUE)
tissue_mapping_path <- resolve_input_path(args[["tissue-mapping"]])
project_rules_path <- resolve_input_path(args$rules)
dna_repair_relations <- trimws(strsplit(
  as.character(args[["dna-repair-relations"]]),
  ",",
  fixed = TRUE
)[[1L]])
dna_repair_relations <- unique(dna_repair_relations[nzchar(dna_repair_relations)])
allowed_dna_repair_relations <- c(
  "is_a", "part_of", "regulates",
  "positively_regulates", "negatively_regulates"
)
assert(
  length(dna_repair_relations) > 0L &&
    all(dna_repair_relations %chin% allowed_dna_repair_relations),
  "--dna-repair-relations must be a nonempty subset of the declared GO relations."
)

input_dir <- file.path(outdir, "prepared_inputs")
proliferation_expression_audit_path <- file.path(
  root,
  "03_analysis_results/06_cross_tissue_integration/05_proliferation_state/tables",
  "02_analysis_unit_expression_and_group_audit.csv"
)
skill_prefix <- if (ignore_pvalue) {
  "seven_pathway_no_pvalue"
} else {
  "seven_pathway"
}
skill_outdir <- file.path(
  outdir,
  if (ignore_pvalue) "skill_output_no_pvalue" else "skill_output"
)
project_catalog_path <- file.path(
  input_dir,
  "07_project_curated_seven_pathway_GO_catalog.csv"
)
project_catalog_provenance_path <- file.path(
  input_dir,
  "07_project_curated_seven_pathway_GO_catalog_provenance.csv"
)
project_catalog_log_path <- file.path(
  input_dir,
  "07_project_curated_seven_pathway_GO_catalog_build_log.txt"
)
run_log_path <- file.path(
  outdir,
  if (ignore_pvalue) {
    "run_command_and_log_no_pvalue.txt"
  } else {
    "run_command_and_log.txt"
  }
)
run_summary_path <- file.path(
  outdir,
  if (ignore_pvalue) "run_summary_no_pvalue.csv" else "run_summary.csv"
)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(skill_outdir, recursive = TRUE, showWarnings = FALSE)

dataset_manifest_path <- file.path(
  root,
  "07_scripts/config/dataset_manifest.csv"
)
tissue_registry_path <- file.path(
  root,
  "01_osdr_database/主分析数据集_OSDR数据库抽取_20260617.csv"
)
gaf_path <- file.path(
  root,
  "08_GO_annotation/mgi_mod_gaf_five_relation_20260822/source/MOUSE-mod.gaf.gz"
)
obo_path <- file.path(
  root,
  "08_GO_annotation/mgi_mod_gaf_five_relation_20260822/source/go-basic.obo"
)
stable_mapping_audit_path <- file.path(
  root,
  "08_GO_annotation/mgi_mod_gaf_five_relation_20260822",
  "01_MGI_to_Ensembl116_multi_route_mapping_audit.csv"
)
annotation_source_manifest_path <- file.path(
  root,
  "08_GO_annotation/mgi_mod_gaf_five_relation_20260822",
  "00_source_manifest.csv"
)
catalog_builder_script_path <- file.path(
  skill_root,
  "scripts/build_seven_pathway_go_catalog.R"
)
display_path <- file.path(
  skill_root,
  "assets/seven_pathway_display.csv"
)
skill_script_path <- file.path(
  skill_root,
  "scripts/plot_go_pathway_linear_matrix.R"
)

required_files <- c(
  dataset_manifest_path,
  tissue_registry_path,
  tissue_mapping_path,
  proliferation_expression_audit_path,
  gaf_path,
  obo_path,
  stable_mapping_audit_path,
  annotation_source_manifest_path,
  project_rules_path,
  catalog_builder_script_path,
  display_path,
  skill_script_path
)
assert(
  all(file.exists(required_files)),
  paste(
    "Missing required file(s):",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
)

normalize_ensembl <- function(values) {
  values <- trimws(as.character(values))
  sub("\\.[0-9]+$", "", values)
}

collapse_values <- function(values) {
  values <- sort(unique(trimws(as.character(values))))
  values <- values[!is.na(values) & nzchar(values)]
  paste(values, collapse = ";")
}

first_display_label <- function(values) {
  values <- sort(unique(trimws(as.character(values))))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values)) values[[1L]] else ""
}

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    args = c("-a", "256", shQuote(path)),
    stdout = TRUE
  )
  sub(paste0("  ", path, "$"), "", output[[1L]])
}

bh_adjust <- function(values) {
  result <- rep(NA_real_, length(values))
  keep <- is.finite(values)
  result[keep] <- p.adjust(values[keep], method = "BH")
  result
}

companion_column <- function(log2fc_column, prefix) {
  assert(
    startsWith(log2fc_column, "Log2fc_"),
    paste("Invalid OSDR log2FC column:", log2fc_column)
  )
  sub("^Log2fc_", paste0(prefix, "_"), log2fc_column)
}

signed_z_from_p <- function(p_value, effect) {
  p_value <- as.numeric(p_value)
  effect <- as.numeric(effect)
  result <- rep(NA_real_, length(p_value))
  keep <- is.finite(p_value) &
    p_value >= 0 &
    p_value <= 1 &
    is.finite(effect)
  bounded_p <- pmin(
    pmax(p_value[keep], .Machine$double.xmin),
    1
  )
  result[keep] <- qnorm(bounded_p / 2, lower.tail = FALSE) *
    sign(effect[keep])
  result
}

read_mgi_gaf <- function(path) {
  gaf_columns <- c(
    "DB",
    "DB_Object_ID",
    "SYMBOL",
    "QUALIFIER",
    "GO_TERM",
    "DB_REFERENCE",
    "EVIDENCE_CODE",
    "WITH_FROM",
    "GO_ASPECT",
    "DB_OBJECT_NAME",
    "DB_OBJECT_SYNONYM",
    "DB_OBJECT_TYPE",
    "TAXON",
    "DATE",
    "ASSIGNED_BY",
    "ANNOTATION_EXTENSION",
    "GENE_PRODUCT_FORM_ID"
  )
  command <- paste(
    "gzip -dc",
    shQuote(path),
    "| awk 'BEGIN{FS=\"\\t\"} $0 !~ /^!/'"
  )
  result <- fread(
    cmd = command,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    col.names = gaf_columns,
    showProgress = FALSE
  )
  assert(ncol(result) == length(gaf_columns), "Unexpected MGI GAF schema.")
  result[, source_row_number := seq_len(.N)]
  result
}

parse_obo_relation_edges <- function(path, allowed_relations) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  starts <- which(lines == "[Term]")
  assert(length(starts) > 0L, "Frozen GO OBO contains no [Term] blocks.")
  ends <- c(starts[-1L] - 1L, length(lines))
  edge_rows <- vector("list", length(starts))
  for (index in seq_along(starts)) {
    block <- lines[starts[[index]]:ends[[index]]]
    id_line <- grep("^id: GO:", block, value = TRUE)
    if (!length(id_line) || any(block == "is_obsolete: true")) next
    child <- sub("^id: ", "", id_line[[1L]])
    term_edges <- list()
    is_a_parents <- sub(
      "^is_a: (GO:[0-9]{7}).*$", "\\1",
      grep("^is_a: GO:", block, value = TRUE)
    )
    is_a_parents <- is_a_parents[grepl("^GO:[0-9]{7}$", is_a_parents)]
    if (length(is_a_parents)) {
      term_edges[[length(term_edges) + 1L]] <- data.table(
        Child = child,
        Parent = is_a_parents,
        Relation = "is_a"
      )
    }
    relationship_lines <- grep("^relationship: [^ ]+ GO:", block, value = TRUE)
    if (length(relationship_lines)) {
      relationship_edges <- data.table(
        Child = child,
        Relation = sub(
          "^relationship: ([^ ]+) GO:.*$", "\\1", relationship_lines
        ),
        Parent = sub(
          "^relationship: [^ ]+ (GO:[0-9]{7}).*$", "\\1", relationship_lines
        )
      )
      relationship_edges <- relationship_edges[
        Relation %chin% allowed_relations & grepl("^GO:[0-9]{7}$", Parent)
      ]
      if (nrow(relationship_edges)) {
        term_edges[[length(term_edges) + 1L]] <- relationship_edges
      }
    }
    term_edges <- rbindlist(term_edges, use.names = TRUE, fill = TRUE)
    if (nrow(term_edges)) {
      edge_rows[[index]] <- term_edges[Relation %chin% allowed_relations]
    }
  }
  unique(rbindlist(edge_rows, use.names = TRUE, fill = TRUE))
}

ontology_descendants <- function(edges, root_go_id) {
  seen <- root_go_id
  frontier <- root_go_id
  repeat {
    children <- unique(edges[Parent %chin% frontier, Child])
    children <- setdiff(children, seen)
    if (!length(children)) break
    seen <- c(seen, children)
    frontier <- children
  }
  seen
}

read_flight_expression_unit <- function(manifest_row, expression_audit) {
  audit_row <- expression_audit[
    analysis_unit_id == manifest_row$analysis_unit_id
  ]
  assert(
    nrow(audit_row) == 1L,
    paste("Missing or duplicated proliferation expression audit for", manifest_row$analysis_unit_id)
  )
  input_path <- file.path(root, manifest_row$input_relative_path)
  assert(
    identical(normalizePath(audit_row$input_relative_path), normalizePath(input_path)),
    paste("Expression audit input does not match manifest for", manifest_row$analysis_unit_id)
  )
  flight_samples <- strsplit(audit_row$flight_samples, ";", fixed = TRUE)[[1L]]
  header <- names(fread(input_path, nrows = 0L, check.names = FALSE))
  assert(
    all(flight_samples %chin% header),
    paste("Flight/uG sample columns are absent for", manifest_row$analysis_unit_id)
  )
  selected_columns <- unique(c(
    "ENSEMBL",
    if ("SYMBOL" %in% header) "SYMBOL" else character(),
    flight_samples
  ))
  raw <- fread(
    input_path,
    select = selected_columns,
    check.names = FALSE,
    showProgress = FALSE
  )
  if (!"SYMBOL" %in% names(raw)) raw[, SYMBOL := ""]
  counts <- as.matrix(raw[, ..flight_samples])
  storage.mode(counts) <- "double"
  valid_expression <- rowSums(!is.finite(counts) | counts <= 0) == 0L
  result <- data.table(
    EnsemblID = normalize_ensembl(raw$ENSEMBL),
    gene_symbol = as.character(raw$SYMBOL),
    FlightExpressionLog2NormalizedCount = rowMeans(log2(counts))
  )
  valid_id <- grepl("^ENSMUSG[0-9]+$", result$EnsemblID)
  result <- result[valid_id & valid_expression]
  assert(
    nrow(result) > 0L && !anyDuplicated(result$EnsemblID),
    paste("Flight/uG expression has no unique stable IDs for", manifest_row$analysis_unit_id)
  )
  result[, `:=`(
    Group = manifest_row$Group,
    accession = manifest_row$accession,
    analysis_unit_id = manifest_row$analysis_unit_id,
    mission_cluster = manifest_row$mission_cluster,
    n_flight_samples = length(flight_samples)
  )]
  list(
    data = result,
    audit = data.table(
      Group = manifest_row$Group,
      accession = manifest_row$accession,
      analysis_unit_id = manifest_row$analysis_unit_id,
      mission_cluster = manifest_row$mission_cluster,
      n_flight_samples = length(flight_samples),
      expression_scale = "log2(GeneLab sample-column normalized count with +1 pseudocount)",
      valid_unique_Ensembl_expression_rows = nrow(result),
      invalid_or_nonpositive_expression_rows = sum(!valid_id | !valid_expression)
    )
  )
}

map_gradient_colours <- function(values, low, high, palette, missing = "#F0F2F4") {
  result <- rep(missing, length(values))
  keep <- is.finite(values)
  if (!any(keep)) return(result)
  if (!is.finite(low) || !is.finite(high) || high <= low) {
    result[keep] <- palette[[length(palette)]]
    return(result)
  }
  scaled <- pmin(1, pmax(0, (values[keep] - low) / (high - low)))
  indices <- pmax(1L, pmin(length(palette), floor(scaled * (length(palette) - 1L)) + 1L))
  result[keep] <- palette[indices]
  result
}

safe_figure_filename <- function(value) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", as.character(value))
  value <- gsub("^_+|_+$", "", value)
  ifelse(nzchar(value), value, "group")
}

render_expression_matrix_outputs <- function(
  skill_outdir,
  skill_prefix,
  tissue_rna,
  expected_pathways
) {
  matrix_path <- file.path(skill_outdir, "gene_pathway_matrix.csv")
  gene_order_path <- file.path(skill_outdir, "gene_order.csv")
  summary_path <- file.path(skill_outdir, "pathway_summary.csv")
  matrix_data <- fread(matrix_path)
  gene_order <- fread(gene_order_path)
  expression_columns <- tissue_rna[, .(
    Group,
    EnsemblID,
    FlightExpressionLog2NormalizedCount,
    tissue_meta_p,
    log2FoldChange,
    padj
  )]
  # Make presentation-only re-rendering idempotent: a prior rendering adds
  # these display columns to gene_pathway_matrix.csv, but they must never be
  # used as a second source during a subsequent merge.
  previous_render_columns <- intersect(
    names(matrix_data),
    c(
      setdiff(names(expression_columns), c("Group", "EnsemblID")),
      "ExpressionScaleLower", "ExpressionScaleUpper", "ExpressionColor",
      "TopPValue", "TopLog2FC", "TopNegLog10P", "TopLog2FCVisible",
      "TopPValueVisible", "TopColourCap", "TopPColourCap", "TopColor",
      "TopPColor"
    )
  )
  if (length(previous_render_columns)) {
    matrix_data[, (previous_render_columns) := NULL]
  }
  matrix_data <- merge(
    matrix_data,
    expression_columns,
    by = c("Group", "EnsemblID"),
    all.x = TRUE,
    sort = FALSE
  )
  assert(
    all(matrix_data[PathwayPresent == TRUE, is.finite(FlightExpressionLog2NormalizedCount)]),
    "A displayed GO member lacks flight/uG normalized-count expression."
  )
  expression_palette <- colorRampPalette(c("#FFFFFF", "#FEE0D2", "#CB181D"))(256L)
  top_palette <- colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(257L)
  group_scale <- matrix_data[
    PathwayPresent == TRUE,
    {
      bounds <- quantile(
        FlightExpressionLog2NormalizedCount,
        probs = c(0.05, 0.95),
        na.rm = TRUE,
        names = FALSE
      )
      .(ExpressionScaleLower = bounds[[1L]], ExpressionScaleUpper = bounds[[2L]])
    },
    by = Group
  ]
  matrix_data <- merge(matrix_data, group_scale, by = "Group", all.x = TRUE, sort = FALSE)
  matrix_data[, ExpressionColor := map_gradient_colours(
    FlightExpressionLog2NormalizedCount,
    ExpressionScaleLower[[1L]],
    ExpressionScaleUpper[[1L]],
    expression_palette
  ), by = Group]
  matrix_data[PathwayPresent == FALSE, ExpressionColor := "#F0F2F4"]
  matrix_data[, `:=`(
    TopPValue = tissue_meta_p,
    TopLog2FC = fifelse(
      is.finite(tissue_meta_p) & tissue_meta_p < 0.05 & is.finite(log2FoldChange),
      log2FoldChange,
      NA_real_
    ),
    TopNegLog10P = fifelse(
      is.finite(tissue_meta_p) & tissue_meta_p > 0 & tissue_meta_p < 0.05,
      -log10(tissue_meta_p),
      NA_real_
    )
  )]
  top_rows <- unique(matrix_data[, .(
    Group,
    EnsemblID,
    DisplayLabel,
    GeneRank,
    EligibleGeneCount,
    TopPValue,
    TopLog2FC,
    TopNegLog10P
  )])
  top_rows[, TopLog2FCVisible := is.finite(TopLog2FC)]
  top_rows[, TopPValueVisible := is.finite(TopNegLog10P)]
  top_rows[, TopColourCap := {
    cap <- quantile(abs(TopLog2FC[TopLog2FCVisible]), 0.95, na.rm = TRUE, names = FALSE)
    if (!is.finite(cap) || cap <= 0) cap <- 1
    cap
  }, by = Group]
  top_rows[, TopColor := {
    values <- TopLog2FC
    cap <- TopColourCap[[1L]]
    shifted <- pmin(cap, pmax(-cap, values))
    map_gradient_colours(shifted, -cap, cap, top_palette)
  }, by = Group]
  top_rows[TopLog2FCVisible == FALSE, TopColor := "#FFFFFF"]
  top_rows[, TopPColourCap := {
    cap <- quantile(TopNegLog10P[TopPValueVisible], 0.95, na.rm = TRUE, names = FALSE)
    if (!is.finite(cap) || cap <= 0) cap <- 1
    cap
  }, by = Group]
  p_palette <- colorRampPalette(c("#FFFFFF", "#FEE391", "#D7301F"))(256L)
  top_rows[, TopPColor := map_gradient_colours(
    pmin(TopNegLog10P, TopPColourCap[[1L]]), 0, TopPColourCap[[1L]], p_palette
  ), by = Group]
  top_rows[TopPValueVisible == FALSE, TopPColor := "#FFFFFF"]
  matrix_data <- merge(
    matrix_data,
    top_rows[, .(
      Group, EnsemblID, TopLog2FCVisible, TopPValueVisible,
      TopColourCap, TopPColourCap, TopColor, TopPColor
    )],
    by = c("Group", "EnsemblID"),
    all.x = TRUE,
    sort = FALSE
  )
  fwrite(matrix_data, matrix_path)

  pathway_summary <- matrix_data[
    PathwayPresent == TRUE,
    .(
      PresentGeneCount = uniqueN(EnsemblID),
      MeanFlightExpressionLog2NormalizedCount = mean(FlightExpressionLog2NormalizedCount),
      MedianFlightExpressionLog2NormalizedCount = median(FlightExpressionLog2NormalizedCount),
      NominalPValueLt005GeneCount = uniqueN(EnsemblID[tissue_meta_p < 0.05])
    ),
    by = .(Group, Pathway, PathwayOrder)
  ]
  pathway_summary[, Pathway := factor(Pathway, levels = expected_pathways)]
  setorder(pathway_summary, Group, PathwayOrder)
  fwrite(pathway_summary, summary_path)

  metadata_path <- file.path(skill_outdir, "run_metadata.csv")
  metadata <- fread(metadata_path)
  custom_metadata <- data.table(
    Key = c(
      "Primary rendering value",
      "Primary rendering colour",
      "P-value annotation row",
      "Top row value",
      "Top row colour",
      "Expression sample selection"
    ),
    Value = c(
      "Mission-equal mean flight/uG log2 normalized count",
      "white to red within each tissue (5th to 95th percentile clipping)",
      "-log10(unadjusted tissue meta P), white to red; shown only for P < 0.05",
      "log2FoldChange",
      "blue negative to white zero to red positive",
      "Flight/uG sample columns reconstructed and audited by proliferation-state workflow"
    )
  )
  metadata <- metadata[!Key %chin% custom_metadata$Key]
  fwrite(rbindlist(list(metadata, custom_metadata)), metadata_path)

  y_levels <- c(rev(expected_pathways), "log2FC\nP < 0.05", "-log10(P)\nP < 0.05")
  for (group_name in unique(matrix_data$Group)) {
    panel <- copy(matrix_data[Group == group_name])
    top_panel <- copy(top_rows[Group == group_name])
    panel[, PlotRow := factor(Pathway, levels = y_levels)]
    top_fc_panel <- copy(top_panel)
    top_fc_panel[, PlotRow := factor("log2FC\nP < 0.05", levels = y_levels)]
    top_p_panel <- copy(top_panel)
    top_p_panel[, PlotRow := factor("-log10(P)\nP < 0.05", levels = y_levels)]
    plot <- ggplot() +
      geom_tile(
        data = panel,
        aes(x = GeneRank, y = PlotRow),
        fill = "#F0F2F4", colour = "white", linewidth = 0.08
      ) +
      geom_tile(
        data = panel[PathwayPresent == TRUE],
        aes(x = GeneRank, y = PlotRow, fill = ExpressionColor),
        colour = "white", linewidth = 0.08
      ) +
      geom_tile(
        data = top_fc_panel,
        aes(x = GeneRank, y = PlotRow, fill = TopColor),
        colour = "white", linewidth = 0.08
      ) +
      geom_tile(
        data = top_p_panel,
        aes(x = GeneRank, y = PlotRow, fill = TopPColor),
        colour = "white", linewidth = 0.08
      ) +
      scale_fill_identity() +
      scale_y_discrete(drop = FALSE) +
      labs(
        title = paste0(group_name, ": seven DNA-repair pathways"),
        subtitle = "Main rows: flight/uG log2(normalized count), white to red; annotation rows: log2FC and -log10(unadjusted meta P), only where P < 0.05",
        x = "Genes ordered only by seven-pathway GO membership score",
        y = NULL,
        caption = "Counts are shown for all GO members. log2FC and P-value tiles are shown only for unadjusted tissue-meta P < 0.05; expression colours are independently scaled within tissue."
      ) +
      theme_minimal(base_family = "Arial", base_size = 9) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8, face = "bold"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5),
        plot.caption = element_text(size = 7.5, colour = "#4B5563"),
        plot.margin = margin(8, 12, 8, 8)
      )
    width <- max(8, min(22, 3.5 + unique(panel$EligibleGeneCount)[[1L]] * 0.045))
    stem <- paste0(skill_prefix, "_matrix_", safe_figure_filename(group_name))
    ggsave(file.path(skill_outdir, paste0(stem, ".png")), plot, width = width, height = 4.8, dpi = 300, bg = "white")
    ggsave(file.path(skill_outdir, paste0(stem, ".pdf")), plot, width = width, height = 4.8, device = cairo_pdf, bg = "white")

    summary_panel <- pathway_summary[Group == group_name]
    summary_plot <- ggplot(summary_panel, aes(x = Pathway, y = MeanFlightExpressionLog2NormalizedCount)) +
      geom_col(fill = "#D73027", width = 0.72) +
      labs(
        title = paste0(group_name, ": pathway-member expression summary"),
        subtitle = "Mean flight/uG log2(normalized count) among GO-member genes; no P-value filtering",
        x = NULL,
        y = "Mean log2(normalized count)"
      ) +
      theme_minimal(base_family = "Arial", base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5),
        axis.text.x = element_text(face = "bold")
      )
    summary_stem <- paste0(skill_prefix, "_summary_", safe_figure_filename(group_name))
    ggsave(file.path(skill_outdir, paste0(summary_stem, ".png")), summary_plot, width = 7, height = 4.3, dpi = 300, bg = "white")
    ggsave(file.path(skill_outdir, paste0(summary_stem, ".pdf")), summary_plot, width = 7, height = 4.3, device = cairo_pdf, bg = "white")
  }
  invisible(list(matrix = matrix_data, summary = pathway_summary))
}

# -------------------------------------------------------------------------
# 1. Load raw OSDR rows by Ensembl ID and aggregate every primary analysis
#    unit into the 26 normalized original-Material-Type tissue groups.
# -------------------------------------------------------------------------

tissue_registry <- fread(tissue_registry_path)
required_registry_columns <- c(
  "Accession",
  "原tissue英文",
  "主分析mission",
  "方向校正sign"
)
assert(
  all(required_registry_columns %in% names(tissue_registry)),
  "The primary tissue registry schema has changed."
)

tissue_mapping <- fread(tissue_mapping_path)
required_mapping_columns <- c(
  "accession",
  "analysis_tissue",
  "original_material_type",
  "included_in_exact_go_analysis"
)
assert(
  all(required_mapping_columns %in% names(tissue_mapping)),
  "The normalized-original-Material-Type tissue mapping schema has changed."
)
tissue_mapping[, included := tolower(trimws(as.character(
  included_in_exact_go_analysis
))) %chin% c("true", "t", "1", "yes", "y")]
tissue_mapping <- unique(tissue_mapping[
  included == TRUE,
  .(
    accession = trimws(as.character(accession)),
    analysis_tissue = trimws(as.character(analysis_tissue)),
    original_material_type = trimws(as.character(original_material_type))
  )
])
assert(
  !anyDuplicated(tissue_mapping$accession) &&
    !anyNA(tissue_mapping$analysis_tissue) &&
    all(nzchar(tissue_mapping$analysis_tissue)),
  "The 26-tissue mapping has duplicated accessions or empty tissue labels."
)
selected_tissues <- sort(unique(tissue_mapping$analysis_tissue))
assert(
  length(selected_tissues) == 26L,
  paste("The current tissue contract requires exactly 26 tissues; found", length(selected_tissues))
)

tissue_registry <- merge(
  tissue_registry,
  tissue_mapping,
  by.x = "Accession",
  by.y = "accession",
  all = FALSE,
  sort = FALSE
)
assert(
  nrow(tissue_registry) == nrow(tissue_mapping),
  "The 26-tissue mapping is not one-to-one with the primary tissue registry."
)

selected_accessions <- unique(
  tissue_registry[
    ,
    .(
      accession = Accession,
      Group = analysis_tissue,
      original_tissue = original_material_type,
      registry_mission = `主分析mission`,
      registry_orientation_multiplier = `方向校正sign`
    )
  ]
)
assert(
  setequal(selected_accessions$Group, selected_tissues),
  "Not all 26 analysis tissues were found in the primary registry."
)
assert(
  !anyDuplicated(selected_accessions[, .(accession, Group)]),
  "The tissue registry contains duplicated accession-tissue assignments."
)

contrast_manifest <- fread(dataset_manifest_path, fill = TRUE)
required_manifest_columns <- c(
  "accession",
  "contrast_id",
  "analysis_role",
  "input_relative_path",
  "raw_log2fc_column",
  "orientation_multiplier",
  "mission_cluster",
  "unified_contrast_label"
)
assert(
  all(required_manifest_columns %in% names(contrast_manifest)),
  "The selected-contrast manifest schema has changed."
)
primary_contrast_manifest <- contrast_manifest[analysis_role == "primary"]
assert(
  setequal(unique(primary_contrast_manifest$accession), selected_accessions$accession),
  "The 26-tissue mapping does not cover exactly the selected primary accessions."
)
contrast_manifest <- merge(
  primary_contrast_manifest,
  selected_accessions,
  by = "accession",
  all = FALSE,
  sort = FALSE
)
contrast_manifest[, contrast_id := trimws(as.character(contrast_id))]
contrast_manifest[is.na(contrast_id), contrast_id := ""]
contrast_manifest[, analysis_unit_id := fifelse(
  nzchar(contrast_id),
  paste(accession, contrast_id, sep = "__"),
  accession
)]
assert(
  nrow(contrast_manifest) == nrow(primary_contrast_manifest) &&
    setequal(unique(contrast_manifest$Group), selected_tissues),
  "The 26-tissue mapping does not resolve every primary analysis unit."
)
assert(
  !anyDuplicated(contrast_manifest$analysis_unit_id),
  "The selected contrast manifest has duplicated analysis-unit IDs."
)
assert(
  all(contrast_manifest$unified_contrast_label == "Space Flight vs Ground"),
  "The selected RNA effects are not uniformly Space Flight vs Ground."
)

load_ensembl_keyed_unit <- function(manifest_row) {
  input_path <- file.path(root, manifest_row$input_relative_path)
  assert(file.exists(input_path), paste("Missing raw DE file:", input_path))

  log2fc_column <- manifest_row$raw_log2fc_column
  stat_column <- companion_column(log2fc_column, "Stat")
  raw_p_column <- companion_column(log2fc_column, "P.value")
  adjusted_p_column <- companion_column(log2fc_column, "Adj.p.value")
  header <- names(fread(input_path, nrows = 0L, showProgress = FALSE))
  required_columns <- c(
    "ENSEMBL",
    log2fc_column,
    raw_p_column,
    adjusted_p_column
  )
  assert(
    all(required_columns %in% header),
    paste(
      manifest_row$analysis_unit_id,
      "lacks raw DE columns:",
      paste(setdiff(required_columns, header), collapse = ", ")
    )
  )
  selected_columns <- unique(c(
    "ENSEMBL",
    if ("SYMBOL" %in% header) "SYMBOL" else character(),
    log2fc_column,
    if (stat_column %in% header) stat_column else character(),
    raw_p_column,
    adjusted_p_column
  ))
  raw <- fread(
    input_path,
    select = selected_columns,
    showProgress = FALSE
  )
  if (!"SYMBOL" %in% names(raw)) raw[, SYMBOL := ""]

  multiplier <- as.integer(manifest_row$orientation_multiplier)
  assert(
    multiplier %in% c(-1L, 1L),
    paste("Invalid orientation multiplier for", manifest_row$analysis_unit_id)
  )
  prepared <- data.table(
    EnsemblID = normalize_ensembl(raw$ENSEMBL),
    gene_symbol = as.character(raw$SYMBOL),
    log2fc_spaceflight_vs_ground =
      as.numeric(raw[[log2fc_column]]) * multiplier,
    raw_p_value = as.numeric(raw[[raw_p_column]]),
    adjusted_p_value = as.numeric(raw[[adjusted_p_column]])
  )
  statistic_source <- if (stat_column %in% names(raw)) {
    prepared[, stat_spaceflight_vs_ground :=
      as.numeric(raw[[stat_column]]) * multiplier]
    "OSDR_Stat"
  } else {
    prepared[, stat_spaceflight_vs_ground :=
      signed_z_from_p(raw_p_value, log2fc_spaceflight_vs_ground)]
    "signed_Z_from_raw_P"
  }
  invalid_id <- !grepl("^ENSMUSG[0-9]+$", prepared$EnsemblID)
  nonfinite_effect <- !is.finite(prepared$log2fc_spaceflight_vs_ground)
  nonfinite_stat <- !is.finite(prepared$stat_spaceflight_vs_ground)
  keep <- !invalid_id & !nonfinite_effect & !nonfinite_stat
  prepared <- prepared[keep]

  assert(
    !anyDuplicated(prepared$EnsemblID),
    paste(
      manifest_row$analysis_unit_id,
      "contains duplicated raw Ensembl Gene IDs."
    )
  )
  assert(
    nrow(prepared) > 0L,
    paste("No usable Ensembl-keyed rows for", manifest_row$analysis_unit_id)
  )

  prepared[, `:=`(
    Group = manifest_row$Group,
    accession = manifest_row$accession,
    analysis_unit_id = manifest_row$analysis_unit_id,
    analysis_role = manifest_row$analysis_role,
    mission_cluster = manifest_row$mission_cluster,
    unified_contrast_label = manifest_row$unified_contrast_label
  )]

  audit <- data.table(
    Group = manifest_row$Group,
    accession = manifest_row$accession,
    analysis_unit_id = manifest_row$analysis_unit_id,
    mission_cluster = manifest_row$mission_cluster,
    original_tissue = manifest_row$original_tissue,
    registry_mission = manifest_row$registry_mission,
    orientation_multiplier = multiplier,
    unified_contrast_label = manifest_row$unified_contrast_label,
    input_relative_path = manifest_row$input_relative_path,
    input_sha256 = sha256_file(input_path),
    statistic_source = statistic_source,
    raw_rows = nrow(raw),
    valid_unique_Ensembl_rows = nrow(prepared),
    invalid_Ensembl_rows = sum(invalid_id),
    nonfinite_effect_rows = sum(nonfinite_effect),
    nonfinite_signed_stat_rows = sum(nonfinite_stat),
    duplicate_valid_Ensembl_rows =
      sum(duplicated(normalize_ensembl(raw$ENSEMBL)[!invalid_id]))
  )
  list(data = prepared, audit = audit)
}

loaded_units <- lapply(seq_len(nrow(contrast_manifest)), function(index) {
  load_ensembl_keyed_unit(contrast_manifest[index])
})
selected_stats <- rbindlist(
  lapply(loaded_units, function(value) value$data),
  use.names = TRUE
)
unit_audit <- rbindlist(
  lapply(loaded_units, function(value) value$audit),
  use.names = TRUE
)

proliferation_expression_audit <- fread(proliferation_expression_audit_path)
required_expression_audit_columns <- c(
  "analysis_unit_id",
  "input_relative_path",
  "flight_samples"
)
assert(
  all(required_expression_audit_columns %in% names(proliferation_expression_audit)),
  "Proliferation expression audit schema has changed."
)
flight_expression_units <- lapply(seq_len(nrow(contrast_manifest)), function(index) {
  read_flight_expression_unit(
    contrast_manifest[index],
    proliferation_expression_audit
  )
})
selected_flight_expression <- rbindlist(
  lapply(flight_expression_units, function(value) value$data),
  use.names = TRUE
)
flight_expression_audit <- rbindlist(
  lapply(flight_expression_units, function(value) value$audit),
  use.names = TRUE
)
assert(
  !anyDuplicated(selected_flight_expression[, .(analysis_unit_id, EnsemblID)]) &&
    all(is.finite(selected_flight_expression$FlightExpressionLog2NormalizedCount)),
  "Flight/uG normalized-count expression is invalid or duplicated."
)

assert(
  all(grepl("^ENSMUSG[0-9]+$", selected_stats$EnsemblID)),
  "Selected RNA rows contain invalid mouse Ensembl Gene IDs."
)
assert(
  all(is.finite(selected_stats$log2fc_spaceflight_vs_ground)) &&
    all(is.finite(selected_stats$stat_spaceflight_vs_ground)),
  "Selected RNA effects or signed statistics contain non-finite values."
)
assert(
  !anyDuplicated(selected_stats[, .(analysis_unit_id, EnsemblID)]),
  "Selected RNA input has duplicated analysis-unit + Ensembl rows."
)

group_order <- data.table(
  Group = selected_tissues,
  GroupOrder = seq_along(selected_tissues)
)

unit_audit <- merge(unit_audit, group_order, by = "Group", sort = FALSE)
setorder(unit_audit, GroupOrder, mission_cluster, analysis_unit_id)
unit_audit[, GroupOrder := NULL]

display_labels <- rbindlist(list(
  selected_stats[, .(EnsemblID, gene_symbol)],
  selected_flight_expression[, .(EnsemblID, gene_symbol)]
))[
  ,
  .(DisplayLabel = first_display_label(gene_symbol)),
  by = EnsemblID
]

# Units from the same mission are averaged before missions are combined.
# This avoids counting same-mission data as independent mission evidence.
mission_gene <- selected_stats[
  ,
  .(
    mission_log2FoldChange =
      median(log2fc_spaceflight_vs_ground),
    mission_signed_stat =
      mean(stat_spaceflight_vs_ground),
    n_analysis_units = uniqueN(analysis_unit_id),
    analysis_units = collapse_values(analysis_unit_id),
    accessions = collapse_values(accession)
  ),
  by = .(Group, EnsemblID, mission_cluster)
]

tissue_rna <- mission_gene[
  ,
  {
    combined_z <- sum(mission_signed_stat) / sqrt(.N)
    .(
      log2FoldChange = median(mission_log2FoldChange),
      tissue_signed_stouffer_z = combined_z,
      tissue_meta_p = 2 * pnorm(-abs(combined_z)),
      n_missions = .N,
      n_analysis_units = sum(n_analysis_units),
      missions = collapse_values(mission_cluster),
      analysis_units = collapse_values(analysis_units),
      accessions = collapse_values(accessions)
    )
  },
  by = .(Group, EnsemblID)
]
tissue_rna[, padj := bh_adjust(tissue_meta_p), by = Group]
tissue_flight_expression <- selected_flight_expression[
  ,
  .(
    mission_flight_expression_log2_normalized_count =
      mean(FlightExpressionLog2NormalizedCount),
    n_analysis_units_expression = uniqueN(analysis_unit_id),
    n_flight_samples = sum(n_flight_samples)
  ),
  by = .(Group, EnsemblID, mission_cluster)
][
  ,
  .(
    FlightExpressionLog2NormalizedCount =
      mean(mission_flight_expression_log2_normalized_count),
    expression_n_missions = .N,
    expression_n_analysis_units = sum(n_analysis_units_expression),
    expression_n_flight_samples = sum(n_flight_samples)
  ),
  by = .(Group, EnsemblID)
]
tissue_rna <- merge(
  tissue_rna,
  tissue_flight_expression,
  by = c("Group", "EnsemblID"),
  all.x = TRUE,
  sort = FALSE
)
tissue_rna <- merge(
  tissue_rna,
  display_labels,
  by = "EnsemblID",
  all.x = TRUE,
  sort = FALSE
)
tissue_rna[is.na(DisplayLabel), DisplayLabel := ""]
tissue_rna <- merge(tissue_rna, group_order, by = "Group", sort = FALSE)
setorder(tissue_rna, GroupOrder, EnsemblID)
tissue_rna[, GroupOrder := NULL]

assert(
  !anyDuplicated(tissue_rna[, .(Group, EnsemblID)]),
  "Prepared RNA input has duplicated Group + EnsemblID rows."
)
assert(
  all(is.finite(tissue_rna$log2FoldChange)) &&
    all(is.finite(tissue_rna$padj)) &&
    all(is.finite(tissue_rna$FlightExpressionLog2NormalizedCount)),
  "Prepared tissue RNA effects, BH-FDR values, or flight/uG expression are non-finite."
)

tissue_meta_audit <- tissue_rna[
  ,
  .(
    genes = .N,
    accessions = uniqueN(unlist(strsplit(accessions, ";", fixed = TRUE))),
    analysis_units =
      uniqueN(unlist(strsplit(analysis_units, ";", fixed = TRUE))),
    missions = uniqueN(unlist(strsplit(missions, ";", fixed = TRUE))),
    median_log2FoldChange = median(log2FoldChange),
    min_log2FoldChange = min(log2FoldChange),
    max_log2FoldChange = max(log2FoldChange),
    median_flight_expression_log2_normalized_count =
      median(FlightExpressionLog2NormalizedCount),
    min_flight_expression_log2_normalized_count =
      min(FlightExpressionLog2NormalizedCount),
    max_flight_expression_log2_normalized_count =
      max(FlightExpressionLog2NormalizedCount),
    padj_le_0_05 = sum(padj <= 0.05),
    positive_padj_le_0_05 =
      sum(padj <= 0.05 & log2FoldChange > 0),
    negative_padj_le_0_05 =
      sum(padj <= 0.05 & log2FoldChange < 0)
  ),
  by = Group
]
tissue_meta_audit <- merge(
  tissue_meta_audit,
  group_order,
  by = "Group",
  sort = FALSE
)
setorder(tissue_meta_audit, GroupOrder)
tissue_meta_audit[, GroupOrder := NULL]

# -------------------------------------------------------------------------
# 2. Build a mouse-specific catalog and direct mouse gene-to-GO input.
#    BP terms must descend from DNA repair (GO:0006281) through the declared
#    relation scope for the selected seven-pathway contract.
#    MF and CC are retained only if already pathway-specific in the skill
#    catalog and directly supported by a positive mouse MGI annotation.
# -------------------------------------------------------------------------

catalog_build_output <- system2(
  file.path(R.home("bin"), "Rscript"),
  args = c(
    shQuote(catalog_builder_script_path),
    "--obo", shQuote(obo_path),
    "--out", shQuote(project_catalog_path),
    "--provenance-out", shQuote(project_catalog_provenance_path),
    "--rules", shQuote(project_rules_path),
    "--source-url", "project_frozen_go_basic_obo"
  ),
  stdout = TRUE,
  stderr = TRUE
)
catalog_build_status <- attr(catalog_build_output, "status")
if (is.null(catalog_build_status)) catalog_build_status <- 0L
writeLines(catalog_build_output, project_catalog_log_path)
assert(
  identical(as.integer(catalog_build_status), 0L),
  paste(
    "The project-curated seven-pathway GO catalog build failed:",
    paste(catalog_build_output, collapse = "\n")
  )
)
full_catalog <- fread(project_catalog_path)
required_catalog_columns <- c(
  "GO_ID",
  "CanonicalGOID",
  "IsAlternativeID",
  "GO_Name",
  "Ontology",
  "Pathway",
  "PathwayOrder",
  "Coefficient",
  "OntologyRelease"
)
assert(
  all(required_catalog_columns %in% names(full_catalog)),
  "The skill GO catalog schema has changed."
)
assert(
  length(unique(full_catalog$OntologyRelease)) == 1L &&
    nzchar(unique(full_catalog$OntologyRelease)[[1L]]),
  "The project-curated GO catalog has no unique ontology release."
)

dna_repair_go_id <- "GO:0006281"
obo_edges <- parse_obo_relation_edges(obo_path, dna_repair_relations)
dna_repair_bp_terms <- ontology_descendants(obo_edges, dna_repair_go_id)
catalog <- full_catalog[
  Ontology %chin% c("CC", "MF") |
    (Ontology == "BP" & CanonicalGOID %chin% dna_repair_bp_terms)
]
assert(
  setequal(
    unique(catalog$Pathway),
    c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
  ),
  "DNA-repair branch filtering removed an entire required pathway."
)

gaf <- read_mgi_gaf(gaf_path)
eligible_gaf <- gaf[
  DB == "MGI" &
    grepl("(^|\\|)(taxon:|NCBITaxon:)10090($|\\|)", TAXON) &
    !grepl("(^|\\|)NOT($|\\|)", QUALIFIER) &
    GO_TERM %chin% catalog$GO_ID
]
eligible_gaf <- unique(
  eligible_gaf[
    ,
    .(
      mgi_id = DB_Object_ID,
      source_symbol = SYMBOL,
      GO_ID = GO_TERM,
      GO_ASPECT,
      QUALIFIER,
      EVIDENCE_CODE,
      DB_REFERENCE,
      source_row_number
    )
  ]
)
assert(nrow(eligible_gaf) > 0L, "No MGI GAF rows matched the skill catalog.")

catalog[, DirectMouseMGIAnnotation := GO_ID %chin% unique(eligible_gaf$GO_ID)]
catalog[, `:=`(
  BPDescendsFromDNARepair = Ontology != "BP" | CanonicalGOID %chin% dna_repair_bp_terms,
  IncludedInMouseCatalog = DirectMouseMGIAnnotation
)]
catalog_scope_audit <- merge(
  full_catalog[, .(GO_ID, CanonicalGOID, Pathway, Ontology)],
  catalog[, .(
    GO_ID,
    CanonicalGOID,
    Pathway,
    BPDescendsFromDNARepair,
    DirectMouseMGIAnnotation,
    IncludedInMouseCatalog
  )],
  by = c("GO_ID", "CanonicalGOID", "Pathway"),
  all.x = TRUE,
  sort = FALSE
)
catalog_scope_audit[is.na(BPDescendsFromDNARepair), BPDescendsFromDNARepair := FALSE]
catalog_scope_audit[is.na(DirectMouseMGIAnnotation), DirectMouseMGIAnnotation := FALSE]
catalog_scope_audit[is.na(IncludedInMouseCatalog), IncludedInMouseCatalog := FALSE]
catalog <- catalog[IncludedInMouseCatalog == TRUE]
catalog[, c("DirectMouseMGIAnnotation", "BPDescendsFromDNARepair", "IncludedInMouseCatalog") := NULL]
assert(
  setequal(
    unique(catalog$Pathway),
    c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
  ),
  "Mouse direct annotations did not support every required pathway."
)

stable_mapping <- fread(stable_mapping_audit_path, showProgress = FALSE)
required_stable_mapping_columns <- c(
  "mgi_id",
  "final_ensembl_gene_id",
  "final_mapping_decision",
  "n_mapping_routes",
  "mapping_routes",
  "candidate_ensembl_gene_ids",
  "manual_decision",
  "evidence"
)
assert(
  all(required_stable_mapping_columns %in% names(stable_mapping)),
  "Five-relation MGI-to-Ensembl mapping audit has an unexpected schema."
)
stable_mapping <- stable_mapping[
  mgi_id %chin% eligible_gaf$mgi_id,
  ..required_stable_mapping_columns
]

mapping_audit <- data.table(mgi_id = sort(unique(eligible_gaf$mgi_id)))
mapping_audit <- merge(
  mapping_audit,
  stable_mapping,
  by = "mgi_id",
  all.x = TRUE,
  sort = FALSE
)
mapping_audit <- merge(
  mapping_audit,
  eligible_gaf[
    ,
    .(
      source_symbol = collapse_values(source_symbol),
      annotated_GO_IDs = uniqueN(GO_ID),
      annotation_rows = .N
    ),
    by = mgi_id
  ],
  by = "mgi_id",
  all.x = TRUE,
  sort = FALSE
)
mapping_audit[, included :=
  !is.na(final_ensembl_gene_id) &
    grepl("^ENSMUSG[0-9]+$", final_ensembl_gene_id)
]
setorder(mapping_audit, mgi_id)

assert(
  !anyNA(mapping_audit$final_mapping_decision),
  paste(
    "Catalog-matched MGI gene(s) are absent from the audited five-relation ",
    "MGI-to-Ensembl table:",
    paste(mapping_audit[is.na(final_mapping_decision), mgi_id], collapse = ", ")
  )
)
assert(
  !anyDuplicated(
    mapping_audit[included == TRUE, .(mgi_id, final_ensembl_gene_id)]
  ),
  "The final MGI-to-Ensembl mapping contains duplicated ID pairs."
)

mapped_gaf <- merge(
  eligible_gaf,
  mapping_audit[
    included == TRUE,
    .(
      mgi_id,
      EnsemblID = final_ensembl_gene_id,
      final_mapping_decision
    )
  ],
  by = "mgi_id",
  all = FALSE,
  sort = FALSE
)
go_names <- catalog[
  ,
  .(GO_Name = first_display_label(GO_Name)),
  by = GO_ID
]
go_input <- unique(mapped_gaf[, .(EnsemblID, GO_ID)])
go_input <- merge(
  go_input,
  go_names,
  by = "GO_ID",
  all.x = TRUE,
  sort = FALSE
)
setcolorder(go_input, c("EnsemblID", "GO_ID", "GO_Name"))
setorder(go_input, EnsemblID, GO_ID)

assert(
  all(grepl("^ENSMUSG[0-9]+$", go_input$EnsemblID)) &&
    all(grepl("^GO:[0-9]{7}$", go_input$GO_ID)),
  "Prepared GO input contains invalid stable identifiers."
)
assert(
  !anyDuplicated(go_input[, .(EnsemblID, GO_ID)]),
  "Prepared GO input contains duplicated Ensembl-GO pairs."
)

mapped_gaf_pathways <- merge(
  unique(mapped_gaf[, .(mgi_id, EnsemblID, GO_ID)]),
  unique(catalog[, .(GO_ID, Pathway)]),
  by = "GO_ID",
  all.x = TRUE,
  allow.cartesian = TRUE,
  sort = FALSE
)
assert(
  setequal(
    unique(mapped_gaf_pathways$Pathway),
    c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
  ),
  "Stable-ID mapping left at least one required repair pathway without a mapped mouse gene."
)
go_coverage_audit <- merge(
  catalog[
    ,
    .(
      catalog_rows = .N,
      canonical_GO_IDs = uniqueN(CanonicalGOID)
    ),
    by = .(Pathway, PathwayOrder)
  ],
  mapped_gaf_pathways[
    ,
    .(
      observed_GO_IDs = uniqueN(GO_ID),
      mapped_MGI_genes = uniqueN(mgi_id),
      mapped_Ensembl_genes = uniqueN(EnsemblID),
      direct_gene_GO_pairs = uniqueN(paste(EnsemblID, GO_ID))
    ),
    by = Pathway
  ],
  by = "Pathway",
  all.x = TRUE,
  sort = FALSE
)
for (column in c(
  "observed_GO_IDs",
  "mapped_MGI_genes",
  "mapped_Ensembl_genes",
  "direct_gene_GO_pairs"
)) {
  set(
    go_coverage_audit,
    i = which(is.na(go_coverage_audit[[column]])),
    j = column,
    value = 0L
  )
}
setorder(go_coverage_audit, PathwayOrder)

# -------------------------------------------------------------------------
# 3. Write prepared inputs, provenance, and invoke the named skill.
# -------------------------------------------------------------------------

rna_input_path <- file.path(
  input_dir,
  "02_selected_tissue_RNA_meta_Ensembl116.tsv.gz"
)
go_input_path <- file.path(
  input_dir,
  "03_mouse_direct_GO_annotations_for_skill_catalog.tsv.gz"
)
mouse_catalog_path <- file.path(
  input_dir,
  "08_mouse_DNA_repair_branch_seven_pathway_catalog.csv"
)

fwrite(
  unit_audit,
  file.path(input_dir, "01_selected_tissue_analysis_units.csv")
)
fwrite(
  tissue_mapping,
  file.path(input_dir, "00_normalized_original_material_26_tissue_mapping.csv")
)
fwrite(tissue_rna, rna_input_path, sep = "\t")
fwrite(go_input, go_input_path, sep = "\t")
fwrite(
  mapping_audit,
  file.path(input_dir, "04_MGI_to_Ensembl116_mapping_audit.csv")
)
fwrite(
  tissue_meta_audit,
  file.path(input_dir, "05_tissue_RNA_meta_audit.csv")
)
fwrite(
  go_coverage_audit,
  file.path(input_dir, "06_GO_catalog_annotation_coverage.csv")
)
fwrite(catalog_scope_audit, file.path(input_dir, "08_mouse_catalog_scope_audit.csv"))
fwrite(catalog, mouse_catalog_path)

provenance <- data.table(
  Item = c(
    "Selected contrast manifest",
    "Selected-tissue registry",
    "Normalized original-Material-Type 26-tissue mapping",
    "MGI MOD-GAF",
    "Frozen GO ontology",
    "Audited five-relation MGI-to-Ensembl mapping",
    "Five-relation annotation source manifest",
    "Project-curated GO seed rules",
    "Project-curated GO catalog",
    "Skill pathway display",
  "Skill plotting script",
  "Proliferation expression sample-membership audit"
  ),
  Path = c(
    dataset_manifest_path,
    tissue_registry_path,
    tissue_mapping_path,
    gaf_path,
    obo_path,
    stable_mapping_audit_path,
    annotation_source_manifest_path,
    project_rules_path,
    project_catalog_path,
    display_path,
  skill_script_path,
  proliferation_expression_audit_path
  )
)
provenance[, SHA256 := vapply(Path, sha256_file, character(1))]
provenance[, Note := c(
  "Raw OSDR DE columns and audited Space Flight vs Ground orientation",
  "Raw OSDR accession metadata retained for traceability",
  "Authoritative 26-tissue grouping; preserves anatomical subparts and distinct tissue types",
  "Mouse direct GO annotations; TAXON taxon:10090 or NCBITaxon:10090; NOT excluded",
  paste0("DNA-repair BP closure uses ", paste(dna_repair_relations, collapse = ", ")),
  "Official multi-route MGI to Entrez to Ensembl release 116 mapping audit",
  "Frozen source hashes and official-ID manual-exception evidence",
  "Seven-pathway seed rules for the selected retained contract",
  "Catalog generated from frozen GO OBO using the project rules",
  "Fixed BER, NER, MMR, FA, HR, AEJ, NHEJ order; primary matrix uses signed response effects",
  "Named skill supplies fixed GO membership/order, signed-effect matrices, and audit tables",
  "Defines audited flight/uG sample columns used for normalized-count expression"
)]
fwrite(
  provenance,
  file.path(input_dir, "07_input_provenance_and_sha256.csv")
)

skill_command_args <- c(
  shQuote(skill_script_path),
  "--rna", shQuote(rna_input_path),
  "--go", shQuote(go_input_path),
  "--outdir", shQuote(skill_outdir),
  "--id-column", "EnsemblID",
  "--go-id-column", "GO_ID",
  "--go-name-column", "GO_Name",
  "--group-column", "Group",
  "--label-column", "DisplayLabel",
  "--effect-column", "log2FoldChange",
  "--formats", "png,pdf",
  "--font", "Arial",
  "--prefix", skill_prefix,
  "--catalog", shQuote(mouse_catalog_path),
  "--display", shQuote(display_path)
)
if (!ignore_pvalue) {
  skill_command_args <- append(
    skill_command_args,
    c("--padj-column", "padj", "--padj-cutoff", "0.05"),
    after = 17L
  )
}

skill_output <- system2(
  file.path(R.home("bin"), "Rscript"),
  args = skill_command_args,
  stdout = TRUE,
  stderr = TRUE
)
skill_status <- attr(skill_output, "status")
if (is.null(skill_status)) skill_status <- 0L
writeLines(
  c(
    paste(
      shQuote(file.path(R.home("bin"), "Rscript")),
      paste(skill_command_args, collapse = " ")
    ),
    "",
    skill_output
  ),
  run_log_path
)
assert(
  identical(as.integer(skill_status), 0L),
  paste(
    "The plot-go-pathway-linear-matrix skill failed:",
    paste(skill_output, collapse = "\n")
  )
)

# Preserve the named skill's stable-ID membership/order audit, then render the
# project matrix with counts, log2FC and P value shown at the same time.
render_expression_matrix_outputs(
  skill_outdir,
  skill_prefix,
  tissue_rna,
  c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
)

run_summary <- data.table(
  Metric = c(
    "Requested tissues",
    "Tissue grouping definition",
    "Tissue grouping mapping",
    "DNA-repair BP relation scope",
    "Primary accessions",
    "Selected analysis units",
    "Selected mission clusters across tissues",
    "Prepared RNA Group-Ensembl rows",
    "Prepared RNA unique Ensembl IDs",
    "Full project-curated catalog rows",
    "Mouse DNA-repair catalog rows",
    "Mouse DNA-repair catalog canonical GO terms",
    "Mouse DNA-repair catalog alternative GO IDs",
    "Catalog-matched positive mouse GAF rows",
    "Catalog-observed direct GO IDs",
    "Catalog GO IDs without a positive direct mouse annotation",
    "Catalog-matched MGI genes",
    "Mapped MGI genes",
    "Prepared unique Ensembl-GO pairs",
    "Primary matrix PNG files",
    "Primary matrix PDF files",
    "Supplementary summary PNG files",
    "Supplementary summary PDF files",
    "Primary matrix value",
    "P-value annotation row",
    "Top log2FC row",
    "Random seed"
  ),
  Value = c(
    length(selected_tissues),
    "26 normalized original-Material-Type groups; anatomical subparts retained",
    tissue_mapping_path,
    paste(dna_repair_relations, collapse = ","),
    uniqueN(contrast_manifest$accession),
    uniqueN(unit_audit$analysis_unit_id),
    uniqueN(paste(unit_audit$Group, unit_audit$mission_cluster)),
    nrow(tissue_rna),
    uniqueN(tissue_rna$EnsemblID),
    nrow(full_catalog),
    nrow(catalog),
    uniqueN(catalog$CanonicalGOID),
    uniqueN(catalog[IsAlternativeID == TRUE, GO_ID]),
    nrow(eligible_gaf),
    uniqueN(eligible_gaf$GO_ID),
    uniqueN(catalog$GO_ID) - uniqueN(eligible_gaf$GO_ID),
    uniqueN(eligible_gaf$mgi_id),
    mapping_audit[included == TRUE, uniqueN(mgi_id)],
    nrow(go_input),
    length(list.files(
      skill_outdir,
      pattern = paste0("^", skill_prefix, "_matrix_.*\\.png$")
    )),
    length(list.files(
      skill_outdir,
      pattern = paste0("^", skill_prefix, "_matrix_.*\\.pdf$")
    )),
    length(list.files(
      skill_outdir,
      pattern = paste0("^", skill_prefix, "_summary_.*\\.png$")
    )),
    length(list.files(
      skill_outdir,
      pattern = paste0("^", skill_prefix, "_summary_.*\\.pdf$")
    )),
    "Mission-equal mean flight/uG log2(normalized count); white to red per tissue",
    "-log10(unadjusted tissue meta P), white to red; shown only for P < 0.05",
    "Shown only for P < 0.05; blue negative to white zero to red positive",
    25L
  )
)
fwrite(run_summary, run_summary_path)
writeLines(
  c(
    "# 七修复通路 26 组织矩阵",
    "",
    paste0("- 组织数：", length(selected_tissues), "。"),
    "- 分组：基于原始 Material Type 的标准化映射；保留解剖亚区和不同组织类型。",
    "- 三个此前被合并的标签现独立保留：Heart / Heart right ventricle、Left lobe of the liver、Spleen-distal。",
    "- RNA 与 GO 的匹配主键均为 Ensembl Gene ID；MGI/Entrez 仅用于来源映射，Symbol 只作显示。",
    paste0("- DNA repair BP 后代关系：", paste(dna_repair_relations, collapse = ", "), "。"),
    "- 主矩阵同时保留成员计数、任务等权 flight/uG log2(normalized count)、log2FC 和 P 值；不将正负方向解释为通路活化或抑制。",
    "- `prepared_inputs/00_normalized_original_material_26_tissue_mapping.csv` 记录每个 accession 的分组。"
  ),
  file.path(outdir, "README.md")
)

message(
  "Created joint counts/log2FC/P-value seven-pathway matrices for ",
  length(selected_tissues),
  " tissues in: ",
  outdir
)
