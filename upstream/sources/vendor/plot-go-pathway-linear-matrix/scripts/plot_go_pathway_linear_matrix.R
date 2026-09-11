#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "scripts/plot_go_pathway_linear_matrix.R",
    mustWork = TRUE
  )
}
skill_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript plot_go_pathway_linear_matrix.R",
      "--rna RNA.tsv --go DIRECT_GO.tsv --outdir OUTPUT_DIR",
      "[--id-column EnsemblID] [--go-id-column GO_ID]",
      "[--go-name-column GO_Name] [--group-column Group]",
      "[--label-column DisplayLabel]",
      "[--effect-column log2FoldChange] [--padj-column padj]",
      "[--padj-cutoff 0.05] [--formats png,pdf,svg]",
      "[--font 'Arial Unicode MS'] [--prefix go_pathway_linear]",
      "[--catalog FILE] [--display FILE]",
      sep = "\n  "
    ),
    "\n"
  )
}

parse_args <- function(arguments) {
  if (!length(arguments) || any(arguments %in% c("-h", "--help"))) {
    usage()
    quit(status = if (length(arguments)) 0L else 2L)
  }
  parsed <- list()
  index <- 1L
  while (index <= length(arguments)) {
    key <- arguments[[index]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    if (index == length(arguments)) {
      stop("Missing value for ", key, call. = FALSE)
    }
    parsed[[sub("^--", "", key)]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  parsed
}

options <- parse_args(commandArgs(trailingOnly = TRUE))
defaults <- list(
  `id-column` = "EnsemblID",
  `go-id-column` = "GO_ID",
  `go-name-column` = "GO_Name",
  `group-column` = "Group",
  `label-column` = "DisplayLabel",
  `effect-column` = "",
  `padj-column` = "",
  `padj-cutoff` = "0.05",
  formats = "png,pdf,svg",
  font = "Arial Unicode MS",
  prefix = "go_pathway_linear",
  catalog = file.path(
    skill_root,
    "assets",
    "seven_pathway_go_term_catalog.csv"
  ),
  display = file.path(
    skill_root,
    "assets",
    "seven_pathway_display.csv"
  )
)
for (name in names(defaults)) {
  if (is.null(options[[name]])) options[[name]] <- defaults[[name]]
}

required_arguments <- c("rna", "go", "outdir")
missing_arguments <- required_arguments[
  !vapply(required_arguments, function(name) {
    !is.null(options[[name]]) && nzchar(options[[name]])
  }, logical(1))
]
if (length(missing_arguments)) {
  stop(
    "Missing required argument(s): ",
    paste(paste0("--", missing_arguments), collapse = ", "),
    call. = FALSE
  )
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

rna_path <- normalizePath(options$rna, mustWork = TRUE)
go_path <- normalizePath(options$go, mustWork = TRUE)
catalog_path <- normalizePath(options$catalog, mustWork = TRUE)
display_path <- normalizePath(options$display, mustWork = TRUE)
outdir <- normalizePath(options$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, mustWork = TRUE)

id_column <- options[["id-column"]]
go_id_column <- options[["go-id-column"]]
go_name_column <- options[["go-name-column"]]
group_column <- options[["group-column"]]
label_column <- options[["label-column"]]
effect_column <- options[["effect-column"]]
padj_column <- options[["padj-column"]]
padj_cutoff <- as.numeric(options[["padj-cutoff"]])
formats <- unique(trimws(strsplit(options$formats, ",", fixed = TRUE)[[1L]]))
font_family <- options$font
file_prefix <- options$prefix

assert(
  is.finite(padj_cutoff) && padj_cutoff >= 0 && padj_cutoff <= 1,
  "--padj-cutoff must be between 0 and 1."
)
assert(
  length(formats) > 0L && all(formats %in% c("png", "pdf", "svg")),
  "--formats may contain only png, pdf and svg."
)

normalize_ensembl <- function(values) {
  values <- trimws(as.character(values))
  values <- sub("\\.[0-9]+$", "", values)
  values[is.na(values)] <- ""
  values
}

safe_filename <- function(value) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", as.character(value))
  value <- gsub("^_+|_+$", "", value)
  ifelse(nzchar(value), value, "group")
}

read_input <- function(file_path) {
  result <- fread(file_path, na.strings = c("", "NA", "NaN", "null"))
  assert(nrow(result) > 0L, paste("Input table is empty:", file_path))
  result
}

rna <- read_input(rna_path)
go_annotations <- read_input(go_path)
display <- fread(display_path)
catalog <- fread(catalog_path)
rna_input_columns <- names(rna)

expected_pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
expected_coefficients <- seq_along(expected_pathways)
setorder(display, PathwayOrder)
assert(
  identical(display$Pathway, expected_pathways) &&
    identical(as.integer(display$PathwayOrder), expected_coefficients) &&
    identical(as.numeric(display$Coefficient), as.numeric(expected_coefficients)) &&
    !anyNA(display$Color),
  paste(
    "The pathway display must be BER, NER, MMR, FA, HR, AEJ, NHEJ",
    "with coefficients 1 through 7."
  )
)
expected_pathway_order <- setNames(expected_coefficients, expected_pathways)
required_catalog_columns <- c(
  "GO_ID",
  "CanonicalGOID",
  "IsAlternativeID",
  "GO_Name",
  "Ontology",
  "Definition",
  "Pathway",
  "PathwayOrder",
  "Coefficient",
  "MatchedSeedGOIDs",
  "ExpansionModes",
  "MatchTypes",
  "RuleLabels",
  "MappingRationale",
  "OntologyRelease",
  "OntologySHA256",
  "RuleSetSHA256"
)
assert(
  all(required_catalog_columns %in% names(catalog)) &&
    setequal(unique(catalog$Pathway), expected_pathways) &&
    all(catalog$PathwayOrder == expected_pathway_order[catalog$Pathway]) &&
    all(catalog$Coefficient == catalog$PathwayOrder) &&
    all(grepl("^GO:[0-9]{7}$", catalog$GO_ID)) &&
    all(grepl("^GO:[0-9]{7}$", catalog$CanonicalGOID)) &&
    !anyDuplicated(catalog[, .(GO_ID, Pathway)]),
  "The complete GO-term catalog violates the fixed seven-pathway contract."
)
assert(id_column %in% names(rna), paste("RNA table lacks", id_column))
assert(
  all(c(id_column, go_id_column) %in% names(go_annotations)),
  paste("GO table must contain", id_column, "and", go_id_column)
)
if (nzchar(effect_column)) {
  assert(
    effect_column %in% names(rna),
    paste("RNA table lacks effect column", effect_column)
  )
}
if (nzchar(padj_column)) {
  assert(
    nzchar(effect_column),
    "--padj-column requires --effect-column."
  )
  assert(
    padj_column %in% names(rna),
    paste("RNA table lacks adjusted-P column", padj_column)
  )
}

rna[, EnsemblID := normalize_ensembl(get(id_column))]
assert(
  all(nzchar(rna$EnsemblID)),
  "RNA table contains missing or empty Ensembl IDs."
)
if (group_column %in% names(rna)) {
  rna[, Group := trimws(as.character(get(group_column)))]
  rna[is.na(Group) | !nzchar(Group), Group := "all"]
} else {
  rna[, Group := "all"]
}
if (label_column %in% names(rna)) {
  rna[, DisplayLabel := as.character(get(label_column))]
  rna[is.na(DisplayLabel), DisplayLabel := ""]
} else {
  rna[, DisplayLabel := ""]
}
if (nzchar(effect_column)) {
  rna[, Effect := as.numeric(get(effect_column))]
  assert(
    all(is.finite(rna$Effect)),
    "The selected RNA effect column contains missing or non-finite values."
  )
} else {
  rna[, Effect := NA_real_]
}
if (nzchar(padj_column)) {
  rna[, Padj := as.numeric(get(padj_column))]
  assert(
    all(is.na(rna$Padj) | (is.finite(rna$Padj) & rna$Padj >= 0 & rna$Padj <= 1)),
    "The selected adjusted-P column contains values outside 0 to 1."
  )
} else {
  rna[, Padj := NA_real_]
}
assert(
  !anyDuplicated(rna[, .(Group, EnsemblID)]),
  paste(
    "RNA table has duplicated Group + EnsemblID rows.",
    "Resolve the biological replicate model before plotting."
  )
)
rna <- rna[
  ,
  .(Group, EnsemblID, DisplayLabel, Effect, Padj)
]

go_annotations[, EnsemblID := normalize_ensembl(get(id_column))]
go_annotations[, GO_ID := trimws(as.character(get(go_id_column)))]
assert(
  all(nzchar(go_annotations$EnsemblID)) &&
    all(grepl("^GO:[0-9]{7}$", go_annotations$GO_ID)),
  "GO table contains invalid Ensembl IDs or GO identifiers."
)
if (go_name_column %in% names(go_annotations)) {
  go_annotations[, GO_Name := as.character(get(go_name_column))]
  go_annotations[is.na(GO_Name), GO_Name := ""]
} else {
  go_annotations[, GO_Name := ""]
}
go_annotations <- unique(
  go_annotations[
    EnsemblID %in% rna$EnsemblID,
    .(EnsemblID, GO_ID, GO_Name)
  ]
)
assert(
  nrow(go_annotations) > 0L,
  "No direct GO annotations matched the RNA Ensembl IDs."
)

observed_terms <- unique(
  go_annotations[, .(GO_ID, InputGOName = GO_Name)]
)
term_decisions <- merge(
  observed_terms,
  catalog[
    ,
    .(
      GO_ID,
      CanonicalGOID,
      IsAlternativeID,
      CatalogGOName = GO_Name,
      Ontology,
      CatalogDefinition = Definition,
      Pathway,
      PathwayOrder,
      Coefficient,
      MatchedSeedGOIDs,
      ExpansionModes,
      MatchTypes,
      RuleLabels,
      MappingRationale,
      OntologyRelease,
      OntologySHA256,
      RuleSetSHA256
    )
  ],
  by = "GO_ID",
  all.x = TRUE,
  allow.cartesian = TRUE,
  sort = FALSE
)
term_decisions[
  is.na(Pathway),
  `:=`(
    CanonicalGOID = "",
    IsAlternativeID = FALSE,
    CatalogGOName = "",
    Ontology = "",
    CatalogDefinition = "",
    Pathway = "Others",
    PathwayOrder = 8L,
    Coefficient = 0L,
    MatchedSeedGOIDs = "",
    ExpansionModes = "",
    MatchTypes = "not_in_complete_seven_pathway_catalog",
    RuleLabels = "Others",
    MappingRationale = paste(
      "Observed direct GO term did not match the frozen complete",
      "seven-pathway ontology catalog."
    ),
    OntologyRelease = unique(catalog$OntologyRelease)[[1L]],
    OntologySHA256 = unique(catalog$OntologySHA256)[[1L]],
    RuleSetSHA256 = unique(catalog$RuleSetSHA256)[[1L]]
  )
]
setorder(term_decisions, PathwayOrder, GO_ID, CanonicalGOID)

gene_term_pathway <- merge(
  go_annotations,
  term_decisions[
    ,
    .(
      GO_ID,
      CanonicalGOID,
      IsAlternativeID,
      CatalogGOName,
      Ontology,
      Pathway,
      PathwayOrder,
      Coefficient,
      MatchedSeedGOIDs,
      ExpansionModes,
      MatchTypes,
      RuleLabels,
      MappingRationale
    )
  ],
  by = "GO_ID",
  all.x = TRUE,
  allow.cartesian = TRUE,
  sort = FALSE
)
gene_term_pathway <- unique(
  gene_term_pathway,
  by = c("EnsemblID", "GO_ID", "Pathway")
)

seven_pathway_terms <- gene_term_pathway[Pathway %in% expected_pathways]
seven_pathway_terms[, TermCountKey := fifelse(
  nzchar(CanonicalGOID),
  CanonicalGOID,
  GO_ID
)]
term_counts <- seven_pathway_terms[
  ,
  .(DirectTermCount = uniqueN(TermCountKey)),
  by = .(EnsemblID, Pathway)
]
term_count_wide <- dcast(
  term_counts,
  EnsemblID ~ Pathway,
  value.var = "DirectTermCount",
  fill = 0L
)
for (pathway in expected_pathways) {
  if (!pathway %in% names(term_count_wide)) {
    set(term_count_wide, j = pathway, value = 0L)
  }
}
setcolorder(term_count_wide, c("EnsemblID", expected_pathways))
term_count_wide[, SevenPathwayDirectTermCount :=
  rowSums(.SD), .SDcols = expected_pathways]
term_count_wide[, WeightedOrderScore :=
  as.numeric(as.matrix(.SD) %*% display$Coefficient),
  .SDcols = expected_pathways
]

excluded_genes <- rna[
  !EnsemblID %in% term_count_wide$EnsemblID,
  .(Group, EnsemblID, DisplayLabel, Reason = "No direct GO term assigned to the seven pathways")
]
gene_order <- merge(
  rna,
  term_count_wide,
  by = "EnsemblID",
  all = FALSE,
  sort = FALSE
)
assert(
  nrow(gene_order) > 0L,
  "No RNA genes had direct GO terms assigned to the seven pathways."
)
setorder(gene_order, Group, WeightedOrderScore, EnsemblID)
gene_order[, GeneRank := seq_len(.N), by = Group]
gene_order[, EligibleGeneCount := .N, by = Group]

gene_grid_base <- gene_order[
  ,
  .(
    Group,
    EnsemblID,
    DisplayLabel,
    Effect,
    Padj,
    WeightedOrderScore,
    SevenPathwayDirectTermCount,
    GeneRank,
    EligibleGeneCount
  )
]
gene_grid_base[, JoinKey := 1L]
display_grid <- copy(display)
display_grid[, JoinKey := 1L]
pathway_grid <- merge(
  gene_grid_base,
  display_grid,
  by = "JoinKey",
  allow.cartesian = TRUE,
  sort = FALSE
)
pathway_grid[, JoinKey := NULL]
matrix_data <- merge(
  pathway_grid,
  term_counts,
  by = c("EnsemblID", "Pathway"),
  all.x = TRUE,
  sort = FALSE
)
matrix_data[is.na(DirectTermCount), DirectTermCount := 0L]
matrix_data[, PathwayPresent := DirectTermCount > 0L]

has_effect <- nzchar(effect_column)
has_padj <- nzchar(padj_column)
fdr_caption <- if (has_padj) {
  paste0(
    "grey = padj > ",
    format(padj_cutoff, scientific = FALSE),
    " or missing"
  )
} else {
  "no padj filter; all present GO members shown"
}
if (has_effect) {
  matrix_data[, Significant :=
    if (has_padj) !is.na(Padj) & Padj <= padj_cutoff else TRUE]
  matrix_data[, RNAState := fifelse(
    !Significant | Effect == 0,
    "zero_or_nonsignificant",
    fifelse(Effect > 0, "positive", "negative")
  )]
  # Preserve the sign explicitly. Positive and negative effects are scaled
  # within separate visual channels instead of collapsing them to abs(effect):
  # filled saturation for positive effects and outline strength for negative
  # effects.
  matrix_data[, `:=`(
    PositiveEffect = fifelse(Effect > 0, Effect, 0),
    NegativeEffect = fifelse(Effect < 0, -Effect, 0)
  )]
  matrix_data[
    ,
    PositiveIntensityCap := {
      eligible <- PositiveEffect[Significant & PositiveEffect > 0]
      if (length(eligible)) {
        as.numeric(quantile(
          eligible,
          probs = 0.95,
          names = FALSE,
          type = 8
        ))
      } else {
        1
      }
    },
    by = Group
  ]
  matrix_data[
    ,
    NegativeIntensityCap := {
      eligible <- NegativeEffect[Significant & NegativeEffect > 0]
      if (length(eligible)) {
        as.numeric(quantile(
          eligible,
          probs = 0.95,
          names = FALSE,
          type = 8
        ))
      } else {
        1
      }
    },
    by = Group
  ]
  matrix_data[
    !is.finite(PositiveIntensityCap) | PositiveIntensityCap <= 0,
    PositiveIntensityCap := 1
  ]
  matrix_data[
    !is.finite(NegativeIntensityCap) | NegativeIntensityCap <= 0,
    NegativeIntensityCap := 1
  ]
  matrix_data[, NormalizedIntensity := fifelse(
    Effect > 0,
    pmin(PositiveEffect / PositiveIntensityCap, 1),
    fifelse(
      Effect < 0,
      pmin(NegativeEffect / NegativeIntensityCap, 1),
      0
    )
  )]
  matrix_data[, DisplayAlpha := 0.25 + 0.75 * NormalizedIntensity]
} else {
  matrix_data[, `:=`(
    Significant = NA,
    RNAState = "not_provided",
    PositiveEffect = NA_real_,
    NegativeEffect = NA_real_,
    PositiveIntensityCap = NA_real_,
    NegativeIntensityCap = NA_real_,
    NormalizedIntensity = NA_real_,
    DisplayAlpha = 1
  )]
}
setorder(matrix_data, Group, PathwayOrder, GeneRank)

group_sizes <- unique(gene_order[, .(Group, EligibleGeneCount)])
pathway_summary <- matrix_data[
  ,
  .(
    PresentGeneCount = sum(PathwayPresent),
    PositiveGeneCount = sum(PathwayPresent & RNAState == "positive"),
    NegativeGeneCount = sum(PathwayPresent & RNAState == "negative"),
    NeutralGeneCount = sum(
      PathwayPresent &
        RNAState %in% c("zero_or_nonsignificant", "not_provided")
    )
  ),
  by = .(Group, Pathway, PathwayOrder, Coefficient, Color)
]
pathway_summary <- merge(
  pathway_summary,
  group_sizes,
  by = "Group",
  all.x = TRUE,
  sort = FALSE
)
pathway_summary[, `:=`(
  PresentFraction = PresentGeneCount / EligibleGeneCount,
  PositiveFraction = PositiveGeneCount / EligibleGeneCount,
  NegativeFraction = NegativeGeneCount / EligibleGeneCount,
  NeutralFraction = NeutralGeneCount / EligibleGeneCount
)]
setorder(pathway_summary, Group, PathwayOrder)

write_output <- function(value, filename) {
  fwrite(value, file.path(outdir, filename))
}

write_output(term_decisions, "go_term_decisions.csv")
write_output(gene_term_pathway, "gene_go_pathway_long.csv")
write_output(gene_order, "gene_order.csv")
write_output(matrix_data, "gene_pathway_matrix.csv")
write_output(pathway_summary, "pathway_summary.csv")
write_output(excluded_genes, "excluded_genes_without_seven_pathway.csv")
write_output(catalog, "seven_pathway_go_catalog_used.csv")

run_metadata <- data.table(
  Key = c(
    "RNA input",
    "Direct GO input",
    "Complete GO-term catalog",
    "GO ontology release",
    "GO ontology SHA256",
    "GO rule-set SHA256",
    "Catalog canonical GO terms",
    "Catalog alternative GO IDs",
    "Pathway display",
    "RNA ID column",
    "GO ID column",
    "Group column",
    "Display-label column",
    "Effect column",
    "Adjusted-P column",
    "Adjusted-P cutoff",
    "Signed-effect encoding",
    "FDR encoding",
    "Pathway order",
    "Ordering formula",
    "RNA rows",
    "RNA unique genes",
    "RNA groups",
    "Matched direct GO pairs",
    "Observed direct GO terms",
    "Eligible seven-pathway genes",
    "Excluded RNA group-gene rows",
    "R version",
    "Random seed"
  ),
  Value = c(
    rna_path,
    go_path,
    catalog_path,
    paste(unique(catalog$OntologyRelease), collapse = ";"),
    paste(unique(catalog$OntologySHA256), collapse = ";"),
    paste(unique(catalog$RuleSetSHA256), collapse = ";"),
    uniqueN(catalog$CanonicalGOID),
    uniqueN(catalog[IsAlternativeID == TRUE, GO_ID]),
    display_path,
    id_column,
    go_id_column,
    if (group_column %in% rna_input_columns) group_column else "all",
    if (label_column %in% rna_input_columns) label_column else "",
    effect_column,
    padj_column,
    if (has_padj) format(padj_cutoff, scientific = FALSE) else "not used",
    paste(
      "positive log2FC = filled pathway colour;",
      "negative log2FC = white fill plus pathway-colour outline;",
      "positive and negative intensity are scaled separately to their",
      "within-group 95th percentiles"
    ),
    if (has_padj) {
      paste0(
        "adjusted P <= ", format(padj_cutoff, scientific = FALSE),
        " uses signed encoding; adjusted P above cutoff or missing is grey"
      )
    } else {
      "No adjusted-P column supplied; all present GO members use signed log2FC"
    },
    paste(expected_pathways, collapse = " -> "),
    paste0(
      "BER_count*1 + NER_count*2 + MMR_count*3 + FA_count*4 + ",
      "HR_count*5 + AEJ_count*6 + NHEJ_count*7"
    ),
    nrow(rna),
    uniqueN(rna$EnsemblID),
    uniqueN(rna$Group),
    nrow(go_annotations),
    uniqueN(go_annotations$GO_ID),
    uniqueN(gene_order$EnsemblID),
    nrow(excluded_genes),
    as.character(getRversion()),
    25L
  )
)
write_output(run_metadata, "run_metadata.csv")

matrix_breaks <- function(n_genes) {
  if (n_genes <= 10L) return(seq_len(n_genes))
  unique(as.integer(round(seq(1, n_genes, length.out = 5L))))
}

save_plot <- function(plot, filename_stem, width, height) {
  for (format in formats) {
    output_path <- file.path(outdir, paste0(filename_stem, ".", format))
    if (format == "png") {
      ggsave(
        output_path,
        plot,
        width = width,
        height = height,
        units = "in",
        dpi = 300,
        device = ragg::agg_png,
        bg = "white"
      )
    } else if (format == "pdf") {
      ggsave(
        output_path,
        plot,
        width = width,
        height = height,
        units = "in",
        device = grDevices::cairo_pdf,
        bg = "white"
      )
    } else {
      ggsave(
        output_path,
        plot,
        width = width,
        height = height,
        units = "in",
        device = grDevices::svg,
        bg = "white"
      )
    }
  }
}

make_matrix_plot <- function(group_name) {
  panel <- copy(matrix_data[Group == group_name])
  n_genes <- unique(panel$EligibleGeneCount)
  assert(length(n_genes) == 1L, paste("Invalid gene count for group", group_name))
  panel[, Pathway := factor(Pathway, levels = rev(expected_pathways))]
  background <- panel
  present <- panel[PathwayPresent == TRUE]

  plot <- ggplot() +
    geom_tile(
      data = background,
      aes(x = GeneRank, y = Pathway),
      fill = "#F0F2F4",
      colour = NA,
      width = 1,
      height = 0.84
    )
  if (!has_effect) {
    plot <- plot +
      geom_tile(
        data = present,
        aes(x = GeneRank, y = Pathway, fill = Color),
        colour = NA,
        width = 1,
        height = 0.84
      )
  } else {
    plot <- plot +
      geom_tile(
        data = present[RNAState == "positive"],
        aes(
          x = GeneRank,
          y = Pathway,
          fill = Color,
          alpha = DisplayAlpha
        ),
        colour = NA,
        width = 1,
        height = 0.84
      ) +
      geom_tile(
        data = present[RNAState == "negative"],
        aes(
          x = GeneRank,
          y = Pathway,
          colour = Color,
          alpha = DisplayAlpha
        ),
        fill = "white",
        linewidth = 0.45,
        width = 0.96,
        height = 0.80
      ) +
      geom_tile(
        data = present[RNAState == "zero_or_nonsignificant"],
        aes(x = GeneRank, y = Pathway),
        fill = "#D1D5DB",
        colour = NA,
        width = 1,
        height = 0.84
      )
  }

  plot +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_alpha_identity() +
    scale_x_continuous(
      limits = c(0.5, n_genes + 0.5),
      breaks = matrix_breaks(n_genes),
      expand = c(0, 0)
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = group_name,
      x = "Gene rank (ascending GO-weighted score; RNA intensity does not affect order)",
      y = NULL,
      caption = paste0(
        "Signed log2FC: filled = up, outlined = down; ",
        fdr_caption,
        ". Colour strength is scaled separately within direction."
      )
    ) +
    theme_minimal(base_family = font_family, base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        colour = "#D1D5DB",
        linewidth = 0.32
      ),
      axis.text.y = element_text(
        face = "bold",
        colour = "#374151",
        size = 10
      ),
      axis.text.x = element_text(colour = "#4B5563", size = 10),
      axis.title.x = element_text(colour = "#374151", size = 9.5),
      plot.title = element_text(
        colour = "#111827",
        face = "bold",
        size = 12
      ),
      plot.caption = element_text(
        colour = "#4B5563",
        size = 8.3,
        hjust = 0
      ),
      plot.margin = margin(7, 12, 7, 10),
      legend.position = "none"
    )
}

make_summary_plot <- function(group_name) {
  panel <- copy(pathway_summary[Group == group_name])
  panel[, Pathway := factor(Pathway, levels = rev(expected_pathways))]
  panel[, PresentLabel := sprintf(
    "%d (%.1f%%)",
    PresentGeneCount,
    100 * PresentFraction
  )]

  ggplot(panel, aes(x = PresentFraction, y = Pathway, fill = Color)) +
    geom_col(width = 0.58) +
    geom_text(
      aes(x = PresentFraction + 0.01, label = PresentLabel),
      family = font_family,
      hjust = 0,
      size = 3.2,
      colour = "#374151"
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      limits = c(0, max(0.12, max(panel$PresentFraction) + 0.18)),
      labels = function(value) paste0(round(100 * value), "%"),
      expand = c(0, 0)
    ) +
    labs(x = "Genes with at least one direct pathway GO term", y = NULL) +
    theme_minimal(base_family = font_family, base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        colour = "#E5E7EB",
        linewidth = 0.32
      ),
      axis.text.y = element_text(
        face = "bold",
        colour = "#374151",
        size = 10
      ),
      axis.text.x = element_text(colour = "#4B5563", size = 10),
      axis.title.x = element_text(colour = "#374151", size = 9.5),
      plot.margin = margin(8, 18, 8, 8),
      legend.position = "none"
    )
}

figure_manifest <- rbindlist(lapply(unique(gene_order$Group), function(group_name) {
  safe_group <- safe_filename(group_name)
  n_genes <- gene_order[Group == group_name, uniqueN(EnsemblID)]
  matrix_width <- min(24, max(8, 4 + n_genes / 24))
  matrix_height <- 4.6
  summary_width <- 7.2
  summary_height <- 4.2

  matrix_stem <- paste(file_prefix, "matrix", safe_group, sep = "_")
  summary_stem <- paste(file_prefix, "summary", safe_group, sep = "_")
  save_plot(
    make_matrix_plot(group_name),
    matrix_stem,
    matrix_width,
    matrix_height
  )
  save_plot(
    make_summary_plot(group_name),
    summary_stem,
    summary_width,
    summary_height
  )
  rbindlist(list(
    data.table(
      Group = group_name,
      FigureType = "linear_matrix",
      Format = formats,
      File = paste0(matrix_stem, ".", formats)
    ),
    data.table(
      Group = group_name,
      FigureType = "pathway_summary",
      Format = formats,
      File = paste0(summary_stem, ".", formats)
    )
  ))
}))
write_output(figure_manifest, "figure_manifest.csv")

message(
  "Created GO pathway linear matrices for ",
  uniqueN(gene_order$Group),
  " group(s), ",
  uniqueN(gene_order$EnsemblID),
  " eligible Ensembl genes, and ",
  uniqueN(go_annotations$GO_ID),
  " observed direct GO terms in: ",
  outdir
)
