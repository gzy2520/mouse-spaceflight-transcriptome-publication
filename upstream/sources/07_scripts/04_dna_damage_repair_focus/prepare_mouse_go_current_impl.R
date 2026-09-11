#!/usr/bin/env Rscript

# Build mouse GO gene sets directly from the MGI MOD-GAF.
#
# Primary annotation key: MGI ID from the frozen MOUSE-mod.gaf.gz snapshot.
# Analysis key: Ensembl Gene ID, release 116 (GRCm39).
#
# Mapping routes:
#   1. Ensembl 116 BioMart MGI cross-reference
#   2. current MGI MRK_ENSEMBL report, restricted to Ensembl 116 universe
#   3. current MGI Entrez report -> Ensembl 116 Entrez cross-reference
#
# Symbol is display-only and is never used to choose an Ensembl candidate.

suppressPackageStartupMessages(library(data.table))

options(stringsAsFactors = FALSE, datatable.allow.cartesian = TRUE)

root <- normalizePath(".", mustWork = TRUE)
env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
annotation_root <- env_or(
  "MOUSE_GO_ANNOTATION_ROOT",
  file.path(root, "08_GO_annotation/mgi_mod_gaf_five_relation_20260822")
)
source_dir <- env_or(
  "MOUSE_GO_SOURCE_DIR",
  file.path(annotation_root, "source")
)
go_relation_scope <- trimws(strsplit(
  env_or(
    "MOUSE_GO_RELATIONS",
    paste(
      c(
        "is_a", "part_of", "regulates",
        "positively_regulates", "negatively_regulates"
      ),
      collapse = ","
    )
  ),
  ",",
  fixed = TRUE
)[[1L]])
valid_go_relations <- c(
  "is_a", "part_of", "regulates",
  "positively_regulates", "negatively_regulates"
)
if (!length(go_relation_scope) ||
    any(!go_relation_scope %chin% valid_go_relations)) {
  stop(
    "MOUSE_GO_RELATIONS must contain only: ",
    paste(valid_go_relations, collapse = ",")
  )
}
output_dir <- file.path(
  annotation_root, "strict_consensus"
)
expanded_dir <- file.path(
  annotation_root, "expanded_unique"
)
final_dir <- file.path(
  annotation_root, "final_manual_curated"
)
manual_decision_file <- env_or(
  "MOUSE_GO_MANUAL_DECISION_FILE",
  file.path(
    annotation_root,
    "manual_review",
    "00_five_relation_multi_route_manual_decisions.csv"
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(expanded_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

gaf_file <- file.path(source_dir, "MOUSE-mod.gaf.gz")
obo_file <- file.path(source_dir, "go-basic.obo")
mgi_ensembl_file <- file.path(source_dir, "MRK_ENSEMBL.rpt")
mgi_entrez_file <- file.path(source_dir, "MGI_EntrezGene.rpt")
biomart_full_file <- file.path(
  source_dir,
  "ensembl116_full_mgi_entrez_mapping.tsv.gz"
)

required_files <- c(
  gaf_file, obo_file, mgi_ensembl_file, mgi_entrez_file,
  biomart_full_file
)
if (!all(file.exists(required_files))) {
  stop(
    "Missing required MOD-GAF source files: ",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
}
if (!file.exists(manual_decision_file)) {
  stop("Missing manual mapping decision file: ", manual_decision_file)
}

term_config <- data.table(
  requested_term_key = c(
    "DNA_damage_response_GO0006974",
    "DNA_damage_signal_transduction_GO0042770",
    "DNA_repair_GO0006281",
    "DNA_damage_tolerance_GO0006301",
    "intrinsic_apoptotic_signaling_GO0008630",
    "cell_cycle_GO0007049",
    "translation_GO0006412",
    "DNA_templated_transcription_GO0006351",
    "DNA_replication_GO0006260",
    "mechanosensation_GO0050954",
    "cytoskeleton_GO0005856",
    "mitochondrion_GO0005739",
    "telomere_region_GO0000781",
    "telomere_maintenance_GO0043247",
    "cell_death_GO0008219"
  ),
  requested_go_id = c(
    "GO:0006974", "GO:0042770", "GO:0006281", "GO:0006301",
    "GO:0008630", "GO:0007049", "GO:0006412", "GO:0006351",
    "GO:0006260", "GO:0050954", "GO:0005856", "GO:0005739",
    "GO:0000781", "GO:0043247", "GO:0008219"
  ),
  requested_term_name = c(
    "DNA damage response", "Signal transduction in response to DNA damage",
    "DNA repair", "DNA damage tolerance", "Intrinsic apoptotic signaling pathway",
    "Cell cycle", "Translation", "DNA-templated transcription", "DNA replication",
    "Sensory perception of mechanical stimulus", "Cytoskeleton",
    "Mitochondrion", "Chromosome, telomeric region",
    "Telomere maintenance in response to DNA damage", "Cell death"
  ),
  ontology_namespace = c(
    rep("biological_process", 10L),
    rep("cellular_component", 3L),
    rep("biological_process", 2L)
  ),
  go_aspect = c(rep("P", 10L), rep("C", 3L), rep("P", 2L))
)

# Compatibility table consumed by the extended cross-tissue runner.  The
# authoritative columns above remain requested_* so that membership records
# retain their explicit root-term provenance.
official_term_config <- term_config[, .(
  term_key = requested_term_key,
  go_id = requested_go_id,
  term_name = requested_term_name,
  ontology_namespace,
  go_aspect,
  official_plot_order = seq_len(.N)
)]

sha256_file <- function(path) {
  unname(system2(
    "shasum",
    c("-a", "256", shQuote(path)),
    stdout = TRUE
  )) |>
    sub(paste0("  ", path, "$"), "", x = _)
}

collapse_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  paste(x, collapse = ";")
}

parse_go_obo <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  term_starts <- which(lines == "[Term]")
  term_ends <- c(term_starts[-1L] - 1L, length(lines))
  terms <- vector("list", length(term_starts))
  alt_rows <- list()

  for (i in seq_along(term_starts)) {
    block <- lines[term_starts[i]:term_ends[i]]
    id_line <- grep("^id: GO:", block, value = TRUE)
    if (!length(id_line)) next
    go_id <- sub("^id: ", "", id_line[1L])
    name_line <- grep("^name: ", block, value = TRUE)
    namespace_line <- grep("^namespace: ", block, value = TRUE)
    obsolete <- any(block == "is_obsolete: true")
    is_a <- if ("is_a" %chin% go_relation_scope) {
      sub(
        "^is_a: (GO:[0-9]+).*$", "\\1",
        grep("^is_a: GO:", block, value = TRUE)
      )
    } else {
      character()
    }
    relationship_scope <- setdiff(go_relation_scope, "is_a")
    relationship_pattern <- if (length(relationship_scope)) {
      paste(relationship_scope, collapse = "|")
    } else {
      ""
    }
    relation_lines <- if (nzchar(relationship_pattern)) {
      grep(
        paste0("^relationship: (", relationship_pattern, ") GO:"),
        block,
        value = TRUE
      )
    } else {
      character()
    }
    relation_parents <- if (length(relation_lines)) {
      sub(
        paste0(
          "^relationship: (", relationship_pattern,
          ") (GO:[0-9]+).*$"
        ),
        "\\2",
        relation_lines
      )
    } else {
      character()
    }
    alt_ids <- sub(
      "^alt_id: ", "",
      grep("^alt_id: GO:", block, value = TRUE)
    )
    terms[[i]] <- data.table(
      go_id = go_id,
      go_name = if (length(name_line)) {
        sub("^name: ", "", name_line[1L])
      } else {
        NA_character_
      },
      namespace = if (length(namespace_line)) {
        sub("^namespace: ", "", namespace_line[1L])
      } else {
        NA_character_
      },
      is_obsolete = obsolete,
      parents = list(unique(c(is_a, relation_parents)))
    )
    if (length(alt_ids)) {
      alt_rows[[length(alt_rows) + 1L]] <- data.table(
        alt_go_id = alt_ids,
        go_id = go_id
      )
    }
  }
  term_table <- rbindlist(terms, fill = TRUE)
  term_table <- term_table[
    !is.na(go_id) & !is_obsolete &
      namespace %chin% c("biological_process", "cellular_component")
  ]
  alt_table <- if (length(alt_rows)) {
    rbindlist(alt_rows)
  } else {
    data.table(alt_go_id = character(), go_id = character())
  }
  list(terms = term_table, alt = alt_table)
}

descendants_of <- function(term_table, root_go_id) {
  edges <- term_table[
    lengths(parents) > 0L,
    .(child = go_id, parent = unlist(parents)),
    by = go_id
  ][, .(child, parent)]
  children_by_parent <- split(edges$child, edges$parent)
  seen <- root_go_id
  frontier <- root_go_id
  while (length(frontier)) {
    children <- unique(unlist(children_by_parent[frontier], use.names = FALSE))
    children <- setdiff(children, seen)
    seen <- unique(c(seen, children))
    frontier <- children
  }
  seen
}

read_gaf <- function(path) {
  gaf_columns <- c(
    "DB", "DB_Object_ID", "SYMBOL", "QUALIFIER", "GO_TERM",
    "DB_REFERENCE", "EVIDENCE_CODE", "WITH_FROM", "GO_ASPECT",
    "DB_OBJECT_NAME", "DB_OBJECT_SYNONYM", "DB_OBJECT_TYPE", "TAXON",
    "DATE", "ASSIGNED_BY", "ANNOTATION_EXTENSION",
    "GENE_PRODUCT_FORM_ID"
  )
  command <- paste(
    "gzip -dc", shQuote(path),
    "| awk 'BEGIN{FS=\"\\t\"} $0 !~ /^!/'"
  )
  gaf <- fread(
    cmd = command,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    col.names = gaf_columns,
    showProgress = FALSE
  )
  gaf[, source_row_number := seq_len(.N)]
  gaf
}

split_ids <- function(x, pattern = "[ |]+") {
  values <- unlist(strsplit(
    as.character(x[!is.na(x) & nzchar(as.character(x))]),
    pattern
  ))
  values[nzchar(values)]
}

build_mapping_table <- function(mgi_ids) {
  biomart_full <- fread(biomart_full_file)
  biomart_full <- biomart_full[
    grepl("^ENSMUSG[0-9]+$", ensembl_gene_id)
  ]
  biomart_mgi <- unique(biomart_full[
    mgi_id %chin% mgi_ids &
      grepl("^ENSMUSG[0-9]+$", ensembl_gene_id),
    .(
      mgi_id,
      ensembl_gene_id,
      entrezgene_id = as.character(entrezgene_id),
      external_gene_name
    )
  ])
  biomart_all <- biomart_full[
    ,
    .(ensembl_gene_id, entrezgene_id, external_gene_name)
  ]
  ensembl116_universe <- unique(biomart_full$ensembl_gene_id)

  mgi_ensembl <- fread(
    mgi_ensembl_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    select = c(1L, 2L, 6L),
    col.names = c("mgi_id", "mgi_symbol", "ensembl_raw")
  )
  mgi_ensembl <- mgi_ensembl[mgi_id %chin% mgi_ids]
  mgi_ensembl <- mgi_ensembl[
    ,
    .(ensembl_gene_id = split_ids(ensembl_raw)),
    by = .(mgi_id, mgi_symbol)
  ][
    ensembl_gene_id %chin% ensembl116_universe &
      grepl("^ENSMUSG[0-9]+$", ensembl_gene_id)
  ]

  mgi_entrez <- fread(
    mgi_entrez_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    select = c(1L, 2L, 3L, 9L),
    col.names = c("mgi_id", "mgi_symbol", "marker_status", "entrezgene_id")
  )
  mgi_entrez <- unique(mgi_entrez[
    mgi_id %chin% mgi_ids &
      marker_status == "O" &
      !is.na(entrezgene_id) &
      nzchar(entrezgene_id),
    .(
      mgi_id,
      mgi_symbol,
      entrezgene_id = as.character(entrezgene_id)
    )
  ])
  entrez_route <- merge(
    mgi_entrez,
    biomart_all[
      !is.na(entrezgene_id),
      .(
        entrezgene_id = as.character(entrezgene_id),
        ensembl_gene_id,
        external_gene_name
      )
    ],
    by = "entrezgene_id",
    allow.cartesian = TRUE
  )

  route_rows <- rbindlist(list(
    biomart_mgi[
      ,
      .(
        mgi_id,
        ensembl_gene_id,
        route = "ensembl116_biomart_mgi"
      )
    ],
    mgi_ensembl[
      ,
      .(
        mgi_id,
        ensembl_gene_id,
        route = "mgi_mrk_ensembl_restricted_to_release116"
      )
    ],
    entrez_route[
      ,
      .(
        mgi_id,
        ensembl_gene_id,
        route = "mgi_entrez_to_ensembl116"
      )
    ]
  ))
  route_rows <- unique(route_rows)

  mapped <- rbindlist(lapply(mgi_ids, function(mgi_id_value) {
    routes <- route_rows[mgi_id == mgi_id_value]
    route_sets <- split(routes$ensembl_gene_id, routes$route)
    route_sets <- lapply(route_sets, unique)
    nonempty_routes <- names(route_sets)[lengths(route_sets) > 0L]
    route_sets <- route_sets[nonempty_routes]
    candidates_union <- unique(unlist(route_sets, use.names = FALSE))
    candidates_intersection <- if (length(route_sets)) {
      Reduce(intersect, route_sets)
    } else {
      character()
    }
    n_routes <- length(route_sets)
    route_support <- if (length(route_sets)) {
      sort(table(unlist(route_sets, use.names = FALSE)), decreasing = TRUE)
    } else {
      integer()
    }
    max_route_support <- if (length(route_support)) {
      as.integer(max(route_support))
    } else {
      0L
    }
    top_supported_candidates <- if (length(route_support)) {
      names(route_support)[route_support == max_route_support]
    } else {
      character()
    }

    accepted_strict <- character()
    accepted_expanded <- character()
    decision <- "unmapped_mgi_id"
    if (
      max_route_support >= 2L &&
        length(top_supported_candidates) == 1L
    ) {
      accepted_strict <- top_supported_candidates
      accepted_expanded <- accepted_strict
      decision <- if (max_route_support == n_routes) {
        "accepted_all_nonempty_routes_consensus"
      } else {
        "accepted_unique_majority_route_consensus"
      }
    } else if (n_routes == 1L && length(candidates_union) == 1L) {
      accepted_expanded <- candidates_union
      decision <- "single_route_unique_expanded_only"
    } else if (length(candidates_union) > 1L) {
      decision <- "ambiguous_without_symbol_resolution"
    }

    data.table(
      mgi_id = mgi_id_value,
      n_mapping_routes = n_routes,
      mapping_routes = paste(sort(nonempty_routes), collapse = ";"),
      candidate_ensembl_gene_count = length(candidates_union),
      candidate_ensembl_gene_ids = paste(
        sort(candidates_union), collapse = ";"
      ),
      consensus_candidate_count = length(candidates_intersection),
      max_route_support = max_route_support,
      top_supported_candidates = paste(
        sort(top_supported_candidates), collapse = ";"
      ),
      strict_ensembl_gene_id = if (length(accepted_strict)) {
        accepted_strict
      } else {
        NA_character_
      },
      expanded_ensembl_gene_id = if (length(accepted_expanded)) {
        accepted_expanded
      } else {
        NA_character_
      },
      mapping_decision = decision
    )
  }), fill = TRUE)

  symbol_lookup <- unique(rbindlist(list(
    biomart_mgi[, .(mgi_id, external_gene_name)],
    mgi_ensembl[, .(mgi_id, external_gene_name = mgi_symbol)],
    mgi_entrez[, .(mgi_id, external_gene_name = mgi_symbol)]
  ))[!is.na(external_gene_name) & nzchar(external_gene_name)])
  symbol_lookup <- symbol_lookup[
    ,
    .(external_gene_name = collapse_values(external_gene_name)),
    by = mgi_id
  ]
  entrez_lookup <- mgi_entrez[
    ,
    .(entrezgene_id = collapse_values(entrezgene_id)),
    by = mgi_id
  ]
  mapped <- merge(mapped, symbol_lookup, by = "mgi_id", all.x = TRUE)
  mapped <- merge(mapped, entrez_lookup, by = "mgi_id", all.x = TRUE)
  mapped
}

write_variant <- function(
  annotation_membership,
  mapping_table,
  variant_dir,
  mapping_id_column,
  variant_name,
  mapping_decision_column = "mapping_decision"
) {
  mapped <- merge(
    annotation_membership,
    mapping_table,
    by = "mgi_id",
    all.x = TRUE,
    sort = FALSE
  )
  mapped[, ensembl_gene_id := get(mapping_id_column)]
  mapped[, mapping_decision := get(mapping_decision_column)]
  mapped[, `:=`(
    gene_product_db = "MGI",
    gene_product_id = mgi_id,
    uniprot_base = NA_character_,
    included_in_gene_set =
      !is.na(ensembl_gene_id) & nzchar(ensembl_gene_id),
    analysis_id_type = "Ensembl Gene ID",
    symbol_used_as_analysis_key = FALSE,
    mapping_variant = variant_name
  )]
  mapped[
    included_in_gene_set == FALSE,
    mapping_decision := fifelse(
      mapping_decision == "single_route_unique_expanded_only" &
        variant_name == "strict_consensus",
      "single_route_unique_excluded_from_strict",
      mapping_decision
    )
  ]

  requested_order <- term_config$requested_term_key
  mapped[, term_order := match(requested_term_key, requested_order)]
  setorder(mapped, term_order, source_row_number, ensembl_gene_id)
  mapped[, term_order := NULL]

  fwrite(
    mapped,
    file.path(variant_dir, "03_mouse_GO_annotations_with_gene_ids.tsv.gz"),
    sep = "\t"
  )

  audit <- mapped[
    ,
    .(
      go_id = unique(requested_go_id),
      term_name = unique(requested_term_name),
      source_file = "08_GO_annotation/mgi_mod_gaf_five_relation_20260822/source/MOUSE-mod.gaf.gz",
      source_sha256 = sha256_file(gaf_file),
      source_rows_all_taxa = NA_integer_,
      source_rows_taxon10090_positive_in_scope_aspect = .N,
      source_rows_uniprotkb = 0L,
      source_rows_other_gene_product_db = .N,
      unique_uniprot_accessions = 0L,
      mapped_uniprot_accessions = 0L,
      unique_mgi_ids = uniqueN(mgi_id),
      mapped_mgi_ids = uniqueN(
        mgi_id[included_in_gene_set == TRUE]
      ),
      unique_ensembl_gene_ids = uniqueN(
        ensembl_gene_id[included_in_gene_set == TRUE]
      ),
      unique_entrez_gene_ids = uniqueN(
        entrezgene_id[
          included_in_gene_set == TRUE &
            !is.na(entrezgene_id) & nzchar(entrezgene_id)
        ]
      ),
      excluded_source_rows = uniqueN(
        source_row_number[included_in_gene_set == FALSE]
      ),
      excluded_mapping_candidate_rows = sum(
        included_in_gene_set == FALSE
      ),
      analysis_id_type = "Ensembl Gene ID",
      symbol_used_as_analysis_key = FALSE,
      ensembl_release = 116L,
      ensembl_assembly = "GRCm39",
      mapping_source = paste(
        "MGI MOUSE-mod GAF; Ensembl116 BioMart MGI;",
        "MGI MRK_ENSEMBL; MGI Entrez to Ensembl116"
      ),
      mapping_variant = variant_name
    ),
    by = .(term_key = requested_term_key)
  ]
  fwrite(
    audit,
    file.path(variant_dir, "04_mouse_GO_ID_mapping_audit.csv"),
    bom = TRUE
  )

  for (key in term_config$requested_term_key) {
    fwrite(
      mapped[requested_term_key == key],
      file.path(
        variant_dir,
        paste0(key, "_with_gene_ids.tsv.gz")
      ),
      sep = "\t"
    )
  }
  audit
}

obo <- parse_go_obo(obo_file)
root_validation <- merge(
  term_config,
  obo$terms[, .(requested_go_id = go_id, obo_namespace = namespace)],
  by = "requested_go_id",
  all.x = TRUE,
  sort = FALSE
)
if (root_validation[, any(is.na(obo_namespace))]) {
  stop(
    "Requested GO root missing from the frozen OBO: ",
    paste(
      root_validation[is.na(obo_namespace), requested_go_id],
      collapse = ", "
    )
  )
}
if (root_validation[, any(ontology_namespace != obo_namespace)]) {
  stop(
    "Requested GO root has an unexpected namespace: ",
    paste(
      root_validation[
        ontology_namespace != obo_namespace,
        paste0(requested_go_id, " (expected ", ontology_namespace,
               ", OBO ", obo_namespace, ")")
      ],
      collapse = "; "
    )
  )
}
fwrite(
  official_term_config,
  file.path(annotation_root, "00_extended_GO_term_config.csv"),
  bom = TRUE
)
descendant_sets <- setNames(
  lapply(term_config$requested_go_id, function(go_id) {
    descendants_of(obo$terms, go_id)
  }),
  term_config$requested_go_id
)

gaf <- read_gaf(gaf_file)
if (ncol(gaf) != 18L) {
  stop("Unexpected GAF column count")
}
gaf[, normalized_go_id := GO_TERM]
if (nrow(obo$alt)) {
  gaf[
    obo$alt,
    on = .(normalized_go_id = alt_go_id),
    normalized_go_id := i.go_id
  ]
}
eligible <- gaf[
  DB == "MGI" &
    grepl("(^|\\|)(taxon:|NCBITaxon:)10090($|\\|)", TAXON) &
    !grepl("(^|\\|)NOT($|\\|)", QUALIFIER)
]

annotation_membership <- rbindlist(lapply(seq_len(nrow(term_config)), function(i) {
  cfg <- term_config[i]
  members <- eligible[
    GO_ASPECT == cfg$go_aspect &
      normalized_go_id %chin% descendant_sets[[cfg$requested_go_id]]
  ]
  members[
    ,
    .(
      mgi_id = DB_Object_ID,
      SYMBOL,
      QUALIFIER,
      `GO TERM` = GO_TERM,
      normalized_go_id,
      GO_NAME = obo$terms$go_name[
        match(normalized_go_id, obo$terms$go_id)
      ],
      EVIDENCE_CODE,
      DB_REFERENCE,
      WITH_FROM,
      TAXON,
      ASSIGNED_BY,
      ANNOTATION_EXTENSION,
      DB_OBJECT_TYPE,
      source_row_number,
      requested_term_key = cfg$requested_term_key,
      requested_go_id = cfg$requested_go_id,
      requested_term_name = cfg$requested_term_name,
      ontology_namespace = cfg$ontology_namespace,
      requested_go_aspect = cfg$go_aspect,
      eligible_non_NOT_mouse_annotation = TRUE
    )
  ]
}), fill = TRUE)

if (!nrow(annotation_membership)) {
  stop("No target-term MGI annotations were found")
}
message(
  "Target annotation rows: ", nrow(annotation_membership),
  "; unique MGI IDs: ", uniqueN(annotation_membership$mgi_id)
)
message(
  paste(
    annotation_membership[
      ,
      .(rows = .N, mgi_ids = uniqueN(mgi_id)),
      by = requested_term_key
    ][
      ,
      paste0(requested_term_key, "=", rows, "/", mgi_ids)
    ],
    collapse = "; "
  )
)

mapping_table <- build_mapping_table(unique(annotation_membership$mgi_id))

manual_decisions <- fread(manual_decision_file)
required_manual_columns <- c(
  "mgi_id", "source_symbol", "manual_decision",
  "manual_ensembl_gene_id", "replacement_mgi_id",
  "evidence", "rationale", "reviewed_on"
)
if (!all(required_manual_columns %chin% names(manual_decisions))) {
  stop("Manual decision file has missing columns")
}
if (manual_decisions[, anyDuplicated(mgi_id)] > 0L) {
  stop("Manual decision file has duplicated MGI IDs")
}
if (!all(
  manual_decisions$manual_decision %chin%
    c("include", "exclude_no_ensembl116_id")
)) {
  stop("Unexpected manual_decision value")
}

ensembl116_universe <- unique(
  fread(
    biomart_full_file,
    select = "ensembl_gene_id"
  )$ensembl_gene_id
)
manual_include_ids <- manual_decisions[
  manual_decision == "include",
  manual_ensembl_gene_id
]
if (!all(manual_include_ids %chin% ensembl116_universe)) {
  stop(
    "Manual included IDs absent from Ensembl 116: ",
    paste(
      setdiff(manual_include_ids, ensembl116_universe),
      collapse = ", "
    )
  )
}

strict_missing_ids <- mapping_table[
  is.na(strict_ensembl_gene_id) | !nzchar(strict_ensembl_gene_id),
  mgi_id
]
if (length(setdiff(strict_missing_ids, manual_decisions$mgi_id))) {
  stop(
    "Manual decision file is missing strict mapping exceptions. ",
    "Only-in-automatic: ",
    paste(setdiff(strict_missing_ids, manual_decisions$mgi_id), collapse = ";")
  )
}
manual_decisions <- manual_decisions[mgi_id %chin% strict_missing_ids]

mapping_table <- merge(
  mapping_table,
  manual_decisions,
  by = "mgi_id",
  all.x = TRUE,
  sort = FALSE
)
mapping_table[, final_ensembl_gene_id := strict_ensembl_gene_id]
mapping_table[
  manual_decision == "include",
  final_ensembl_gene_id := manual_ensembl_gene_id
]
mapping_table[, final_mapping_decision := mapping_decision]
mapping_table[
  manual_decision == "include",
  final_mapping_decision := "accepted_manual_official_id_chain"
]
mapping_table[
  manual_decision == "exclude_no_ensembl116_id",
  final_mapping_decision := "excluded_manual_no_ensembl116_id"
]

fwrite(
  mapping_table,
  file.path(annotation_root, "01_MGI_to_Ensembl116_multi_route_mapping_audit.csv"),
  bom = TRUE
)

manual_candidate_summary <- annotation_membership[
  mgi_id %chin% strict_missing_ids,
  .(
    gaf_symbols = collapse_values(SYMBOL),
    object_types = collapse_values(DB_OBJECT_TYPE),
    terms = collapse_values(requested_term_name),
    go_ids = collapse_values(`GO TERM`),
    evidence_codes = collapse_values(EVIDENCE_CODE),
    assigned_by = collapse_values(ASSIGNED_BY),
    n_annotation_rows = .N
  ),
  by = mgi_id
]
manual_candidate_summary <- merge(
  manual_candidate_summary,
  mapping_table,
  by = "mgi_id",
  all.x = TRUE,
  sort = FALSE
)
fwrite(
  manual_candidate_summary,
  file.path(annotation_root, "03_manual_review_candidates_unmapped_MGI.csv"),
  bom = TRUE
)

strict_audit <- write_variant(
  annotation_membership,
  mapping_table,
  output_dir,
  "strict_ensembl_gene_id",
  "strict_consensus"
)
expanded_audit <- write_variant(
  annotation_membership,
  mapping_table,
  expanded_dir,
  "expanded_ensembl_gene_id",
  "expanded_unique"
)
final_audit <- write_variant(
  annotation_membership,
  mapping_table,
  final_dir,
  "final_ensembl_gene_id",
  "final_manual_curated",
  "final_mapping_decision"
)

manifest_files <- c(required_files, manual_decision_file)
source_manifest <- data.table(
  file = manifest_files,
  sha256 = vapply(manifest_files, sha256_file, character(1)),
  downloaded_or_frozen_on = "2026-08-22",
  role = c(
    "MGI MOUSE-mod GAF snapshot retrieved 2026-08-22",
    paste0(
      "GO ontology used for descendant expansion: ",
      paste(go_relation_scope, collapse = ",")
    ),
    "MGI to Ensembl current cross-reference",
    "MGI to Entrez current cross-reference",
    "Full Ensembl release 116 BioMart MGI/Entrez mapping",
    "Manual official-ID review decisions for strict mapping exceptions"
  )
)
fwrite(
  source_manifest,
  file.path(annotation_root, "00_source_manifest.csv"),
  bom = TRUE
)

comparison <- merge(
  strict_audit[
    ,
    .(
      term_key,
      strict_mapped_mgi_ids = mapped_mgi_ids,
      strict_unique_ensembl_gene_ids = unique_ensembl_gene_ids
    )
  ],
  expanded_audit[
    ,
    .(
      term_key,
      expanded_mapped_mgi_ids = mapped_mgi_ids,
      expanded_unique_ensembl_gene_ids = unique_ensembl_gene_ids
    )
  ],
  by = "term_key"
)
comparison <- merge(
  comparison,
  final_audit[
    ,
    .(
      term_key,
      final_mapped_mgi_ids = mapped_mgi_ids,
      final_unique_ensembl_gene_ids = unique_ensembl_gene_ids
    )
  ],
  by = "term_key"
)
fwrite(
  comparison,
  file.path(annotation_root, "02_strict_vs_expanded_gene_set_sizes.csv"),
  bom = TRUE
)

message(
  "Prepared MGI MOUSE-mod mappings (",
  paste(go_relation_scope, collapse = ","),
  "): ",
  annotation_root
)
