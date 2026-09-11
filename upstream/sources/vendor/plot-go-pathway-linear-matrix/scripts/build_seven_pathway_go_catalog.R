#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "scripts/build_seven_pathway_go_catalog.R",
    mustWork = TRUE
  )
}
skill_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript build_seven_pathway_go_catalog.R",
      "--obo /path/to/go-basic.obo",
      "--out /path/to/seven_pathway_go_term_catalog.csv",
      "[--provenance-out /path/to/go_catalog_provenance.csv]",
      "[--rules /path/to/go_term_pathway_seed_rules.csv]",
      "[--source-url https://purl.obolibrary.org/obo/go/go-basic.obo]",
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

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

options <- parse_args(commandArgs(trailingOnly = TRUE))
assert(
  !is.null(options$obo) && nzchar(options$obo),
  "--obo is required."
)
assert(
  !is.null(options$out) && nzchar(options$out),
  "--out is required."
)

obo_path <- normalizePath(options$obo, mustWork = TRUE)
rules_path <- normalizePath(
  if (!is.null(options$rules)) {
    options$rules
  } else {
    file.path(skill_root, "assets", "go_term_pathway_seed_rules.csv")
  },
  mustWork = TRUE
)
output_path <- normalizePath(options$out, mustWork = FALSE)
provenance_path <- normalizePath(
  if (!is.null(options[["provenance-out"]])) {
    options[["provenance-out"]]
  } else {
    file.path(dirname(output_path), "go_catalog_provenance.csv")
  },
  mustWork = FALSE
)
source_url <- if (!is.null(options[["source-url"]])) {
  options[["source-url"]]
} else {
  "https://purl.obolibrary.org/obo/go/go-basic.obo"
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(provenance_path), recursive = TRUE, showWarnings = FALSE)

obo_lines <- readLines(obo_path, warn = FALSE, encoding = "UTF-8")
assert(length(obo_lines) > 0L, "The GO OBO file is empty.")
term_starts <- which(obo_lines == "[Term]")
assert(length(term_starts) > 0L, "No [Term] blocks were found in the OBO file.")
term_ends <- c(term_starts[-1L] - 1L, length(obo_lines))

header_lines <- obo_lines[seq_len(term_starts[[1L]] - 1L)]
release_line <- header_lines[
  startsWith(header_lines, "data-version: ")
]
ontology_release <- if (length(release_line)) {
  sub("^data-version: ", "", release_line[[1L]])
} else {
  "not_declared"
}

extract_first <- function(block, prefix, default = "") {
  values <- block[startsWith(block, prefix)]
  if (!length(values)) return(default)
  sub(paste0("^", prefix), "", values[[1L]])
}

parse_term_block <- function(block) {
  go_id <- extract_first(block, "id: ")
  if (!grepl("^GO:[0-9]{7}$", go_id)) return(NULL)

  definition_line <- extract_first(block, "def: ")
  definition <- if (nzchar(definition_line)) {
    definition_line <- sub('^"', "", definition_line)
    sub('" \\[.*$', "", definition_line)
  } else {
    ""
  }
  alternative_ids <- sub(
    "^alt_id: ",
    "",
    block[startsWith(block, "alt_id: ")]
  )

  term <- data.table(
    CanonicalGOID = go_id,
    GO_Name = extract_first(block, "name: "),
    Namespace = extract_first(block, "namespace: "),
    Definition = definition,
    Obsolete = any(block == "is_obsolete: true"),
    AlternativeGOIDs = paste(alternative_ids, collapse = ";")
  )

  is_a_lines <- block[startsWith(block, "is_a: GO:")]
  relationship_lines <- block[startsWith(block, "relationship: ")]
  edges <- list()
  if (length(is_a_lines)) {
    edges[[length(edges) + 1L]] <- data.table(
      Child = go_id,
      Relation = "is_a",
      Parent = sub(
        "^is_a: (GO:[0-9]{7}).*$",
        "\\1",
        is_a_lines
      )
    )
  }
  if (length(relationship_lines)) {
    edges[[length(edges) + 1L]] <- data.table(
      Child = go_id,
      Relation = sub(
        "^relationship: ([^ ]+) .*$",
        "\\1",
        relationship_lines
      ),
      Parent = sub(
        "^relationship: [^ ]+ (GO:[0-9]{7}).*$",
        "\\1",
        relationship_lines
      )
    )
  }

  list(
    term = term,
    edges = rbindlist(edges, use.names = TRUE, fill = TRUE)
  )
}

parsed_blocks <- lapply(seq_along(term_starts), function(index) {
  parse_term_block(obo_lines[term_starts[[index]]:term_ends[[index]]])
})
parsed_blocks <- Filter(Negate(is.null), parsed_blocks)
terms <- rbindlist(lapply(parsed_blocks, `[[`, "term"), use.names = TRUE)
edges <- rbindlist(
  lapply(parsed_blocks, `[[`, "edges"),
  use.names = TRUE,
  fill = TRUE
)

terms <- terms[Obsolete == FALSE]
terms[, Ontology := fifelse(
  Namespace == "biological_process",
  "BP",
  fifelse(
    Namespace == "cellular_component",
    "CC",
    fifelse(Namespace == "molecular_function", "MF", "")
  )
)]
assert(
  all(nzchar(terms$Ontology)) &&
    !anyDuplicated(terms$CanonicalGOID),
  "Active GO terms contain an unsupported namespace or duplicated ID."
)
active_ids <- terms$CanonicalGOID
edges <- unique(
  edges[
    Child %in% active_ids &
      Parent %in% active_ids &
      grepl("^GO:[0-9]{7}$", Parent)
  ]
)

rules <- fread(rules_path)
required_rule_columns <- c(
  "Pathway",
  "PathwayOrder",
  "SeedGOID",
  "ExpansionMode",
  "RuleLabel",
  "Rationale"
)
expected_pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
expected_orders <- setNames(seq_along(expected_pathways), expected_pathways)
allowed_modes <- c("exact", "is_a", "structural", "pathway_related")
assert(
  all(required_rule_columns %in% names(rules)) &&
    setequal(unique(rules$Pathway), expected_pathways) &&
    all(rules$PathwayOrder == expected_orders[rules$Pathway]) &&
    all(rules$ExpansionMode %in% allowed_modes) &&
    !anyDuplicated(rules[, .(Pathway, SeedGOID)]),
  "The seven-pathway seed-rule table violates its fixed schema."
)
assert(
  all(rules$SeedGOID %in% active_ids),
  paste(
    "Seed terms absent from the active GO release:",
    paste(setdiff(rules$SeedGOID, active_ids), collapse = ", ")
  )
)

relations_for_mode <- list(
  exact = character(),
  is_a = "is_a",
  structural = c("is_a", "part_of"),
  pathway_related = c(
    "is_a",
    "part_of",
    "regulates",
    "positively_regulates",
    "negatively_regulates"
  )
)

ontology_closure <- function(seed_go_id, expansion_mode) {
  relations <- relations_for_mode[[expansion_mode]]
  if (!length(relations)) return(seed_go_id)

  seen <- seed_go_id
  frontier <- seed_go_id
  repeat {
    new_ids <- unique(
      edges[
        Parent %in% frontier &
          Relation %in% relations,
        Child
      ]
    )
    new_ids <- setdiff(new_ids, seen)
    if (!length(new_ids)) break
    seen <- c(seen, new_ids)
    frontier <- new_ids
  }
  seen
}

rule_hits <- rbindlist(
  lapply(seq_len(nrow(rules)), function(index) {
    rule <- rules[index]
    candidate_ids <- ontology_closure(
      rule$SeedGOID,
      rule$ExpansionMode
    )
    data.table(
      CanonicalGOID = candidate_ids,
      Pathway = rule$Pathway,
      PathwayOrder = as.integer(rule$PathwayOrder),
      Coefficient = as.integer(rule$PathwayOrder),
      MatchedSeedGOID = rule$SeedGOID,
      ExpansionMode = rule$ExpansionMode,
      MatchType = fifelse(
        candidate_ids == rule$SeedGOID,
        "seed_exact",
        "ontology_related"
      ),
      RuleLabel = rule$RuleLabel,
      Rationale = rule$Rationale
    )
  }),
  use.names = TRUE
)
rule_hits <- unique(
  rule_hits,
  by = c("CanonicalGOID", "Pathway", "MatchedSeedGOID")
)

canonical_catalog <- rule_hits[
  ,
  .(
    PathwayOrder = first(PathwayOrder),
    Coefficient = first(Coefficient),
    MatchedSeedGOIDs = paste(
      sort(unique(MatchedSeedGOID)),
      collapse = ";"
    ),
    ExpansionModes = paste(
      sort(unique(ExpansionMode)),
      collapse = ";"
    ),
    MatchTypes = paste(sort(unique(MatchType)), collapse = ";"),
    RuleLabels = paste(sort(unique(RuleLabel)), collapse = "; "),
    MappingRationale = paste(
      sort(unique(Rationale)),
      collapse = " | "
    )
  ),
  by = .(CanonicalGOID, Pathway)
]
canonical_catalog <- merge(
  canonical_catalog,
  terms[
    ,
    .(
      CanonicalGOID,
      GO_Name,
      Ontology,
      Definition,
      AlternativeGOIDs
    )
  ],
  by = "CanonicalGOID",
  all.x = TRUE,
  sort = FALSE
)
assert(
  !anyNA(canonical_catalog$GO_Name),
  "At least one catalog term could not be joined to the active GO terms."
)

canonical_rows <- copy(canonical_catalog)
canonical_rows[, `:=`(
  GO_ID = CanonicalGOID,
  IsAlternativeID = FALSE
)]

alternative_index <- terms[nzchar(AlternativeGOIDs)]
alternative_index <- alternative_index[
  ,
  .(GO_ID = unlist(strsplit(AlternativeGOIDs, ";", fixed = TRUE))),
  by = CanonicalGOID
]
alternative_rows <- merge(
  alternative_index,
  canonical_catalog,
  by = "CanonicalGOID",
  all = FALSE,
  sort = FALSE
)
alternative_rows[, IsAlternativeID := TRUE]

catalog <- rbindlist(
  list(canonical_rows, alternative_rows),
  use.names = TRUE,
  fill = TRUE
)
catalog[, `:=`(
  OntologyRelease = ontology_release,
  OntologySHA256 = digest(
    file = obo_path,
    algo = "sha256",
    serialize = FALSE
  ),
  RuleSetSHA256 = digest(
    file = rules_path,
    algo = "sha256",
    serialize = FALSE
  )
)]
catalog[, AlternativeGOIDs := NULL]
setcolorder(
  catalog,
  c(
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
)
setorder(catalog, PathwayOrder, CanonicalGOID, IsAlternativeID, GO_ID)
assert(
  all(catalog$GO_ID %in% c(active_ids, alternative_index$GO_ID)) &&
    all(catalog$PathwayOrder == expected_orders[catalog$Pathway]) &&
    all(catalog$Coefficient == catalog$PathwayOrder) &&
    !anyDuplicated(catalog[, .(GO_ID, Pathway)]),
  "The generated GO catalog violates the seven-pathway contract."
)

fwrite(catalog, output_path)
catalog_sha256 <- digest(
  file = output_path,
  algo = "sha256",
  serialize = FALSE
)

pathway_counts <- canonical_catalog[
  ,
  .(
    CanonicalTermCount = uniqueN(CanonicalGOID),
    BP = uniqueN(CanonicalGOID[Ontology == "BP"]),
    CC = uniqueN(CanonicalGOID[Ontology == "CC"]),
    MF = uniqueN(CanonicalGOID[Ontology == "MF"])
  ),
  by = .(Pathway, PathwayOrder)
][order(PathwayOrder)]

provenance <- rbindlist(
  list(
    data.table(
      RecordType = "metadata",
      Key = c(
        "GO OBO input",
        "GO OBO source URL",
        "GO ontology release",
        "GO OBO SHA256",
        "Seed-rule input",
        "Seed-rule SHA256",
        "Catalog output",
        "Catalog SHA256",
        "Active canonical GO terms in OBO",
        "Active ontology edges",
        "Mapped canonical GO terms",
        "Mapped alternative GO IDs",
        "GO ID-pathway rows",
        "R version",
        "data.table version",
        "digest version"
      ),
      Value = c(
        obo_path,
        source_url,
        ontology_release,
        unique(catalog$OntologySHA256),
        rules_path,
        unique(catalog$RuleSetSHA256),
        output_path,
        catalog_sha256,
        nrow(terms),
        nrow(edges),
        uniqueN(catalog$CanonicalGOID),
        uniqueN(catalog[IsAlternativeID == TRUE, GO_ID]),
        nrow(catalog),
        as.character(getRversion()),
        as.character(packageVersion("data.table")),
        as.character(packageVersion("digest"))
      ),
      Pathway = "",
      PathwayOrder = NA_integer_,
      CanonicalTermCount = NA_integer_,
      BP = NA_integer_,
      CC = NA_integer_,
      MF = NA_integer_
    ),
    pathway_counts[
      ,
      .(
        RecordType = "pathway_count",
        Key = "Canonical term count",
        Value = as.character(CanonicalTermCount),
        Pathway,
        PathwayOrder,
        CanonicalTermCount,
        BP,
        CC,
        MF
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)
fwrite(provenance, provenance_path)

message(
  "Created a seven-pathway GO catalog with ",
  uniqueN(catalog$CanonicalGOID),
  " canonical terms, ",
  uniqueN(catalog[IsAlternativeID == TRUE, GO_ID]),
  " alternative IDs, and ",
  nrow(catalog),
  " GO ID-pathway rows from ",
  ontology_release,
  ": ",
  output_path
)
