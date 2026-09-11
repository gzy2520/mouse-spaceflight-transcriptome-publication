#!/usr/bin/env Rscript

# Current mouse GO annotation entry point.  It keeps the prior dual-sensing
# tree untouched and writes a dated parallel tree composed only of explicit
# GO root terms (no derived DDR-minus-repair set).

root <- normalizePath(".", mustWork = TRUE)
Sys.setenv(
  MOUSE_GO_ANNOTATION_ROOT = file.path(
    root,
    "08_GO_annotation/mgi_mod_gaf_five_relation_concrete_terms_20260824"
  ),
  MOUSE_GO_SOURCE_DIR = file.path(
    root,
    "08_GO_annotation/mgi_mod_gaf_five_relation_20260822/source"
  ),
  MOUSE_GO_MANUAL_DECISION_FILE = file.path(
    root,
    "08_GO_annotation/mgi_mod_gaf_five_relation_20260822/manual_review",
    "00_five_relation_multi_route_manual_decisions.csv"
  )
)
source(file.path(
  root,
  "07_scripts/04_dna_damage_repair_focus",
  "prepare_mouse_go_current_impl.R"
), chdir = FALSE)
