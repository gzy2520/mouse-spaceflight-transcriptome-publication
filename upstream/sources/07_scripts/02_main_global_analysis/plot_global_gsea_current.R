#!/usr/bin/env Rscript

# Descriptive global-GSEA figures for the current stable-ID rerun.
# These figures deliberately do not pool heterogeneous tissues into a common
# inferential meta-analysis. Ensembl IDs were used upstream; labels below are
# MSigDB gene-set names only.

suppressPackageStartupMessages({
  library(data.table)
})
set.seed(25)

root <- normalizePath(".", mustWork = TRUE)
global_dir <- file.path(
  root, "03_analysis_results/10_full_stable_id_rerun_20260823/02_global_gsea"
)
term_path <- file.path(global_dir, "03_global_descriptive_term_summary.csv.gz")
mission_path <- file.path(global_dir, "02_mission_descriptive_GSEA_summary.csv.gz")
stopifnot(file.exists(term_path), file.exists(mission_path))

figure_dir <- file.path(global_dir, "figures")
table_dir <- file.path(global_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

term <- fread(term_path)
mission <- fread(mission_path)
expected_libraries <- c(
  "MSigDB_2026.1.Mm_Hallmark",
  "MSigDB_2026.1.Mm_GO_BP",
  "MSigDB_2026.1.Mm_Reactome",
  "MSigDB_2026.1.Mm_WikiPathways"
)
stopifnot(all(expected_libraries %chin% unique(term$library)))

library_label <- c(
  "MSigDB_2026.1.Mm_Hallmark" = "Hallmark",
  "MSigDB_2026.1.Mm_GO_BP" = "GO BP",
  "MSigDB_2026.1.Mm_Reactome" = "Reactome",
  "MSigDB_2026.1.Mm_WikiPathways" = "WikiPathways"
)
clean_term_label <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  x
}

term[, `:=`(
  library_display = factor(library_label[library], levels = unname(library_label)),
  term_display = clean_term_label(gene_set_name)
)]
stopifnot(all(is.finite(term$mean_of_mission_mean_NES)))

# One figure contains one MSigDB database only. Each figure selects six terms
# at each NES extreme from that database using descriptive global mean NES.
top_n <- 6L
file_stem <- c(
  "MSigDB_2026.1.Mm_Hallmark" = "01_Hallmark_global_top_terms",
  "MSigDB_2026.1.Mm_GO_BP" = "02_GO_BP_global_top_terms",
  "MSigDB_2026.1.Mm_Reactome" = "03_Reactome_global_top_terms",
  "MSigDB_2026.1.Mm_WikiPathways" = "04_WikiPathways_global_top_terms"
)
selection <- rbindlist(lapply(expected_libraries, function(one_library) {
  # Require full mission coverage for a cross-mission display ranking. Terms
  # absent after fgsea size filtering in one or more missions stay in tables.
  x <- term[library == one_library & n_missions == 13L]
  stopifnot(nrow(x) >= 2L * top_n)
  selected <- unique(rbindlist(list(
    head(x[order(-mean_of_mission_mean_NES, gene_set_name)], top_n),
    head(x[order(mean_of_mission_mean_NES, gene_set_name)], top_n)
  )), by = "gene_set_name")
  selected[, `:=`(
    direction = fifelse(mean_of_mission_mean_NES >= 0, "Positive NES", "Negative NES"),
    figure_file_stem = file_stem[[one_library]]
  )]
  selected[, term_plot := factor(
    term_display,
    levels = selected[order(mean_of_mission_mean_NES), term_display]
  )]
  selected
}), use.names = TRUE, fill = TRUE)
fwrite(
  selection[, .(
    library, figure_file_stem, gene_set_name, term_display, n_missions,
    mean_of_mission_mean_NES, median_of_mission_mean_NES,
    positive_missions, negative_missions, direction
  )],
  file.path(table_dir, "01_global_top_term_figure_selection.csv")
)

message("Selected the 48 manuscript global enrichment terms.")
