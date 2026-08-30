#!/usr/bin/env Rscript

# Conservative visual refresh. This workflow copies the frozen publication
# tables, renders the same 16 figures with FIGURE_STYLE=A, and writes to a new
# directory so the approved results/ tree is never overwritten.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate workflow/run_style_A.R", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else file.path(root, "results_style_A")
if (!grepl("^/", out)) out <- file.path(root, out)
if (file.exists(out)) stop("Output already exists; choose a new directory: ", out, call. = FALSE)
dir.create(file.path(out, "Main"), recursive = TRUE)
dir.create(file.path(out, "Suppl"), recursive = TRUE)
dir.create(file.path(out, "tables"), recursive = TRUE)
out <- normalizePath(out, mustWork = TRUE)

table_sources <- list.files(file.path(root, "results/tables"), full.names = TRUE)
if (length(table_sources) != 27L) stop("Expected 27 frozen publication tables", call. = FALSE)
copied <- file.copy(table_sources, file.path(out, "tables"), overwrite = FALSE, copy.date = TRUE)
if (!all(copied)) stop("Failed to copy publication tables", call. = FALSE)

Sys.setenv(
  PROJECT_ROOT = root,
  PUBLICATION_OUTPUT_DIR = out,
  FIGURE_STYLE = "A"
)
old_wd <- setwd(root)
on.exit(setwd(old_wd), add = TRUE)
r_scripts <- file.path("R/figures", c(
  "01_go_overview.R", "02_hallmark_gsea.R", "03_kidney_thymus_venn.R",
  "04_qsmooth_tree_heatmaps.R", "05_meta_log2fc_heatmaps.R",
  "06_go_upsets.R", "07_seven_pathway_counts.R"
))
for (script in r_scripts) {
  status <- system2("Rscript", script)
  if (!identical(status, 0L)) stop("Figure script failed: ", basename(script), call. = FALSE)
}

python <- Sys.getenv("PYTHON", unset = "python3")
status <- system2(python, "python/plot_common_direction_heatmaps.py")
if (!identical(status, 0L)) stop("Python heatmap script failed", call. = FALSE)

status <- system2("Rscript", c("tests/validate_style_A_contract.R", shQuote(out)))
if (!identical(status, 0L)) stop("Style A validation failed", call. = FALSE)
message("Style A outputs rendered in: ", out)
