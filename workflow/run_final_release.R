#!/usr/bin/env Rscript

# One-command build of the final manuscript figures from frozen publication
# inputs. The destination must be new: approved outputs are never overwritten.

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate workflow/run_final_release.R", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else file.path(root, "release")
if (!grepl("^/", out)) out <- file.path(root, out)
if (file.exists(out)) stop("Output already exists; choose a new directory: ", out, call. = FALSE)

for (directory in c("Main", "Suppl", "tables", "provenance")) {
  dir.create(file.path(out, directory), recursive = TRUE, showWarnings = FALSE)
}
out <- normalizePath(out, mustWork = TRUE)

sha256_file <- function(path) {
  answer <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  if (!length(answer)) stop("Unable to hash: ", path, call. = FALSE)
  sub("[[:space:]].*$", "", answer[[1L]])
}

input_manifest <- fread(file.path(root, "provenance/publication_input_sha256.csv"))
input_paths <- file.path(root, "data/publication_input", input_manifest$relative_path)
if (!all(file.exists(input_paths))) stop("A frozen publication input is missing", call. = FALSE)
input_hashes <- vapply(input_paths, sha256_file, character(1L))
if (!identical(unname(input_hashes), unname(input_manifest$sha256))) {
  changed <- input_manifest$relative_path[input_hashes != input_manifest$sha256]
  stop("Frozen publication input SHA-256 mismatch: ", paste(changed, collapse = "; "), call. = FALSE)
}

# The final manuscript table set predates the two later Fig. 6 metadata-audit
# exports. Keep exactly the 27 selected manuscript tables.
excluded_tables <- c(
  "Fig_6_mouse_metadata_grouped.csv",
  "Fig_6_mouse_sample_metadata_complete_audit.csv"
)
table_sources <- list.files(file.path(root, "results/tables"), pattern = "[.]csv$", full.names = TRUE)
table_sources <- sort(table_sources[!basename(table_sources) %in% excluded_tables])
if (length(table_sources) != 27L) stop("Expected 27 frozen final-result tables", call. = FALSE)
reference_manifest <- fread(file.path(root, "provenance/reference_sha256.csv"))
table_relative_paths <- file.path("tables", basename(table_sources))
frozen_table_hashes <- reference_manifest$sha256[
  match(table_relative_paths, reference_manifest$relative_path)
]
if (anyNA(frozen_table_hashes)) {
  stop("The frozen reference manifest does not cover all 27 final tables", call. = FALSE)
}
source_table_hashes <- vapply(table_sources, sha256_file, character(1L))
if (!identical(unname(source_table_hashes), unname(frozen_table_hashes))) {
  stop("A source manuscript table drifted from its frozen SHA-256", call. = FALSE)
}
copied <- file.copy(table_sources, file.path(out, "tables"), overwrite = FALSE, copy.date = TRUE)
if (!all(copied)) stop("Failed to copy frozen final-result tables", call. = FALSE)

table_manifest <- data.table(
  relative_path = table_relative_paths,
  sha256 = frozen_table_hashes
)
fwrite(table_manifest, file.path(out, "provenance/publication_table_sha256.csv"))
copied_sources <- file.copy(
  file.path(root, "provenance/final_figure_sources.csv"),
  file.path(out, "provenance/final_figure_sources.csv"),
  overwrite = FALSE
)
copied_palette <- file.copy(
  file.path(root, "provenance/final_palette_contract.csv"),
  file.path(out, "provenance/final_palette_contract.csv"),
  overwrite = FALSE
)
if (!copied_sources || !copied_palette) stop("Failed to copy final provenance contracts", call. = FALSE)

Sys.setenv(
  PROJECT_ROOT = root,
  PUBLICATION_OUTPUT_DIR = out,
  FIGURE_FONT = Sys.getenv("FIGURE_FONT", unset = "Arial"),
  FIGURE1_LABELS = "none",
  FINAL_RELEASE_PALETTE_AUDIT = "1"
)
old_wd <- setwd(root)
on.exit(setwd(old_wd), add = TRUE)

run_r <- function(script) {
  status <- system2("Rscript", script)
  if (!identical(status, 0L)) stop("Figure script failed: ", script, call. = FALSE)
}

# Fig. S1 retains the committed publication renderer. All remaining figures use
# the isolated final-figure scripts below.
run_r("R/figures/02_hallmark_gsea.R")

for (script in file.path("R/final_figures", c(
  "02_fig1b_metadata_bubble.R",
  "03_fig1c_fig2_go.R",
  "04_upset_membership.R",
  "05_expression_tree_heatmaps.R",
  "06_meta_dotplots.R",
  "07_figS5_counts.R"
))) {
  run_r(script)
}

run_r("R/figures/10_pathway_gene_tissue_boxplots_linear.R")

status <- system2("bash", c("workflow/run_log2fc_per_tissue_pathway_spearman_20260902.sh", shQuote(out), "name_only"))
if (!identical(status, 0L)) stop("Sample-level Spearman renderer failed", call. = FALSE)
spearman_manifest_path <- file.path(out, "provenance/figure_manifest.csv")
spearman_manifest <- fread(spearman_manifest_path)
spearman_manifest[, path := fifelse(kind == "overall", file.path("Main", paste0("Fig_2.", tools::file_ext(path))),
  fifelse(kind == "merged_per_tissue", file.path("Suppl", paste0("Fig_S6_log2fc_per_tissue_merged.", tools::file_ext(path))),
    file.path("Suppl/Fig_S6_per_tissue", basename(path))))]
if (!all(file.exists(file.path(out, spearman_manifest$path)))) stop("Spearman figure manifest has a broken path")
fwrite(spearman_manifest, spearman_manifest_path)

python <- Sys.getenv("PYTHON", unset = "python3")
status <- system2(python, "python/render_final_common_direction.py")
if (!identical(status, 0L)) stop("Final Python heatmap renderer failed", call. = FALSE)

# Record every frozen input, including later sample-level revisions.
all_inputs <- sort(list.files(file.path(root, "data/publication_input"), recursive = TRUE, full.names = TRUE))
fwrite(data.table(relative_path = substring(all_inputs, nchar(root) + 2L),
                  sha256 = vapply(all_inputs, sha256_file, character(1L))),
       file.path(out, "provenance/all_input_sha256.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "provenance/sessionInfo.txt"))
writeLines(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), file.path(out, "provenance/source_revision.txt"))
file.copy(file.path(root, "README.md"), file.path(out, "README.md"))
file.copy(file.path(root, "results/tables/Fig_6_mouse_sample_metadata_complete_audit.csv"), file.path(out, "provenance"))

status <- system2("Rscript", c("tests/validate_final_release.R", shQuote(out), "full"))
if (!identical(status, 0L)) stop("Final-release validation failed", call. = FALSE)
status <- system2(python, c("workflow/build_review_gallery.py", shQuote(out)))
if (!identical(status, 0L)) stop("Figure review gallery failed", call. = FALSE)
output_files <- sort(list.files(out, recursive = TRUE, full.names = TRUE))
fwrite(data.table(relative_path = substring(output_files, nchar(out) + 2L),
                  sha256 = vapply(output_files, sha256_file, character(1L))),
       file.path(out, "provenance/output_sha256.csv"))
message("FINAL_RELEASE_BUILD_PASS: ", out, " (Fig. 1a intentionally excluded)")
