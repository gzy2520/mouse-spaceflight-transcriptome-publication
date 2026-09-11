#!/usr/bin/env Rscript
# Render the approved manuscript from publication inputs into a new destination.
suppressPackageStartupMessages(library(data.table))
set.seed(25)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."))
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "reproduced_results/final"
if (!grepl("^/", out)) out <- file.path(root, out)
if (file.exists(out)) stop("Choose a new output directory: ", out)
for (d in c("Main", "Suppl", "tables", "provenance")) dir.create(file.path(out, d), recursive = TRUE)
out <- normalizePath(out)
setwd(root)
Sys.setenv(PROJECT_ROOT = root, PUBLICATION_OUTPUT_DIR = out,
           FIGURE_FONT = Sys.getenv("FIGURE_FONT", "Arial"), FIGURE1_LABELS = "none",
           FINAL_RELEASE_PALETTE_AUDIT = "1")
run <- function(program, args) {
  code <- system2(program, vapply(args, shQuote, character(1)))
  if (code != 0L) stop("Failed: ", program, " ", paste(args, collapse = " "))
}
python <- Sys.getenv("PYTHON", "python3")
tables <- list.files("results/tables", pattern = "[.]csv$", full.names = TRUE)
tables <- tables[!basename(tables) %in% c("Fig_6_mouse_metadata_grouped.csv", "Fig_6_mouse_sample_metadata_complete_audit.csv")]
file.copy(tables, file.path(out, "tables"))
file.copy("provenance/final_figure_sources.csv", file.path(out, "provenance"))
run("Rscript", "R/figures/02_hallmark_gsea.R")
for (s in c("02_fig1b_metadata_bubble.R", "03_fig1c_fig2_go.R", "04_upset_membership.R",
            "05_expression_tree_heatmaps.R", "06_meta_dotplots.R", "07_figS5_counts.R")) {
  run("Rscript", file.path("R/final_figures", s))
}
run("Rscript", "R/figures/09_pathway_gene_tissue_boxplots_log.R")
run("bash", c("workflow/run_log2fc_per_tissue_pathway_spearman_20260902.sh", out, "name_only"))
m <- fread(file.path(out, "provenance/figure_manifest.csv"))
m[, path := fifelse(kind == "overall", file.path("Main", paste0("Fig_2.", tools::file_ext(path))),
  fifelse(kind == "merged_per_tissue", file.path("Suppl", paste0("Fig_S6_log2fc_per_tissue_merged.", tools::file_ext(path))),
    file.path("Suppl/Fig_S6_per_tissue", basename(path))))]
fwrite(m, file.path(out, "provenance/figure_manifest.csv"))
run(python, "python/render_final_common_direction.py")
file.copy("results/tables/Fig_6_mouse_sample_metadata_complete_audit.csv", file.path(out, "provenance"))
run(python, c("workflow/build_all_publication_supplementary_tables.py", "--output-dir",
              file.path(out, "Supplementary_Tables"), "--figure-dir", out))
writeLines(capture.output(sessionInfo()), file.path(out, "provenance/sessionInfo.txt"))
writeLines(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), file.path(out, "provenance/source_revision.txt"))
writeLines(c("# Reproduced manuscript figures and tables", "",
  "Main/ and Suppl/: 51 PNG/PDF pairs. Supplementary_Tables/: eight workbooks.",
  "The supplied two-way ANOVA workbook is preserved verbatim; its model is not recalculated.",
  "See provenance/source_revision.txt and provenance/sessionInfo.txt for the build context.",
  "Reproduction instructions: https://github.com/gzy2520/mouse-spaceflight-transcriptome-publication"),
  file.path(out, "README.md"))
run(python, c("workflow/build_review_gallery.py", out))
message("Rebuilt manuscript figures and tables: ", out)
