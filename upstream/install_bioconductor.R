#!/usr/bin/env Rscript
# Optional isolated dependency setup for the original YARN/qsmooth algorithm.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
lib <- args[[1L]]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(lib), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib)
BiocManager::install(c("Biobase", "yarn"), version = "3.20", lib = lib, update = FALSE, ask = FALSE, Ncpus = 4)
stopifnot(requireNamespace("Biobase", quietly = TRUE), requireNamespace("yarn", quietly = TRUE))
writeLines(capture.output(sessionInfo()), file.path(lib, "upstream_install_sessionInfo.txt"))
