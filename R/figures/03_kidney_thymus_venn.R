#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(ggplot2)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/figures/style_helpers.R"))
input_dir <- file.path(root, "data/publication_input/venn")

spec <- list(
  kidney_five_datasets = list(
    title = "GO:0006974 repertoire | Kidney: five OSDR datasets",
    output = "Suppl/Fig_S2a.png"
  ),
  thymus_four_comparisons = list(
    title = "GO:0006974 repertoire | Thymus: four requested comparisons",
    output = "Main/Fig_3a.png"
  )
)

read_sets <- function(family) {
  membership <- read.csv(
    file.path(input_dir, paste0("04_membership_matrix_", family, ".csv")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  sizes <- read.csv(
    file.path(input_dir, paste0("07_set_sizes_", family, ".csv")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  stopifnot(all(sizes$matrix_column %in% names(membership)))
  sets <- lapply(sizes$matrix_column, function(column) {
    unique(as.character(membership$ensembl_id[membership[[column]] == 1]))
  })
  names(sets) <- as.character(sizes$set_name)
  sets
}

for (family in names(spec)) {
  sets <- read_sets(family)
  n_sets <- length(sets)
  stopifnot(n_sets >= 2L, n_sets <= 7L)
  union_n <- length(unique(unlist(sets, use.names = FALSE)))
  all_n <- length(Reduce(intersect, sets))
  subtitle <- paste0(
    "n = 893 GO:0006974 Ensembl IDs | dataset presence: any included sample with source value - 1 > 0",
    " | all-set intersection = ", all_n, " | union = ", union_n
  )
  plot <- ggVennDiagram(
    sets, label = "count", label_alpha = 0,
    label_size = if (n_sets > 5L) 3.4 else 4.2,
    set_size = if (n_sets > 5L) 4.0 else 5.0
  ) +
    scale_fill_gradient(
      low = if (figure_style_a) "#F3F5F7" else "#F7FBFF",
      high = FIGURE_COLOURS$red, name = "Exact region\ncount"
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = spec[[family]]$title, subtitle = subtitle,
      caption = "Geometry is schematic; exact counts are in tables/06_exact_region_counts_*.csv"
    ) +
    if (figure_style_a) {
      theme_void(base_size = 12) +
        theme(
          plot.title = style_title(13, hjust = 0.5),
          plot.subtitle = style_subtitle(8.7, hjust = 0.5),
          plot.caption = style_caption(7.2, hjust = 0.5),
          legend.position = "right", legend.title = element_text(size = 8.2),
          legend.text = element_text(size = 7.5),
          plot.margin = margin(16, 72, 16, 80)
        )
    } else {
      theme_minimal(base_size = 13) +
        theme(
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "#374151"),
          plot.caption = element_text(size = 8, colour = "#4B5563"),
          legend.position = "right", panel.grid = element_blank(),
          plot.margin = margin(20, 100, 20, 130)
        )
    }
  ggsave(
    file.path(out, spec[[family]]$output), plot,
    width = 10, height = 8.5, units = "in", dpi = if (figure_style_a) 600 else 300,
    bg = "white", limitsize = FALSE
  )
}
