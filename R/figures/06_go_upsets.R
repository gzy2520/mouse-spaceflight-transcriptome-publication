#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
input_dir <- file.path(root, "results/tables")

spec <- list(
  up_tissues = list(
    title = "GO:0006974 dataset-union repertoire | NES > 1 tissues (Kidney excluded)",
    output = "Suppl/Fig_S3a.png"
  ),
  down_tissues = list(
    title = "GO:0006974 dataset-union repertoire | NES < -1 tissues (Lung and Thymus excluded)",
    output = "Suppl/Fig_S4a.png"
  )
)

render_one <- function(family) {
  membership <- fread(file.path(input_dir, paste0("04_membership_matrix_", family, ".csv")))
  sizes <- fread(file.path(input_dir, paste0("07_set_sizes_", family, ".csv")))
  exact <- fread(file.path(input_dir, paste0("06_exact_region_counts_", family, ".csv")))
  flag_cols <- sizes$matrix_column
  stopifnot(all(flag_cols %chin% names(membership)), uniqueN(membership$ensembl_id) == nrow(membership))

  membership[, pattern := do.call(paste0, lapply(.SD, as.integer)), .SDcols = flag_cols]
  zero_pattern <- paste(rep("0", length(flag_cols)), collapse = "")
  intersections <- membership[pattern != zero_pattern, .(n_genes = .N), by = pattern]
  intersections <- intersections[order(-n_genes, pattern)]
  exact_nonzero <- exact[
    set_membership_pattern != 0 & n_genes > 0,
    .(pattern = sprintf(paste0("%0", length(flag_cols), "d"), as.integer(set_membership_pattern)), n_genes)
  ]
  check <- merge(intersections, exact_nonzero, by = "pattern", suffixes = c("_matrix", "_table"), all = TRUE)
  stopifnot(nrow(check) == nrow(intersections), all(check$n_genes_matrix == check$n_genes_table))

  intersections[, intersection_id := seq_len(.N)]
  intersections[, intersection_label := sprintf("%02d", intersection_id)]
  bits <- do.call(rbind, strsplit(intersections$pattern, split = "", fixed = TRUE))
  set_names <- as.character(sizes$set_name)
  display_set_names <- ifelse(
    set_names == "Heart / Heart right ventricle", "Heart right ventricle", set_names
  )
  set_levels <- rev(display_set_names)
  matrix_long <- rbindlist(lapply(seq_len(nrow(intersections)), function(i) {
    data.table(
      intersection_id = intersections$intersection_id[[i]],
      intersection_label = intersections$intersection_label[[i]],
      set_name = factor(display_set_names, levels = set_levels),
      included = as.integer(bits[i, ])
    )
  }))
  matrix_long[, intersection_label := factor(intersection_label, levels = intersections$intersection_label)]

  set_plot_data <- copy(sizes)
  set_plot_data[, set_name := factor(display_set_names, levels = set_levels)]
  p_sets <- ggplot(set_plot_data, aes(y = set_name, x = n_genes)) +
    geom_col(width = 0.72, fill = "#F1C4BB", colour = "#C97A70", linewidth = 0.25) +
    geom_text(aes(label = n_genes), hjust = -0.12, size = 3.0, colour = "#7A2E2A") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.17))) +
    labs(x = "Set size", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9.2, colour = "#374151"),
      axis.text.x = element_text(size = 8.5, colour = "#4B5563"),
      axis.title.x = element_text(size = 9.5), plot.margin = margin(4, 7, 4, 4)
    )
  p_intersections <- ggplot(
    intersections,
    aes(x = factor(intersection_label, levels = intersections$intersection_label), y = n_genes)
  ) +
    geom_col(width = 0.76, fill = "#B94C4A", colour = "#7A2E2A", linewidth = 0.2) +
    geom_text(aes(label = n_genes), vjust = -0.22, size = 2.8, colour = "#374151") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16)), labels = comma) +
    labs(x = NULL, y = "Intersection size") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8.5, colour = "#4B5563"),
      axis.title.y = element_text(size = 9.5), plot.margin = margin(4, 7, 0, 4)
    )
  p_matrix <- ggplot(matrix_long, aes(x = intersection_label, y = set_name)) +
    geom_point(aes(alpha = included), size = 2.5, colour = "#374151") +
    geom_segment(
      data = matrix_long[included == 1L, .(
        ymin = min(as.integer(set_name)), ymax = max(as.integer(set_name))
      ), by = intersection_label],
      aes(x = intersection_label, xend = intersection_label, y = ymin, yend = ymax),
      inherit.aes = FALSE, colour = "#374151", linewidth = 0.65
    ) +
    scale_alpha_continuous(limits = c(0, 1), range = c(0.12, 1), guide = "none") +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Set intersection", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
      panel.grid.major.y = element_line(colour = "#E5E7EB", linewidth = 0.35),
      axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, colour = "#4B5563"),
      axis.text.y = element_text(size = 9.2, colour = "#374151"),
      axis.title.x = element_text(size = 9.5), plot.margin = margin(0, 7, 4, 4)
    )

  union_n <- sum(membership$pattern != zero_pattern)
  all_set_n <- intersections[pattern == paste(rep("1", length(flag_cols)), collapse = ""), n_genes]
  if (!length(all_set_n)) all_set_n <- 0L
  subtitle <- paste0(
    "n = 893 GO:0006974 Ensembl IDs | dataset presence: any included sample with source value - 1 > 0",
    " | all-set intersection = ", all_set_n, " | union = ", union_n
  )
  plot <- p_sets + (p_intersections / p_matrix + plot_layout(heights = c(1.22, 1))) +
    plot_layout(widths = c(1.05, 3.9)) +
    plot_annotation(
      title = spec[[family]]$title, subtitle = subtitle,
      caption = "Exact Ensembl membership counts are unchanged; geometry is an UpSet display. Set-size bars use a pale warm colour."
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, colour = "#374151"),
      plot.caption = element_text(size = 8.2, colour = "#4B5563", hjust = 0.5),
      plot.margin = margin(8, 12, 8, 12)
    )
  ggsave(
    file.path(out, spec[[family]]$output), plot,
    width = 14, height = 8.6, units = "in", dpi = 320,
    bg = "white", limitsize = FALSE
  )
}

invisible(lapply(names(spec), render_one))
