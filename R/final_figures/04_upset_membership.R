#!/usr/bin/env Rscript

# Reproducible UpSet displays from frozen Ensembl-ID membership matrices.
# Fig 3a uses a half-height layout; supplementary panels retain a taller page.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(25)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/final_figures/helpers.R"))

specs <- list(
  list(
    figure = "Fig_3a", output = "Main/Fig_3a",
    membership = "data/publication_input/venn/04_membership_matrix_thymus_four_comparisons.csv",
    sizes = "data/publication_input/venn/07_set_sizes_thymus_four_comparisons.csv",
    title = "GO:0006974 repertoire | thymus comparison sets", compact = TRUE
  ),
  list(
    figure = "Fig_S2a", output = "Suppl/Fig_S2a",
    membership = "data/publication_input/venn/04_membership_matrix_kidney_five_datasets.csv",
    sizes = "data/publication_input/venn/07_set_sizes_kidney_five_datasets.csv",
    title = "GO:0006974 repertoire | kidney dataset sets", compact = FALSE
  ),
  list(
    figure = "Fig_S3a", output = "Suppl/Fig_S3a",
    membership = "results/tables/04_membership_matrix_up_tissues.csv",
    sizes = "results/tables/07_set_sizes_up_tissues.csv",
    title = "GO:0006974 repertoire | positive-NES tissue sets", compact = FALSE
  ),
  list(
    figure = "Fig_S4a", output = "Suppl/Fig_S4a",
    membership = "results/tables/04_membership_matrix_down_tissues.csv",
    sizes = "results/tables/07_set_sizes_down_tissues.csv",
    title = "GO:0006974 repertoire | negative-NES tissue sets", compact = FALSE
  )
)
upset_palette <- c(intersection = "#C65A5A", panel = "#F4F6F8")

render_upset <- function(spec) {
  membership <- fread(file.path(root, spec$membership))
  sizes <- fread(file.path(root, spec$sizes))
  flag_cols <- as.character(sizes$matrix_column)
  assert_final(all(flag_cols %chin% names(membership)), paste(spec$figure, "membership columns missing"))
  membership_values <- unlist(membership[, ..flag_cols], use.names = FALSE)
  assert_final(all(membership_values %in% c(0, 1)), paste(spec$figure, "membership is not binary"))
  assert_final(!anyDuplicated(membership$ensembl_id), paste(spec$figure, "Ensembl IDs are duplicated"))

  membership[, pattern := do.call(paste0, lapply(.SD, as.integer)), .SDcols = flag_cols]
  zero_pattern <- paste(rep("0", length(flag_cols)), collapse = "")
  intersections <- membership[pattern != zero_pattern, .(n_genes = .N), by = pattern]
  setorder(intersections, -n_genes, pattern)
  intersections[, `:=`(intersection_id = .I, intersection_label = sprintf("%02d", .I))]
  core_pattern <- paste(rep("1", length(flag_cols)), collapse = "")
  intersections[, is_core := pattern == core_pattern]

  display_sets <- short_tissue_name_final(as.character(sizes$set_name))
  set_levels <- rev(display_sets)
  bits <- do.call(rbind, strsplit(intersections$pattern, "", fixed = TRUE))
  matrix_long <- rbindlist(lapply(seq_len(nrow(intersections)), function(i) {
    data.table(
      intersection_id = intersections$intersection_id[[i]],
      intersection_label = intersections$intersection_label[[i]],
      set_name = display_sets,
      set_y = match(display_sets, set_levels),
      included = as.integer(bits[i, ])
    )
  }))
  matrix_long[, intersection_label := factor(intersection_label, levels = intersections$intersection_label)]
  segment_data <- matrix_long[included == 1L,
                              .(ymin = min(set_y), ymax = max(set_y)), by = intersection_label]

  set_plot_data <- copy(sizes)
  set_plot_data[, `:=`(display_set = factor(display_sets, levels = set_levels))]
  p_sets <- ggplot(set_plot_data, aes(y = display_set, x = n_genes)) +
    geom_col(width = 0.70, fill = upset_palette[["panel"]],
             colour = FINAL_COLOURS$ink, linewidth = 0.42) +
    geom_text(aes(label = n_genes), hjust = -0.15,
              size = if (spec$compact) 4.70 else 3.2, colour = FINAL_COLOURS$ink) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
    labs(x = "Set size", y = NULL) +
    theme_final(if (spec$compact) 11.5 else 10.2) +
    theme(
      axis.text.y = element_text(size = if (spec$compact) 14.5 else 9.0),
      axis.text.x = element_text(size = if (spec$compact) 12.5 else 8.3),
      plot.margin = margin(2, 8, 2, 2)
    )

  p_intersections <- ggplot(
    intersections,
    aes(factor(intersection_label, levels = intersection_label), n_genes, fill = is_core)
  ) +
    geom_col(width = 0.74, colour = FINAL_COLOURS$white, linewidth = 0.22) +
    geom_text(aes(label = n_genes), vjust = -0.20,
              size = if (spec$compact) 4.50 else 2.75, colour = FINAL_COLOURS$ink) +
    scale_fill_manual(values = c(`FALSE` = upset_palette[["intersection"]],
                                 `TRUE` = upset_palette[["intersection"]]), guide = "none") +
    scale_y_continuous(
      trans = pseudo_log_trans(sigma = 1), labels = comma,
      breaks = c(1, 10, 100, 1000),
      expand = expansion(mult = c(0, 0.16))
    ) +
    labs(x = NULL, y = "Intersection size") +
    theme_final(if (spec$compact) 11.5 else 10.0) +
    theme(
      panel.grid.major.y = element_line(colour = FINAL_COLOURS$grid, linewidth = 0.30),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = if (spec$compact) 12.0 else 8.0),
      axis.title.y = element_text(size = if (spec$compact) 12.5 else 8.8,
                                  margin = margin(r = 4)),
      plot.margin = margin(2, 8, 0, 2)
    )

  p_matrix <- ggplot(matrix_long, aes(intersection_label, set_y)) +
    geom_segment(
      data = segment_data,
      aes(x = intersection_label, xend = intersection_label, y = ymin, yend = ymax),
      inherit.aes = FALSE, colour = FINAL_COLOURS$ink,
      linewidth = if (spec$compact) 0.78 else 0.62
    ) +
    geom_point(aes(alpha = included), size = if (spec$compact) 3.8 else 2.45,
               colour = FINAL_COLOURS$ink) +
    scale_alpha_continuous(limits = c(0, 1), range = c(0.10, 1), guide = "none") +
    scale_y_continuous(breaks = seq_along(set_levels), labels = set_levels,
                       expand = expansion(add = 0.35)) +
    labs(x = "Set intersection", y = NULL) +
    theme_final(if (spec$compact) 11.5 else 10.0) +
    theme(
      panel.grid.major.y = element_line(colour = FINAL_COLOURS$grid, linewidth = 0.30),
      axis.text.x = element_text(size = if (spec$compact) 10.5 else 6.7,
                                 angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = if (spec$compact) 14.5 else 8.8),
      axis.title.x = element_text(size = if (spec$compact) 13.5 else 8.8),
      plot.margin = margin(0, 8, 2, 2)
    )

  union_n <- sum(membership$pattern != zero_pattern)
  core_n <- intersections[pattern == core_pattern, n_genes]
  if (!length(core_n)) core_n <- 0L
  subtitle <- paste0(
    length(flag_cols), " sets | exact Ensembl membership | all-set intersection = ",
    core_n, " | union = ", union_n
  )
  plot <- p_sets + (p_intersections / p_matrix +
                      plot_layout(heights = if (spec$compact) c(0.78, 1.02) else c(1.20, 1.00))) +
    plot_layout(widths = if (spec$compact) c(1.22, 4.65) else c(1.05, 4.15)) +
    plot_annotation(
      title = spec$title,
      subtitle = subtitle,
      caption = "Exact Ensembl membership and intersection ordering are retained; the all-set core is reported in the subtitle and audit.",
      theme = theme(
        plot.title = element_text(family = FINAL_FONT, face = "bold",
                                  size = if (spec$compact) 22 else 14.5,
                                  hjust = 0.5, colour = FINAL_COLOURS$ink),
        plot.subtitle = element_text(family = FINAL_FONT,
                                     size = if (spec$compact) 14.5 else 9.2,
                                     hjust = 0.5, colour = FINAL_COLOURS$muted),
        plot.caption = element_text(family = FINAL_FONT,
                                    size = if (spec$compact) 11.5 else 7.5,
                                    hjust = 0.5, colour = FINAL_COLOURS$muted),
        plot.margin = margin(6, 10, 6, 10)
      )
    )

  if (spec$compact) {
    save_final_figure(plot, out, spec$output, width = 15.5, height = 4.2, dpi = 500)
  } else {
    save_final_figure(plot, out, spec$output, width = 13.5, height = 8.5, dpi = 450)
  }
  data.table(
    figure = spec$figure, n_sets = length(flag_cols), n_membership_rows = nrow(membership),
    union_n = union_n, core_n = core_n, n_nonzero_intersections = nrow(intersections),
    membership_file = spec$membership, sizes_file = spec$sizes
  )
}

audit <- rbindlist(lapply(specs, render_upset), fill = TRUE)
dir.create(file.path(out, "provenance"), recursive = TRUE, showWarnings = FALSE)
fwrite(audit, file.path(out, "provenance/upset_render_audit.csv"))
for (figure_id in audit$figure) append_final_palette_audit(out, figure_id, upset_palette)
message("FINAL_UPSET_PASS: ", paste(audit$figure, collapse = ", "))
