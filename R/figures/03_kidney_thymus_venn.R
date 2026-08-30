#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(scales)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/figures/style_helpers.R"))
if (figure_style_c) source(file.path(root, "R/figures/style_C_helpers.R"))
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

  if (figure_style_c) {
    # C-story set panels use an UpSet display.  The membership matrix is read
    # directly from the frozen publication input; only the geometry changes.
    membership <- fread(file.path(input_dir, paste0("04_membership_matrix_", family, ".csv")))
    sizes <- fread(file.path(input_dir, paste0("07_set_sizes_", family, ".csv")))
    flag_cols <- sizes$matrix_column
    assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
    assert(all(flag_cols %chin% names(membership)), "Set columns are missing")
    membership[, pattern := do.call(paste0, lapply(.SD, as.integer)), .SDcols = flag_cols]
    zero_pattern <- paste(rep("0", n_sets), collapse = "")
    intersections <- membership[pattern != zero_pattern, .(n_genes = .N), by = pattern]
    intersections <- intersections[order(-n_genes, pattern)]
    intersections[, `:=`(
      intersection_id = seq_len(.N),
      intersection_label = sprintf("%02d", seq_len(.N))
    )]
    set_names <- as.character(sizes$set_name)
    display_set_names <- set_names
    display_set_names[display_set_names == "Heart / Heart right ventricle"] <- "Heart right ventricle"
    set_levels <- rev(display_set_names)
    bits <- do.call(rbind, strsplit(intersections$pattern, split = "", fixed = TRUE))
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
    p_sets_c <- ggplot(set_plot_data, aes(y = set_name, x = n_genes)) +
      geom_col(width = 0.70, fill = STYLE_C_COLOURS$panel, colour = STYLE_C_COLOURS$blue, linewidth = 0.35) +
      geom_text(aes(label = n_genes), hjust = -0.14, size = 3.0, colour = STYLE_C_COLOURS$ink) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(x = "Set size", y = NULL) +
      theme_style_C(9.3) +
      theme(
        panel.grid = element_blank(), axis.text.y = element_text(size = 8.8),
        axis.text.x = element_text(size = 8), plot.margin = margin(4, 8, 4, 4)
      )
    core_pattern <- paste(rep("1", n_sets), collapse = "")
    intersections[, is_core := pattern == core_pattern]
    p_intersections_c <- ggplot(
      intersections,
      aes(x = factor(intersection_label, levels = intersections$intersection_label), y = n_genes, fill = is_core)
    ) +
      geom_col(width = 0.74, colour = STYLE_C_COLOURS$white, linewidth = 0.20) +
      geom_text(aes(label = n_genes), vjust = -0.20, size = 2.65, colour = STYLE_C_COLOURS$ink) +
      scale_fill_manual(values = c(`FALSE` = STYLE_C_COLOURS$blue, `TRUE` = STYLE_C_COLOURS$gold), guide = "none") +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(sigma = 1),
        labels = scales::comma,
        expand = expansion(mult = c(0, 0.14))
      ) +
      labs(x = NULL, y = "Intersection size (pseudo-log scale)") +
      theme_style_C(9.1) +
      theme(
        panel.grid.major.y = element_line(colour = STYLE_C_COLOURS$grid, linewidth = 0.30),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8.8),
        plot.margin = margin(4, 8, 0, 4)
      )
    p_matrix_c <- ggplot(matrix_long, aes(intersection_label, set_name)) +
      geom_point(aes(alpha = included), size = 2.35, colour = STYLE_C_COLOURS$ink) +
      geom_segment(
        data = matrix_long[included == 1L, .(ymin = min(as.integer(set_name)), ymax = max(as.integer(set_name))), by = intersection_label],
        aes(x = intersection_label, xend = intersection_label, y = ymin, yend = ymax),
        inherit.aes = FALSE, colour = STYLE_C_COLOURS$ink, linewidth = 0.60
      ) +
      scale_alpha_continuous(limits = c(0, 1), range = c(0.12, 1), guide = "none") +
      labs(x = "Set intersection", y = NULL) +
      theme_style_C(9.1) +
      theme(
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = STYLE_C_COLOURS$grid, linewidth = 0.30),
        axis.text.x = element_text(size = 6.6, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8.8), axis.title.x = element_text(size = 8.8),
        plot.margin = margin(0, 8, 4, 4)
      )
    subtitle_c <- paste0(
      "GO:0006974 Ensembl IDs | union = ", union_n,
      " | all-set intersection = ", all_n,
      " | gold marks the all-set core"
    )
    plot_c <- p_sets_c + (p_intersections_c / p_matrix_c + plot_layout(heights = c(1.20, 1))) +
      plot_layout(widths = c(1.05, 4.15)) +
      plot_annotation(
        title = if (family == "thymus_four_comparisons") {
          "Thymus: shared and comparison-specific GO:0006974 membership"
        } else {
          "Kidney: shared and dataset-specific GO:0006974 membership"
        },
        subtitle = subtitle_c,
        caption = "Exact Ensembl membership counts are unchanged; the UpSet geometry makes the shared core and minority intersections readable."
      ) &
      theme(
        plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5, colour = STYLE_C_COLOURS$ink),
        plot.subtitle = element_text(size = 8.8, hjust = 0.5, colour = STYLE_C_COLOURS$muted),
        plot.caption = element_text(size = 7.2, hjust = 0.5, colour = STYLE_C_COLOURS$muted),
        plot.margin = margin(8, 12, 8, 12)
      )
    stem_c <- sub("[.]png$", "", spec[[family]]$output)
    save_style_C(plot_c, out, stem_c, 13.5, 8.5)
    next
  }

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
