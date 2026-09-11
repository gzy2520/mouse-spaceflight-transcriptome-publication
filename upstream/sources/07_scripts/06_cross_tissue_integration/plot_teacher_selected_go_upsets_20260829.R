#!/usr/bin/env Rscript

# Re-render the two cross-tissue GO:0006974 repertoire figures in one
# consistent UpSet layout.  The membership matrices and exact region tables
# are copied from the frozen 2026-08-26 dataset-level build; this script only
# changes the geometry and the display palette.  Ensembl IDs remain the only
# membership key.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- normalizePath(Sys.getenv("GRADUATE_DESIGN_ROOT", unset = getwd()), mustWork = TRUE)
source_dir <- file.path(
  root, "03_analysis_results/19_ddr_go0006974_dataset_level_venn_20260826"
)
out_dir <- Sys.getenv(
  "DDR_UPSET_OUTPUT_DIR",
  unset = file.path(root, "03_analysis_results/24_teacher_selected_figure_release_20260829/go_upsets")
)
if (!grepl("^/", out_dir)) out_dir <- file.path(root, out_dir)
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)

families <- c("up_tissues", "down_tissues")
titles <- c(
  up_tissues = "GO:0006974 dataset-union repertoire | NES > 1 tissues (Kidney excluded)",
  down_tissues = "GO:0006974 dataset-union repertoire | NES < -1 tissues (Lung and Thymus excluded)"
)
stems <- c(
  up_tissues = "Fig_S3a_NES_gt1_Kidney_excluded_dataset_UpSet_pale_warm",
  down_tissues = "Fig_S4a_NES_lt_minus1_Lung_Thymus_excluded_dataset_UpSet_pale_warm"
)

make_one <- function(family) {
  matrix_path <- file.path(source_dir, "tables", paste0("04_membership_matrix_", family, ".csv"))
  sizes_path <- file.path(source_dir, "tables", paste0("07_set_sizes_", family, ".csv"))
  exact_path <- file.path(source_dir, "tables", paste0("06_exact_region_counts_", family, ".csv"))
  assert(file.exists(matrix_path) && file.exists(sizes_path) && file.exists(exact_path),
         paste("Missing frozen GO membership input for", family))

  membership <- fread(matrix_path, showProgress = FALSE)
  sizes <- fread(sizes_path, showProgress = FALSE)
  exact <- fread(exact_path, showProgress = FALSE)
  flag_cols <- sizes$matrix_column
  assert(all(flag_cols %chin% names(membership)), paste("Membership columns do not match", family))
  assert(uniqueN(membership$ensembl_id) == nrow(membership), paste("Duplicate Ensembl IDs in", family))
  assert(nrow(sizes) >= 2L, paste("Need at least two sets for", family))

  # Recompute the displayed intersection counts directly from the matrix and
  # verify them against the frozen exact-region table.  This guards against a
  # style-only rerender accidentally changing the set definition.
  membership[, pattern := do.call(paste0, lapply(.SD, as.integer)), .SDcols = flag_cols]
  zero_pattern <- paste(rep("0", length(flag_cols)), collapse = "")
  intersections <- membership[pattern != zero_pattern, .(n_genes = .N), by = pattern]
  intersections <- intersections[order(-n_genes, pattern)]
  exact_nonzero <- exact[set_membership_pattern != 0 & n_genes > 0,
                         .(pattern = sprintf(paste0("%0", length(flag_cols), "d"), as.integer(set_membership_pattern)),
                           n_genes)]
  exact_nonzero[, pattern := as.character(pattern)]
  check <- merge(intersections, exact_nonzero, by = "pattern", suffixes = c("_matrix", "_table"), all = TRUE)
  assert(nrow(check) == nrow(intersections) && all(check$n_genes_matrix == check$n_genes_table),
         paste("Matrix/exact-region count mismatch for", family))

  intersections[, intersection_id := seq_len(.N)]
  intersections[, intersection_label := sprintf("%02d", intersection_id)]
  bits <- do.call(rbind, strsplit(intersections$pattern, split = "", fixed = TRUE))
  set_names <- as.character(sizes$set_name)
  # Keep the source set name in the audit tables, but use the concise display
  # label requested for the figures.
  display_set_names <- ifelse(
    set_names == "Heart / Heart right ventricle",
    "Heart right ventricle",
    set_names
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

  # Horizontal set-size bars use a pale warm tone to replace the previous
  # grey block.  The exact count labels remain visible on the bars.
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
      axis.title.x = element_text(size = 9.5),
      plot.margin = margin(4, 7, 4, 4)
    )

  p_intersections <- ggplot(intersections, aes(x = factor(intersection_label, levels = intersections$intersection_label), y = n_genes)) +
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

  # Every set has a faint row guide; only the members of an intersection are
  # filled.  Connecting lines are drawn per intersection to form a standard
  # UpSet matrix without relying on a package that is not installed locally.
  p_matrix <- ggplot(matrix_long, aes(x = intersection_label, y = set_name)) +
    geom_point(aes(alpha = included), size = 2.5, colour = "#374151") +
    geom_segment(
      data = matrix_long[included == 1L,
                         .(ymin = min(as.integer(set_name)), ymax = max(as.integer(set_name))),
                         by = intersection_label],
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

  # Each matrix row is already one unique Ensembl ID, so nonzero patterns
  # directly give the union size; the all-ones pattern gives the exact
  # all-set intersection.
  union_n <- sum(membership$pattern != zero_pattern)
  all_set_n <- intersections[pattern == paste(rep("1", length(flag_cols)), collapse = ""), n_genes]
  if (!length(all_set_n)) all_set_n <- 0L
  subtitle <- paste0(
    "n = 893 GO:0006974 Ensembl IDs | dataset presence: any included sample with source value - 1 > 0",
    " | all-set intersection = ", all_set_n, " | union = ", union_n
  )
  main <- p_sets + (p_intersections / p_matrix + plot_layout(heights = c(1.22, 1))) +
    plot_layout(widths = c(1.05, 3.9)) +
    plot_annotation(
      title = unname(titles[[family]]), subtitle = subtitle,
      caption = "Exact Ensembl membership counts are unchanged; geometry is an UpSet display. Set-size bars use a pale warm colour."
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, colour = "#374151"),
      plot.caption = element_text(size = 8.2, colour = "#4B5563", hjust = 0.5),
      plot.margin = margin(8, 12, 8, 12)
    )

  stem <- unname(stems[[family]])
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  ggsave(png_path, main, width = 14, height = 8.6, units = "in", dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, main, width = 14, height = 8.6, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)

  # Retain self-contained style/input audit copies beside the figures.
  file.copy(matrix_path, file.path(table_dir, basename(matrix_path)), overwrite = TRUE)
  file.copy(sizes_path, file.path(table_dir, basename(sizes_path)), overwrite = TRUE)
  file.copy(exact_path, file.path(table_dir, basename(exact_path)), overwrite = TRUE)
  fwrite(intersections[, .(family = family, intersection_id, membership_pattern = pattern, n_genes)],
         file.path(table_dir, paste0("intersection_counts_", family, ".csv")), bom = TRUE)
  data.table(
    family = family, n_sets = length(flag_cols), n_union = union_n,
    n_all_set_intersection = all_set_n, n_intersections_displayed = nrow(intersections),
    set_bar_fill = "#F1C4BB", geometry = "UpSet",
    png = file.path("figures", basename(png_path)), pdf = file.path("figures", basename(pdf_path))
  )
}

audit <- rbindlist(lapply(families, make_one), fill = TRUE)
fwrite(audit, file.path(table_dir, "upset_render_audit.csv"), bom = TRUE)
cat("Wrote", nrow(audit), "teacher-selected GO UpSet figures to", figure_dir, "\n")
