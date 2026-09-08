#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/figures/style_C_helpers.R"))
selected <- fread(file.path(
  root, "data/publication_input/hallmark/01_global_top_term_figure_selection.csv"
))[library == "MSigDB_2026.1.Mm_Hallmark"]
stopifnot(nrow(selected) == 12L, all(selected$n_missions == 13L))
selected[, term_plot := factor(
  term_display,
  levels = selected[order(mean_of_mission_mean_NES), term_display]
)]

{
  selected_c <- copy(selected)
  selected_c[, direction := factor(direction, levels = c("Negative NES", "Positive NES"))]
  fig_s1_c <- ggplot(selected_c, aes(mean_of_mission_mean_NES, term_plot, colour = direction)) +
    geom_segment(
      aes(x = 0, xend = mean_of_mission_mean_NES, yend = term_plot),
      linewidth = 1.15, lineend = "round", colour = STYLE_C_COLOURS$grey_dark
    ) +
    geom_point(size = 3.4) +
    geom_text(
      aes(
        label = sprintf("%+.2f", mean_of_mission_mean_NES),
        hjust = ifelse(mean_of_mission_mean_NES >= 0, -0.35, 1.35)
      ), size = 2.9, colour = STYLE_C_COLOURS$ink
    ) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = STYLE_C_COLOURS$ink) +
    scale_colour_manual(
      values = c("Negative NES" = STYLE_C_COLOURS$blue, "Positive NES" = STYLE_C_COLOURS$red),
      name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.16))) +
    labs(
      title = "Hallmark terms with the strongest global NES",
      subtitle = "Descriptive mean across 13 missions; six terms in each direction",
      x = "Mean mission-level NES", y = NULL,
      caption = "This is a descriptive summary across heterogeneous tissues, not a pooled effect estimate."
    ) +
    theme_style_C(9.4) +
    theme(axis.text.y = element_text(size = 8), plot.margin = margin(8, 14, 8, 8))
  save_style_C(fig_s1_c, out, "Suppl/Fig_S1", 10.8, 7.0)
}

{
  source(file.path(root, "R/final_figures/helpers.R"))
  append_final_palette_audit(out, "Fig_S1", c(
    negative = STYLE_C_COLOURS$blue,
    positive = STYLE_C_COLOURS$red,
    segment = STYLE_C_COLOURS$grey_dark
  ))
}
