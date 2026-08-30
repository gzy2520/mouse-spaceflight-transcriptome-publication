#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/figures/style_helpers.R"))
selected <- fread(file.path(
  root, "data/publication_input/hallmark/01_global_top_term_figure_selection.csv"
))[library == "MSigDB_2026.1.Mm_Hallmark"]
stopifnot(nrow(selected) == 12L, all(selected$n_missions == 13L))
selected[, term_plot := factor(
  term_display,
  levels = selected[order(mean_of_mission_mean_NES), term_display]
)]

plot <- ggplot(
  selected,
  aes(x = mean_of_mission_mean_NES, y = term_plot, fill = direction)
) +
  geom_col(width = if (figure_style_a) 0.70 else 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#4B5563") +
  scale_fill_manual(values = c("Positive NES" = FIGURE_COLOURS$red, "Negative NES" = FIGURE_COLOURS$blue)) +
  labs(
    title = if (figure_style_a) "Hallmark: strongest global GSEA terms" else "Hallmark: descriptively strongest global GSEA terms",
    subtitle = if (figure_style_a) "Six positive and six negative terms; complete 13-mission coverage" else "Six positive and six negative terms with complete 13-mission coverage; mission-equal mean NES",
    x = "Mean mission-level NES (positive = higher in flight/uG)", y = NULL, fill = NULL,
    caption = "Descriptive summary across heterogeneous tissues; not a pooled effect estimate or a significance test."
  ) +
  if (figure_style_a) {
    theme_minimal(base_family = "Arial", base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8.5, colour = FIGURE_COLOURS$ink), legend.position = "top",
        plot.title = style_title(13), plot.subtitle = style_subtitle(8.5), plot.caption = style_caption(7.3),
        plot.margin = margin(8, 16, 8, 8)
      )
  } else {
    theme_minimal(base_family = "Arial", base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8), legend.position = "top",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9),
        plot.caption = element_text(size = 8, colour = "#4B5563"),
        plot.margin = margin(8, 16, 8, 8)
      )
  }
ggsave(
  file.path(out, "Suppl/Fig_S1.png"), plot,
  width = 10, height = 7.2, dpi = if (figure_style_a) 600 else 320, bg = "white"
)
