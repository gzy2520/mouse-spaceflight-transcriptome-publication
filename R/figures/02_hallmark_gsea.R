#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
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
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#4B5563") +
  scale_fill_manual(values = c("Positive NES" = "#B2182B", "Negative NES" = "#2166AC")) +
  labs(
    title = "Hallmark: descriptively strongest global GSEA terms",
    subtitle = "Six positive and six negative terms with complete 13-mission coverage; mission-equal mean NES",
    x = "Mean mission-level NES (positive = higher in flight/uG)", y = NULL, fill = NULL,
    caption = "Descriptive summary across heterogeneous tissues; not a pooled effect estimate or a significance test."
  ) +
  theme_minimal(base_family = "Arial", base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8), legend.position = "top",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 8, colour = "#4B5563"),
    plot.margin = margin(8, 16, 8, 8)
  )
ggsave(
  file.path(out, "Suppl/Fig_S1.png"), plot,
  width = 10, height = 7.2, dpi = 320, bg = "white"
)
