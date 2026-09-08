#!/usr/bin/env Rscript

# Shared helpers for the C-story publication variant.  These helpers only
# affect display: all values, filters, identifiers and ordering decisions are
# supplied by the frozen publication_input/ tables.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(grid)
})

STYLE_C_COLOURS <- list(
  ink = "#263442",
  muted = "#667585",
  grid = "#E3E8ED",
  panel = "#F4F6F8",
  blue = "#3B6FB6",
  blue_dark = "#234E7A",
  red = "#C65A5A",
  red_dark = "#8E3E43",
  teal = "#3E8F8C",
  gold = "#D49A36",
  grey = "#D6DCE2",
  grey_dark = "#9AA6B2",
  white = "#FFFFFF"
)

STYLE_C_FONT <- "DejaVu Sans"

theme_style_C <- function(base_size = 9, base_family = STYLE_C_FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = STYLE_C_COLOURS$ink, linewidth = 0.35),
      axis.ticks = element_line(colour = STYLE_C_COLOURS$muted, linewidth = 0.3),
      axis.text = element_text(colour = STYLE_C_COLOURS$ink),
      axis.title = element_text(colour = STYLE_C_COLOURS$ink, face = "bold"),
      plot.title = element_text(
        colour = STYLE_C_COLOURS$ink, face = "bold", size = rel(1.35),
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        colour = STYLE_C_COLOURS$muted, size = rel(0.92), margin = margin(b = 6)
      ),
      plot.caption = element_text(
        colour = STYLE_C_COLOURS$muted, size = rel(0.76), hjust = 0
      ),
      legend.position = "top",
      legend.title = element_text(colour = STYLE_C_COLOURS$ink, face = "bold"),
      legend.text = element_text(colour = STYLE_C_COLOURS$ink),
      legend.key = element_blank(),
      strip.background = element_rect(fill = STYLE_C_COLOURS$panel, colour = NA),
      strip.text = element_text(colour = STYLE_C_COLOURS$ink, face = "bold"),
      plot.background = element_rect(fill = STYLE_C_COLOURS$white, colour = NA),
      panel.background = element_rect(fill = STYLE_C_COLOURS$white, colour = NA)
    )
}

save_style_C <- function(plot, out_dir, relative_stem, width, height, dpi = 600) {
  png_path <- file.path(out_dir, paste0(relative_stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(relative_stem, ".pdf"))
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    png_path, plot, width = width, height = height, units = "in", dpi = dpi,
    bg = "white", limitsize = FALSE
  )
  ggsave(
    pdf_path, plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
  )
  invisible(c(png = png_path, pdf = pdf_path))
}

message("C-story figure helpers loaded")
