#!/usr/bin/env Rscript

# Display-only helpers for the final manuscript figures. Analysis inputs,
# stable identifiers, filters and numeric matrices remain in publication_input.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

FINAL_FONT <- {
  requested <- Sys.getenv("FIGURE_FONT", unset = "Arial")
  if (requested %in% c("Arial", "Helvetica", "DejaVu Sans")) requested else "Arial"
}

FINAL_COLOURS <- list(
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
  white = "#FFFFFF"
)

# Compatibility aliases used by the retained approved plotting blocks.
FIGURE_COLOURS <- c(FINAL_COLOURS, list(
  low = FINAL_COLOURS$blue,
  mid = "#F7F7F5",
  high = FINAL_COLOURS$red,
  ground = FINAL_COLOURS$blue,
  flight = FINAL_COLOURS$red
))

style_title <- function(size = 14, hjust = 0, family = FINAL_FONT) {
  element_text(face = "bold", size = size, hjust = hjust,
               colour = FINAL_COLOURS$ink, family = family)
}

style_subtitle <- function(size = 9, hjust = 0, family = FINAL_FONT) {
  element_text(size = size, hjust = hjust, colour = FINAL_COLOURS$muted, family = family)
}

style_caption <- function(size = 7.5, hjust = 0, family = FINAL_FONT) {
  element_text(size = size, hjust = hjust, colour = FINAL_COLOURS$muted, family = family)
}

assert_final <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

append_final_palette_audit <- function(out_dir, figure_id, palette) {
  assert_final(length(palette) > 0L && !is.null(names(palette)),
               paste("Palette roles are missing for", figure_id))
  roles <- names(palette)
  palette <- toupper(as.character(palette))
  assert_final(all(grepl("^#[0-9A-F]{6}$", palette)),
               paste("Invalid palette value for", figure_id))
  audit_path <- file.path(out_dir, "provenance/figure_palette_audit.csv")
  dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
  current <- if (file.exists(audit_path)) data.table::fread(audit_path) else
    data.table::data.table(figure = character(), role = character(), hex = character())
  current <- current[!(figure == figure_id & role %in% roles)]
  added <- data.table::data.table(figure = figure_id, role = roles, hex = palette)
  answer <- data.table::rbindlist(list(current, added), use.names = TRUE)
  data.table::setorder(answer, figure, role)
  data.table::fwrite(answer, audit_path)
  invisible(answer)
}

theme_final <- function(base_size = 12, base_family = FINAL_FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = FINAL_COLOURS$ink, linewidth = 0.42),
      axis.ticks = element_line(colour = FINAL_COLOURS$muted, linewidth = 0.36),
      axis.text = element_text(colour = FINAL_COLOURS$ink),
      axis.title = element_text(colour = FINAL_COLOURS$ink, face = "bold"),
      plot.title = element_text(
        colour = FINAL_COLOURS$ink, face = "bold", size = rel(1.45),
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        colour = FINAL_COLOURS$muted, size = rel(0.98), margin = margin(b = 8)
      ),
      plot.caption = element_text(
        colour = FINAL_COLOURS$muted, size = rel(0.78), hjust = 0,
        margin = margin(t = 7)
      ),
      legend.title = element_text(colour = FINAL_COLOURS$ink, face = "bold"),
      legend.text = element_text(colour = FINAL_COLOURS$ink),
      legend.key = element_blank(),
      strip.background = element_rect(fill = FINAL_COLOURS$panel, colour = NA),
      strip.text = element_text(colour = FINAL_COLOURS$ink, face = "bold"),
      plot.background = element_rect(fill = FINAL_COLOURS$white, colour = NA),
      panel.background = element_rect(fill = FINAL_COLOURS$white, colour = NA)
    )
}

save_final_figure <- function(plot, out_dir, relative_stem, width, height, dpi = 400, svg = FALSE) {
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
  if (isTRUE(svg)) {
    ggsave(
      file.path(out_dir, paste0(relative_stem, ".svg")), plot,
      width = width, height = height, units = "in", device = grDevices::svg,
      bg = "white", limitsize = FALSE
    )
  }
  invisible(c(png = png_path, pdf = pdf_path))
}

final_go_display_name <- function(term_name) {
  answer <- as.character(term_name)
  answer[answer == "Intrinsic apoptotic signaling pathway"] <- "Intrinsic apoptotic signaling"
  answer[answer == "Chromosome, telomeric region"] <- "Telomeric region"
  answer
}

balanced_wrap <- function(label, width = 18L, max_lines = 2L) {
  words <- strsplit(gsub("\\s+", " ", trimws(label)), " ")[[1L]]
  if (length(words) <= 1L) return(label)
  if (max_lines == 2L) {
    candidates <- lapply(seq_len(length(words) - 1L), function(k) {
      c(paste(words[seq_len(k)], collapse = " "),
        paste(words[(k + 1L):length(words)], collapse = " "))
    })
    scores <- vapply(candidates, function(x) {
      overflow <- sum(pmax(nchar(x) - width, 0))
      overflow * 100 + max(nchar(x)) + abs(diff(nchar(x))) / 10
    }, numeric(1L))
    return(paste(candidates[[which.min(scores)]], collapse = "\n"))
  }
  paste(strwrap(label, width = width), collapse = "\n")
}

final_go_axis_label <- function(term_name, go_id, width = 18L, max_lines = 2L) {
  vapply(seq_along(term_name), function(i) {
    paste0(
      balanced_wrap(final_go_display_name(term_name[[i]]), width, max_lines = max_lines),
      "\n", go_id[[i]]
    )
  }, character(1L))
}

short_tissue_name_final <- function(x) {
  x <- as.character(x)
  x[x == "Heart / Heart right ventricle"] <- "Heart right ventricle"
  x
}
