# Shared publication style switch. The default keeps the approved renderer
# byte-for-byte compatible; FIGURE_STYLE=A enables the conservative visual
# refresh used by workflow/run_style_A.R.
figure_style_a <- identical(Sys.getenv("FIGURE_STYLE", unset = "approved"), "A")

FIGURE_COLOURS <- list(
  low = if (figure_style_a) "#3B6FB6" else "#2166AC",
  mid = if (figure_style_a) "#F7F7F5" else "#FFFFFF",
  high = if (figure_style_a) "#C65A5A" else "#B2182B",
  blue = if (figure_style_a) "#3B6FB6" else "#2166AC",
  red = if (figure_style_a) "#C65A5A" else "#B2182B",
  ink = if (figure_style_a) "#27364A" else "#2F2F2F",
  muted = if (figure_style_a) "#667085" else "#4B5563",
  grid = if (figure_style_a) "#E5E7EB" else "#E5E7EB",
  panel = if (figure_style_a) "#F3F5F7" else "#F0F0F0",
  ground = if (figure_style_a) "#3B6FB6" else "#0072B2",
  flight = if (figure_style_a) "#C65A5A" else "#D55E00"
)

style_title <- function(size = 14, hjust = 0, family = "Arial") {
  element_text(face = "bold", size = size, hjust = hjust, colour = FIGURE_COLOURS$ink, family = family)
}

style_subtitle <- function(size = 9, hjust = 0, family = "Arial") {
  element_text(size = size, hjust = hjust, colour = FIGURE_COLOURS$muted, family = family)
}

style_caption <- function(size = 7.5, hjust = 0, family = "Arial") {
  element_text(size = size, hjust = hjust, colour = FIGURE_COLOURS$muted, family = family)
}
