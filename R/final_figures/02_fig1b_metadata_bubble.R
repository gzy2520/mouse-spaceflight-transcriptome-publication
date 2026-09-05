#!/usr/bin/env Rscript

# Final mouse-level experimental-design overview (Fig 1b).
# The frozen CSV was rebuilt directly from the 48 NASA OSDR ISA metadata ZIPs.
# No network call is required when reproducing the released output; the refresh
# script and source-file hashes are retained with the publication inputs.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

requested_root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (nzchar(requested_root)) {
  repo_root <- normalizePath(requested_root, mustWork = TRUE)
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script_arg)) stop("Cannot locate 02_fig1b_metadata_bubble.R", call. = FALSE)
  script_path <- sub("^--file=", "", script_arg[[1L]])
  # Some R front ends encode spaces in --file paths as ~+~.  Decode only when
  # the supplied path is absent and the decoded local path exists.
  decoded_path <- gsub("~[+]~", " ", script_path)
  if (!file.exists(script_path) && file.exists(decoded_path)) script_path <- decoded_path
  script_path <- normalizePath(script_path, mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
}
out <- Sys.getenv("PUBLICATION_OUTPUT_DIR", unset = file.path(repo_root, "results"))
if (!grepl("^/", out)) out <- file.path(repo_root, out)
out <- normalizePath(out, mustWork = FALSE)
dir.create(file.path(out, "Main"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "provenance"), recursive = TRUE, showWarnings = FALSE)
source(file.path(repo_root, "R/final_figures/helpers.R"))

input_path <- file.path(repo_root, "data/publication_input/mouse_metadata/03_mouse_level_metadata.csv")
go_path <- file.path(repo_root, "data/publication_input/go/04_tissue_statistics_concrete_terms_and_context.csv")
if (!file.exists(input_path)) stop("Missing frozen mouse metadata input: ", input_path, call. = FALSE)
if (!file.exists(go_path)) stop("Missing frozen GO tissue table: ", go_path, call. = FALSE)

mouse <- fread(input_path, na.strings = c("", "NA"), encoding = "UTF-8")
go <- fread(go_path, na.strings = c("", "NA"), encoding = "UTF-8")

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert(nrow(mouse) == 761L, "The frozen mouse metadata table must contain 761 sample rows")
assert(!anyDuplicated(mouse$sample_column), "Duplicate source sample columns in mouse metadata")
assert(uniqueN(mouse$accession) == 48L, "The refreshed metadata must contain 48 OSDR accessions")
assert(all(mouse$sex %chin% c("Female", "Male")), "Unexpected or missing sex metadata")
assert(all(!is.na(mouse$age_order_weeks)), "Unexpected or missing age metadata")
assert(all(!is.na(mouse$mission_cluster)), "Unexpected or missing canonical OSDR mission")
assert(all(mouse$isa_status_matches_analysis), "OSDR sample status conflicts with the analysis scope")

# The x-axis order is inherited from the first (DNA damage response) block of
# the cross-tissue GO table rather than alphabetized or re-derived here.
tissue_rows <- go[
  term_key == "DNA_damage_response_GO0006974" & term_order == 1,
  .(Group = analysis_tissue, tissue_order = seq_len(.N))
]
assert(nrow(tissue_rows) == 26L && !anyDuplicated(tissue_rows$Group),
       "The GO table does not provide the expected 26-tissue order")

mouse <- merge(mouse, tissue_rows, by = "Group", all.x = TRUE, sort = FALSE)
assert(!anyNA(mouse$tissue_order), "Mouse metadata contains a tissue outside the GO order")

# Chronological mission order keeps the paired female/male rows readable.  The
# explicit vector also makes the order stable if the input table is re-sorted.
mission_order <- c(
  "SpaceX-4", "SpaceX-8", "SpaceX-9", "SpaceX-11", "SpaceX-12",
  "SpaceX-14", "SpaceX-15", "SpaceX-16", "NG-11", "SpaceX-18",
  "SpaceX-21", "SpaceX-27"
)
observed_missions <- unique(mouse$mission_cluster)
assert(setequal(observed_missions, mission_order),
       "Mission set changed; update the explicit display order after auditing it")
sex_order <- c("Female", "Male")

row_grid <- data.table(
  mission_cluster = rep(mission_order, each = 2L),
  sex = rep(sex_order, times = length(mission_order))
)
row_grid[, row := .I]
row_grid[, sex_label := fifelse(sex == "Female", "♀ Female", "♂ Male")]

# One point is one mission x tissue x sex x age stratum.  Counts use unique
# accession-scoped ISA Source Name values.  This collapses technical replicates
# without assuming that a short source label is globally unique across OSDR
# accessions. Flight and ground mouse counts are retained for auditability.
mouse[, mouse_key := paste(accession, mouse_id, sep = "::")]
grouped <- mouse[
  , .(
    n_samples = uniqueN(sample_column),
    n_mice = uniqueN(mouse_key),
    n_technical_replicates = sum(sample_is_technical_replicate),
    n_flight = uniqueN(mouse_key[sample_status == "Flight"]),
    n_ground = uniqueN(mouse_key[sample_status == "Ground"]),
    n_accessions = uniqueN(accession),
    accessions = paste(sort(unique(accession)), collapse = ";"),
    n_analysis_units = uniqueN(analysis_unit_id),
    analysis_units = paste(sort(unique(analysis_unit_id)), collapse = ";"),
    source_age_values = paste(sort(unique(paste0(age_value_raw, " {", age_unit_raw, "}"))), collapse = ";"),
    age_source_fields = paste(sort(unique(age_source_field)), collapse = ";"),
    age_is_range = any(age_is_range)
  ),
  by = .(mission_cluster, Group, sex, age_label, age_order_weeks, tissue_order)
]
setorder(grouped, mission_cluster, tissue_order, sex, age_order_weeks, age_label)
grouped[, age_rank := seq_len(.N), by = .(mission_cluster, Group, sex)]
grouped[, n_age_levels := .N, by = .(mission_cluster, Group, sex)]
grouped[, row := row_grid$row[match(paste(mission_cluster, sex), paste(row_grid$mission_cluster, row_grid$sex))]]
grouped[, x := tissue_order + fifelse(
  n_age_levels > 1,
  (age_rank - (n_age_levels + 1) / 2) * 0.28,
  0
)]
assert(!anyNA(grouped$row), "Failed to assign a mission-sex row")
assert(!anyDuplicated(grouped[, .(mission_cluster, Group, sex, age_label)]),
       "Duplicate plot strata after aggregation")

mission_centres <- row_grid[, .(row = mean(row)), by = mission_cluster]
mission_centres[, band := (seq_len(.N) %% 2L) + 1L]
mission_bands <- row_grid[, .(
  xmin = 0.75, xmax = 26.5, ymin = min(row) - 0.5, ymax = max(row) + 0.5
), by = mission_cluster]
mission_bands <- merge(mission_bands, mission_centres[, .(mission_cluster, band)], by = "mission_cluster")

age_levels <- unique(grouped[, .(age_label, age_order_weeks)])
setorder(age_levels, age_order_weeks, age_label)
age_levels <- age_levels$age_label
# Reference Fig. 1b age palette supplied for the final release.  These exact
# chronological values retain the intended yellow -> green -> blue -> purple
# progression while keeping age as a discrete source variable.
age_palette <- c(
  "#FDE333", "#D4E02D", "#A6DA42", "#73D25B", "#25C771",
  "#00BA82", "#00AC8E", "#009B95", "#008A98", "#007796",
  "#006290", "#1E4D85", "#3C3777", "#471D67", "#4B0055"
)
assert(length(age_levels) <= length(age_palette), "Add an audited discrete colour before adding another age level")
age_palette <- age_palette[seq_along(age_levels)]
names(age_palette) <- age_levels
grouped[, age_colour := unname(age_palette[age_label])]
assert(!anyNA(grouped$age_colour), "Every age stratum must have a fixed colour")
append_final_palette_audit(
  out, "Fig_1b", setNames(unname(age_palette), sprintf("age_%02d", seq_along(age_palette)))
)
fwrite(
  grouped[, .(
    mission_cluster, Group, tissue_order, sex, age_label, age_order_weeks,
    age_is_range, source_age_values, age_source_fields, age_colour, n_age_levels,
    n_samples, n_mice, n_technical_replicates, n_flight, n_ground,
    n_accessions, accessions, n_analysis_units, analysis_units
  )],
  file.path(out, "provenance/Fig_1b_mouse_metadata_grouped.csv")
)

wrap_label <- function(value, width = 16L) {
  vapply(strwrap(value, width = width, simplify = FALSE), paste, collapse = "\n", FUN.VALUE = character(1L))
}
tissue_labels <- wrap_label(tissue_rows$Group, width = 15L)
names(tissue_labels) <- as.character(tissue_rows$tissue_order)

max_mice <- max(grouped$n_mice)
size_breaks <- pretty(c(1, max_mice), n = 4)
size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= max_mice]

ink <- "#263442"
muted <- "#667585"
grid_col <- "#E4EAF0"
band_cols <- c("#FFFFFF", "#F6F8FA")
mission_bands[, band_colour := band_cols[band]]

fig <- ggplot() +
  geom_rect(
    data = mission_bands,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = mission_bands$band_colour, inherit.aes = FALSE, colour = NA
  ) +
  geom_vline(xintercept = seq(0.5, 26.5, by = 1), colour = "#EEF2F5", linewidth = 0.25) +
  geom_hline(
    yintercept = seq(2.5, max(row_grid$row) - 0.5, by = 2),
    colour = grid_col, linewidth = 0.55
  ) +
  geom_vline(xintercept = 0.72, colour = "#B8C4CF", linewidth = 0.55) +
  geom_point(
    data = grouped,
    aes(x = x, y = row, fill = age_label, size = n_mice),
    shape = 21, colour = "white", stroke = 0.45, alpha = 0.96
  ) +
  geom_text(
    data = row_grid,
    aes(x = 0.55, y = row, label = sex_label),
    hjust = 1, size = 5.35, colour = ink, family = "DejaVu Sans"
  ) +
  geom_text(
    data = mission_centres,
    aes(x = -1.05, y = row, label = mission_cluster),
    hjust = 1, size = 5.55, fontface = "bold", colour = ink, family = "DejaVu Sans"
  ) +
  scale_x_continuous(
    breaks = tissue_rows$tissue_order,
    labels = tissue_labels,
    limits = c(-2.8, 26.6),
    expand = c(0, 0)
  ) +
  scale_y_reverse(
    breaks = row_grid$row,
    labels = NULL,
    limits = c(max(row_grid$row) + 0.55, 0.45),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = age_palette,
    breaks = age_levels,
    drop = FALSE,
    name = "Age reported by OSDR (discrete)"
  ) +
  scale_size_area(
    max_size = 14,
    breaks = size_breaks,
    limits = c(0, max_mice),
    name = "Unique mice"
  ) +
  guides(
    fill = guide_legend(
      order = 1, nrow = 2, byrow = TRUE,
      override.aes = list(shape = 21, size = 5.5, colour = "white", alpha = 1)
    ),
    size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#718096", colour = "white"))
  ) +
  labs(
    title = "Mouse age and cohort size across missions and tissues",
    subtitle = "NASA OSDR ISA metadata refreshed for 48 datasets. Each row is a sex-specific mission stratum; point area shows accession-scoped unique mice across selected flight and ground samples.",
    caption = "Tissue order follows the cross-tissue DNA-damage-response GO table. Ages are discrete source values; ranges are retained only where the downloaded OSDR ISA table reports a range. Multiple age strata within one cell are horizontally offset.",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = "DejaVu Sans", base_size = 16.0) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.text.x = element_text(colour = ink, size = 15.0, angle = 55, hjust = 1, vjust = 1, lineheight = 0.95),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(colour = ink, face = "bold", size = 27, hjust = 0, margin = margin(b = 5)),
    plot.subtitle = element_text(colour = muted, size = 16.0, hjust = 0, margin = margin(b = 10)),
    plot.caption = element_text(colour = muted, size = 13.0, hjust = 0, lineheight = 1.1, margin = margin(t = 10)),
    legend.position = "top",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(colour = ink, face = "bold", size = 22),
    legend.text = element_text(colour = ink, size = 20),
    legend.key = element_blank(),
    legend.margin = margin(0, 0, 2, 0),
    plot.margin = margin(10, 14, 12, 14)
  )

png_path <- file.path(out, "Main/Fig_1b.png")
pdf_path <- file.path(out, "Main/Fig_1b.pdf")
ggsave(png_path, fig, width = 27, height = 17, units = "in", dpi = 400, bg = "white", limitsize = FALSE)
ggsave(pdf_path, fig, width = 27, height = 17, units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)

message("FINAL_FIG1B_PASS: ", png_path, " and ", pdf_path)
