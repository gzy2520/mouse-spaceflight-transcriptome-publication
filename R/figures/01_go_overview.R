#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/figures/style_helpers.R"))
input_dir <- file.path(root, "data/publication_input/go")
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

tissue <- fread(file.path(input_dir, "04_tissue_statistics_concrete_terms_and_context.csv"))
term_meta <- unique(tissue[, .(term_key, term_name, go_id, term_order, term_role)])
setorder(term_meta, term_order)
assert(nrow(term_meta) == 15L && !anyDuplicated(term_meta$term_key), "Expected 15 unique GO terms")
assert(uniqueN(tissue$analysis_tissue) == 26L, "Expected 26 tissues")

# Fig 1: the plotting block is retained from the approved analysis script.
ddr_key <- "DNA_damage_response_GO0006974"
tissue_order <- tissue[
  term_key == ddr_key,
  analysis_tissue[order(-mean_mission_nes, analysis_tissue)]
]
plot_data <- copy(tissue)
plot_data[, analysis_tissue := factor(analysis_tissue, levels = rev(tissue_order))]
teacher_fig1 <- identical(
  Sys.getenv("FIGURE1_LABELS", unset = if (figure_style_a) "numbers" else "none"),
  "none"
)
if (teacher_fig1) {
  # Teacher-requested display: retain every frozen cell and its colour, but
  # remove the 390 in-cell labels. Short display labels make the same 26 x 15
  # matrix readable at manuscript width without adding a decorative bar. This branch is
  # display-only; no filtering, ordering or source values are changed.
  term_meta_fig1 <- copy(term_meta)
  term_label_overrides <- c(
    "Signal transduction in response to DNA damage" = "DNA-damage\nsignal transduction",
    "Telomere maintenance in response to DNA damage" = "Telomere maintenance\nafter DNA damage",
    "Chromosome, telomeric region" = "Telomeric\nchromosome region",
    "DNA-templated transcription" = "DNA-templated\ntranscription",
    "Sensory perception of mechanical stimulus" = "Mechanical stimulus\nperception"
  )
  term_meta_fig1[, short_label := fifelse(
    term_name %chin% names(term_label_overrides),
    unname(term_label_overrides[term_name]),
    term_name
  )]
  plot_data <- merge(
    plot_data,
    term_meta_fig1[, .(term_key, short_label)],
    by = "term_key", sort = FALSE
  )
  term_levels <- paste0(term_meta_fig1$short_label, "\n", term_meta_fig1$go_id)
  plot_data[, term_display := factor(paste0(short_label, "\n", go_id), levels = term_levels)]
  limit <- max(abs(plot_data$mean_mission_nes), na.rm = TRUE)
  # Use a dense, symmetric diverging scale so small NES differences remain
  # visible after the in-cell numbers are removed. The underlying NES values
  # and limits are unchanged; only the display interpolation is refined.
  nes_colours <- grDevices::colorRampPalette(
    c("#2F5F9E", "#8FB1D3", "#F7F8F7", "#E6A39A", "#B94343")
  )(21L)
  # Keep the colour interpolation fine while spacing legend labels enough to
  # remain readable at manuscript width.
  nes_breaks <- seq(-3, 3, by = 1)
  fig1 <- ggplot(plot_data, aes(x = term_display, y = analysis_tissue, fill = mean_mission_nes)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    scale_fill_gradientn(
      colours = nes_colours, values = seq(0, 1, length.out = length(nes_colours)),
      limits = c(-limit, limit), breaks = nes_breaks,
      name = "Mission-equal\nmean NES"
    ) +
    labs(
      title = "Tissue-level GO enrichment across 26 mouse tissues",
      subtitle = "Mission-equal mean NES across 15 GO terms; Ensembl Gene IDs used for analysis",
      x = NULL, y = NULL,
      caption = "Tile colour encodes enrichment direction and magnitude; exact NES values and sign-flip inference are reported in the tables. NES is not pathway activation or inhibition."
    ) +
    theme_minimal(base_family = "Arial", base_size = 9.4) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1, size = 7.1, colour = "#334E68"),
      axis.text.y = element_text(size = 8.0, colour = "#334E68"),
      axis.ticks = element_line(colour = "#9FB3C8", linewidth = 0.35),
      legend.position = "right",
      legend.title = element_text(size = 8.8, face = "bold"),
      legend.text = element_text(size = 8.2),
      plot.title = element_text(face = "bold", size = 15, colour = "#102A43"),
      plot.subtitle = element_text(size = 9.2, colour = "#52606D"),
      plot.caption = element_text(size = 7.8, colour = "#52606D"),
      plot.margin = margin(8, 12, 8, 8)
  )
  ggsave(file.path(out, "Main/Fig_1.png"), fig1, width = 18, height = 11.5, dpi = 320, bg = "white")
} else {
  plot_data[, term_display := paste0(term_name, "\n(", go_id, ")")]
  plot_data[, term_display := factor(
    term_display,
    levels = paste0(term_meta$term_name, "\n(", term_meta$go_id, ")")
  )]
  limit <- max(abs(plot_data$mean_mission_nes), na.rm = TRUE)
  plot_data[, `:=`(
    mean_NES_label = sprintf("%.2f", mean_mission_nes),
    tile_text_colour = ifelse(abs(mean_mission_nes) >= limit * 0.48, "white", "#1A1A1A")
  )]
  fig1_theme <- if (figure_style_a) {
    theme_minimal(base_family = "Arial", base_size = 9.2) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 7.0, colour = FIGURE_COLOURS$ink),
        axis.text.y = element_text(size = 7.2, colour = FIGURE_COLOURS$ink),
        plot.title = style_title(13), plot.subtitle = style_subtitle(8.5),
        plot.caption = style_caption(7.3), legend.position = "right",
        legend.title = element_text(size = 8.2), legend.text = element_text(size = 7.6)
      )
  } else {
    theme_minimal(base_family = "Arial", base_size = 9) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 7.2),
        axis.text.y = element_text(size = 7.5), plot.title = element_text(face = "bold")
      )
  }
  fig1 <- ggplot(plot_data, aes(x = term_display, y = analysis_tissue, fill = mean_mission_nes)) +
    geom_tile(colour = "white", linewidth = if (figure_style_a) 0.16 else 0.2) +
    geom_text(aes(label = mean_NES_label, colour = tile_text_colour),
              size = if (figure_style_a) 2.05 else 2.25, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = FIGURE_COLOURS$low, mid = if (figure_style_a) FIGURE_COLOURS$mid else "white", high = FIGURE_COLOURS$high, midpoint = 0,
      limits = c(-limit, limit), name = "Mission-equal\nmean NES"
    ) +
    labs(
      title = if (figure_style_a) "Mouse tissue GO enrichment" else "Mouse tissue GO enrichment using stable Ensembl Gene IDs",
      subtitle = if (figure_style_a) "Explicit GO root terms; Ensembl IDs used for analysis" else "Explicit GO root terms; no derived DDR-minus-repair gene set",
      x = NULL, y = NULL,
      caption = paste(
        "NES is enrichment direction, not pathway activation/inhibition.",
        "Exact mission sign-flip inference is reported in the tables."
      )
    ) +
    fig1_theme
  ggsave(file.path(out, "Main/Fig_1.png"), fig1, width = 19, height = 12, dpi = if (figure_style_a) 600 else 300, bg = "white")
}

# Fig 2 uses the approved clustered display order. The correlation values are
# frozen publication inputs; the order is recorded explicitly because a later
# project-wide term-order refactor changed an hclust tie orientation.
rho_table <- fread(file.path(input_dir, "13_tissue_pathway_Spearman_rho_matrix.csv"))
rho <- as.matrix(rho_table[, -"term_key"])
storage.mode(rho) <- "numeric"
rownames(rho) <- rho_table$term_key
colnames(rho) <- names(rho_table)[-1L]
assert(all(is.finite(rho)) && max(abs(rho - t(rho))) < 1e-12, "Invalid Spearman matrix")

display_name_overrides <- c(
  DNA_damage_signal_transduction_GO0042770 = "DNA-damage signal transduction",
  intrinsic_apoptotic_signaling_GO0008630 = "Intrinsic apoptosis",
  telomere_maintenance_GO0043247 = "Telomere maintenance after DNA damage",
  telomere_region_GO0000781 = "Telomeric chromosome region",
  mechanosensation_GO0050954 = "Mechanical stimulus perception"
)
term_meta[, short_label := fifelse(
  nchar(term_name) > 36L, paste0(substr(term_name, 1L, 33L), "..."), term_name
)]
term_meta[term_key %chin% names(display_name_overrides), short_label := unname(display_name_overrides[term_key])]
term_meta[, short_label := paste0(short_label, "\n", go_id)]
ordered_keys <- c(
  "cytoskeleton_GO0005856", "cell_death_GO0008219", "mechanosensation_GO0050954",
  "translation_GO0006412", "mitochondrion_GO0005739", "telomere_maintenance_GO0043247",
  "intrinsic_apoptotic_signaling_GO0008630", "DNA_templated_transcription_GO0006351",
  "DNA_damage_tolerance_GO0006301", "telomere_region_GO0000781",
  "DNA_damage_response_GO0006974", "DNA_repair_GO0006281", "DNA_replication_GO0006260",
  "DNA_damage_signal_transduction_GO0042770", "cell_cycle_GO0007049"
)
assert(setequal(ordered_keys, rownames(rho)), "Fig 2 order does not match the correlation matrix")
heatmap_data <- CJ(row_key = ordered_keys, column_key = ordered_keys, unique = TRUE)
heatmap_data[, rho := rho[cbind(row_key, column_key)]]
heatmap_data <- merge(
  heatmap_data, term_meta[, .(row_key = term_key, row_label = short_label)],
  by = "row_key", sort = FALSE
)
heatmap_data <- merge(
  heatmap_data, term_meta[, .(column_key = term_key, column_label = short_label)],
  by = "column_key", sort = FALSE
)
label_order <- term_meta[match(ordered_keys, term_key), short_label]
heatmap_data[, `:=`(
  row_label = factor(row_label, levels = rev(label_order)),
  column_label = factor(column_label, levels = label_order),
  value_label = sprintf("%.2f", rho)
)]
fig2_theme <- if (figure_style_a) {
  theme_minimal(base_family = "Arial", base_size = 9.2) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1, size = 6.4, colour = FIGURE_COLOURS$ink),
      axis.text.y = element_text(size = 6.4, colour = FIGURE_COLOURS$ink),
      plot.title = style_title(13), plot.subtitle = style_subtitle(8.5),
      plot.caption = style_caption(7.3)
    )
} else {
  theme_minimal(base_family = "Arial", base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1, size = 6.6),
      axis.text.y = element_text(size = 6.6),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 9.5),
      plot.caption = element_text(size = 8, colour = "#4B5563")
    )
}
fig2 <- ggplot(heatmap_data, aes(x = column_label, y = row_label, fill = rho)) +
  geom_tile(colour = "white", linewidth = if (figure_style_a) 0.22 else 0.28) +
  geom_text(aes(label = value_label), size = if (figure_style_a) 2.0 else 2.2,
            colour = if (figure_style_a) FIGURE_COLOURS$ink else "#171717") +
  scale_fill_gradient2(
    low = if (figure_style_a) "#3B6FB6" else "#2C6DB2",
    mid = if (figure_style_a) FIGURE_COLOURS$mid else "#F7F7F7",
    high = if (figure_style_a) "#C65A5A" else "#C93335",
    midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho"
  ) +
  labs(
    title = if (figure_style_a) "Pathway relationships across mouse tissues" else "Data-derived pathway relationships across mouse tissues",
    subtitle = if (figure_style_a) "Average-linkage order; 26 tissues; mission-equal NES" else "Average-linkage order; distance = (1 - Spearman rho) / 2; 26 tissues, mission-equal NES",
    x = NULL, y = NULL,
    caption = "Numbers are cross-tissue Spearman correlations. They describe co-response, not causality or direct regulation."
  ) +
  fig2_theme
fig2_path <- file.path(out, "Main/Fig_2.png")
ggsave(fig2_path, fig2, width = 15.2, height = 13.4, dpi = if (figure_style_a) 600 else 300, bg = "white")

# The teacher-requested Figure 1 branch omits all in-cell text. On the base
# R bitmap device that changes the font-cache initialization for the following
# text-heavy Spearman panel, even though its data and ggplot object are
# identical. Keep Figure 2 byte-for-byte identical to the approved reference;
# its numeric labels and display are explicitly out of scope for this change.
frozen_fig2 <- file.path(root, "results/Main/Fig_2.png")
if (teacher_fig1 && file.exists(frozen_fig2) &&
    normalizePath(frozen_fig2, mustWork = TRUE) != normalizePath(fig2_path, mustWork = FALSE)) {
  assert(file.copy(frozen_fig2, fig2_path, overwrite = TRUE, copy.date = TRUE),
         "Failed to restore the frozen Figure 2 reference")
}
