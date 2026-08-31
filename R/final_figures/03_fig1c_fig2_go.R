#!/usr/bin/env Rscript

# Final display of the frozen 26 x 15 GO matrix and its 15 x 15 cross-tissue
# Spearman matrix. Only labels, typography, legend geometry and output size are
# changed; all numeric values and the approved order remain unchanged.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

set.seed(25)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/final_figures/helpers.R"))
input_dir <- file.path(root, "data/publication_input/go")

tissue <- fread(file.path(input_dir, "04_tissue_statistics_concrete_terms_and_context.csv"))
term_meta <- unique(tissue[, .(term_key, term_name, go_id, term_order, term_role)])
setorder(term_meta, term_order)
assert_final(nrow(term_meta) == 15L && !anyDuplicated(term_meta$term_key), "Expected 15 GO terms")
assert_final(nrow(tissue) == 390L && uniqueN(tissue$analysis_tissue) == 26L,
             "Expected the complete 26 x 15 GO matrix")

ddr_key <- "DNA_damage_response_GO0006974"
tissue_order <- tissue[term_key == ddr_key,
                       analysis_tissue[order(-mean_mission_nes, analysis_tissue)]]
term_meta[, display_name := final_go_display_name(term_name)]
term_meta[, axis_label := final_go_axis_label(display_name, go_id, width = 14L, max_lines = 3L)]
assert_final(term_meta[term_order == 4L, display_name] == "Intrinsic apoptotic signaling",
             "Fourth GO label was not updated")
assert_final(term_meta[term_order == 7L, display_name] == "Telomeric region",
             "Seventh GO label was not updated")

plot_data <- merge(
  copy(tissue), term_meta[, .(term_key, axis_label)], by = "term_key", sort = FALSE
)
plot_data[, analysis_tissue := factor(analysis_tissue, levels = rev(tissue_order))]
plot_data[, term_display := factor(axis_label, levels = term_meta$axis_label)]
nes_limit <- max(abs(plot_data$mean_mission_nes), na.rm = TRUE)
fig1c_palette <- c(low = "#2F5F9E", mid = "#F7F8F7", high = "#B94343")
nes_colours <- grDevices::colorRampPalette(c(
  fig1c_palette[["low"]], "#8FB1D3", fig1c_palette[["mid"]],
  "#E6A39A", fig1c_palette[["high"]]
))(21L)

fig1c <- ggplot(plot_data, aes(term_display, analysis_tissue, fill = mean_mission_nes)) +
  geom_tile(colour = "white", linewidth = 0.24) +
  scale_fill_gradientn(
    colours = nes_colours,
    values = seq(0, 1, length.out = length(nes_colours)),
    limits = c(-nes_limit, nes_limit), breaks = seq(-3, 3, 1),
    name = "Mission-equal mean NES",
    guide = guide_colorbar(
      title.position = "right",
      title.theme = element_text(
        angle = 90, hjust = 0.5, vjust = 0.5, size = 16.0,
        face = "bold", family = FINAL_FONT, colour = FINAL_COLOURS$ink
      ),
      barheight = unit(0.77, "npc"), barwidth = unit(0.55, "cm"),
      ticks.colour = FINAL_COLOURS$white, frame.colour = FINAL_COLOURS$muted
    )
  ) +
  labs(
    title = "Tissue-level GO enrichment across 26 mouse tissues",
    subtitle = "Mission-equal mean NES across 15 GO terms; Ensembl Gene IDs used for analysis",
    x = NULL, y = NULL,
    caption = "Tile colour encodes enrichment direction and magnitude; exact NES values and sign-flip inference are reported in the tables. NES is not pathway activation or inhibition."
  ) +
  theme_final(12.2) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 15.0,
                               lineheight = 0.90, margin = margin(t = 7)),
    axis.text.y = element_text(size = 15.0),
    axis.ticks = element_line(colour = "#9FB3C8", linewidth = 0.42),
    legend.position = "right",
    legend.text = element_text(size = 14.0),
    plot.title = element_text(size = 26),
    plot.subtitle = element_text(size = 16.0),
    plot.caption = element_text(size = 12.8),
    plot.margin = margin(10, 16, 10, 10)
  )

save_final_figure(fig1c, out, "Main/Fig_1c", width = 25.0, height = 16.5, dpi = 400)

rho_table <- fread(file.path(input_dir, "13_tissue_pathway_Spearman_rho_matrix.csv"))
rho <- as.matrix(rho_table[, -1L, with = FALSE])
storage.mode(rho) <- "numeric"
rownames(rho) <- rho_table$term_key
colnames(rho) <- names(rho_table)[-1L]
assert_final(all(is.finite(rho)) && max(abs(rho - t(rho))) < 1e-12,
             "Invalid Spearman matrix")

ordered_keys <- c(
  "cytoskeleton_GO0005856", "cell_death_GO0008219", "mechanosensation_GO0050954",
  "translation_GO0006412", "mitochondrion_GO0005739", "telomere_maintenance_GO0043247",
  "intrinsic_apoptotic_signaling_GO0008630", "DNA_templated_transcription_GO0006351",
  "DNA_damage_tolerance_GO0006301", "telomere_region_GO0000781",
  "DNA_damage_response_GO0006974", "DNA_repair_GO0006281", "DNA_replication_GO0006260",
  "DNA_damage_signal_transduction_GO0042770", "cell_cycle_GO0007049"
)
assert_final(setequal(ordered_keys, rownames(rho)), "Fig 2 approved order does not match the matrix")
labels <- term_meta[match(ordered_keys, term_key), axis_label]

rho_long <- rbindlist(lapply(seq_along(ordered_keys), function(i) {
  data.table(
    i = i,
    j = seq_along(ordered_keys),
    rho = rho[ordered_keys[[i]], ordered_keys]
  )
}))[i >= j]
rho_long[, `:=`(
  x = factor(j, levels = seq_along(ordered_keys), labels = labels),
  y = factor(i, levels = rev(seq_along(ordered_keys)), labels = rev(labels)),
  value_label = sprintf("%.2f", rho)
)]
assert_final(nrow(rho_long) == 120L, "Fig 2 lower triangle must contain 120 cells")

fig2_palette <- c(low = "#3B6FB6", mid = "#FFFFFF", high = "#C65A5A")
fig2 <- ggplot(rho_long, aes(x, y, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.30) +
  geom_text(aes(label = value_label), size = 4.05, colour = FINAL_COLOURS$ink,
            family = FINAL_FONT) +
  coord_fixed() +
  scale_fill_gradient2(
    low = fig2_palette[["low"]], mid = fig2_palette[["mid"]], high = fig2_palette[["high"]],
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho",
    guide = guide_colorbar(
      title.position = "right",
      title.theme = element_text(
        angle = 90, hjust = 0.5, vjust = 0.5, size = 15.5,
        face = "bold", family = FINAL_FONT, colour = FINAL_COLOURS$ink
      ),
      barheight = unit(0.55, "npc"), barwidth = unit(0.55, "cm")
    )
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Cross-tissue relationships among GO response terms",
    subtitle = "Lower triangle shown once; the approved average-linkage display order is retained",
    x = NULL, y = NULL,
    caption = "Spearman rho describes co-response across tissues, not causality or direct regulation."
  ) +
  theme_final(12.0) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 13.8,
                               lineheight = 0.88, margin = margin(t = 8)),
    axis.text.y = element_text(size = 14.5, lineheight = 0.90),
    legend.position = "right",
    legend.text = element_text(size = 13.5),
    plot.title = element_text(size = 26),
    plot.subtitle = element_text(size = 15.5),
    plot.caption = element_text(size = 12.5),
    plot.margin = margin(10, 18, 10, 10)
  )

save_final_figure(fig2, out, "Main/Fig_2", width = 24.0, height = 20.0, dpi = 400)

dir.create(file.path(out, "provenance"), recursive = TRUE, showWarnings = FALSE)
append_final_palette_audit(out, "Fig_1c", fig1c_palette)
append_final_palette_audit(out, "Fig_2", fig2_palette)
fwrite(term_meta[, .(term_order, term_key, source_term_name = term_name,
                     display_name, go_id, axis_label)],
       file.path(out, "provenance/GO_display_label_audit.csv"))
fwrite(data.table(figure = c("Fig_1c", "Fig_2"), n_cells = c(390L, 120L),
                  numeric_source = c("04_tissue_statistics_concrete_terms_and_context.csv",
                                     "13_tissue_pathway_Spearman_rho_matrix.csv")),
       file.path(out, "provenance/GO_figure_scope_audit.csv"))
message("FINAL_GO_FIGURES_PASS: Fig 1c = 390 cells; Fig 2 = 120 cells")
