#!/usr/bin/env Rscript

# Final dot-matrix display of the frozen tissue meta-log2FC matrices.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(grid)
})

set.seed(25)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/final_figures/helpers.R"))
meta_palette <- c(low = "#3B6FB6", mid = "#FFFFFF", high = "#C65A5A")

vertical_meta_guide <- function() {
  guide_colorbar(
    title.position = "right",
    title.theme = element_text(
      angle = 90, hjust = 0.5, vjust = 0.5, size = 15.5,
      face = "bold", family = FINAL_FONT, colour = FINAL_COLOURS$ink
    ),
    barheight = unit(0.58, "npc"), barwidth = unit(0.52, "cm")
  )
}

render_meta_dot <- function(data, pathways, output_stem, title, subtitle,
                            p_col, fc_col, tissue_col, component_col,
                            pathway_col, tissue_order_col = NULL,
                            colour_limit = NULL, width = 22.0, height = 11.5) {
  plot_data <- copy(data)
  if (is.null(tissue_order_col)) {
    tissue_levels <- rev(sort(unique(as.character(plot_data[[tissue_col]]))))
  } else {
    tissue_levels <- rev(plot_data[order(get(tissue_order_col)), unique(as.character(get(tissue_col)))])
  }
  component_levels <- plot_data[order(get(pathway_col), ComponentOrder),
                                unique(as.character(get(component_col)))]
  plot_data[, `:=`(
    TissueFinal = factor(short_tissue_name_final(as.character(get(tissue_col))),
                         levels = short_tissue_name_final(tissue_levels)),
    ComponentFinal = factor(as.character(get(component_col)), levels = component_levels),
    PathwayFinal = factor(as.character(get(pathway_col)), levels = pathways),
    Significant = is.finite(get(p_col)) & get(p_col) < 0.05 & is.finite(get(fc_col)),
    Effect = as.numeric(get(fc_col))
  )]
  sig <- plot_data[Significant == TRUE]
  if (is.null(colour_limit)) {
    colour_limit <- quantile(abs(plot_data$Effect), 0.95, na.rm = TRUE, names = FALSE)
  }
  assert_final(nrow(sig) > 0L, paste(output_stem, "has no significant cells"))

  fig <- ggplot(plot_data, aes(ComponentFinal, TissueFinal)) +
    geom_tile(fill = FINAL_COLOURS$white, colour = FINAL_COLOURS$grid, linewidth = 0.18) +
    geom_point(
      data = sig, aes(colour = Effect, size = abs(Effect)), alpha = 0.96
    ) +
    facet_grid(. ~ PathwayFinal, scales = "free_x", space = "free_x", switch = "x") +
    scale_colour_gradient2(
      low = meta_palette[["low"]], mid = meta_palette[["mid"]], high = meta_palette[["high"]],
      midpoint = 0, limits = c(-colour_limit, colour_limit), oob = squish,
      name = "Meta log2FC (P < 0.05)", guide = vertical_meta_guide()
    ) +
    scale_size_continuous(
      range = c(4.0, 12.0), name = "|log2FC|",
      breaks = pretty(c(0, max(abs(sig$Effect))), n = 3),
      guide = guide_legend(
        title.position = "top", title.hjust = 0.5,
        keyheight = unit(1.5, "lines"),
        keywidth = unit(1.5, "lines")
      )
    ) +
    labs(
      title = title, subtitle = subtitle, x = NULL, y = NULL,
      caption = "Meta-log2FC is the median of mission-level median log2FC values; exact values and P values are retained in the publication tables."
    ) +
    theme_final(12.0) +
    theme(
      axis.text.x = element_text(angle = 52, hjust = 1, vjust = 1, size = 12.8),
      axis.text.y = element_text(size = 13.8),
      strip.placement = "outside",
      strip.background = element_rect(fill = FINAL_COLOURS$white, colour = NA),
      strip.text.x.bottom = element_text(size = 19.5, face = "bold",
                                         margin = margin(t = 5, b = 4)),
      panel.spacing.x = unit(0.52, "lines"),
      panel.background = element_rect(fill = FINAL_COLOURS$white, colour = NA),
      plot.background = element_rect(fill = FINAL_COLOURS$white, colour = NA),
      legend.position = "right",
      legend.title = element_text(size = 15.5, face = "bold"),
      legend.text = element_text(size = 13.5),
      plot.title = element_text(size = 25),
      plot.subtitle = element_text(size = 15.5),
      plot.caption = element_text(size = 12.0),
      plot.margin = margin(8, 12, 9, 8)
    )
  save_final_figure(fig, out, output_stem, width = width, height = height, dpi = 400)
  data.table(
    figure = basename(output_stem), n_rows = nrow(plot_data),
    n_tissues = uniqueN(plot_data$TissueFinal),
    n_components = uniqueN(plot_data$ComponentFinal), n_significant_dots = nrow(sig),
    colour_limit = colour_limit
  )
}

dsb <- fread(file.path(root, "data/publication_input/meta/08_essential_components_counts_log2FC_tissue_matrix.csv"))
assert_final(nrow(dsb) == 26L * 29L && uniqueN(dsb$EnsemblID) == 29L,
             "Invalid DSB meta matrix")
dsb_audit <- render_meta_dot(
  dsb, c("NHEJ", "HR", "HR–A-EJ", "A-EJ"), "Main/Fig_4e",
  "Double-strand break response across tissues",
  "Dots mark nominal tissue-meta P < 0.05; colour is direction and size is effect magnitude",
  p_col = "tissue_meta_p", fc_col = "Log2FCDisplayed", tissue_col = "Group",
  component_col = "OfficialMouseSymbol", pathway_col = "PathwayDisplay",
  width = 31.0, height = 13.0
)

ssb <- fread(file.path(root, "results/tables/01_SSB_tissue_meta_log2FC_pvalue_matrix_without_Neil2.csv"))
run_audit <- fread(file.path(root, "results/tables/04_run_audit.csv"))
assert_final(nrow(ssb) == 26L * 45L && uniqueN(ssb$EnsemblID) == 45L,
             "Invalid Neil2-excluded SSB meta matrix")
ssb_audit <- render_meta_dot(
  ssb, c("BER", "NER", "MMR", "FA"), "Main/Fig_5f",
  "Single-strand break response across tissues",
  "45 components; Neil2 excluded; dots mark nominal tissue-meta P < 0.05",
  p_col = "tissue_meta_p", fc_col = "Log2FCDisplayed", tissue_col = "Group",
  component_col = "DisplayComponent", pathway_col = "PathwayDisplay",
  tissue_order_col = "TissueOrder",
  colour_limit = as.numeric(run_audit$common_symmetric_fc_limit[[1L]]),
  width = 31.0, height = 13.0
)

audit <- rbindlist(list(dsb_audit, ssb_audit), fill = TRUE)
dir.create(file.path(out, "provenance"), recursive = TRUE, showWarnings = FALSE)
fwrite(audit, file.path(out, "provenance/meta_dotplot_scope_audit.csv"))
append_final_palette_audit(out, "Fig_4e", meta_palette)
append_final_palette_audit(out, "Fig_5f", meta_palette)
message("FINAL_META_DOTPLOT_PASS: Fig 4e and Fig 5f")
