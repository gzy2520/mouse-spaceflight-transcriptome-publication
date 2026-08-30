#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/figures/style_helpers.R"))
if (figure_style_c) source(file.path(root, "R/figures/style_C_helpers.R"))
counts <- fread(file.path(
  root, "results/tables/02_seven_pathway_member_gene_counts_per_source_sample_long.csv"
))

pathway_levels <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
pathway_shapes <- c(BER = 24L, NER = 25L, MMR = 23L, FA = 0L, HR = 8L, AEJ = 3L, NHEJ = 4L)
pathway_legend_labels <- c(
  BER = "BER  Base excision", NER = "NER  Nucleotide excision",
  MMR = "MMR  Mismatch repair", FA = "FA  Fanconi anemia",
  HR = "HR  Homologous recombination", AEJ = "A-EJ  Alternative end joining",
  NHEJ = "NHEJ  Non-homologous end joining"
)
tissue_levels <- c(
  "Adrenal gland", "Bone marrow", "Cecum", "Cerebellum", "Colon", "Dorsal skin",
  "Extensor digitorum longus", "Eye", "Femoral lateral skin", "Femoral skin",
  "Gastrocnemius", "Heart", "Heart / Heart right ventricle", "Kidney",
  "Left lobe of the liver", "Liver", "Lung", "Mammary gland", "Optic nerve",
  "Quadriceps femoris", "Retina", "Soleus", "Spleen", "Spleen-distal", "Thymus",
  "Tibialis anterior"
)
tissue_display <- setNames(tissue_levels, tissue_levels)
tissue_display[["Heart / Heart right ventricle"]] <- "Heart right ventricle"
stopifnot(uniqueN(counts$Group) == 26L, uniqueN(counts$Pathway) == 7L)
sample_colours <- c(Flight = FIGURE_COLOURS$flight, Ground = FIGURE_COLOURS$ground)

plot_data <- counts[, .(
  Tissue = factor(TissueLabel, levels = unname(tissue_display[tissue_levels])),
  TissueOrder,
  Pathway = factor(as.character(Pathway), levels = pathway_levels),
  sample_status, sample_column, n_member_genes
)]
plot_data[, sample_status := factor(sample_status, levels = c("Ground", "Flight"))]
setorder(plot_data, TissueOrder, sample_status, sample_column, Pathway)
max_y <- ceiling((max(plot_data$n_member_genes) + 12) / 10) * 10

if (figure_style_c) {
  # C-story S5: use pathway facets so the seven catalog sizes are readable
  # without 26 tiny tissue panels. Every source-sample point is retained.
  catalog <- fread(file.path(root, "results/tables/08_Table_S3_pathway_repertoire_sizes.csv"))
  catalog_map <- setNames(catalog$MemberGeneCount, catalog$Pathway)
  plot_data_c <- copy(plot_data)
  plot_data_c[, TissueC := factor(as.character(Tissue), levels = rev(unname(tissue_display[tissue_levels])))]
  plot_data_c[, PathwayLabel := factor(
    paste0(as.character(Pathway), "  (catalog n = ", unname(catalog_map[as.character(Pathway)]), ")"),
    levels = paste0(pathway_levels, "  (catalog n = ", unname(catalog_map[pathway_levels]), ")")
  )]
  median_c <- plot_data_c[, .(median_count = median(n_member_genes)), by = .(TissueC, PathwayLabel)]
  fig_s5_c <- ggplot(plot_data_c, aes(x = n_member_genes, y = TissueC, group = interaction(TissueC, sample_status))) +
    geom_boxplot(
      aes(fill = sample_status), width = 0.58, outlier.shape = NA,
      position = position_dodge(width = 0.68), colour = STYLE_C_COLOURS$ink, linewidth = 0.30
    ) +
    geom_point(
      aes(colour = sample_status),
      position = position_jitterdodge(jitter.height = 0.14, jitter.width = 0, dodge.width = 0.68, seed = 25),
      size = 1.25, alpha = 0.65, show.legend = FALSE
    ) +
    geom_point(
      data = median_c, aes(x = median_count, y = TissueC), inherit.aes = FALSE,
      shape = 23, size = 2.0, fill = STYLE_C_COLOURS$gold,
      colour = STYLE_C_COLOURS$ink, stroke = 0.30
    ) +
    facet_wrap(~PathwayLabel, ncol = 2, scales = "free_x") +
    scale_fill_manual(values = c(Flight = STYLE_C_COLOURS$red, Ground = STYLE_C_COLOURS$blue), name = "Sample status") +
    scale_colour_manual(values = c(Flight = STYLE_C_COLOURS$red, Ground = STYLE_C_COLOURS$blue), guide = "none") +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.10))) +
    labs(
      title = "Seven DNA-repair pathway repertoires across 26 tissues",
      subtitle = "Every point is one source sample; boxplots show status-specific distributions and diamonds show tissue medians",
      x = "Member genes detected", y = NULL,
      caption = "Presence is finite source value − 1 > 0; Ensembl Gene IDs are counted once per pathway. Counts describe repertoire presence, not pathway activity or expression magnitude."
    ) +
    theme_style_C(8.6) +
    theme(
      axis.text.y = element_text(size = 6.3),
      axis.text.x = element_text(size = 7.4),
      strip.background = element_rect(fill = STYLE_C_COLOURS$panel, colour = NA),
      strip.text = element_text(size = 9.0, face = "bold"),
      panel.grid.major.x = element_line(colour = STYLE_C_COLOURS$grid, linewidth = 0.28),
      legend.position = "top", legend.title = element_text(size = 8.5), legend.text = element_text(size = 8.0),
      plot.margin = margin(8, 10, 8, 8)
    )
  save_style_C(fig_s5_c, out, "Suppl/Fig_S5", 15.5, 14.0)
} else {

make_tissue_plot <- function(idx) {
  panel <- plot_data[TissueOrder == idx]
  panel_summary <- panel[, .(median_count = median(n_member_genes)), by = Pathway]
  panel_summary[, label_y := pmin(median_count + 6.5, max_y - 4)]
  tissue_name <- as.character(panel$Tissue[1L])
  show_y_axis <- ((idx - 1L) %% 4L) == 0L
  show_y_title <- idx == 1L
  ggplot(panel, aes(x = Pathway, y = n_member_genes)) +
    geom_boxplot(
      aes(group = Pathway), width = 0.88, outlier.shape = NA,
      fill = if (figure_style_a) "#EEF2F7" else "#E2E8F0",
      colour = if (figure_style_a) "#52606D" else "#475569", linewidth = 0.75
    ) +
    geom_point(
      aes(shape = Pathway, colour = sample_status),
      position = position_jitter(width = 0.39, height = 0, seed = 25),
      size = 2.75, alpha = 0.94, stroke = 0.85, show.legend = FALSE
    ) +
    geom_text(
      data = panel_summary, aes(x = Pathway, y = label_y, label = round(median_count)),
      inherit.aes = FALSE, size = 4.10,
      colour = if (figure_style_a) FIGURE_COLOURS$red else "#7A2E2A", fontface = "bold"
    ) +
    scale_shape_manual(values = pathway_shapes, limits = pathway_levels, drop = FALSE) +
    scale_colour_manual(values = sample_colours, guide = "none") +
    scale_x_discrete(limits = pathway_levels, labels = rep("", length(pathway_levels)), drop = FALSE) +
    scale_y_continuous(
      limits = c(0, max_y), breaks = pretty_breaks(n = 4),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = tissue_name, x = NULL, y = if (show_y_title) "Member genes" else NULL) +
    if (figure_style_a) {
      theme_minimal(base_size = 8.2, base_family = "Arial") +
        theme(
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = FIGURE_COLOURS$grid, linewidth = 0.25),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = if (show_y_axis) element_text(size = 7.5, colour = FIGURE_COLOURS$ink) else element_blank(),
          axis.ticks.y = if (show_y_axis) element_line(colour = FIGURE_COLOURS$muted, linewidth = 0.35) else element_blank(),
          axis.title.y = if (show_y_title) element_text(size = 9.5, colour = FIGURE_COLOURS$ink) else element_blank(),
          plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5, colour = FIGURE_COLOURS$ink, margin = margin(2, 1, 2, 1)),
          plot.background = element_rect(fill = "white", colour = FIGURE_COLOURS$ink, linewidth = 0.85),
          panel.border = element_rect(fill = NA, colour = "#CBD5E1", linewidth = 0.45), plot.margin = margin(0.5, 0.5, 0.5, 0.5)
        )
    } else {
      theme_minimal(base_size = 8, base_family = "Arial Unicode MS") +
        theme(
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#E5E7EB", linewidth = 0.25),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = if (show_y_axis) element_text(size = 7.5, colour = "#374151") else element_blank(),
          axis.ticks.y = if (show_y_axis) element_line(colour = "#6B7280", linewidth = 0.35) else element_blank(),
          axis.title.y = if (show_y_title) element_text(size = 9.5, colour = "#374151") else element_blank(),
          plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5, colour = "#1F2937", margin = margin(2, 1, 2, 1)),
          plot.background = element_rect(fill = "white", colour = "#374151", linewidth = 1.25),
          panel.border = element_rect(fill = NA, colour = "#CBD5E1", linewidth = 0.55), plot.margin = margin(0.5, 0.5, 0.5, 0.5)
        )
    }
}

legend_data <- data.table(
  Pathway = factor(pathway_levels, levels = pathway_levels),
  y = rev(seq_along(pathway_levels)),
  label = unname(pathway_legend_labels[pathway_levels])
)
status_data <- data.table(
  status = factor(c("Flight", "Ground"), levels = c("Flight", "Ground")),
  y = c(1.05, 0.62), label = c("Flight sample", "Ground sample")
)
legend_plot <- ggplot() +
  geom_point(
    data = legend_data, aes(x = 0.20, y = y, shape = Pathway),
    size = 5.2, colour = if (figure_style_a) FIGURE_COLOURS$ink else "#374151", fill = NA, show.legend = FALSE
  ) +
  geom_text(data = legend_data, aes(x = 0.48, y = y, label = label),
            hjust = 0, size = 3.6, colour = if (figure_style_a) FIGURE_COLOURS$ink else "#374151") +
  geom_point(data = status_data, aes(x = 2.22, y = y, colour = status),
             shape = 16, size = 4.5, show.legend = FALSE) +
  geom_text(data = status_data, aes(x = 2.38, y = y, label = label),
            hjust = 0, size = 3.4, colour = if (figure_style_a) FIGURE_COLOURS$ink else "#374151") +
  annotate("text", x = 0.04, y = 8.28, label = "Pathway symbols", hjust = 0,
           size = 4.7, fontface = "bold", colour = if (figure_style_a) FIGURE_COLOURS$ink else "#111827") +
  annotate("text", x = 0.04, y = 7.88, label = "Symbol = pathway", hjust = 0,
           size = 3.2, colour = if (figure_style_a) FIGURE_COLOURS$muted else "#6B7280") +
  annotate("text", x = 2.07, y = 1.53, label = "Colour = sample status", hjust = 0,
           size = 3.2, fontface = "bold", colour = if (figure_style_a) FIGURE_COLOURS$ink else "#111827") +
  scale_shape_manual(values = pathway_shapes, limits = pathway_levels, drop = FALSE) +
  scale_colour_manual(values = sample_colours) +
  scale_x_continuous(limits = c(0, 3.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.25, 8.65), expand = c(0, 0)) +
  theme_void(base_size = 8, base_family = if (figure_style_a) "Arial" else "Arial Unicode MS") +
  theme(
    plot.background = element_rect(fill = "white", colour = if (figure_style_a) "#CBD5E1" else "#9CA3AF", linewidth = 0.55),
    plot.margin = margin(5, 5, 5, 5)
  )

tissue_plots <- setNames(lapply(seq_along(tissue_levels), make_tissue_plot), LETTERS[seq_along(tissue_levels)])
plot_list <- c(tissue_plots, list(a = legend_plot))
layout_design <- paste(c("ABCD", "EFGH", "IJKL", "MNOP", "QRST", "UVWX", "YZaa"), collapse = "\n")
plot <- wrap_plots(plot_list, design = layout_design) +
  plot_annotation(
    title = if (figure_style_a) "Seven DNA-repair pathway member-gene counts" else "Seven DNA-repair pathway member-gene counts per source sample",
    subtitle = paste0(
      if (figure_style_a) "26 mouse tissues; one point per source sample; symbols encode pathway and colour encodes sample status.\n" else "26 mouse tissues in a 4 × 7 layout; each tissue is framed, and the lower-right two cells contain the symbol legend.\n",
      "Each point is one source sample column; symbol = pathway, colour = sample status; median labels are shown.\n",
      "Presence = finite source value − 1 > 0; Ensembl Gene IDs are counted once per pathway.\n",
      "Teacher-text/Table S3 catalog sizes: BER 54 | NER 53 | MMR 27 | FA 39 | HR 135 | A-EJ 9 | NHEJ 18."
    ),
    caption = "Counts describe repertoire presence, not pathway activity or expression magnitude. Source column labels, including technical-replicate markers, are retained in the audit table; no imputation was performed.",
    theme = if (figure_style_a) {
      theme(
        plot.title = style_title(18, hjust = 0), plot.subtitle = style_subtitle(10.2, hjust = 0),
        plot.caption = style_caption(8.2, hjust = 0), plot.margin = margin(7, 10, 6, 10)
      )
    } else {
      theme(
        plot.title = element_text(face = "bold", size = 20.5, hjust = 0, colour = "#111827"),
        plot.subtitle = element_text(size = 11.5, colour = "#374151", hjust = 0),
        plot.caption = element_text(size = 9.2, colour = "#4B5563", hjust = 0),
        plot.margin = margin(7, 10, 6, 10)
      )
    }
  )

ggsave(file.path(out, "Suppl/Fig_S5.png"), plot, width = 18.5, height = 26.5,
       units = "in", dpi = if (figure_style_a) 600 else 320, bg = "white", limitsize = FALSE)

# Cairo embeds the current wall-clock time in a PDF even when the page content
# is identical. Keep the approved vector export byte-for-byte; the PNG above is
# freshly rendered and provides the deterministic figure reproduction check.
pdf_source <- file.path(root, "results/Suppl/Fig_S5.pdf")
copied <- file.copy(pdf_source, file.path(out, "Suppl/Fig_S5.pdf"), overwrite = FALSE, copy.date = TRUE)
if (!copied) stop("Failed to copy the immutable approved Fig S5 PDF", call. = FALSE)
}
