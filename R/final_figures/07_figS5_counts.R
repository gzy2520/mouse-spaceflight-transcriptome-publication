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
source(file.path(root, "R/final_figures/helpers.R"))
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
sample_colours <- c(Flight = "#D55E00", Ground = "#0072B2")

plot_data <- counts[, .(
  Tissue = factor(TissueLabel, levels = unname(tissue_display[tissue_levels])),
  TissueOrder,
  Pathway = factor(as.character(Pathway), levels = pathway_levels),
  sample_status, sample_column, n_member_genes
)]
plot_data[, sample_status := factor(sample_status, levels = c("Ground", "Flight"))]
setorder(plot_data, TissueOrder, sample_status, sample_column, Pathway)
max_y <- ceiling((max(plot_data$n_member_genes) + 12) / 10) * 10

make_tissue_plot <- function(idx) {
  panel <- plot_data[TissueOrder == idx]
  panel_summary <- panel[, .(median_count = median(n_member_genes)), by = Pathway]
  panel_summary[, label_y := pmin(median_count + 6.5, max_y - 4)]
  tissue_name <- balanced_wrap(as.character(panel$Tissue[1L]), width = 22L, max_lines = 2L)
  show_y_axis <- ((idx - 1L) %% 4L) == 0L
  show_y_title <- idx == 1L
  ggplot(panel, aes(x = Pathway, y = n_member_genes)) +
    geom_boxplot(
      aes(group = Pathway), width = 0.88, outlier.shape = NA,
      fill = "#E2E8F0",
      colour = "#475569", linewidth = 0.75
    ) +
    geom_point(
      aes(shape = Pathway, colour = sample_status),
      position = position_jitter(width = 0.39, height = 0, seed = 25),
      size = 2.75, alpha = 0.94, stroke = 0.85, show.legend = FALSE
    ) +
    geom_text(
      data = panel_summary, aes(x = Pathway, y = label_y, label = round(median_count)),
      inherit.aes = FALSE, size = 5.60,
      colour = "#7A2E2A", fontface = "bold"
    ) +
    scale_shape_manual(values = pathway_shapes, limits = pathway_levels, drop = FALSE) +
    scale_colour_manual(values = sample_colours, guide = "none") +
    scale_x_discrete(limits = pathway_levels, labels = rep("", length(pathway_levels)), drop = FALSE) +
    scale_y_continuous(
      limits = c(0, max_y), breaks = pretty_breaks(n = 4),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = tissue_name, x = NULL, y = if (show_y_title) "Member genes" else NULL) +
    {
      theme_minimal(base_size = 11, base_family = FINAL_FONT) +
        theme(
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#E5E7EB", linewidth = 0.25),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = if (show_y_axis) element_text(size = 10.5, colour = "#374151") else element_blank(),
          axis.ticks.y = if (show_y_axis) element_line(colour = "#6B7280", linewidth = 0.35) else element_blank(),
          axis.title.y = if (show_y_title) element_text(size = 13.0, colour = "#374151") else element_blank(),
          plot.title = element_text(face = "bold", size = 23.0, hjust = 0.5, lineheight = 0.86, colour = "#1F2937", margin = margin(3, 1, 2, 1)),
          plot.background = element_rect(fill = "white", colour = "#374151", linewidth = 1.25),
          panel.border = element_rect(fill = NA, colour = "#CBD5E1", linewidth = 0.55), plot.margin = margin(0.5, 0.5, 0.5, 0.5)
        )
    }
}

legend_data <- data.table(
  Pathway = factor(pathway_levels, levels = pathway_levels),
  x_symbol = c(rep(0.20, 4L), rep(4.15, 3L)),
  x_text = c(rep(0.52, 4L), rep(4.47, 3L)),
  y = c(5.00, 4.00, 3.00, 2.00, 5.00, 4.00, 3.00),
  label = unname(pathway_legend_labels[pathway_levels])
)
status_data <- data.table(
  status = factor(c("Flight", "Ground"), levels = c("Flight", "Ground")),
  x_point = c(2.55, 5.65), x_text = c(2.88, 5.98),
  y = c(0.72, 0.72), label = c("Flight sample", "Ground sample")
)
legend_plot <- ggplot() +
  geom_point(
    data = legend_data, aes(x = x_symbol, y = y, shape = Pathway),
    size = 7.8, colour = "#374151", fill = NA, show.legend = FALSE
  ) +
  geom_text(data = legend_data, aes(x = x_text, y = y, label = label),
            hjust = 0, size = 6.0, colour = "#374151") +
  geom_point(data = status_data, aes(x = x_point, y = y, colour = status),
             shape = 16, size = 7.0, show.legend = FALSE) +
  geom_text(data = status_data, aes(x = x_text, y = y, label = label),
            hjust = 0, size = 5.8, colour = "#374151") +
  annotate("text", x = 0.04, y = 6.05, label = "Pathway symbols", hjust = 0,
           size = 7.2, fontface = "bold", colour = "#111827") +
  annotate("text", x = 0.10, y = 0.72, label = "Sample status", hjust = 0,
           size = 6.4, fontface = "bold", colour = "#111827") +
  scale_shape_manual(values = pathway_shapes, limits = pathway_levels, drop = FALSE) +
  scale_colour_manual(values = sample_colours) +
  scale_x_continuous(limits = c(0, 9.40), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.25, 6.45), expand = c(0, 0)) +
  theme_void(base_size = 10, base_family = FINAL_FONT) +
  theme(
    plot.background = element_rect(fill = "white", colour = "#9CA3AF", linewidth = 0.55),
    plot.margin = margin(5, 5, 5, 5)
  )

tissue_plots <- setNames(lapply(seq_along(tissue_levels), make_tissue_plot), LETTERS[seq_along(tissue_levels)])
plot_list <- c(tissue_plots, list(a = legend_plot))
layout_design <- paste(c("ABCD", "EFGH", "IJKL", "MNOP", "QRST", "UVWX", "YZaa"), collapse = "\n")
plot <- wrap_plots(plot_list, design = layout_design) +
  plot_annotation(
    title = "Seven DNA-repair pathway member-gene counts per source sample",
    subtitle = paste0(
      "26 mouse tissues in a 4 × 7 layout; each tissue is framed, and the lower-right two cells form a dedicated legend panel.\n",
      "Each point is one source sample column; symbol = pathway, colour = sample status; median labels are shown.\n",
      "Presence = finite source value − 1 > 0; Ensembl Gene IDs are counted once per pathway.\n",
      "Teacher-text/Table S3 catalog sizes: BER 54 | NER 53 | MMR 27 | FA 39 | HR 135 | A-EJ 9 | NHEJ 18."
    ),
    caption = "Counts describe repertoire presence, not pathway activity or expression magnitude. Source column labels, including technical-replicate markers, are retained in the audit table; no imputation was performed.",
    theme = {
      theme(
        plot.title = element_text(face = "bold", size = 26, hjust = 0, colour = "#111827"),
        plot.subtitle = element_text(size = 14.5, colour = "#374151", hjust = 0),
        plot.caption = element_text(size = 11.5, colour = "#4B5563", hjust = 0),
        plot.margin = margin(7, 10, 6, 10)
      )
    }
  )

save_final_figure(plot, out, "Suppl/Fig_S5", width = 20.0, height = 28.0, dpi = 400)
dir.create(file.path(out, "provenance"), recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(
  n_rows = nrow(counts), n_tissues = uniqueN(counts$Group),
  n_pathways = uniqueN(counts$Pathway), n_samples = uniqueN(counts$sample_column), seed = 25L
), file.path(out, "provenance/Fig_S5_scope_audit.csv"))
append_final_palette_audit(out, "Fig_S5", sample_colours)
message("FINAL_FIGS5_PASS: 26 tissues, 7 pathways, seed 25")
