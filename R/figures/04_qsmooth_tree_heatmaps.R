#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)
source(file.path(root, "R/figures/style_helpers.R"))
if (figure_style_c) source(file.path(root, "R/figures/style_C_helpers.R"))
input_dir <- file.path(root, "data/publication_input/qsmooth")
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

component_audit <- fread(file.path(input_dir, "01_filtered_component_stable_id_mapping_audit.csv"))
tissue_matrix <- fread(file.path(input_dir, "05_component_tissue_yarn_qsmooth_log2_matrix.csv"))
setorder(component_audit, PlotComponentOrder)
setorder(tissue_matrix, TissueOrder, PlotComponentOrder)
assert(nrow(component_audit) == 74L, "Expected 74 Neil2-excluded components")
assert(nrow(tissue_matrix) == 26L * 74L && all(is.finite(tissue_matrix$TissueYARNNormalizedLog2)),
       "Invalid 26 x 74 YARN matrix")

tissue_labels <- c(
  "Adrenal gland" = "肾上腺", "Bone marrow" = "骨髓", "Cecum" = "盲肠", "Cerebellum" = "小脑",
  "Colon" = "结肠", "Dorsal skin" = "背部皮肤", "Extensor digitorum longus" = "趾长伸肌", "Eye" = "眼",
  "Femoral lateral skin" = "股外侧皮肤", "Femoral skin" = "股部皮肤", "Gastrocnemius" = "腓肠肌",
  "Heart" = "心脏", "Heart / Heart right ventricle" = "右心室", "Kidney" = "肾脏",
  "Left lobe of the liver" = "肝左叶", "Liver" = "肝脏", "Lung" = "肺", "Mammary gland" = "乳腺",
  "Optic nerve" = "视神经", "Quadriceps femoris" = "股四头肌", "Retina" = "视网膜", "Soleus" = "比目鱼肌",
  "Spleen" = "脾脏", "Spleen-distal" = "脾脏（远端）", "Thymus" = "胸腺", "Tibialis anterior" = "胫骨前肌"
)
tissue_universe <- tissue_matrix[order(TissueOrder), unique(Group)]

spearman_cluster <- function(plot_data, figure_name) {
  profile <- dcast(
    plot_data[, .(Group, OfficialMouseSymbol, TissueYARNNormalizedLog2)],
    Group ~ OfficialMouseSymbol, value.var = "TissueYARNNormalizedLog2"
  )
  profile <- profile[match(tissue_universe, Group)]
  component_columns <- setdiff(names(profile), "Group")
  matrix <- as.matrix(profile[, ..component_columns])
  storage.mode(matrix) <- "numeric"
  rownames(matrix) <- profile$Group
  rho <- suppressWarnings(cor(t(matrix), method = "spearman", use = "pairwise.complete.obs"))
  rho[rho > 1] <- 1
  rho[rho < -1] <- -1
  distance <- 1 - rho
  distance <- (distance + t(distance)) / 2
  diag(distance) <- 0
  hc <- hclust(as.dist(distance), method = "average")
  list(
    figure = figure_name,
    component_n = length(component_columns),
    method = "average-linkage hierarchical clustering",
    distance_definition = "1 - Spearman rho across displayed components",
    hc = hc,
    ordered_tissues = rownames(matrix)[hc$order]
  )
}

dsb_pathways <- c("NHEJ", "HR", "HR–A-EJ", "A-EJ")
ssb_pathways <- c("BER", "NER", "MMR", "FA")
dsb <- spearman_cluster(tissue_matrix[PlotPathwayDisplay %chin% dsb_pathways], "DSB")
ssb <- spearman_cluster(tissue_matrix[PlotPathwayDisplay %chin% ssb_pathways], "SSB")
yarn_min <- min(tissue_matrix$TissueYARNNormalizedLog2)
yarn_max <- max(tissue_matrix$TissueYARNNormalizedLog2)

build_tree_segments <- function(hc) {
  n <- length(hc$labels)
  ordered <- hc$labels[hc$order]
  leaf_y <- setNames(seq_len(n), ordered)
  node_x <- numeric(n - 1L)
  node_y <- numeric(n - 1L)
  segments <- vector("list", n - 1L)
  child_coord <- function(code) {
    if (code < 0L) c(x = 0, y = unname(leaf_y[hc$labels[-code]])) else c(x = node_x[[code]], y = node_y[[code]])
  }
  for (i in seq_len(n - 1L)) {
    left <- child_coord(hc$merge[i, 1L])
    right <- child_coord(hc$merge[i, 2L])
    parent_x <- hc$height[[i]]
    parent_y <- mean(c(left[["y"]], right[["y"]]))
    node_x[[i]] <- parent_x
    node_y[[i]] <- parent_y
    segments[[i]] <- rbind(
      data.table(x = left[["x"]], xend = parent_x, y = left[["y"]], yend = left[["y"]]),
      data.table(x = right[["x"]], xend = parent_x, y = right[["y"]], yend = right[["y"]]),
      data.table(x = parent_x, xend = parent_x, y = min(left[["y"]], right[["y"]]), yend = max(left[["y"]], right[["y"]]))
    )
  }
  list(segments = rbindlist(segments), ordered_tissues = ordered)
}

display_tissue_name <- function(x) {
  x[x == "Heart / Heart right ventricle"] <- "Heart right ventricle"
  x
}

make_tree_plot <- function(cluster) {
  tree <- build_tree_segments(cluster$hc)
  max_height <- max(cluster$hc$height)
  tree$segments[, `:=`(x = max_height - x, xend = max_height - xend)]
  ggplot(tree$segments) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = if (figure_style_a) 0.52 else 0.62,
                 colour = if (figure_style_a) FIGURE_COLOURS$ink else "#334E68", lineend = "round") +
    scale_x_continuous(limits = c(0, max_height), breaks = NULL, name = NULL, expand = c(0, 0)) +
    scale_y_reverse(limits = c(26.7, 0.3), breaks = NULL, expand = c(0, 0)) +
    if (figure_style_a) {
      theme_void(base_size = 9.2, base_family = "Arial") +
        theme(plot.margin = margin(2, 0, 10, 0))
    } else {
      theme_minimal(base_size = 9.2, base_family = "Arial Unicode MS") +
        theme(
          panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.margin = margin(2, 0, 13, 0)
        )
    }
}

make_label_plot <- function(ordered_tissues) {
  label_values <- if (figure_style_c) {
    display_tissue_name(ordered_tissues)
  } else if (figure_style_a) {
    display_tissue_name(ordered_tissues)
  } else {
    paste0(display_tissue_name(ordered_tissues), "\n", unname(tissue_labels[ordered_tissues]))
  }
  labels <- data.table(
    y = seq_along(ordered_tissues),
    Label = label_values
  )
  ggplot(labels, aes(x = 0, y = y, label = Label)) +
    geom_text(hjust = 0, size = if (figure_style_a) 3.25 else 3.05, lineheight = 0.82,
              family = if (figure_style_a) "Arial" else "Arial Unicode MS",
              colour = if (figure_style_a) FIGURE_COLOURS$ink else "#263238") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_reverse(limits = c(26.7, 0.3), breaks = NULL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = if (figure_style_a) "Arial" else "Arial Unicode MS") +
    theme(plot.margin = margin(2, 0, if (figure_style_a) 10 else 13, 0))
}

make_heatmap_plot <- function(plot_data, ordered_tissues) {
  plot_data <- copy(plot_data)
  plot_data[, DisplayComponent := factor(OfficialMouseSymbol, levels = component_audit$OfficialMouseSymbol)]
  plot_data[, PlotPathwayDisplay := factor(
    PlotPathwayDisplay, levels = c("NHEJ", "HR", "HR–A-EJ", "A-EJ", "BER", "NER", "MMR", "FA")
  )]
  plot_data[, TissueKey := factor(Group, levels = rev(ordered_tissues))]
  plot_data[, CellTextColour := fifelse(
    TissueYARNNormalizedLog2 <= yarn_min + 0.16 * (yarn_max - yarn_min) |
      TissueYARNNormalizedLog2 >= yarn_min + 0.84 * (yarn_max - yarn_min),
    "white", "#1F2937"
  )]
  if (figure_style_c) {
    return(
      ggplot(plot_data, aes(x = DisplayComponent, y = TissueKey, fill = TissueYARNNormalizedLog2)) +
        geom_tile(colour = STYLE_C_COLOURS$white, linewidth = 0.16) +
        facet_grid(. ~ PlotPathwayDisplay, scales = "free_x", space = "free_x", switch = "x", drop = TRUE) +
        scale_x_discrete(expand = expansion(add = 0)) +
        scale_y_discrete(drop = FALSE, labels = function(x) rep("", length(x))) +
        scale_style_C_expression(c(yarn_min, yarn_max), "YARN qsmooth log2 expression") +
        labs(x = NULL, y = NULL) +
        theme_style_C(8.7) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 58, hjust = 1, vjust = 1, size = 6.2),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.background = element_rect(fill = STYLE_C_COLOURS$panel, colour = NA),
          strip.text.x.bottom = element_text(face = "bold", size = 9.0),
          strip.placement = "outside", panel.spacing.x = unit(0.42, "lines"),
          legend.position = "top", legend.justification = "center",
          legend.title = element_text(size = 8.5, face = "bold"),
          legend.text = element_text(size = 7.5),
          plot.margin = margin(2, 0, 10, 0)
        )
    )
  }
  ggplot(plot_data, aes(x = DisplayComponent, y = TissueKey, fill = TissueYARNNormalizedLog2)) +
    geom_tile(colour = "white", linewidth = if (figure_style_a) 0.12 else 0.16) +
    geom_text(aes(label = sprintf("%.1f", TissueYARNNormalizedLog2), colour = CellTextColour),
              size = if (figure_style_a) 1.35 else 1.55, na.rm = TRUE, show.legend = FALSE) +
    scale_colour_identity() +
    facet_grid(. ~ PlotPathwayDisplay, scales = "free_x", space = "free_x", switch = "x", drop = TRUE) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(drop = FALSE, labels = function(x) rep("", length(x))) +
    scale_fill_gradientn(
      colours = if (figure_style_a) {
        c("#3B6FB6", "#74A9CF", "#F7F7F5", "#E9A66F", "#C65A5A")
      } else {
        c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")
      },
      values = scales::rescale(seq(yarn_min, yarn_max, length.out = 5L), from = c(yarn_min, yarn_max)),
      limits = c(yarn_min, yarn_max), oob = scales::squish,
      name = "YARN qsmooth\nlog2 expression",
      breaks = pretty(c(yarn_min, yarn_max), n = 4),
      labels = function(x) formatC(x, format = "fg", digits = 3),
      guide = guide_colorbar(
        title.position = "top", title.hjust = 0.5, direction = "horizontal",
        barwidth = unit(15, "cm"), barheight = unit(0.38, "cm"),
        ticks.colour = "#52606D", frame.colour = "#8A9BA8"
      )
    ) +
    if (figure_style_a) {
      theme_minimal(base_size = 8.8, base_family = "Arial") +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 58, hjust = 1, vjust = 1, size = 6.5, colour = FIGURE_COLOURS$ink),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
          strip.background = element_rect(fill = FIGURE_COLOURS$panel, colour = NA),
          strip.text.x.bottom = element_text(face = "bold", size = 9.0, colour = FIGURE_COLOURS$ink),
          strip.placement = "outside", panel.spacing.x = unit(0.36, "lines"),
          legend.position = "top", legend.justification = "center",
          legend.title = element_text(size = 8.5, face = "bold"), legend.text = element_text(size = 7.6),
          plot.margin = margin(2, 0, 10, 0)
        )
    } else {
      theme_minimal(base_size = 8.8, base_family = "Arial Unicode MS") +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.3, colour = "#3D3D3D"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
          strip.background = element_rect(fill = "#F0F0F0", colour = NA),
          strip.text.x.bottom = element_text(face = "bold", size = 9.2),
          strip.placement = "outside", panel.spacing.x = unit(0.42, "lines"),
          legend.position = "top", legend.justification = "center",
          legend.title = element_text(size = 8.7, face = "bold"), legend.text = element_text(size = 7.8),
          plot.margin = margin(2, 0, 13, 0)
        )
    }
}

render_combined <- function(cluster, pathways, title, subtitle, output_file) {
  selected <- tissue_matrix[PlotPathwayDisplay %chin% pathways]
  component_n <- component_audit[PlottedPathwayDisplay %chin% pathways, .N]
  assert(nrow(selected) == 26L * component_n, paste("Unexpected", cluster$figure, "cells"))
  if (figure_style_c) {
    c_title <- if (cluster$figure == "DSB") {
      "Double-strand break repair components across tissues"
    } else {
      "Single-strand break repair components across tissues"
    }
    c_subtitle <- if (cluster$figure == "DSB") {
      "29 components across 26 tissues; rows ordered by tissue-level Spearman clustering"
    } else {
      "45 components across 26 tissues; Neil2 excluded; rows ordered by tissue-level Spearman clustering"
    }
    combined_c <- make_tree_plot(cluster) + make_heatmap_plot(selected, cluster$ordered_tissues) +
      make_label_plot(cluster$ordered_tissues) +
      plot_layout(widths = c(3.6, 31.0, 8.2), guides = "collect") +
      plot_annotation(
        title = c_title,
        subtitle = c_subtitle,
        caption = "Left: average-linkage tree from tissue Spearman profiles. Centre: frozen YARN qsmooth log2 expression matrix. Right: tissue labels."
      ) &
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13.5, colour = STYLE_C_COLOURS$ink),
        plot.subtitle = element_text(hjust = 0.5, size = 8.6, colour = STYLE_C_COLOURS$muted),
        plot.caption = element_text(hjust = 0, size = 7.2, colour = STYLE_C_COLOURS$muted),
        plot.margin = margin(8, 8, 8, 8)
      ) & theme(legend.position = "top")
    stem_c <- sub("[.]png$", "", output_file)
    save_style_C(combined_c, out, stem_c, 25.0, 9.2)
    return(invisible(NULL))
  }
  combined <- make_tree_plot(cluster) + make_heatmap_plot(selected, cluster$ordered_tissues) +
    make_label_plot(cluster$ordered_tissues) +
    plot_layout(widths = if (figure_style_a) c(3.6, 28, 6.0) else c(4.4, 28, 7.0), guides = "collect") +
    plot_annotation(
      title = title, subtitle = subtitle,
      caption = if (figure_style_a) {
        "Left: tissue Spearman tree; centre: YARN qsmooth log2 expression; right: tissue labels. Values use the complete shared Ensembl-gene fit and are unchanged."
      } else {
        "左：从组织表达谱重算的 Spearman 行聚类树，根在左、叶端贴合对应组织行；中：YARN qsmooth log2 expression 热图；右：组织名称。颜色图例置于顶部。YARN 数值沿用完整共同 Ensembl 基因全集拟合结果，未因删除 Neil2 而改变。"
      },
      theme = if (figure_style_a) {
        theme(
          plot.title = style_title(13, hjust = 0.5, family = "Arial"),
          plot.subtitle = style_subtitle(8.5, hjust = 0.5, family = "Arial"),
          plot.caption = style_caption(7.2, hjust = 0, family = "Arial"), plot.margin = margin(7, 8, 7, 8)
        )
      } else {
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14, family = "Arial Unicode MS"),
          plot.subtitle = element_text(hjust = 0.5, size = 8.8, family = "Arial Unicode MS"),
          plot.caption = element_text(hjust = 0, size = 7.3, colour = "#52606D", family = "Arial Unicode MS"),
          plot.margin = margin(8, 8, 8, 8)
        )
      }
    ) & theme(legend.position = "top")
  ggsave(file.path(out, output_file), combined,
         width = if (figure_style_a) 36 else 39, height = if (figure_style_a) 12.5 else 13,
         dpi = if (figure_style_a) 600 else 320, bg = "white", limitsize = FALSE)
}

render_combined(
  dsb, dsb_pathways,
  if (figure_style_a) "DSB components: tissue tree and YARN expression" else "DSB essential components: tissue Spearman tree + YARN qsmooth heatmap",
  if (figure_style_a) "26 mouse tissues; 29 DSB components; rows ordered by the DSB Spearman tree" else "26 mouse tissues; 29 DSB components; Exo1/Dna2 displayed in HR; rows ordered by the recomputed DSB tree",
  "Main/Fig_4a.png"
)
render_combined(
  ssb, ssb_pathways,
  if (figure_style_a) "SSB components: tissue tree and YARN expression" else "SSB essential components: tissue Spearman tree + YARN qsmooth heatmap",
  if (figure_style_a) "26 mouse tissues; 45 SSB components; Neil2 excluded from BER; rows ordered by the SSB Spearman tree" else "26 mouse tissues; 45 SSB components after excluding Neil2 from BER; rows ordered by the recomputed SSB tree",
  "Main/Fig_5a.png"
)
