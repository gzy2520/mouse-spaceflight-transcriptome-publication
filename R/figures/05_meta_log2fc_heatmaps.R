#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(25)
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out <- normalizePath(Sys.getenv("PUBLICATION_OUTPUT_DIR"), mustWork = TRUE)

# Fig 4b: DSB components. This is the approved plotting block applied to the
# frozen 26 x 29 publication matrix.
dsb <- fread(file.path(root, "data/publication_input/meta/08_essential_components_counts_log2FC_tissue_matrix.csv"))
stopifnot(nrow(dsb) == 26L * 29L, uniqueN(dsb$EnsemblID) == 29L)
# fread correctly recognizes the displayed values as numeric, but that drops
# the leading plus sign and fixed two-decimal formatting used in the approved
# labels. Reconstruct the display-only strings from the frozen numeric fields.
dsb[, Log2FCLabel := fifelse(
  NominalPDisplayed,
  sprintf("%+.2f", log2FoldChange),
  ""
)]
tissue_labels <- c(
  "Adrenal gland" = "肾上腺", "Bone marrow" = "骨髓", "Cecum" = "盲肠", "Cerebellum" = "小脑",
  "Colon" = "结肠", "Dorsal skin" = "背部皮肤", "Extensor digitorum longus" = "趾长伸肌", "Eye" = "眼",
  "Femoral lateral skin" = "股外侧皮肤", "Femoral skin" = "股部皮肤", "Gastrocnemius" = "腓肠肌",
  "Heart" = "心脏", "Heart / Heart right ventricle" = "心脏 / 右心室", "Kidney" = "肾脏",
  "Left lobe of the liver" = "肝左叶", "Liver" = "肝脏", "Lung" = "肺", "Mammary gland" = "乳腺",
  "Optic nerve" = "视神经", "Quadriceps femoris" = "股四头肌", "Retina" = "视网膜", "Soleus" = "比目鱼肌",
  "Spleen" = "脾脏", "Spleen-distal" = "脾脏（远端）", "Thymus" = "胸腺", "Tibialis anterior" = "胫骨前肌"
)
tissue_order <- sort(unique(dsb$Group))
tissue_display_labels <- setNames(paste0(tissue_order, "\n", tissue_labels[tissue_order]), tissue_order)
component_levels <- dsb[order(ComponentOrder), unique(OfficialMouseSymbol)]
pathway_levels <- dsb[order(ComponentOrder), unique(PathwayDisplay)]
dsb[, `:=`(
  TissueLabel = factor(
    paste0(Group, "\n", tissue_labels[Group]),
    levels = unname(tissue_display_labels[rev(tissue_order)])
  ),
  DisplayComponent = factor(OfficialMouseSymbol, levels = component_levels),
  PathwayDisplay = factor(PathwayDisplay, levels = pathway_levels)
)]
dsb_limit <- quantile(abs(dsb$Log2FCDisplayed), probs = 0.95, na.rm = TRUE, names = FALSE)
dsb_theme <- theme_minimal(base_size = 8.5, base_family = "Arial Unicode MS") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 6.7, lineheight = 0.85), axis.ticks = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", colour = NA),
    strip.text.x.bottom = element_text(face = "bold", size = 9),
    strip.placement = "outside", panel.spacing.x = grid::unit(0.45, "lines"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5),
    plot.caption = element_text(size = 7.2, hjust = 0, colour = "#444444"),
    legend.position = "right", plot.margin = margin(6, 10, 6, 7)
  )
fig4b <- ggplot(dsb, aes(x = DisplayComponent, y = TissueLabel, fill = Log2FCDisplayed)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  geom_text(aes(label = Log2FCLabel), size = 2.25, colour = "#2F2F2F") +
  facet_grid(. ~ PathwayDisplay, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0,
    limits = c(-dsb_limit, dsb_limit), oob = squish,
    na.value = "#EEF1F4", name = "Meta log2FC\n(P < 0.05)"
  ) +
  labs(
    title = "Flight response: tissue meta-log2FC",
    subtitle = "Only cells with unadjusted tissue-meta P < 0.05 are coloured and labelled; blue = lower in flight, red = higher in flight.",
    x = NULL, y = NULL,
    caption = "Meta-log2FC is the median of mission-level median log2FC values. Nominal tissue-meta P is from mission-level signed Stouffer aggregation."
  ) + dsb_theme
ggsave(file.path(out, "Main/Fig_4b.png"), fig4b, width = 20, height = 10.5, dpi = 320, bg = "white")

# Fig 5b: SSB components after removing Neil2. The symmetric limit is retained
# from the original 75-component source so the approved colour mapping is exact.
ssb <- fread(file.path(root, "results/tables/01_SSB_tissue_meta_log2FC_pvalue_matrix_without_Neil2.csv"))
run_audit <- fread(file.path(root, "results/tables/04_run_audit.csv"))
stopifnot(nrow(ssb) == 26L * 45L, uniqueN(ssb$EnsemblID) == 45L)
ssb[, Log2FCLabel := fifelse(
  P_lt_0_05,
  sprintf("%+.2f", log2FoldChange),
  ""
)]
ssb_tissues <- ssb[order(TissueOrder), unique(TissueLabel)]
ssb_components <- ssb[order(PathwayOrder, ComponentOrder), unique(DisplayComponent)]
ssb[, `:=`(
  TissueLabel = factor(as.character(TissueLabel), levels = rev(ssb_tissues)),
  DisplayComponent = factor(as.character(DisplayComponent), levels = ssb_components),
  PathwayDisplay = factor(as.character(PathwayDisplay), levels = c("BER", "NER", "MMR", "FA"))
)]
ssb_limit <- as.numeric(run_audit$common_symmetric_fc_limit[[1L]])
ssb_theme <- theme_minimal(base_size = 8.8, base_family = "Arial Unicode MS") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.5, colour = "#3D3D3D"),
    axis.text.y = element_text(size = 7.1, lineheight = 0.84, colour = "#3D3D3D"),
    axis.ticks = element_blank(), strip.background = element_rect(fill = "#F0F0F0", colour = NA),
    strip.text.x.bottom = element_text(face = "bold", size = 9.5), strip.placement = "outside",
    panel.spacing.x = grid::unit(0.42, "lines"), plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 7.5, hjust = 0, colour = "#444444"),
    legend.position = "right", plot.margin = margin(8, 12, 10, 8)
  )
fig5b <- ggplot(ssb, aes(x = DisplayComponent, y = TissueLabel, fill = Log2FCDisplayed)) +
  geom_tile(colour = "white", linewidth = 0.16) +
  geom_text(aes(label = Log2FCLabel), size = 1.9, colour = "#2F2F2F", family = "Arial Unicode MS") +
  facet_grid(. ~ PathwayDisplay, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0,
    limits = c(-ssb_limit, ssb_limit), oob = squish,
    na.value = "#EEF1F4", name = "Meta log2FC\n(P < 0.05)"
  ) +
  labs(
    title = "SSB essential components: tissue meta-log2FC (P < 0.05)",
    subtitle = "26 mouse tissues; 45 SSB components; Neil2 removed to match Fig 5a; blue = lower in flight, red = higher in flight",
    x = NULL, y = NULL,
    caption = paste(
      "Only unadjusted tissue_meta_p < 0.05 cells are coloured and labelled; grey cells are not displayed as significant.",
      "log2FC is the median of mission-level median log2FC values; P is from mission-level signed Stouffer aggregation."
    )
  ) + ssb_theme
ggsave(file.path(out, "Main/Fig_5b.png"), fig5b, width = 28, height = 12.5,
       units = "in", dpi = 320, bg = "white", limitsize = FALSE)
