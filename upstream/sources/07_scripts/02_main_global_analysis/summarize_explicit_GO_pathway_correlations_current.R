#!/usr/bin/env Rscript

# Data-derived association, clustering and two-dimensional mapping of the
# explicit mouse GO pathways. The primary observation unit is a standardised
# tissue. Registered analysis units are retained as a sensitivity analysis.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})
set.seed(25)

root <- normalizePath(".", mustWork = TRUE)
run_root <- Sys.getenv(
  "STABLE_ID_RERUN_ROOT",
  unset = file.path(root, "03_analysis_results/10_full_stable_id_rerun_20260823")
)
go_dir <- Sys.getenv(
  "STABLE_GO_OUTPUT_DIR",
  unset = file.path(run_root, "03_mouse_go_concrete_terms_20260824")
)
table_dir <- file.path(go_dir, "tables")
figure_dir <- file.path(go_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

tissue_file <- file.path(
  table_dir, "04_tissue_statistics_concrete_terms_and_context.csv"
)
unit_file <- file.path(
  table_dir, "02_dataset_GSEA_concrete_terms_and_context.csv.gz"
)
stopifnot(file.exists(tissue_file), file.exists(unit_file))

tissue <- fread(tissue_file)
unit <- fread(unit_file)
term_meta <- unique(tissue[, .(
  term_key, term_name, go_id, term_order, term_role
)])
order_audit_file <- file.path(table_dir, "00_teacher_requested_term_order.csv")
main_order_file <- file.path(
  run_root,
  "03_mouse_go_concrete_terms_20260824/tables/00_teacher_requested_term_order.csv"
)
project_order_file <- file.path(
  root, "07_scripts/config/mouse_full_go_term_display_order.csv"
)
if (!file.exists(order_audit_file)) {
  order_audit_file <- if (file.exists(main_order_file)) {
    main_order_file
  } else {
    project_order_file
  }
}
if (file.exists(order_audit_file)) {
  order_audit <- fread(order_audit_file)
  required_order_columns <- c(
    "original_term_order", "display_term_order", "term_key"
  )
  if (!all(required_order_columns %chin% names(order_audit)) ||
      nrow(order_audit) != 15L ||
      anyDuplicated(order_audit$term_key) ||
      anyDuplicated(order_audit$display_term_order) ||
      !identical(sort(as.integer(order_audit$display_term_order)), seq_len(15L)) ||
      !identical(sort(as.integer(order_audit$original_term_order)), seq_len(15L))) {
    stop("GO term order audit must contain a one-to-one permutation of the 15 terms.")
  }
  term_meta <- merge(
    term_meta,
    order_audit[, .(term_key, original_term_order, display_term_order)],
    by = "term_key", all.x = TRUE, sort = FALSE
  )
} else {
  stop("No validated GO term display-order contract found: ", order_audit_file)
}
stopifnot(
  nrow(term_meta) == 15L,
  !anyDuplicated(term_meta$term_key),
  identical(sort(unique(unit$term_key)), sort(term_meta$term_key))
)
teacher_order <- term_meta[order(display_term_order), term_key]
stopifnot(length(teacher_order) == 15L)
term_meta <- term_meta[match(teacher_order, term_key)]
term_meta[, term_order := display_term_order]
fwrite(
  term_meta[, .(
    original_term_order, display_term_order, term_key, term_name, go_id, term_role
  )],
  file.path(table_dir, "00_teacher_requested_term_order.csv"),
  bom = TRUE
)
term_meta[, term_label := paste0(term_name, "\n", go_id)]
term_meta[, short_label := fifelse(
  nchar(term_name) > 36L,
  paste0(substr(term_name, 1L, 33L), "..."),
  term_name
)]
display_name_overrides <- c(
  DNA_damage_signal_transduction_GO0042770 =
    "DNA-damage signal transduction",
  intrinsic_apoptotic_signaling_GO0008630 =
    "Intrinsic apoptosis",
  telomere_maintenance_GO0043247 =
    "Telomere maintenance after DNA damage",
  telomere_region_GO0000781 =
    "Telomeric chromosome region",
  mechanosensation_GO0050954 =
    "Mechanical stimulus perception"
)
term_meta[
  term_key %chin% names(display_name_overrides),
  short_label := unname(display_name_overrides[term_key])
]
term_meta[, short_label := paste0(short_label, "\n", go_id)]

make_wide <- function(data, id_columns, value_column, expected_rows) {
  pathway_columns <- term_meta$term_key
  formula <- as.formula(paste(
    paste(id_columns, collapse = " + "), "~ term_key"
  ))
  wide <- dcast(data, formula, value.var = value_column)
  if (nrow(wide) != expected_rows ||
      any(!term_meta$term_key %chin% names(wide)) ||
      any(!is.finite(as.matrix(wide[, ..pathway_columns])))) {
    stop("Pathway NES matrix is incomplete or contains non-finite values")
  }
  setcolorder(wide, c(id_columns, term_meta$term_key))
  wide
}

tissue_wide <- make_wide(
  tissue,
  c("analysis_tissue", "tissue_definition"),
  "mean_mission_nes",
  uniqueN(tissue$analysis_tissue)
)
unit_wide <- make_wide(
  unit,
  c(
    "analysis_unit_id", "accession", "mission_cluster",
    "analysis_tissue", "tissue_definition"
  ),
  "NES",
  uniqueN(unit$analysis_unit_id)
)
fwrite(
  tissue_wide,
  file.path(table_dir, "10_tissue_mission_equal_NES_pathway_matrix.csv"),
  bom = TRUE
)
fwrite(
  unit_wide,
  file.path(table_dir, "11_analysis_unit_NES_pathway_matrix.csv"),
  bom = TRUE
)

pairwise_associations <- function(wide, grain, n_label) {
  rows <- vector("list", choose(nrow(term_meta), 2L))
  index <- 0L
  for (i in seq_len(nrow(term_meta) - 1L)) {
    for (j in (i + 1L):nrow(term_meta)) {
      x <- wide[[term_meta$term_key[i]]]
      y <- wide[[term_meta$term_key[j]]]
      pearson <- cor.test(x, y, method = "pearson")
      spearman <- suppressWarnings(
        cor.test(x, y, method = "spearman", exact = FALSE)
      )
      index <- index + 1L
      rows[[index]] <- data.table(
        analysis_grain = grain,
        pathway_1_key = term_meta$term_key[i],
        pathway_1_name = term_meta$term_name[i],
        pathway_1_go_id = term_meta$go_id[i],
        pathway_2_key = term_meta$term_key[j],
        pathway_2_name = term_meta$term_name[j],
        pathway_2_go_id = term_meta$go_id[j],
        n_observations = nrow(wide),
        observation_label = n_label,
        pearson_r = unname(pearson$estimate),
        pearson_p_two_sided = pearson$p.value,
        spearman_rho = unname(spearman$estimate),
        spearman_p_two_sided = spearman$p.value
      )
    }
  }
  result <- rbindlist(rows)
  result[, `:=`(
    pearson_FDR_BH_all_pairs = p.adjust(
      pearson_p_two_sided, method = "BH"
    ),
    spearman_FDR_BH_all_pairs = p.adjust(
      spearman_p_two_sided, method = "BH"
    ),
    n_pathway_pairs_in_FDR_family = nrow(result),
    analysis_id_type = "Ensembl Gene ID",
    symbol_used_as_analysis_key = FALSE,
    interpretation_limit = paste(
      "Association of spaceflight-response NES profiles;",
      "not direct regulation, pathway activity, energy use or causality."
    )
  )]
  result
}

pairwise <- rbindlist(list(
  pairwise_associations(
    tissue_wide, "standardised_tissue_mission_equal",
    "26 standardised tissues"
  ),
  pairwise_associations(
    unit_wide, "registered_analysis_unit",
    "59 registered analysis units"
  )
))
pairwise[, absolute_spearman_rho := abs(spearman_rho)]
setorder(pairwise, analysis_grain, -absolute_spearman_rho)
fwrite(
  pairwise,
  file.path(
    table_dir,
    "12_pairwise_pathway_NES_correlations_tissue_and_analysis_unit.csv"
  ),
  bom = TRUE
)

tissue_pathway_columns <- term_meta$term_key
tissue_cor <- cor(
  as.matrix(tissue_wide[, ..tissue_pathway_columns]),
  method = "spearman"
)
stopifnot(
  all(is.finite(tissue_cor)),
  max(abs(tissue_cor - t(tissue_cor))) < 1e-12
)
fwrite(
  data.table(term_key = rownames(tissue_cor), tissue_cor),
  file.path(table_dir, "13_tissue_pathway_Spearman_rho_matrix.csv"),
  bom = TRUE
)

distance_matrix <- as.dist((1 - tissue_cor) / 2)
cluster <- hclust(distance_matrix, method = "average")
cluster_order <- data.table(
  cluster_position = seq_along(cluster$order),
  term_key = cluster$labels[cluster$order]
)
cluster_order <- merge(
  cluster_order, term_meta,
  by = "term_key", all.x = TRUE, sort = FALSE
)
setorder(cluster_order, cluster_position)
fwrite(
  cluster_order[, !c("term_label", "short_label")],
  file.path(table_dir, "14_tissue_Spearman_hierarchical_cluster_order.csv"),
  bom = TRUE
)

nearest <- rbindlist(lapply(seq_len(nrow(term_meta)), function(i) {
  key <- term_meta$term_key[i]
  candidates <- data.table(
    term_key = key,
    neighbor_key = setdiff(term_meta$term_key, key)
  )
  candidates[, spearman_rho := tissue_cor[key, neighbor_key]]
  candidates <- candidates[order(-spearman_rho)]
  candidates[, neighbor_rank := seq_len(.N)]
  candidates <- candidates[neighbor_rank <= 3L]
  candidates <- merge(
    candidates,
    term_meta[, .(term_key, term_name, go_id)],
    by = "term_key", all.x = TRUE, sort = FALSE
  )
  candidates <- merge(
    candidates,
    term_meta[, .(
      neighbor_key = term_key,
      neighbor_name = term_name,
      neighbor_go_id = go_id
    )],
    by = "neighbor_key", all.x = TRUE, sort = FALSE
  )
  candidates
}))
nearest[, term_display_order := match(term_key, term_meta$term_key)]
setorder(nearest, term_display_order, neighbor_rank)
fwrite(
  nearest,
  file.path(table_dir, "15_top3_pathway_neighbors_tissue_Spearman.csv"),
  bom = TRUE
)

key_pairs <- data.table(
  comparison = c(
    "Translation vs mitochondrion",
    "Cell death vs mechanosensation",
    "Cell death vs cytoskeleton",
    "Cell death vs intrinsic apoptosis"
  ),
  key_1 = c(
    "translation_GO0006412",
    "cell_death_GO0008219",
    "cell_death_GO0008219",
    "cell_death_GO0008219"
  ),
  key_2 = c(
    "mitochondrion_GO0005739",
    "mechanosensation_GO0050954",
    "cytoskeleton_GO0005856",
    "intrinsic_apoptotic_signaling_GO0008630"
  )
)
extract_key_pair <- function(one_pair) {
  result <- pairwise[
    (pathway_1_key == one_pair$key_1 &
       pathway_2_key == one_pair$key_2) |
      (pathway_1_key == one_pair$key_2 &
       pathway_2_key == one_pair$key_1)
  ]
  result[, comparison := one_pair$comparison]
  result
}
key_results <- rbindlist(lapply(
  seq_len(nrow(key_pairs)),
  function(i) extract_key_pair(key_pairs[i])
))
setcolorder(key_results, c(
  "comparison", setdiff(names(key_results), "comparison")
))
fwrite(
  key_results,
  file.path(table_dir, "16_teacher_requested_key_pathway_associations.csv"),
  bom = TRUE
)

ordered_keys <- cluster_order$term_key
heatmap_data <- CJ(
  row_key = ordered_keys,
  column_key = ordered_keys,
  unique = TRUE
)
heatmap_data[, rho := tissue_cor[cbind(row_key, column_key)]]
heatmap_data <- merge(
  heatmap_data,
  term_meta[, .(row_key = term_key, row_label = short_label)],
  by = "row_key", sort = FALSE
)
heatmap_data <- merge(
  heatmap_data,
  term_meta[, .(column_key = term_key, column_label = short_label)],
  by = "column_key", sort = FALSE
)
label_order <- term_meta[match(ordered_keys, term_key), short_label]
heatmap_data[, `:=`(
  row_label = factor(row_label, levels = rev(label_order)),
  column_label = factor(column_label, levels = label_order),
  value_label = sprintf("%.2f", rho)
)]

heatmap_plot <- ggplot(
  heatmap_data,
  aes(x = column_label, y = row_label, fill = rho)
) +
  geom_tile(colour = "white", linewidth = 0.28) +
  geom_text(aes(label = value_label), size = 2.2, colour = "#171717") +
  scale_fill_gradient2(
    low = "#2C6DB2", mid = "#F7F7F7", high = "#C93335",
    midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho"
  ) +
  labs(
    title = "Data-derived pathway relationships across mouse tissues",
    subtitle = paste(
      "Average-linkage order; distance = (1 - Spearman rho) / 2;",
      "26 tissues, mission-equal NES"
    ),
    x = NULL, y = NULL,
    caption = paste(
      "Numbers are cross-tissue Spearman correlations.",
      "They describe co-response, not causality or direct regulation."
    )
  ) +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 48, hjust = 1, vjust = 1, size = 6.6
    ),
    axis.text.y = element_text(size = 6.6),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9.5),
    plot.caption = element_text(size = 8, colour = "#4B5563")
  )
ggsave(
  file.path(
    figure_dir,
    "09_data_derived_pathway_Spearman_clustered_heatmap.png"
  ),
  heatmap_plot, width = 15.2, height = 13.4, dpi = 300, bg = "white"
)
ggsave(
  file.path(
    figure_dir,
    "09_data_derived_pathway_Spearman_clustered_heatmap.pdf"
  ),
  heatmap_plot, width = 15.2, height = 13.4,
  device = cairo_pdf, bg = "white"
)

cluster_labels <- setNames(term_meta$short_label, term_meta$term_key)
cluster_for_plot <- cluster
cluster_for_plot$labels <- cluster_labels[cluster$labels]
save_dendrogram <- function(path, device) {
  if (device == "png") {
    png(path, width = 4800, height = 2800, res = 300)
  } else {
    cairo_pdf(path, width = 16, height = 9.3)
  }
  par(
    mar = c(12, 5, 4.5, 2),
    family = "Arial",
    las = 2,
    cex = 0.78
  )
  plot(
    cluster_for_plot,
    main = "Data-derived hierarchical clustering of mouse pathways",
    sub = paste(
      "26-tissue mission-equal NES; average linkage;",
      "distance = (1 - Spearman rho) / 2"
    ),
    xlab = "",
    ylab = "Correlation distance",
    hang = -1
  )
  mtext(
    "Cluster proximity is descriptive co-response, not a causal edge.",
    side = 1, line = 10.5, cex = 0.8, col = "#4B5563", las = 0
  )
  dev.off()
}
save_dendrogram(
  file.path(
    figure_dir,
    "10_data_derived_pathway_hierarchical_clustering.png"
  ),
  "png"
)
save_dendrogram(
  file.path(
    figure_dir,
    "10_data_derived_pathway_hierarchical_clustering.pdf"
  ),
  "pdf"
)

mds <- cmdscale(distance_matrix, k = 2L, eig = TRUE, add = TRUE)
mds_data <- data.table(
  term_key = rownames(mds$points),
  MDS1 = mds$points[, 1],
  MDS2 = mds$points[, 2]
)
mds_data <- merge(
  mds_data, term_meta,
  by = "term_key", all.x = TRUE, sort = FALSE
)
mds_data[, ontology_domain := fifelse(
  go_id %chin% c("GO:0000781", "GO:0005856", "GO:0005739"),
  "Cellular component",
  "Biological process"
)]
fwrite(
  mds_data[, !c("term_label", "short_label")],
  file.path(table_dir, "17_tissue_Spearman_MDS_coordinates.csv"),
  bom = TRUE
)
mds_plot <- ggplot(
  mds_data,
  aes(x = MDS1, y = MDS2, colour = ontology_domain)
) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#D9D9D9") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#D9D9D9") +
  geom_point(size = 4.2, alpha = 0.95) +
  geom_text_repel(
    aes(label = paste0(term_name, "\n", go_id)),
    seed = 25,
    size = 3.1,
    box.padding = 0.55,
    point.padding = 0.35,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.colour = "#9CA3AF"
  ) +
  scale_colour_manual(values = c(
    "Biological process" = "#2E6E9E",
    "Cellular component" = "#D47A2A"
  )) +
  coord_equal() +
  labs(
    title = "Two-dimensional map of data-derived pathway relationships",
    subtitle = paste(
      "Classical MDS of Spearman correlation distance;",
      "nearby points have similar 26-tissue NES profiles"
    ),
    x = "MDS dimension 1",
    y = "MDS dimension 2",
    colour = NULL,
    caption = paste(
      "This deterministic map is used instead of UMAP because only 15",
      "pathways are available; relative distances are descriptive."
    )
  ) +
  theme_minimal(base_family = "Arial", base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "top",
    plot.caption = element_text(size = 8, colour = "#4B5563")
  )
ggsave(
  file.path(figure_dir, "11_data_derived_pathway_MDS_map.png"),
  mds_plot, width = 13, height = 10, dpi = 300, bg = "white"
)
ggsave(
  file.path(figure_dir, "11_data_derived_pathway_MDS_map.pdf"),
  mds_plot, width = 13, height = 10, device = cairo_pdf, bg = "white"
)

translation_tissue <- key_results[
  comparison == "Translation vs mitochondrion" &
    analysis_grain == "standardised_tissue_mission_equal"
]
death_mechanical_tissue <- key_results[
  comparison == "Cell death vs mechanosensation" &
    analysis_grain == "standardised_tissue_mission_equal"
]
death_cytoskeleton_tissue <- key_results[
  comparison == "Cell death vs cytoskeleton" &
    analysis_grain == "standardised_tissue_mission_equal"
]
death_intrinsic_tissue <- key_results[
  comparison == "Cell death vs intrinsic apoptosis" &
    analysis_grain == "standardised_tissue_mission_equal"
]
death_intrinsic_unit <- key_results[
  comparison == "Cell death vs intrinsic apoptosis" &
    analysis_grain == "registered_analysis_unit"
]
writeLines(c(
  "# 基于本项目数据的通路关系、聚类与二维图",
  "",
  "## 计算对象",
  "",
  paste0(
    "- 主分析：", nrow(tissue_wide),
    " 个标准化组织的任务等权平均 NES；每条通路是一个跨组织响应轮廓。"
  ),
  paste0(
    "- 敏感性分析：", nrow(unit_wide),
    " 个注册分析单元的 NES；保留同一组织内的不同实验/分层。"
  ),
  paste0(
    "- 共 ", choose(nrow(term_meta), 2L),
    " 对通路；每个分析层级分别对全部通路对做 BH-FDR。"
  ),
  "- 分析键始终是 Ensembl Gene ID；Symbol 只用于显示。",
  "- 展示顺序以 tables/00_teacher_requested_term_order.csv 的 display_term_order 为准；原始序号保留。",
  "",
  "## 关系与聚类方法",
  "",
  "- 关联系数：Spearman rho 为主，Pearson r 为辅。",
  "- 聚类距离：(1-rho)/2；average linkage；不人为指定簇数。",
  "- 二维图：对同一相关距离做 classical MDS。只有 15 个通路时，MDS 比依赖邻居数和随机初始化的 UMAP 更稳定、可复现。",
  "",
  "## 老师提出的四个关系",
  "",
  sprintf(
    "- Translation–mitochondrion：组织层 rho=%.3f，P=%.3g，FDR=%.3g；显示两者跨组织响应高度同向。",
    translation_tissue$spearman_rho,
    translation_tissue$spearman_p_two_sided,
    translation_tissue$spearman_FDR_BH_all_pairs
  ),
  sprintf(
    "- Cell death–mechanosensation：组织层 rho=%.3f，P=%.3g，FDR=%.3g。",
    death_mechanical_tissue$spearman_rho,
    death_mechanical_tissue$spearman_p_two_sided,
    death_mechanical_tissue$spearman_FDR_BH_all_pairs
  ),
  sprintf(
    "- Cell death–cytoskeleton：组织层 rho=%.3f，P=%.3g，FDR=%.3g。",
    death_cytoskeleton_tissue$spearman_rho,
    death_cytoskeleton_tissue$spearman_p_two_sided,
    death_cytoskeleton_tissue$spearman_FDR_BH_all_pairs
  ),
  sprintf(
    "- Cell death–intrinsic apoptosis：组织层 rho=%.3f，P=%.3g，FDR=%.3g；但分析单元层 rho=%.3f，P=%.3g，FDR=%.3g，因此“关系不明显”不是跨层级稳健结论。",
    death_intrinsic_tissue$spearman_rho,
    death_intrinsic_tissue$spearman_p_two_sided,
    death_intrinsic_tissue$spearman_FDR_BH_all_pairs,
    death_intrinsic_unit$spearman_rho,
    death_intrinsic_unit$spearman_p_two_sided,
    death_intrinsic_unit$spearman_FDR_BH_all_pairs
  ),
  "",
  "## 可解释边界",
  "",
  "- Translation 与 mitochondrion 的同向共变可作为“能量代谢与翻译响应耦合”的数据线索，但 NES 相关本身不能证明翻译消耗 ATP 是该相关的原因。",
  "- Cell death 与机械感知/细胞骨架的相关支持机械环境相关死亡机制的候选解释，但当前 spaceflight vs ground bulk RNA 数据同时混合微重力、辐射及其他飞行因素；NES 也不是死亡细胞数。因此不能据此定量断言“微重力造成的细胞死亡大于辐照造成的死亡”。",
  "- GO 集合有基因重叠，相关还可能来自共享基因或共享组织构成；因此所有图均描述共响应，不画成直接调控网络。",
  "",
  "## 输出",
  "",
  "- `figures/09_data_derived_pathway_Spearman_clustered_heatmap.png`：带数值的相关热图和聚类顺序。",
  "- `figures/10_data_derived_pathway_hierarchical_clustering.png`：层次聚类树。",
  "- `figures/11_data_derived_pathway_MDS_map.png`：二维关系图。",
  "- `tables/12_pairwise_pathway_NES_correlations_tissue_and_analysis_unit.csv`：完整系数、P 值及 FDR。",
  "- `tables/16_teacher_requested_key_pathway_associations.csv`：老师点名的四组关系。"
), file.path(go_dir, "README_data_derived_pathway_relationships.md"))

message("Completed data-derived pathway relationship analysis: ", go_dir)
