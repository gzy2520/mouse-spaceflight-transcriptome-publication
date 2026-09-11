#!/usr/bin/env Rscript

# Descriptive global-GSEA figures for the current stable-ID rerun.
# These figures deliberately do not pool heterogeneous tissues into a common
# inferential meta-analysis. Ensembl IDs were used upstream; labels below are
# MSigDB gene-set names only.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
set.seed(25)

root <- normalizePath(".", mustWork = TRUE)
global_dir <- file.path(
  root, "03_analysis_results/10_full_stable_id_rerun_20260823/02_global_gsea"
)
term_path <- file.path(global_dir, "03_global_descriptive_term_summary.csv.gz")
mission_path <- file.path(global_dir, "02_mission_descriptive_GSEA_summary.csv.gz")
stopifnot(file.exists(term_path), file.exists(mission_path))

figure_dir <- file.path(global_dir, "figures")
table_dir <- file.path(global_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

term <- fread(term_path)
mission <- fread(mission_path)
expected_libraries <- c(
  "MSigDB_2026.1.Mm_Hallmark",
  "MSigDB_2026.1.Mm_GO_BP",
  "MSigDB_2026.1.Mm_Reactome",
  "MSigDB_2026.1.Mm_WikiPathways"
)
stopifnot(all(expected_libraries %chin% unique(term$library)))

library_label <- c(
  "MSigDB_2026.1.Mm_Hallmark" = "Hallmark",
  "MSigDB_2026.1.Mm_GO_BP" = "GO BP",
  "MSigDB_2026.1.Mm_Reactome" = "Reactome",
  "MSigDB_2026.1.Mm_WikiPathways" = "WikiPathways"
)
clean_term_label <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  x
}

term[, `:=`(
  library_display = factor(library_label[library], levels = unname(library_label)),
  term_display = clean_term_label(gene_set_name)
)]
stopifnot(all(is.finite(term$mean_of_mission_mean_NES)))

# One figure contains one MSigDB database only. Each figure selects six terms
# at each NES extreme from that database using descriptive global mean NES.
top_n <- 6L
file_stem <- c(
  "MSigDB_2026.1.Mm_Hallmark" = "01_Hallmark_global_top_terms",
  "MSigDB_2026.1.Mm_GO_BP" = "02_GO_BP_global_top_terms",
  "MSigDB_2026.1.Mm_Reactome" = "03_Reactome_global_top_terms",
  "MSigDB_2026.1.Mm_WikiPathways" = "04_WikiPathways_global_top_terms"
)
selection <- rbindlist(lapply(expected_libraries, function(one_library) {
  # Require full mission coverage for a cross-mission display ranking. Terms
  # absent after fgsea size filtering in one or more missions stay in tables.
  x <- term[library == one_library & n_missions == 13L]
  stopifnot(nrow(x) >= 2L * top_n)
  selected <- unique(rbindlist(list(
    head(x[order(-mean_of_mission_mean_NES, gene_set_name)], top_n),
    head(x[order(mean_of_mission_mean_NES, gene_set_name)], top_n)
  )), by = "gene_set_name")
  selected[, `:=`(
    direction = fifelse(mean_of_mission_mean_NES >= 0, "Positive NES", "Negative NES"),
    figure_file_stem = file_stem[[one_library]]
  )]
  selected[, term_plot := factor(
    term_display,
    levels = selected[order(mean_of_mission_mean_NES), term_display]
  )]
  plot <- ggplot(selected, aes(x = mean_of_mission_mean_NES, y = term_plot, fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#4B5563") +
    scale_fill_manual(values = c("Positive NES" = "#B2182B", "Negative NES" = "#2166AC")) +
    labs(
      title = paste0(library_label[[one_library]], ": descriptively strongest global GSEA terms"),
      subtitle = "Six positive and six negative terms with complete 13-mission coverage; mission-equal mean NES",
      x = "Mean mission-level NES (positive = higher in flight/uG)", y = NULL, fill = NULL,
      caption = "Descriptive summary across heterogeneous tissues; not a pooled effect estimate or a significance test."
    ) +
    theme_minimal(base_family = "Arial", base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8), legend.position = "top",
      plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 8, colour = "#4B5563"),
      plot.margin = margin(8, 16, 8, 8)
    )
  stem <- file.path(figure_dir, file_stem[[one_library]])
  ggsave(paste0(stem, ".png"), plot, width = 10, height = 7.2, dpi = 320, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = 10, height = 7.2, device = cairo_pdf, bg = "white")
  selected
}), use.names = TRUE, fill = TRUE)
fwrite(
  selection[, .(
    library, figure_file_stem, gene_set_name, term_display, n_missions,
    mean_of_mission_mean_NES, median_of_mission_mean_NES,
    positive_missions, negative_missions, direction
  )],
  file.path(table_dir, "01_global_top_term_figure_selection.csv")
)

# Dataset-resolved figures follow the same one-figure/one-database rule.
# A dataset means one registered analysis unit, not a Symbol-collapsed or
# accession-merged result. Each panel uses the library-specific fgsea BH FDR.
dataset_path <- file.path(global_dir, "01_dataset_GSEA_mouse_MSigDB_Ensembl.csv.gz")
dataset <- fread(dataset_path)
dataset <- dataset[library %chin% expected_libraries]
dataset_metadata <- unique(dataset[, .(
  analysis_unit_id, accession, mission_cluster, biological_experiment_id,
  analysis_tissue
)])
stopifnot(
  nrow(dataset_metadata) == 59L,
  !anyDuplicated(dataset_metadata$analysis_unit_id),
  uniqueN(dataset[, .(analysis_unit_id, library)]) == 59L * 4L
)
dataset[, term_display := clean_term_label(gene_set_name)]
dataset_figure_dir <- file.path(figure_dir, "by_analysis_unit")
dir.create(dataset_figure_dir, recursive = TRUE, showWarnings = FALSE)
dataset_top_n <- 10L
dataset_selection <- rbindlist(lapply(seq_len(nrow(dataset_metadata)), function(i) {
  meta <- dataset_metadata[i]
  unit <- dataset[analysis_unit_id == meta$analysis_unit_id & padj < 0.05]
  selected_unit <- rbindlist(lapply(expected_libraries, function(one_library) {
    x <- unit[library == one_library]
    stopifnot(nrow(x) > 0L)
    selected <- unique(rbindlist(list(
      head(x[order(-NES, gene_set_name)], dataset_top_n),
      head(x[order(NES, gene_set_name)], dataset_top_n)
    )), by = "gene_set_name")
    selected[, `:=`(
      direction = fifelse(NES >= 0, "Positive NES", "Negative NES"),
      term_plot = factor(term_display,
        levels = selected[order(NES), term_display]
      )
    )]
    unit_dir <- file.path(dataset_figure_dir, meta$analysis_unit_id)
    dir.create(unit_dir, recursive = TRUE, showWarnings = FALSE)
    stem <- file.path(unit_dir, paste0(file_stem[[one_library]], "_FDR005"))
    plot <- ggplot(selected, aes(x = NES, y = term_plot, fill = direction)) +
      geom_col(width = 0.72) +
      geom_vline(xintercept = 0, linewidth = 0.35, colour = "#4B5563") +
      scale_fill_manual(values = c("Positive NES" = "#B2182B", "Negative NES" = "#2166AC")) +
      labs(
        title = paste0(meta$analysis_unit_id, " | ", library_label[[one_library]]),
        subtitle = paste0(
          meta$analysis_tissue, " | ", meta$mission_cluster,
          " | top positive and negative terms with library-specific FDR < 0.05"
        ),
        x = "NES (positive = higher in flight/uG)", y = NULL, fill = NULL,
        caption = paste0(
          "Registered analysis unit ", meta$analysis_unit_id,
          "; FDR is fgsea BH-adjusted within this dataset and database."
        )
      ) +
      theme_minimal(base_family = "Arial", base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7), legend.position = "top",
        plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 9),
        plot.caption = element_text(size = 7.5, colour = "#4B5563"),
        plot.margin = margin(8, 16, 8, 8)
      )
    ggsave(paste0(stem, ".png"), plot, width = 10, height = 7.2, dpi = 320, bg = "white")
    ggsave(paste0(stem, ".pdf"), plot, width = 10, height = 7.2, device = cairo_pdf, bg = "white")
    selected[, `:=`(
      figure_relative_path = file.path(
        "figures", "by_analysis_unit", meta$analysis_unit_id,
        paste0(file_stem[[one_library]], "_FDR005.png")
      ),
      accession = meta$accession,
      mission_cluster = meta$mission_cluster,
      biological_experiment_id = meta$biological_experiment_id,
      analysis_tissue = meta$analysis_tissue
    )]
    selected
  }), use.names = TRUE, fill = TRUE)
  selected_unit
}), use.names = TRUE, fill = TRUE)
fwrite(
  dataset_selection[, .(
    analysis_unit_id, accession, mission_cluster, biological_experiment_id,
    analysis_tissue, library, figure_relative_path, gene_set_name, term_display,
    NES, pval, padj, size, direction
  )],
  file.path(table_dir, "02_dataset_library_FDR005_figure_selection.csv")
)
figure_manifest <- CJ(
  analysis_unit_id = dataset_metadata$analysis_unit_id,
  library = expected_libraries,
  unique = TRUE
)
figure_manifest <- merge(figure_manifest, dataset_metadata,
  by = "analysis_unit_id", all.x = TRUE, sort = FALSE
)
figure_manifest[, `:=`(
  png = file.path(
    "figures", "by_analysis_unit", analysis_unit_id,
    paste0(file_stem[library], "_FDR005.png")
  ),
  pdf = file.path(
    "figures", "by_analysis_unit", analysis_unit_id,
    paste0(file_stem[library], "_FDR005.pdf")
  )
)]
stopifnot(
  nrow(figure_manifest) == 236L,
  all(file.exists(file.path(global_dir, figure_manifest$png))),
  all(file.exists(file.path(global_dir, figure_manifest$pdf)))
)
fwrite(figure_manifest, file.path(table_dir, "03_dataset_library_figure_manifest.csv"))
accession_index <- figure_manifest[, .(
  n_registered_analysis_units = uniqueN(analysis_unit_id),
  analysis_unit_ids = paste(sort(unique(analysis_unit_id)), collapse = ";"),
  mission_clusters = paste(sort(unique(mission_cluster)), collapse = ";"),
  biological_experiment_ids = paste(sort(unique(biological_experiment_id)), collapse = ";"),
  analysis_tissues = paste(sort(unique(analysis_tissue)), collapse = ";"),
  figure_directories = paste(
    sort(unique(dirname(png))), collapse = ";"
  )
), by = accession]
stopifnot(nrow(accession_index) == 48L)
fwrite(accession_index, file.path(table_dir, "04_OSDR_accession_to_analysis_unit_figure_index.csv"))

readme <- c(
  "# 全局GSEA图形说明",
  "",
  "本目录图形基于稳定ID重跑的四个小鼠原生MSigDB库和59个主分析单元。上游GSEA主键为Ensembl Gene ID；图中术语名仅为显示标签。",
  "",
  "- 严格一图一数据库：Hallmark、GO BP、Reactome和WikiPathways各有一张图；每张图仅从在全部13个mission cluster中都有结果的术语中，按任务等权平均NES分别展示六个正向和六个负向描述性最强术语。",
  "- `figures/by_analysis_unit/`：59个注册分析单元各自按四个数据库分别出图，共236张PNG和236张PDF。主图层必须按分析单元保留，因为7个OSDR accession含2–4个预设分层比较，不能为凑48个 accession 而合并不同对比并重算FDR。每张图只显示该分析单元、该数据库内FDR<0.05的术语，并按NES两端各取最多10条。",
  "- `tables/04_OSDR_accession_to_analysis_unit_figure_index.csv`：48个原始OSDR accession的导航表，列出各自对应的分析单元、任务/组织信息与图形目录；用于从48个数据来源定位59个统计上独立的比较。",
  "",
  "跨组织四图为描述性概览，不能把平均NES或条形长度解释为共同组织效应、显著性或通路激活/抑制。逐分析单元图的FDR只适用于各自数据库内的fgsea多重校正。"
)
writeLines(readme, file.path(global_dir, "README.md"))
message("Created global GSEA figures in: ", figure_dir)
