#!/usr/bin/env Rscript

# Build the teacher-requested all-tissue seven-pathway member-gene count
# figure.  Counts are calculated independently for each source sample column
# (flight and ground are retained as a sample-status annotation), rather than
# from a tissue-pooled expression value.  Ensembl Gene ID is the only key.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})
set.seed(25)

root <- normalizePath(Sys.getenv("GRADUATE_DESIGN_ROOT", unset = getwd()), mustWork = TRUE)
source_root <- file.path(root, "03_analysis_results/10_full_stable_id_rerun_20260823/04_seven_pathway_response_26_tissues_20260824")
catalog_root <- Sys.getenv(
  "SEVEN_PATHWAY_CATALOG_ROOT",
  unset = file.path(root, "03_analysis_results/10_full_stable_id_rerun_20260823/04_seven_pathway_modified_contract_updated_20260824")
)
if (!grepl("^/", catalog_root)) catalog_root <- file.path(root, catalog_root)
out_dir <- Sys.getenv(
  "SEVEN_PATHWAY_COUNT_OUTPUT_DIR",
  unset = file.path(root, "03_analysis_results/24_teacher_selected_figure_release_20260829/seven_pathway_counts")
)
if (!grepl("^/", out_dir)) out_dir <- file.path(root, out_dir)
figure_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
input_dir <- file.path(out_dir, "input")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
normalize_ensembl <- function(x) sub("\\.[0-9]+$", "", toupper(trimws(as.character(x))))
split_samples <- function(x) {
  if (is.na(x) || !nzchar(trimws(as.character(x)))) return(character())
  z <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  z[nzchar(z)]
}
collapse_max <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}
sha256_file <- function(path) {
  ans <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", ans[[1L]])
}

pathway_levels <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
pathway_shapes <- c(
  BER = 24L,   # hollow triangle up
  NER = 25L,   # hollow triangle down
  MMR = 23L,   # hollow diamond
  FA = 0L,     # hollow square
  HR = 8L,     # star
  AEJ = 3L,    # plus
  NHEJ = 4L    # x
)
pathway_legend_labels <- c(
  BER = "BER  Base excision",
  NER = "NER  Nucleotide excision",
  MMR = "MMR  Mismatch repair",
  FA = "FA  Fanconi anemia",
  HR = "HR  Homologous recombination",
  AEJ = "A-EJ  Alternative end joining",
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
tissue_display <- c(
  "Adrenal gland" = "Adrenal gland", "Bone marrow" = "Bone marrow", "Cecum" = "Cecum",
  "Cerebellum" = "Cerebellum", "Colon" = "Colon", "Dorsal skin" = "Dorsal skin",
  "Extensor digitorum longus" = "Extensor digitorum longus", "Eye" = "Eye",
  "Femoral lateral skin" = "Femoral lateral skin", "Femoral skin" = "Femoral skin",
  "Gastrocnemius" = "Gastrocnemius", "Heart" = "Heart",
  "Heart / Heart right ventricle" = "Heart right ventricle", "Kidney" = "Kidney",
  "Left lobe of the liver" = "Left lobe of the liver", "Liver" = "Liver", "Lung" = "Lung",
  "Mammary gland" = "Mammary gland", "Optic nerve" = "Optic nerve",
  "Quadriceps femoris" = "Quadriceps femoris", "Retina" = "Retina", "Soleus" = "Soleus",
  "Spleen" = "Spleen", "Spleen-distal" = "Spleen-distal", "Thymus" = "Thymus",
  "Tibialis anterior" = "Tibialis anterior"
)

go_path <- file.path(catalog_root, "skill_output/gene_go_pathway_long.csv")
unit_path <- file.path(source_root, "prepared_inputs/01_selected_tissue_analysis_units.csv")
sample_audit_path <- file.path(
  root, "03_analysis_results/06_cross_tissue_integration/05_proliferation_state/tables",
  "02_analysis_unit_expression_and_group_audit.csv"
)
assert(file.exists(go_path) && file.exists(unit_path) && file.exists(sample_audit_path),
       "Missing seven-pathway or sample-membership input")

go <- fread(go_path, showProgress = FALSE)
assert(all(c("EnsemblID", "Pathway") %chin% names(go)), "GO member table schema changed")
go[, EnsemblID := normalize_ensembl(EnsemblID)]
go <- unique(go[Pathway %chin% pathway_levels, .(EnsemblID, Pathway)])
pathway_sets <- lapply(pathway_levels, function(p) unique(go[Pathway == p, EnsemblID]))
names(pathway_sets) <- pathway_levels
pathway_denominators <- data.table(
  Pathway = pathway_levels,
  n_unique_Ensembl_genes_in_GO_set = vapply(pathway_sets, length, integer(1L))
)
assert(all(pathway_denominators$n_unique_Ensembl_genes_in_GO_set > 0L), "An empty pathway gene set was found")
expected_text_catalog <- c(BER = 54L, NER = 53L, MMR = 27L, FA = 39L, HR = 135L, AEJ = 9L, NHEJ = 18L)
assert(setequal(pathway_denominators$Pathway, names(expected_text_catalog)), "Unexpected pathway catalog levels")
observed_text_catalog <- setNames(pathway_denominators$n_unique_Ensembl_genes_in_GO_set, pathway_denominators$Pathway)
assert(all(observed_text_catalog[names(expected_text_catalog)] == expected_text_catalog),
       paste0("The selected GO catalog does not match the teacher-text member counts: ",
              paste(names(observed_text_catalog), observed_text_catalog, sep = "=", collapse = ", ")))

units <- fread(unit_path, showProgress = FALSE)
audit <- fread(sample_audit_path, showProgress = FALSE)
required_units <- c("Group", "analysis_unit_id", "accession", "input_relative_path")
assert(all(required_units %chin% names(units)), "Selected tissue-unit manifest schema changed")
required_audit <- c("analysis_unit_id", "accession", "input_relative_path", "flight_samples", "control_samples",
                    "n_flight_samples", "n_control_samples")
assert(all(required_audit %chin% names(audit)), "Sample-membership audit schema changed")
units <- merge(
  units[, ..required_units],
  audit[, ..required_audit],
  by = c("analysis_unit_id", "accession", "input_relative_path"),
  all.x = TRUE, sort = FALSE
)
assert(!anyNA(units$Group) && nrow(units) == 59L && uniqueN(units$Group) == 26L,
       "Expected 59 analysis units across 26 tissues")
assert(setequal(unique(units$Group), tissue_levels), "Tissue list is not the fixed 26-tissue contract")
assert(!anyNA(units$flight_samples) && !anyNA(units$control_samples), "Missing sample lists")

target_ids <- unique(go$EnsemblID)
rows <- list()
audit_rows <- list()
source_groups <- split(units, units$input_relative_path)

for (source_path in names(source_groups)) {
  part <- source_groups[[source_path]]
  full_path <- file.path(root, source_path)
  assert(file.exists(full_path), paste("Missing source expression file:", full_path))
  header <- names(fread(full_path, nrows = 0L, showProgress = FALSE, check.names = FALSE))
  id_candidates <- c("ENSEMBL", "ENSEMBL_ID", "GENE_ID", "GENEID")
  id_col <- header[match(id_candidates, toupper(trimws(header)))][1L]
  assert(!is.na(id_col), paste("No Ensembl ID column in", source_path))
  all_samples <- unique(c(
    unlist(lapply(part$flight_samples, split_samples), use.names = FALSE),
    unlist(lapply(part$control_samples, split_samples), use.names = FALSE)
  ))
  assert(length(all_samples) > 0L, paste("No source samples for", source_path))
  missing_samples <- setdiff(all_samples, header)
  assert(!length(missing_samples), paste("Missing sample columns in", source_path, paste(head(missing_samples), collapse = ", ")))
  frame <- fread(full_path, select = c(id_col, all_samples), showProgress = FALSE, check.names = FALSE)
  setnames(frame, id_col, "EnsemblID")
  frame[, EnsemblID := normalize_ensembl(EnsemblID)]
  frame <- frame[EnsemblID %chin% target_ids]
  n_duplicate_rows <- frame[, sum(duplicated(EnsemblID))]
  if (n_duplicate_rows > 0L) {
    frame <- frame[, lapply(.SD, collapse_max), by = EnsemblID]
  }
  frame <- frame[!duplicated(EnsemblID)]
  setkey(frame, EnsemblID)
  found_ids <- sum(target_ids %chin% frame$EnsemblID)

  for (i in seq_len(nrow(part))) {
    u <- part[i]
    for (status in c("Flight", "Ground")) {
      sample_col <- if (status == "Flight") "flight_samples" else "control_samples"
      samples <- split_samples(u[[sample_col]])
      if (!length(samples)) next
      for (sample_name in samples) {
        one <- data.table(
          Group = as.character(u$Group), analysis_unit_id = as.character(u$analysis_unit_id),
          accession = as.character(u$accession), input_relative_path = as.character(u$input_relative_path),
          sample_status = status, sample_column = sample_name,
          sample_is_technical_replicate = grepl("techrep|technical", sample_name, ignore.case = TRUE)
        )
        for (pathway in pathway_levels) {
          ids <- pathway_sets[[pathway]]
          values <- as.numeric(frame[match(ids, EnsemblID), ..sample_name][[1L]])
          valid <- is.finite(values)
          present <- valid & (values - 1 > 0)
          one[[paste0("n_", pathway)]] <- sum(present)
        }
        rows[[length(rows) + 1L]] <- one
      }
    }
  }
  audit_rows[[length(audit_rows) + 1L]] <- data.table(
    input_relative_path = source_path, input_sha256 = sha256_file(full_path),
    n_analysis_units = nrow(part), n_source_samples = length(all_samples),
    n_target_Ensembl_genes = length(target_ids), n_target_genes_found = found_ids,
    n_duplicate_target_rows_collapsed = n_duplicate_rows,
    presence_rule = "finite source value - 1 > 0; unique Ensembl IDs; no imputation",
    go_catalog_source = normalizePath(catalog_root, mustWork = TRUE),
    go_catalog_sha256 = sha256_file(go_path)
  )
}

counts_wide <- rbindlist(rows, fill = TRUE)
assert(nrow(counts_wide) == 761L, paste("Expected 761 source sample rows; observed", nrow(counts_wide)))
assert(uniqueN(counts_wide[, .(analysis_unit_id, sample_column, sample_status)]) == nrow(counts_wide),
       "Duplicate source sample identifiers across analysis units")
assert(all(vapply(pathway_levels, function(p) all(is.finite(counts_wide[[paste0("n_", p)]])), logical(1L))),
       "Non-finite per-sample counts")

fwrite(counts_wide, file.path(table_dir, "01_seven_pathway_member_gene_counts_per_source_sample_wide.csv"), bom = TRUE)
counts_long <- melt(
  counts_wide,
  id.vars = c("Group", "analysis_unit_id", "accession", "input_relative_path", "sample_status", "sample_column", "sample_is_technical_replicate"),
  measure.vars = paste0("n_", pathway_levels),
  variable.name = "Pathway", value.name = "n_member_genes"
)
counts_long[, Pathway := sub("^n_", "", Pathway)]
counts_long[, Pathway := factor(Pathway, levels = pathway_levels)]
counts_long[, TissueOrder := match(Group, tissue_levels)]
counts_long[, TissueLabel := unname(tissue_display[Group])]
setorder(counts_long, TissueOrder, sample_status, sample_column, Pathway)
fwrite(counts_long, file.path(table_dir, "02_seven_pathway_member_gene_counts_per_source_sample_long.csv"), bom = TRUE)

summary_status <- counts_long[, .(
  n_samples = .N, n_technical_replicate_columns = sum(sample_is_technical_replicate),
  min = as.numeric(min(n_member_genes)), q1 = as.numeric(quantile(n_member_genes, 0.25)),
  median = as.numeric(median(n_member_genes)), mean = as.numeric(mean(n_member_genes)),
  q3 = as.numeric(quantile(n_member_genes, 0.75)), max = as.numeric(max(n_member_genes))
), by = .(Group, TissueLabel, TissueOrder, sample_status, Pathway)]
summary_all <- counts_long[, .(
  n_samples = .N, n_flight_samples = sum(sample_status == "Flight"), n_ground_samples = sum(sample_status == "Ground"),
  min = as.numeric(min(n_member_genes)), q1 = as.numeric(quantile(n_member_genes, 0.25)),
  median = as.numeric(median(n_member_genes)), mean = as.numeric(mean(n_member_genes)),
  q3 = as.numeric(quantile(n_member_genes, 0.75)), max = as.numeric(max(n_member_genes))
), by = .(Group, TissueLabel, TissueOrder, Pathway)]
summary_all <- merge(summary_all, pathway_denominators, by = "Pathway", all.x = TRUE, sort = FALSE)
setorder(summary_status, TissueOrder, sample_status, Pathway)
setorder(summary_all, TissueOrder, Pathway)
fwrite(summary_status, file.path(table_dir, "03_seven_pathway_member_gene_count_summary_by_status.csv"), bom = TRUE)
fwrite(summary_all, file.path(table_dir, "04_seven_pathway_member_gene_count_summary_all_statuses.csv"), bom = TRUE)
fwrite(pathway_denominators, file.path(table_dir, "05_pathway_gene_set_denominators.csv"), bom = TRUE)
fwrite(rbindlist(audit_rows, fill = TRUE), file.path(table_dir, "06_source_expression_read_audit.csv"), bom = TRUE)

# Table S3 is the compact repertoire-size table cited in the teacher's paragraph.
# The per-tissue/per-source-sample distributions remain in tables 03/04.
table_s3_catalog <- data.table(
  Pathway = pathway_levels,
  PathwayOrder = seq_along(pathway_levels),
  MemberGeneCount = observed_text_catalog[pathway_levels],
  DescriptionInTeacherText = c(
    "BER (>50 genes)", "NER (>50 genes)", "MMR (~30 genes)", "FA (~40 genes)",
    "HR (>130 genes)", "A-EJ (~10 genes)", "NHEJ (~20 genes)"
  )
)
fwrite(table_s3_catalog, file.path(table_dir, "08_Table_S3_pathway_repertoire_sizes.csv"), bom = TRUE)
table_s3_tissue <- copy(summary_all)
setcolorder(table_s3_tissue, c(
  "Group", "TissueLabel", "TissueOrder", "Pathway", "n_unique_Ensembl_genes_in_GO_set",
  "n_samples", "n_flight_samples", "n_ground_samples", "min", "q1", "median", "mean", "q3", "max"
))
fwrite(table_s3_tissue, file.path(table_dir, "09_Table_S3_per_tissue_sample_summary.csv"), bom = TRUE)

plot_data <- counts_long[, .(
  Tissue = factor(TissueLabel, levels = unname(tissue_display[tissue_levels])),
  TissueOrder,
  Pathway = factor(as.character(Pathway), levels = pathway_levels),
  sample_status, sample_column, n_member_genes
)]
# Draw Ground first and Flight second so the biologically focal Flight symbols
# remain visible when source-sample values overlap exactly.
plot_data[, sample_status := factor(sample_status, levels = c("Ground", "Flight"))]
setorder(plot_data, TissueOrder, sample_status, sample_column, Pathway)
max_y <- ceiling((max(plot_data$n_member_genes) + 12) / 10) * 10

# Each tissue is a framed mini-panel.  The x positions remain fixed, but their
# pathway names are intentionally hidden; the pathway is encoded by the point
# shape and decoded in the dedicated two-cell legend at the lower right.
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
      fill = "#E2E8F0", colour = "#475569", linewidth = 0.75
    ) +
    geom_point(
      aes(shape = Pathway, colour = sample_status),
      position = position_jitter(width = 0.39, height = 0, seed = 25),
      size = 2.75, alpha = 0.94, stroke = 0.85, show.legend = FALSE
    ) +
    geom_text(
      data = panel_summary,
      aes(x = Pathway, y = label_y, label = round(median_count)),
      inherit.aes = FALSE, size = 4.10, colour = "#7A2E2A", fontface = "bold"
    ) +
    scale_shape_manual(values = pathway_shapes, limits = pathway_levels, drop = FALSE) +
    scale_colour_manual(
      values = c(Flight = "#D55E00", Ground = "#0072B2"), guide = "none"
    ) +
    scale_x_discrete(
      limits = pathway_levels, labels = rep("", length(pathway_levels)), drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, max_y), breaks = pretty_breaks(n = 4),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = tissue_name, x = NULL, y = if (show_y_title) "Member genes" else NULL) +
    theme_minimal(base_size = 8, base_family = "Arial Unicode MS") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#E5E7EB", linewidth = 0.25),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = if (show_y_axis) element_text(size = 7.5, colour = "#374151") else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line(colour = "#6B7280", linewidth = 0.35) else element_blank(),
      axis.title.y = if (show_y_title) element_text(size = 9.5, colour = "#374151") else element_blank(),
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5, colour = "#1F2937",
                                margin = margin(2, 1, 2, 1)),
      plot.background = element_rect(fill = "white", colour = "#374151", linewidth = 1.25),
      panel.border = element_rect(fill = NA, colour = "#CBD5E1", linewidth = 0.55),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5)
    )
}

legend_data <- data.table(
  Pathway = factor(pathway_levels, levels = pathway_levels),
  y = rev(seq_along(pathway_levels)),
  label = unname(pathway_legend_labels[pathway_levels])
)
status_data <- data.table(
  status = factor(c("Flight", "Ground"), levels = c("Flight", "Ground")),
  y = c(1.05, 0.62),
  label = c("Flight sample", "Ground sample")
)
legend_plot <- ggplot() +
  geom_point(
    data = legend_data, aes(x = 0.20, y = y, shape = Pathway),
    size = 5.2, colour = "#374151", fill = NA, show.legend = FALSE
  ) +
  geom_text(
    data = legend_data, aes(x = 0.48, y = y, label = label),
    hjust = 0, size = 3.6, colour = "#374151"
  ) +
  geom_point(
    data = status_data, aes(x = 2.22, y = y, colour = status),
    shape = 16, size = 4.5, show.legend = FALSE
  ) +
  geom_text(
    data = status_data, aes(x = 2.38, y = y, label = label),
    hjust = 0, size = 3.4, colour = "#374151"
  ) +
  annotate("text", x = 0.04, y = 8.28, label = "Pathway symbols", hjust = 0,
           size = 4.7, fontface = "bold", colour = "#111827") +
  annotate("text", x = 0.04, y = 7.88, label = "Symbol = pathway", hjust = 0,
           size = 3.2, colour = "#6B7280") +
           annotate("text", x = 2.07, y = 1.53, label = "Colour = sample status", hjust = 0,
           size = 3.2, fontface = "bold", colour = "#111827") +
  scale_shape_manual(values = pathway_shapes, limits = pathway_levels, drop = FALSE) +
  scale_colour_manual(values = c(Flight = "#D55E00", Ground = "#0072B2")) +
  scale_x_continuous(limits = c(0, 3.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.25, 8.65), expand = c(0, 0)) +
  theme_void(base_size = 8, base_family = "Arial Unicode MS") +
  theme(
    plot.background = element_rect(fill = "white", colour = "#9CA3AF", linewidth = 0.55),
    plot.margin = margin(5, 5, 5, 5)
  )

tissue_plots <- setNames(
  lapply(seq_along(tissue_levels), make_tissue_plot), LETTERS[seq_along(tissue_levels)]
)
plot_list <- c(tissue_plots, list(a = legend_plot))
layout_design <- paste(
  c("ABCD", "EFGH", "IJKL", "MNOP", "QRST", "UVWX", "YZaa"),
  collapse = "\n"
)
plot <- wrap_plots(plot_list, design = layout_design) +
  plot_annotation(
    title = "Seven DNA-repair pathway member-gene counts per source sample",
    subtitle = paste0(
      "26 mouse tissues in a 4 × 7 layout; each tissue is framed, and the lower-right two cells contain the symbol legend.\n",
      "Each point is one source sample column; symbol = pathway, colour = sample status; median labels are shown.\n",
      "Presence = finite source value − 1 > 0; Ensembl Gene IDs are counted once per pathway.\n",
      "Teacher-text/Table S3 catalog sizes: BER 54 | NER 53 | MMR 27 | FA 39 | HR 135 | A-EJ 9 | NHEJ 18."
    ),
    caption = paste(
      "Counts describe repertoire presence, not pathway activity or expression magnitude. Source column labels, including technical-replicate markers, are retained in the audit table; no imputation was performed."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 20.5, hjust = 0, colour = "#111827"),
      plot.subtitle = element_text(size = 11.5, colour = "#374151", hjust = 0),
      plot.caption = element_text(size = 9.2, colour = "#4B5563", hjust = 0),
      plot.margin = margin(7, 10, 6, 10)
    )
  )

png_path <- file.path(figure_dir, "Fig_S5_seven_pathway_member_gene_counts_per_source_sample_26_tissues_4x7_symbols.png")
pdf_path <- file.path(figure_dir, "Fig_S5_seven_pathway_member_gene_counts_per_source_sample_26_tissues_4x7_symbols.pdf")
ggsave(png_path, plot, width = 18.5, height = 26.5, units = "in", dpi = 320, bg = "white", limitsize = FALSE)
ggsave(pdf_path, plot, width = 18.5, height = 26.5, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)

# Keep the teacher-selected release aliases synchronized with the new layout.
release_root <- Sys.getenv(
  "TEACHER_RELEASE_ROOT",
  unset = file.path(root, "03_analysis_results/24_teacher_selected_figure_release_20260829")
)
release_suppl_dir <- file.path(release_root, "Suppl")
if (dir.exists(release_suppl_dir)) {
  file.copy(png_path, file.path(release_suppl_dir, "Fig_S5_4x7_symbols.png"), overwrite = TRUE)
  file.copy(pdf_path, file.path(release_suppl_dir, "Fig_S5_4x7_symbols.pdf"), overwrite = TRUE)
  file.copy(png_path, file.path(release_suppl_dir, "Fig_S5.png"), overwrite = TRUE)
  file.copy(pdf_path, file.path(release_suppl_dir, "Fig_S5.pdf"), overwrite = TRUE)
}

run_audit <- data.table(
  tissue_count = uniqueN(counts_long$Group), analysis_unit_count = uniqueN(counts_long$analysis_unit_id),
  source_sample_count = uniqueN(counts_long[, .(analysis_unit_id, sample_column, sample_status)]),
  flight_sample_count = uniqueN(counts_long[sample_status == "Flight", .(analysis_unit_id, sample_column)]),
  ground_sample_count = uniqueN(counts_long[sample_status == "Ground", .(analysis_unit_id, sample_column)]),
  pathway_count = uniqueN(counts_long$Pathway), presence_rule = "finite source value - 1 > 0",
  id_key = "Ensembl Gene ID", grid_layout = "4 columns x 7 rows; 26 tissue panels + 2-cell symbol legend",
  pathway_marker_encoding = "shape; BER triangle-up, NER triangle-down, MMR diamond, FA square, HR star, AEJ plus, NHEJ x",
  sample_status_encoding = "outline colour (hollow symbols); Ground drawn first (blue #0072B2), Flight drawn last (orange #D55E00)",
  png = file.path("figures", basename(png_path)), pdf = file.path("figures", basename(pdf_path))
)
fwrite(run_audit, file.path(table_dir, "07_run_audit.csv"), bom = TRUE)

readme_lines <- c(
  "# 26 组织七条 DNA 修复通路：按单个来源样品统计成员基因数",
  "",
  "## 统计口径",
  "",
  "- 组织范围为当前 26 个组织、59 个 analysis unit；图采用 4 × 7 宫格，每格一个组织，右下两格为通路符号和样品状态图例。",
  "- 每个点对应一个 GeneLab/OSDR 来源表达表中的一个样品列；点的形状表示通路，Flight 和 Ground 用颜色表示样品状态，不把样品列假定为互相独立的生物学重复。带有 `techrep` 的列在表中保留标记。",
  "- 图形符号为 BER 空心上三角、NER 空心下三角、MMR 空心菱形、FA 空心方形、HR 星号、A-EJ 加号、NHEJ 叉号；轮廓颜色表示样品状态：Flight 使用橙色，Ground 使用蓝色。绘图时 Ground 先画、Flight 后画，以避免 Ground 遮挡 Flight。",
  "- 对每个样品、每个通路，统计 GO 七通路成员中满足“有限 source value − 1 > 0”的唯一 Ensembl Gene ID 数；不使用 Symbol 连接，不做缺失值填补。",
  "- GO 成员集合使用老师文字描述/示例下方成员数柱图对应的冻结七通路目录；该目录的成员基数为 BER=54、NER=53、MMR=27、FA=39、HR=135、AEJ=9、NHEJ=18。一个基因在同一通路内只计数一次。",
  "- 为保持 26 组织要求，以上文字图目录不变，仅将其应用到当前 26 组织的来源样品列；因此图既复现文字图的通路规模口径，又覆盖全部 26 个组织。",
  "- 这是成员基因存在/检测数量的描述性统计，不等于通路活性、表达总量或 DNA 修复能力。",
  "",
  "## 输出",
  "",
  "- `figures/Fig_S5_seven_pathway_member_gene_counts_per_source_sample_26_tissues_4x7_symbols.png/pdf`：26 组织 4 × 7 宫格图；每格显示 7 通路的单样品符号、箱线和中位数数字，右下两格为图例。",
  "- `tables/01_*_wide.csv`：每个来源样品一行、七通路计数列。",
  "- `tables/02_*_long.csv`：绘图用长表。",
  "- `tables/03_*_by_status.csv`、`04_*_all_statuses.csv`：分状态和合并状态汇总。",
  "- `tables/08_Table_S3_pathway_repertoire_sizes.csv`：与老师文字中 Fig S5/Table S3 对应的七通路目录基数。",
  "- `tables/09_Table_S3_per_tissue_sample_summary.csv`：26 个组织的逐样品计数汇总。",
  "- `tables/06_source_expression_read_audit.csv`：每个输入表达文件的读取、稳定 ID 和文字图 GO 目录哈希审计。",
  "- `tables/07_run_audit.csv`：本次运行规模、符号编码与图文件审计。",
  "- `../archive_5x6_point_20260829/`：修改前的 5 × 6 点图、脚本和审计表归档。"
)
writeLines(readme_lines, file.path(out_dir, "README.md"), useBytes = TRUE)
cat("Wrote per-source seven-pathway counts and 4x7 symbol figure to", out_dir, "\n")
