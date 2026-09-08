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

counts <- fread(file.path(root, "results/tables/02_seven_pathway_member_gene_counts_per_source_sample_long.csv"))

pathway_levels <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")

pathway_shapes <- c(BER = 24L, NER = 25L, MMR = 23L, FA = 0L, HR = 8L, AEJ = 3L, NHEJ = 4L)

pathway_legend_labels <- c(BER = "BER  Base excision", NER = "NER  Nucleotide excision", MMR = "MMR  Mismatch repair",
    FA = "FA  Fanconi anemia", HR = "HR  Homologous recombination", AEJ = "A-EJ  Alternative end joining",
    NHEJ = "NHEJ  Non-homologous end joining")

tissue_levels <- c("Adrenal gland", "Bone marrow", "Cecum", "Cerebellum", "Colon", "Dorsal skin", "Extensor digitorum longus",
    "Eye", "Femoral lateral skin", "Femoral skin", "Gastrocnemius", "Heart", "Heart / Heart right ventricle",
    "Kidney", "Left lobe of the liver", "Liver", "Lung", "Mammary gland", "Optic nerve", "Quadriceps femoris",
    "Retina", "Soleus", "Spleen", "Spleen-distal", "Thymus", "Tibialis anterior")

tissue_display <- setNames(tissue_levels, tissue_levels)

tissue_display[["Heart / Heart right ventricle"]] <- "Heart right ventricle"

stopifnot(uniqueN(counts$Group) == 26L, uniqueN(counts$Pathway) == 7L)

sample_colours <- c(Flight = "#D55E00", Ground = "#0072B2")

plot_data <- counts[, .(Tissue = factor(TissueLabel, levels = unname(tissue_display[tissue_levels])),
    TissueOrder, Pathway = factor(as.character(Pathway), levels = pathway_levels), sample_status, sample_column,
    n_member_genes)]

plot_data[, `:=`(sample_status, factor(sample_status, levels = c("Ground", "Flight")))]

setorder(plot_data, TissueOrder, sample_status, sample_column, Pathway)

max_y <- ceiling((max(plot_data$n_member_genes) + 12)/10) * 10

Error in x[[i]] : subscript out of bounds
Calls: cat ... simplify -> simplify -> simplify -> simplify -> simplify
Execution halted
