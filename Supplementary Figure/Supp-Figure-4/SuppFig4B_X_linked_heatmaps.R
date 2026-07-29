{
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "SuppFig4B")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# X-linked genes ----------------------------------------------------------
{
  # Snapshot of the Ensembl release 115 query used for the v14 figure:
  # chromosome_name = "X", attribute = "hgnc_symbol".
  gene_list_x_linked <- read.csv(
    file.path(source_dir, "X_linked_HGNC_genes_Ensembl115.csv"),
    stringsAsFactors = FALSE
  )$Gene
  gene_list_x_linked <- unique(gene_list_x_linked)
}

# Generate one sex-specific heatmap --------------------------------------
make_x_linked_heatmap <- function(
  wt1_file,
  wt2_file,
  ko1_file,
  ko2_file,
  sample_names,
  sex_label,
  output_file
) {
  wt1_counts <- read.table(
    file.path(source_dir, wt1_file),
    header = TRUE
  )
  wt2_counts <- read.table(
    file.path(source_dir, wt2_file),
    header = TRUE
  )
  ko1_counts <- read.table(
    file.path(source_dir, ko1_file),
    header = TRUE
  )
  ko2_counts <- read.table(
    file.path(source_dir, ko2_file),
    header = TRUE
  )

  stopifnot(
    identical(wt1_counts$Geneid, wt2_counts$Geneid),
    identical(wt1_counts$Geneid, ko1_counts$Geneid),
    identical(wt1_counts$Geneid, ko2_counts$Geneid)
  )

  counts_data <- data.frame(
    wt1_counts$counts,
    wt2_counts$counts,
    ko1_counts$counts,
    ko2_counts$counts
  )
  colnames(counts_data) <- sample_names
  rownames(counts_data) <- wt1_counts$Geneid

  col_data <- data.frame(
    condition = c("WT", "WT", "KO", "KO")
  )
  rownames(col_data) <- sample_names

  dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = col_data,
    design = ~ condition
  )
  dds$condition <- relevel(dds$condition, ref = "WT")
  dds <- DESeq(dds)

  dds_vst <- vst(dds, blind = FALSE)
  vst_counts <- assay(dds_vst)

  gene_list_x_linked_clean <- gene_list_x_linked[
    gene_list_x_linked %in% rownames(vst_counts)
  ]
  counts_x_linked <- vst_counts[gene_list_x_linked_clean, ]
  z_matrix_x_linked <- t(scale(t(counts_x_linked)))

  z_matrix_x_linked_clean <- z_matrix_x_linked[
    apply(z_matrix_x_linked, 1, function(x) all(is.finite(x))),
    ,
    drop = FALSE
  ]

  annotation_x_linked <- data.frame(
    Condition = colData(dds)$condition
  )
  rownames(annotation_x_linked) <- colnames(vst_counts)

  annotation_colors <- list(
    Condition = c(
      WT = "black",
      KO = "lightgrey"
    )
  )

  heatmap_plot <- pheatmap(
    z_matrix_x_linked_clean,
    annotation_col = annotation_x_linked,
    annotation_colors = annotation_colors,
    color = colorRampPalette(
      c(
        "#123366",
        "#C9DCEA",
        "#F7F5F2",
        "#E1A080",
        "#7D1B28"
      )
    )(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = paste("Heatmap of X-linked Genes in", sex_label),
    silent = TRUE
  )

  ggsave(
    file.path(figure_dir, output_file),
    plot = heatmap_plot,
    width = 10,
    height = 10,
    units = "in",
    dpi = 600,
    limitsize = FALSE,
    compression = "lzw"
  )

  invisible(heatmap_plot)
}

# Supplementary Figure 4B ------------------------------------------------
female_heatmap <- make_x_linked_heatmap(
  wt1_file = "AFF3_F_WT_1_mm_counts_revised.txt",
  wt2_file = "AFF3_F_WT_2_mm_counts_revised.txt",
  ko1_file = "AFF3_F_KO_1_mm_counts_revised.txt",
  ko2_file = "AFF3_F_KO_2_mm_counts_revised.txt",
  sample_names = c(
    "AFF3_F_WT_1",
    "AFF3_F_WT_2",
    "AFF3_F_KO_1",
    "AFF3_F_KO_2"
  ),
  sex_label = "Female",
  output_file = "heatmap_x_linked_F.tiff"
)

male_heatmap <- make_x_linked_heatmap(
  wt1_file = "AFF3_M_WT_1_mm_counts_revised.txt",
  wt2_file = "AFF3_M_WT_2_mm_counts_revised.txt",
  ko1_file = "AFF3_M_KO_1_1_mm_counts_revised.txt",
  ko2_file = "AFF3_M_KO_1_2_mm_counts_revised.txt",
  sample_names = c(
    "AFF3_M_WT_1",
    "AFF3_M_WT_2",
    "AFF3_M_KO_1",
    "AFF3_M_KO_2"
  ),
  sex_label = "Male",
  output_file = "heatmap_x_linked_M.tiff"
)
