{
  library(DESeq2)
  library(pheatmap)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- "generated_figures/Fig4H"

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Thresholds --------------------------------------------------------------
{
  padj_cutoff <- 0.05
  log2FC_cutoff <- 0.5
}

# Read count tables -------------------------------------------------------
{
  old_WT1_counts <- read.table(
    file.path(source_dir, "AFF3_F_WT_1_mm_counts_revised.txt"),
    header = TRUE
  )
  old_WT2_counts <- read.table(
    file.path(source_dir, "AFF3_F_WT_2_mm_counts_revised.txt"),
    header = TRUE
  )
  old_KO1_counts <- read.table(
    file.path(source_dir, "AFF3_F_KO_1_mm_counts_revised.txt"),
    header = TRUE
  )
  old_KO2_counts <- read.table(
    file.path(source_dir, "AFF3_F_KO_2_mm_counts_revised.txt"),
    header = TRUE
  )
  old_DKO1_counts <- read.table(
    file.path(source_dir, "AFF3_F_KO_Xt_KO_A3_1_mm_counts_revised.txt"),
    header = TRUE
  )
  old_DKO2_counts <- read.table(
    file.path(source_dir, "AFF3_F_KO_Xt_KO_A3_2_mm_counts_revised.txt"),
    header = TRUE
  )

  new_WT1_counts <- read.table(
    file.path(source_dir, "WT_1_mm_counts_revised.txt"),
    header = TRUE
  )
  new_WT2_counts <- read.table(
    file.path(source_dir, "WT_2_mm_counts_revised.txt"),
    header = TRUE
  )
  new_KO1_counts <- read.table(
    file.path(source_dir, "KO_1_mm_counts_revised.txt"),
    header = TRUE
  )
  new_KO2_counts <- read.table(
    file.path(source_dir, "KO_2_mm_counts_revised.txt"),
    header = TRUE
  )
  new_KOX1_1_counts <- read.table(
    file.path(source_dir, "KOX1_1_mm_counts_revised.txt"),
    header = TRUE
  )
  new_KOX1_2_counts <- read.table(
    file.path(source_dir, "KOX1_2_mm_counts_revised.txt"),
    header = TRUE
  )

  XIST_binding_genes <- read.csv(
    file.path(source_dir, "XIST_binding_933_genes.csv"),
    check.names = FALSE
  )$Gene
  XIST_binding_genes <- unique(XIST_binding_genes)
}

# Helper functions for rescued-gene identification -----------------------
{
  run_deseq <- function(countData, colData) {
    dds <- DESeqDataSetFromMatrix(
      countData = countData,
      colData = colData,
      design = ~ condition
    )
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    dds <- DESeq(dds)
    dds
  }

  get_up_genes <- function(res) {
    rownames(res)[
      !is.na(res$log2FoldChange) &
        !is.na(res$padj) &
        res$padj < padj_cutoff &
        res$log2FoldChange > log2FC_cutoff
    ]
  }

  get_down_genes <- function(res) {
    rownames(res)[
      !is.na(res$log2FoldChange) &
        !is.na(res$padj) &
        res$padj < padj_cutoff &
        res$log2FoldChange < -log2FC_cutoff
    ]
  }
}

# Old experiment: identify genes rescued by AFF3/XIST DKO ----------------
{
  old_rescue_countData <- data.frame(
    old_WT_1 = old_WT1_counts$counts,
    old_WT_2 = old_WT2_counts$counts,
    old_KO_1 = old_KO1_counts$counts,
    old_KO_2 = old_KO2_counts$counts,
    old_DKO_1 = old_DKO1_counts$counts,
    old_DKO_2 = old_DKO2_counts$counts
  )
  rownames(old_rescue_countData) <- old_WT1_counts$Geneid

  old_rescue_colData <- data.frame(
    condition = factor(
      c("WT", "WT", "KO", "KO", "DKO", "DKO"),
      levels = c("WT", "KO", "DKO")
    )
  )
  rownames(old_rescue_colData) <- colnames(old_rescue_countData)

  old_rescue_dds <- run_deseq(old_rescue_countData, old_rescue_colData)
  old_KO_vs_WT <- results(
    old_rescue_dds,
    contrast = c("condition", "KO", "WT")
  )
  old_DKO_vs_KO <- results(
    old_rescue_dds,
    contrast = c("condition", "DKO", "KO")
  )

  old_KO_down <- get_down_genes(old_KO_vs_WT)
  old_DKO_up <- get_up_genes(old_DKO_vs_KO)
  old_rescued_genes <- intersect(old_KO_down, old_DKO_up)
}

# New experiment: identify genes rescued by X1 treatment -----------------
{
  new_rescue_countData <- data.frame(
    new_WT_1 = new_WT1_counts$counts,
    new_WT_2 = new_WT2_counts$counts,
    new_KO_1 = new_KO1_counts$counts,
    new_KO_2 = new_KO2_counts$counts,
    new_KOX1_1 = new_KOX1_1_counts$counts,
    new_KOX1_2 = new_KOX1_2_counts$counts
  )
  rownames(new_rescue_countData) <- new_WT1_counts$Geneid

  new_rescue_colData <- data.frame(
    condition = factor(
      c("WT", "WT", "KO", "KO", "KOX1", "KOX1"),
      levels = c("WT", "KO", "KOX1")
    )
  )
  rownames(new_rescue_colData) <- colnames(new_rescue_countData)

  new_rescue_dds <- run_deseq(new_rescue_countData, new_rescue_colData)
  new_KO_vs_WT <- results(
    new_rescue_dds,
    contrast = c("condition", "KO", "WT")
  )
  new_KOX1_vs_KO <- results(
    new_rescue_dds,
    contrast = c("condition", "KOX1", "KO")
  )

  new_KO_down <- get_down_genes(new_KO_vs_WT)
  new_KOX1_up <- get_up_genes(new_KOX1_vs_KO)
  new_rescued_genes <- intersect(new_KO_down, new_KOX1_up)
}

# Shared rescued and XIST-binding genes ----------------------------------
{
  shared_rescued_genes <- intersect(old_rescued_genes, new_rescued_genes)
  shared_rescued_XIST_binding_genes <- intersect(
    shared_rescued_genes,
    XIST_binding_genes
  )
}

# pheatmap-style plotting -------------------------------------------------
{
  plot_three_group_heatmap <- function(counts_Data, col_Data, sample_labels,
                                       figure_file) {
    dds <- DESeqDataSetFromMatrix(
      countData = counts_Data,
      colData = col_Data,
      design = ~ condition
    )
    dds <- DESeq(dds)
    dds_vst <- vst(dds, blind = FALSE)
    vst_counts <- assay(dds_vst)

    missing_genes <- setdiff(
      shared_rescued_XIST_binding_genes,
      rownames(vst_counts)
    )
    if (length(missing_genes) > 0) {
      warning(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
    }

    gene_list_present <- shared_rescued_XIST_binding_genes[
      shared_rescued_XIST_binding_genes %in% rownames(vst_counts)
    ]
    vst_counts_filtered <- vst_counts[gene_list_present, ]
    z_matrix_filtered <- t(scale(t(vst_counts_filtered)))
    z_matrix_filtered <- z_matrix_filtered[complete.cases(z_matrix_filtered), ]
    colnames(z_matrix_filtered) <- sample_labels

    cols_new <- rev(RColorBrewer::brewer.pal(11, "PuOr"))

    tiff(
      figure_file,
      width = 8,
      height = 7,
      units = "in",
      res = 600,
      compression = "lzw",
      type = "cairo"
    )

    pheatmap(
      z_matrix_filtered,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      color = cols_new,
      main = "",
      fontsize = 16,
      fontsize_row = 18,
      fontsize_col = 18,
      border_color = "grey60",
      scale = "none",
      cellwidth = 30,
      cellheight = 30,
      legend = TRUE,
      angle_col = "45"
    )

    dev.off()
  }
}

# Old WT / KO / DKO heatmap ----------------------------------------------
{
  old_heatmap_countData <- data.frame(
    AFF3_F_WT_1 = old_WT1_counts$counts,
    AFF3_F_WT_2 = old_WT2_counts$counts,
    AFF3_F_KO_1 = old_KO1_counts$counts,
    AFF3_F_KO_2 = old_KO2_counts$counts,
    AFF3_F_DKO_1 = old_DKO1_counts$counts,
    AFF3_F_DKO_2 = old_DKO2_counts$counts
  )
  rownames(old_heatmap_countData) <- old_WT1_counts$Geneid

  old_heatmap_colData <- data.frame(
    condition = factor(
      c("WT", "WT", "KO", "KO", "DKO", "DKO"),
      levels = c("WT", "KO", "DKO")
    )
  )
  rownames(old_heatmap_colData) <- colnames(old_heatmap_countData)

  plot_three_group_heatmap(
    counts_Data = old_heatmap_countData,
    col_Data = old_heatmap_colData,
    sample_labels = c("WT 1", "WT 2", "KO 1", "KO 2", "DKO 1", "DKO 2"),
    figure_file = file.path(
      figure_dir,
      "Fig4H_old_WT_KO_DKO_pheatmap.tiff"
    )
  )
}

# New WT / KO / KO+X1 heatmap --------------------------------------------
{
  new_heatmap_countData <- data.frame(
    AFF3_F_WT_1 = new_WT1_counts$counts,
    AFF3_F_WT_2 = new_WT2_counts$counts,
    AFF3_F_KO_1 = new_KO1_counts$counts,
    AFF3_F_KO_2 = new_KO2_counts$counts,
    AFF3_F_DKO_1 = new_KOX1_1_counts$counts,
    AFF3_F_DKO_2 = new_KOX1_2_counts$counts
  )
  rownames(new_heatmap_countData) <- new_WT1_counts$Geneid

  new_heatmap_colData <- data.frame(
    condition = factor(
      c("WT", "WT", "KO", "KO", "DKO", "DKO"),
      levels = c("WT", "KO", "DKO")
    )
  )
  rownames(new_heatmap_colData) <- colnames(new_heatmap_countData)

  plot_three_group_heatmap(
    counts_Data = new_heatmap_countData,
    col_Data = new_heatmap_colData,
    sample_labels = c(
      "WT 1", "WT 2", "KO 1", "KO 2", "KO+X1 1", "KO+X1 2"
    ),
    figure_file = file.path(
      figure_dir,
      "Fig4H_new_WT_KO_KOX1_pheatmap.tiff"
    )
  )
}
