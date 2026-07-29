{
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  fig9c_dir <- file.path("generated_figures", "SuppFig9C")
  fig9d_dir <- file.path("generated_figures", "SuppFig9D")

  dir.create(fig9c_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig9d_dir, recursive = TRUE, showWarnings = FALSE)
}

# Thresholds --------------------------------------------------------------
{
  padj_cutoff <- 0.05
  log2FC_cutoff <- 0.5
}

# Read count tables and XIST-binding genes -------------------------------
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

# Old experiment: AFF3/XIST DKO rescue -----------------------------------
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

  old_rescue_dds <- run_deseq(
    old_rescue_countData,
    old_rescue_colData
  )
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

# New experiment: X1 rescue ----------------------------------------------
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

  new_rescue_dds <- run_deseq(
    new_rescue_countData,
    new_rescue_colData
  )
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
  shared_rescued_genes <- intersect(
    old_rescued_genes,
    new_rescued_genes
  )
  shared_rescued_XIST_binding_genes <- intersect(
    shared_rescued_genes,
    XIST_binding_genes
  )
}

# Supplementary Figure 9C: all rescued-gene heatmaps ---------------------
{
  plot_rescued_heatmap <- function(
      counts_Data,
      col_Data,
      gene_list,
      sample_labels,
      figure_file) {
    dds <- DESeqDataSetFromMatrix(
      countData = counts_Data,
      colData = col_Data,
      design = ~ condition
    )
    dds <- DESeq(dds)
    dds_vst <- vst(dds, blind = FALSE)
    vst_counts <- assay(dds_vst)

    missing_genes <- setdiff(gene_list, rownames(vst_counts))
    if (length(missing_genes) > 0) {
      warning(
        paste("Missing genes:", paste(missing_genes, collapse = ", "))
      )
    }

    gene_list_present <- gene_list[
      gene_list %in% rownames(vst_counts)
    ]
    vst_counts_filtered <- vst_counts[gene_list_present, ]
    z_matrix_filtered <- t(scale(t(vst_counts_filtered)))
    z_matrix_filtered <- z_matrix_filtered[
      complete.cases(z_matrix_filtered),
    ]
    colnames(z_matrix_filtered) <- sample_labels

    cols_new <- rev(RColorBrewer::brewer.pal(11, "PuOr"))

    tiff(
      figure_file,
      width = 4.8,
      height = 12,
      units = "in",
      res = 600,
      compression = "lzw",
      type = "cairo"
    )

    pheatmap(
      z_matrix_filtered,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = FALSE,
      show_colnames = TRUE,
      color = cols_new,
      main = "",
      fontsize = 16,
      fontsize_col = 18,
      border_color = NA,
      scale = "none",
      legend = TRUE,
      angle_col = "45"
    )

    dev.off()
  }

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

  plot_rescued_heatmap(
    counts_Data = old_heatmap_countData,
    col_Data = old_heatmap_colData,
    gene_list = old_rescued_genes,
    sample_labels = c(
      "WT 1", "WT 2", "KO 1", "KO 2", "DKO 1", "DKO 2"
    ),
    figure_file = file.path(
      fig9c_dir,
      "SuppFig9C_old_DKO_rescued_genes_heatmap.tiff"
    )
  )

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

  plot_rescued_heatmap(
    counts_Data = new_heatmap_countData,
    col_Data = new_heatmap_colData,
    gene_list = new_rescued_genes,
    sample_labels = c(
      "WT 1", "WT 2", "KO 1", "KO 2", "KO+X1 1", "KO+X1 2"
    ),
    figure_file = file.path(
      fig9c_dir,
      "SuppFig9C_new_KOX1_rescued_genes_heatmap.tiff"
    )
  )
}

# Supplementary Figure 9D: Venn-diagram functions ------------------------
{
  circle_df <- function(x0, y0, r, set_name, n = 600) {
    theta <- seq(0, 2 * pi, length.out = n)
    data.frame(
      x = x0 + r * cos(theta),
      y = y0 + r * sin(theta),
      set = set_name
    )
  }

  count_label <- function(n, label) {
    paste0(format(n, big.mark = ","), "\n", label)
  }

  make_two_set_venn <- function(
      set_a,
      set_b,
      label_a,
      label_b,
      title,
      subtitle,
      figure_file,
      left_fill = "#4C78A8",
      right_fill = "#F58518",
      left_center = -0.62,
      right_center = 0.62,
      left_radius = 1,
      right_radius = 1,
      region_labels = c("Old only", "Overlap", "New only")) {
    a_only <- setdiff(set_a, set_b)
    b_only <- setdiff(set_b, set_a)
    overlap <- intersect(set_a, set_b)

    circle_data <- bind_rows(
      circle_df(left_center, 0, left_radius, "left"),
      circle_df(right_center, 0, right_radius, "right")
    )

    max_radius <- max(left_radius, right_radius)
    x_min <- min(
      left_center - left_radius,
      right_center - right_radius
    ) - 0.35
    x_max <- max(
      left_center + left_radius,
      right_center + right_radius
    ) + 0.35

    region_label_data <- data.frame(
      x = c(
        left_center - left_radius * 0.43,
        mean(c(left_center, right_center)),
        right_center + right_radius * 0.43
      ),
      y = c(-0.04, -0.04, -0.04),
      label = c(
        count_label(length(a_only), region_labels[1]),
        count_label(length(overlap), region_labels[2]),
        count_label(length(b_only), region_labels[3])
      )
    )

    set_label_data <- data.frame(
      x = c(left_center, right_center),
      y = c(max_radius + 0.31, max_radius + 0.31),
      label = c(
        paste0(
          label_a,
          "\n",
          "n = ",
          format(length(set_a), big.mark = ",")
        ),
        paste0(
          label_b,
          "\n",
          "n = ",
          format(length(set_b), big.mark = ",")
        )
      )
    )

    venn_plot <- ggplot() +
      geom_polygon(
        data = circle_data,
        aes(x = x, y = y, group = set, fill = set),
        color = "#252525",
        linewidth = 1.1,
        alpha = 0.48
      ) +
      geom_text(
        data = region_label_data,
        aes(x = x, y = y, label = label),
        size = 5.8,
        fontface = "bold",
        lineheight = 0.9,
        color = "#111111"
      ) +
      geom_text(
        data = set_label_data,
        aes(x = x, y = y, label = label),
        size = 4.1,
        fontface = "bold",
        lineheight = 0.95,
        color = "#222222"
      ) +
      scale_fill_manual(
        values = c(left = left_fill, right = right_fill),
        guide = "none"
      ) +
      coord_equal(
        xlim = c(x_min, x_max),
        ylim = c(
          -max_radius - 0.22,
          max_radius + 0.62
        ),
        clip = "off"
      ) +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 14, 10, 14)
      )

    ggsave(
      figure_file,
      plot = venn_plot,
      device = "tiff",
      dpi = 600,
      width = 7.2,
      height = 4.8,
      units = "in",
      compression = "lzw"
    )

    invisible(NULL)
  }
}

# Supplementary Figure 9D: rescue overlaps -------------------------------
{
  invisible(make_two_set_venn(
    set_a = old_rescued_genes,
    set_b = new_rescued_genes,
    label_a = "Old rescue\nWT-KO-DKO",
    label_b = "New rescue\nWT-KO-KO+X1",
    title = "Old vs new rescued genes",
    subtitle = "Old DKO rescue compared with new KO+X1 rescue",
    figure_file = file.path(
      fig9d_dir,
      "SuppFig9D_old_vs_new_rescue_gene_venn.tiff"
    ),
    left_fill = "#5B8FF9",
    right_fill = "#F6BD16",
    region_labels = c("Old only", "Shared", "New only")
  ))

  invisible(make_two_set_venn(
    set_a = shared_rescued_genes,
    set_b = XIST_binding_genes,
    label_a = "Shared rescued\n68-gene list",
    label_b = "XIST-binding\nfull list",
    title = "Shared rescued genes vs XIST-binding list",
    subtitle = paste(
      "68 old/new shared rescued genes compared with",
      "the full XIST-binding list"
    ),
    figure_file = file.path(
      fig9d_dir,
      "SuppFig9D_shared_rescued_68_vs_XIST_binding_venn.tiff"
    ),
    left_fill = "#6DC8B5",
    right_fill = "#B279A2",
    left_center = -0.62,
    right_center = 0.62,
    left_radius = 1,
    right_radius = 1,
    region_labels = c(
      "Shared\nonly",
      "Overlap",
      "XIST-binding\nonly"
    )
  ))

}
