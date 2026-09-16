{
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "Fig5D")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the WT and dominant-negative count tables -------------------------
{
  WT_1_counts <- read.table(
    file.path(source_dir, "Fig5D_WT_1_mm_counts_revised.txt"),
    header = TRUE
  )
  WT_2_counts <- read.table(
    file.path(source_dir, "Fig5D_WT_2_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_1_1_counts <- read.table(
    file.path(source_dir, "Fig5D_DN_1_1_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_1_2_counts <- read.table(
    file.path(source_dir, "Fig5D_DN_1_2_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_2_1_counts <- read.table(
    file.path(source_dir, "Fig5D_DN_2_1_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_2_2_counts <- read.table(
    file.path(source_dir, "Fig5D_DN_2_2_mm_counts_revised.txt"),
    header = TRUE
  )

  count_tables <- list(
    WT_1_counts,
    WT_2_counts,
    DN_1_1_counts,
    DN_1_2_counts,
    DN_2_1_counts,
    DN_2_2_counts
  )
  stopifnot(all(vapply(
    count_tables[-1],
    function(x) identical(x$Geneid, count_tables[[1]]$Geneid),
    logical(1)
  )))
}

# DESeq2 analysis used for the WT-versus-DN comparison ------------------
{
  counts_data <- data.frame(
    AFF3_WT_1 = WT_1_counts$counts,
    AFF3_WT_2 = WT_2_counts$counts,
    AFF3_DN_1_1 = DN_1_1_counts$counts,
    AFF3_DN_1_2 = DN_1_2_counts$counts,
    AFF3_DN_2_1 = DN_2_1_counts$counts,
    AFF3_DN_2_2 = DN_2_2_counts$counts
  )
  rownames(counts_data) <- WT_1_counts$Geneid

  col_data <- data.frame(
    condition = c("WT", "WT", "DN", "DN", "DN", "DN")
  )
  rownames(col_data) <- colnames(counts_data)

  dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", "DN", "WT"))
  res <- as.data.frame(res)
  res$gene <- rownames(res)

  norm_counts <- counts(dds, normalized = TRUE)
  res$WT_mean <- rowMeans(
    norm_counts[, c("AFF3_WT_1", "AFF3_WT_2")]
  ) + 1
  res$DN_mean <- rowMeans(
    norm_counts[, c(
      "AFF3_DN_1_1",
      "AFF3_DN_1_2",
      "AFF3_DN_2_1",
      "AFF3_DN_2_2"
    )]
  ) + 1
  res$WT_log2 <- log2(res$WT_mean)
  res$DN_log2 <- log2(res$DN_mean)

  res$regulation <- "Not DEG"
  res$regulation[
    !is.na(res$padj) &
      res$padj < 0.05 &
      res$log2FoldChange > 0.5
  ] <- "Up-regulated"
  res$regulation[
    !is.na(res$padj) &
      res$padj < 0.05 &
      res$log2FoldChange < -0.5
  ] <- "Down-regulated"

  xist_point <- res %>%
    filter(gene == "XIST")
  stopifnot(nrow(xist_point) == 1)
}

# Figure 5C --------------------------------------------------------------
{
  fig6c_plot <- ggplot(res, aes(x = DN_log2, y = WT_log2)) +
    geom_point(
      data = subset(res, regulation == "Not DEG"),
      color = "grey75",
      alpha = 0.45,
      size = 0.35
    ) +
    geom_point(
      data = subset(res, regulation == "Up-regulated"),
      color = "#85312B",
      alpha = 0.75,
      size = 0.55
    ) +
    geom_point(
      data = subset(res, regulation == "Down-regulated"),
      color = "#373671",
      alpha = 0.75,
      size = 0.55
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.35
    ) +
    geom_point(
      data = xist_point,
      shape = 21,
      fill = "#F2C14E",
      color = "black",
      stroke = 0.45,
      size = 2.6
    ) +
    labs(
      x = "DN expression (log2)",
      y = "WT expression (log2)"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7, color = "black"),
      plot.margin = margin(4, 8, 4, 4)
    ) +
    coord_fixed(ratio = 1, clip = "off")

  ggsave(
    file.path(figure_dir, "Fig5D_AFF3_WT_vs_DN_scatterplot.tiff"),
    plot = fig6c_plot,
    device = "tiff",
    dpi = 600,
    width = 3.3,
    height = 3.3,
    units = "in",
    compression = "lzw"
  )
}
