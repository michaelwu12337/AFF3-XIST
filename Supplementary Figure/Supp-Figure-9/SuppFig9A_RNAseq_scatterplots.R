{
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "SuppFig9A")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Count tables ------------------------------------------------------------
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
  new_WTX1_1_counts <- read.table(
    file.path(source_dir, "WT+X1_1_mm_counts_revised.txt"),
    header = TRUE
  )
  new_WTX1_2_counts <- read.table(
    file.path(source_dir, "WT+X1_2_mm_counts_revised.txt"),
    header = TRUE
  )
}

# DESeq2 scatterplot function --------------------------------------------
{
  plot_RNAseq_scatter <- function(
      reference_1,
      reference_2,
      comparison_1,
      comparison_2,
      reference_condition,
      comparison_condition,
      plot_title,
      x_axis_title,
      y_axis_title,
      output_name) {
    stopifnot(
      identical(reference_1$Geneid, reference_2$Geneid),
      identical(reference_1$Geneid, comparison_1$Geneid),
      identical(reference_1$Geneid, comparison_2$Geneid)
    )

    count_data <- data.frame(
      reference_1 = reference_1$counts,
      reference_2 = reference_2$counts,
      comparison_1 = comparison_1$counts,
      comparison_2 = comparison_2$counts
    )
    rownames(count_data) <- reference_1$Geneid

    col_data <- data.frame(
      condition = factor(
        c(
          reference_condition,
          reference_condition,
          comparison_condition,
          comparison_condition
        ),
        levels = c(reference_condition, comparison_condition)
      )
    )
    rownames(col_data) <- colnames(count_data)

    dds <- DESeqDataSetFromMatrix(
      countData = count_data,
      colData = col_data,
      design = ~ condition
    )
    dds <- DESeq(dds)

    res <- results(
      dds,
      contrast = c(
        "condition",
        comparison_condition,
        reference_condition
      )
    )
    res <- as.data.frame(res)
    res$gene <- rownames(res)

    norm_counts <- counts(dds, normalized = TRUE)
    res$reference_mean <- rowMeans(
      norm_counts[, c("reference_1", "reference_2")]
    ) + 1
    res$comparison_mean <- rowMeans(
      norm_counts[, c("comparison_1", "comparison_2")]
    ) + 1
    res$reference_log2 <- log2(res$reference_mean)
    res$comparison_log2 <- log2(res$comparison_mean)

    res$regulation <- "Not DEG"
    res$regulation[
      res$padj < 0.05 & res$log2FoldChange > 0.5
    ] <- "Up-regulated"
    res$regulation[
      res$padj < 0.05 & res$log2FoldChange < -0.5
    ] <- "Down-regulated"

    non_deg_data <- res[res$regulation == "Not DEG", ]
    median_expr <- median(
      non_deg_data$reference_log2 + non_deg_data$comparison_log2,
      na.rm = TRUE
    )

    res$non_deg_group <- "Not DEG - Low"
    res$non_deg_group[
      res$regulation == "Not DEG" &
        (res$reference_log2 + res$comparison_log2) > median_expr
    ] <- "Not DEG - High"
    res$non_deg_group[
      res$regulation != "Not DEG"
    ] <- res$regulation[res$regulation != "Not DEG"]

    genes_to_label <- c()
    res$label <- ifelse(
      res$gene %in% genes_to_label,
      res$gene,
      ""
    )

    scatter_plot <- ggplot(
      res,
      aes(x = comparison_log2, y = reference_log2)
    ) +
      geom_point(
        data = subset(
          res,
          regulation == "Not DEG" &
            non_deg_group == "Not DEG - Low"
        ),
        color = "gray70",
        alpha = 0.6,
        size = 1
      ) +
      geom_point(
        data = subset(
          res,
          regulation == "Not DEG" &
            non_deg_group == "Not DEG - High"
        ),
        color = "gray70",
        alpha = 0.6,
        size = 1
      ) +
      geom_point(
        data = subset(res, regulation == "Up-regulated"),
        color = "#85312B",
        alpha = 0.8,
        size = 1.5
      ) +
      geom_point(
        data = subset(res, regulation == "Down-regulated"),
        color = "#373671",
        alpha = 0.8,
        size = 1.5
      ) +
      geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed",
        color = "blue",
        alpha = 0.5
      ) +
      geom_text_repel(
        aes(label = label),
        size = 10,
        box.padding = 0.5,
        point.padding = 0.2,
        max.overlaps = 20,
        segment.color = "black",
        segment.size = 0.3
      ) +
      theme_minimal() +
      labs(
        title = plot_title,
        x = x_axis_title,
        y = y_axis_title
      ) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 50
        ),
        axis.title = element_text(size = 50),
        axis.text = element_text(size = 50)
      ) +
      coord_fixed(ratio = 1)

    ggsave(
      file.path(figure_dir, paste0(output_name, ".tiff")),
      plot = scatter_plot,
      device = "tiff",
      dpi = 900,
      width = 9,
      height = 9,
      units = "in",
      compression = "lzw"
    )
  }
}

# Supplementary Figure 9A -------------------------------------------------
{
  plot_RNAseq_scatter(
    reference_1 = new_WT1_counts,
    reference_2 = new_WT2_counts,
    comparison_1 = new_WTX1_1_counts,
    comparison_2 = new_WTX1_2_counts,
    reference_condition = "WT",
    comparison_condition = "WT+X1",
    plot_title = "WT vs WT+X1",
    x_axis_title = "WT+X1 (Log2 values)",
    y_axis_title = "WT (Log2 values)",
    output_name = "SuppFig9A_WT_vs_WT_X1"
  )

  plot_RNAseq_scatter(
    reference_1 = new_WT1_counts,
    reference_2 = new_WT2_counts,
    comparison_1 = new_KO1_counts,
    comparison_2 = new_KO2_counts,
    reference_condition = "WT",
    comparison_condition = "KO",
    plot_title = "WT vs KO",
    x_axis_title = "KO (Log2 values)",
    y_axis_title = "WT (Log2 values)",
    output_name = "SuppFig9A_new_WT_vs_KO"
  )

  plot_RNAseq_scatter(
    reference_1 = new_KO1_counts,
    reference_2 = new_KO2_counts,
    comparison_1 = new_KOX1_1_counts,
    comparison_2 = new_KOX1_2_counts,
    reference_condition = "KO",
    comparison_condition = "KO+X1",
    plot_title = "KO vs KO+X1",
    x_axis_title = "KO+X1 (Log2 values)",
    y_axis_title = "KO (Log2 values)",
    output_name = "SuppFig9A_KO_vs_KO_X1"
  )

  plot_RNAseq_scatter(
    reference_1 = old_KO1_counts,
    reference_2 = old_KO2_counts,
    comparison_1 = old_DKO1_counts,
    comparison_2 = old_DKO2_counts,
    reference_condition = "KO",
    comparison_condition = "DKO",
    plot_title = "DKO vs KO",
    x_axis_title = "DKO (Log2 values)",
    y_axis_title = "KO (Log2 values)",
    output_name = "SuppFig9A_DKO_vs_KO"
  )

  plot_RNAseq_scatter(
    reference_1 = old_WT1_counts,
    reference_2 = old_WT2_counts,
    comparison_1 = old_KO1_counts,
    comparison_2 = old_KO2_counts,
    reference_condition = "WT",
    comparison_condition = "KO",
    plot_title = "KO vs WT",
    x_axis_title = "KO (Log2 values)",
    y_axis_title = "WT (Log2 values)",
    output_name = "SuppFig9A_old_KO_vs_WT"
  )
}
