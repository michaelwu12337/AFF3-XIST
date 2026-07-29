{
  library(Seurat)
  library(ggplot2)
  library(patchwork)
}

# Input and output directories -------------------------------------------
{
  source_file <- file.path(
    "..",
    "..",
    "Figure-5",
    "source_files",
    "AFF3_integrated_final_v3.rds"
  )
  figure_dir <- "generated_figures"

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the integrated brain-organoid single-cell object ------------------
{
  integrated_labelv3 <- readRDS(source_file)
  DefaultAssay(integrated_labelv3) <- "RNA"
}

# Supplementary Figure 6: pluripotency and XIST-regulator expression ----
{
  gene_limits <- c(
    "TSIX" = 2,
    "JPX" = 4,
    "ZFP42" = 2.5,
    "YY1" = 3
  )

  plot_feature <- function(obj, genotype_label, gene, upper_limit) {
    FeaturePlot(
      subset(obj, subset = sample == genotype_label),
      features = gene,
      reduction = "umap",
      order = TRUE,
      pt.size = 0.25,
      keep.scale = "all"
    ) +
      scale_color_gradient(
        low = "lightgrey",
        high = "#4B2EFF",
        limits = c(0, upper_limit)
      ) +
      theme_classic() +
      theme(
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 11),
        axis.line = element_line(linewidth = 1.2),
        axis.ticks = element_line(linewidth = 1.2)
      ) +
      labs(
        x = "UMAP 1",
        y = "UMAP 2",
        title = genotype_label
      )
  }

  for (gene in names(gene_limits)) {
    upper_limit <- unname(gene_limits[gene])

    p_WT <- plot_feature(integrated_labelv3, "WT", gene, upper_limit)
    p_KO <- plot_feature(integrated_labelv3, "KO", gene, upper_limit)
    p_DKO <- plot_feature(integrated_labelv3, "DKO", gene, upper_limit)

    p_combined <- (p_WT | p_KO | p_DKO) +
      plot_annotation(title = paste(gene, "expression"))

    ggsave(
      file.path(figure_dir, paste0("supp_", gene, ".tiff")),
      plot = p_combined,
      width = 9,
      height = 3,
      units = "in",
      dpi = 600,
      limitsize = FALSE,
      compression = "lzw"
    )
  }
}
