{
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(stringr)
}

# Input and output directories -------------------------------------------
{
  source_file <- file.path(
    "source_files",
    "AFF3_integrated_final_v3.rds"
  )
  fig6c_dir <- file.path("generated_figures", "Fig6C")
  fig6d_dir <- file.path("generated_figures", "Fig6D")

  dir.create(fig6c_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig6d_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the integrated brain-organoid single-cell object ------------------
{
  AFF3_integrated_final <- readRDS(source_file)
  DefaultAssay(AFF3_integrated_final) <- "RNA"
}

# Cell-type annotations used in the v46 UMAP -----------------------------
{
  cluster_to_celltype <- c(
    "0" = "Neural_crest_or_non_ectodermal",
    "1" = "Meningeal_ChoroidPlexus_like",
    "2" = "Immature_deep_layer_excitatory_neuron",
    "3" = "Maturing_excitatory_neuron",
    "4" = "Cycling_progenitor_G2M",
    "5" = "Immature_deep_layer_excitatory_neuron",
    "6" = "Cycling_progenitor_G2M",
    "7" = "Epithelial_like_non_neural_ectoderm",
    "8" = "Nephron_like_kidney_lineage",
    "9" = "Mesenchymal_or_Vascular_SMC",
    "10" = "Forebrain_progenitor_early_neuron",
    "11" = "Endothelial",
    "12" = "Radial_glia_apical_progenitor",
    "13" = "Epithelial_like_non_neural_ectoderm",
    "14" = "Neural_crest_melanocyte_lineage",
    "15" = "York_sac_endoderm",
    "16" = "Maturing_excitatory_neuron",
    "17" = "Cycling_progenitor_G2M",
    "18" = "Meso_endoderm",
    "19" = "Radial_glia_apical_progenitor",
    "20" = "Meningeal_mixed_signature",
    "21" = "Myeloid_hematopoietic",
    "22" = "Peripheral_sensory_neuron",
    "23" = "Nephron_like_kidney_lineage"
  )

  integrated_labelv3 <- AFF3_integrated_final

  cell_types <- data.frame(
    row.names = rownames(integrated_labelv3@meta.data),
    cell_type = cluster_to_celltype[
      as.character(integrated_labelv3@meta.data$seurat_clusters)
    ]
  )

  integrated_labelv3 <- AddMetaData(
    object = integrated_labelv3,
    metadata = cell_types$cell_type,
    col.name = "cell_type"
  )
}

# Figure 6C: annotated brain-organoid UMAP -------------------------------
{
  cluster_colors <- c(
    "Neural_crest_or_non_ectodermal" = "#C5926D",
    "Meningeal_ChoroidPlexus_like" = "#F2A464",
    "Immature_deep_layer_excitatory_neuron" = "#85BC70",
    "Maturing_excitatory_neuron" = "#D1A740",
    "Cycling_progenitor_G2M" = "#BE6329",
    "Epithelial_like_non_neural_ectoderm" = "#E7E193",
    "Nephron_like_kidney_lineage" = "#6D2A5C",
    "Mesenchymal_or_Vascular_SMC" = "#7095B4",
    "Forebrain_progenitor_early_neuron" = "#C5352E",
    "Endothelial" = "#D17EBF",
    "Radial_glia_apical_progenitor" = "#E45750",
    "Neural_crest_melanocyte_lineage" = "#4F181C",
    "York_sac_endoderm" = "#3F517D",
    "Meso_endoderm" = "#123366",
    "Meningeal_mixed_signature" = "#E45750",
    "Myeloid_hematopoietic" = "#C9DCEA",
    "Peripheral_sensory_neuron" = "#7D1B28"
  )

  point_size <- 0.5

  fig6c_plot <- DimPlot(
    integrated_labelv3,
    group.by = "cell_type",
    reduction = "umap",
    pt.size = point_size,
    label = FALSE
  ) +
    scale_color_manual(
      values = cluster_colors,
      labels = function(x) {
        str_wrap(str_replace_all(x, "_", " "), width = 25)
      }
    ) +
    theme_classic() +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      plot.title = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(x = NULL, y = NULL, title = NULL)

  ggsave(
    file.path(fig6c_dir, "Fig6C_brain_organoid_UMAP.tiff"),
    plot = fig6c_plot,
    width = 9,
    height = 8,
    units = "in",
    dpi = 600,
    limitsize = FALSE,
    compression = "lzw"
  )
}

# Figure 6D: SFRP2 expression in WT, KO, and DKO ------------------------
{
  gene <- "SFRP2"

  plot_feature <- function(obj, genotype_label) {
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
        limits = c(0, 5.5)
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

  p_WT <- plot_feature(integrated_labelv3, "WT")
  p_KO <- plot_feature(integrated_labelv3, "KO")
  p_DKO <- plot_feature(integrated_labelv3, "DKO")

  clean_theme <- theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    panel.background = element_rect(fill = "#F7F7F7", color = NA),
    plot.background = element_rect(fill = "#F7F7F7", color = NA),
    legend.background = element_rect(fill = "#F7F7F7", color = NA)
  )

  p_WT <- p_WT + clean_theme
  p_KO <- p_KO + clean_theme
  p_DKO <- p_DKO + clean_theme

  ggsave(
    file.path(fig6d_dir, "Fig6D_SFRP2_WT.tiff"),
    plot = p_WT,
    width = 3,
    height = 2.5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  ggsave(
    file.path(fig6d_dir, "Fig6D_SFRP2_KO.tiff"),
    plot = p_KO,
    width = 3,
    height = 2.5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  ggsave(
    file.path(fig6d_dir, "Fig6D_SFRP2_DKO.tiff"),
    plot = p_DKO,
    width = 3,
    height = 2.5,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}
