{
  setwd("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Figure 5 new/UMAP")
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(stringr)
  library(patchwork)
}

# Load your data
AFF3_integrated_final <- readRDS("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Figure 5 new/UMAP/AFF3_integrated_final_v3.rds")

# Set default assay to RNA
DefaultAssay(AFF3_integrated_final) <- "RNA"

# Define cluster-specific markers based on previous predictions (without pluripotency markers)
cluster_markers <- list(
  "0" = c("SLITRK6", "ALX3", "TWIST1", "SNAI1"),
  "1" = c("TFPI2", "STC1", "NPY", "PTN"),
  "2" = c("SLC24A2", "SOX5", "EBF1", "CNTN5"),
  "3" = c("GRIA1", "NTRK2", "SST", "CHRM2"),
  "4" = c("H2AC11", "H3C2", "TOP2A", "MKI67"),
  "5" = c("DISC1FP1", "DDAH1", "CNTN5", "CHSY3"),
  "6" = c("TOP2A", "MKI67", "AURKA", "ASPM"),
  "7" = c("CLDN4", "CLDN7", "TP63", "GABRP"),
  "8" = c("SLC34A2", "PODXL", "LRP2", "NPHS1"),
  "9" = c("GREM2", "WNT5A", "TEK", "TGFBR3"),
  "10" = c("PAX6", "SOX2", "EMX2", "OTX2"),
  "11" = c("PECAM1", "CDH5", "CLDN5", "FLT4"),
  "12" = c("SOX2", "PAX6", "HES5"),
  "13" = c("ITGB6", "CLDN6", "KRT7", "GABRP"),
  "14" = c("SOX10", "PMEL", "ERBB3", "FOXD3"),
  "15" = c("AFP", "ALB", "APOA1", "HNF4A"),
  "16" = c("TBR1", "MYT1L", "EOMES", "STMN4"),
  "17" = c("CENPE", "ASPM", "BUB1", "KIF4A"),
  "18" = c("MIXL1", "CDX2", "HAND1", "APLNR"),
  "19" = c("SOX2", "HES5", "PAX6", "OTX2"),
  "20" = c("MGP", "SST", "DES", "PTN"),
  "21" = c("LYZ", "MPO", "AZU1", "GATA1"),
  "22" = c("POU4F1", "TLX3", "PRPH", "SCG2"),
  "23" = c("SERPINE3", "TMEM151A", "SLC34A2", "PODXL")
)

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

# Create a data frame with cell barcodes and cell types
integrated_labelv3 <- AFF3_integrated_final

# Map cluster numbers to cell types
cell_types <- data.frame(
  row.names = rownames(integrated_labelv3@meta.data),
  cell_type = cluster_to_celltype[as.character(integrated_labelv3@meta.data$seurat_clusters)]
)

# Add cell type annotations to the Seurat object
integrated_labelv3 <- AddMetaData(
  object = integrated_labelv3,
  metadata = cell_types$cell_type,
  col.name = "cell_type"
)

# Define cluster names and their assigned colors
cluster_colors <- c(
  "Neural_crest_or_non_ectodermal" = "#C5926D",          #NC
  "Meningeal_ChoroidPlexus_like" = "#F2A464",            #MCP
  "Immature_deep_layer_excitatory_neuron" = "#85BC70",   #IEN
  "Maturing_excitatory_neuron" = "#D1A740",              #EN
  "Cycling_progenitor_G2M" = "#BE6329",                  #CycP
  "Epithelial_like_non_neural_ectoderm" = "#E7E193",     #Epi
  "Nephron_like_kidney_lineage" = "#6D2A5C",             #Kid
  "Mesenchymal_or_Vascular_SMC" = "#7095B4",             #SMC
  "Forebrain_progenitor_early_neuron" = "#C5352E",       #FP
  "Endothelial" = "#D17EBF",                             #EC
  "Radial_glia_apical_progenitor" = "#E45750",           #RGP
  "Neural_crest_melanocyte_lineage" = "#4F181C",         #NCM
  "York_sac_endoderm" = "#3F517D",                       #YSE
  "Meso_endoderm" = "#123366",                           #ME
  "Meningeal_mixed_signature" = "#E45750",               #MM
  "Myeloid_hematopoietic" = "#C9DCEA",                   #BP
  "Peripheral_sensory_neuron" = "#7D1B28"                #PSN
)

# ---- manual knobs ----
{
  legend_title_size <- 16
  legend_text_size  <- 16
  legend_key_lines  <- 1.2
  legend_dot_size   <- 5
  
  axis_title_size   <- 24
  axis_text_size    <- 16
  
  border_width      <- 1.5
  axis_line_width   <- 1.2
  tick_line_width   <- 1.2
  tick_length_cm    <- 0.25
  
  point_size        <- 0.5
}

DimPlot(
  integrated_labelv3,
  group.by = "cell_type",
  reduction = "umap",
  pt.size = point_size,
  label = FALSE   # no superimposed labels
  ) +
  scale_color_manual(
    values = cluster_colors,
    labels = function(x) str_wrap(str_replace_all(x, "_", " "), width = 25)
  ) +
  theme_classic() +
  theme(
    # tight bounding box
    panel.border = element_rect(color = "black", 
                                fill = NA, 
                                linewidth = border_width),
    
    # axis labels + ticks
    axis.title = element_text(size = axis_title_size, face = "bold"),
    axis.text = element_text(size = axis_text_size),
    axis.line = element_line(linewidth = axis_line_width),
    axis.ticks = element_line(linewidth = tick_line_width),
    axis.ticks.length = unit(tick_length_cm, "cm"),
    
    # legend
    legend.text = element_text(size = legend_text_size),
    legend.title = element_text(size = legend_title_size, face = "bold"),
    legend.key.size = unit(legend_key_lines, "lines"),
    legend.spacing.y = unit(1, "lines"),
    
    # 2) remove title
    plot.title = element_blank(),
    
    # margins
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(x = "UMAP 1", 
       y = "UMAP 2", 
       title = NULL) +   # also ensures title is gone
  guides(color = guide_legend(
    override.aes = list(size = legend_dot_size)
    )
  )

DimPlot(
  integrated_labelv3,
  group.by = "cell_type",
  reduction = "umap",
  pt.size = point_size,
  label = FALSE
) +
  scale_color_manual(
    values = cluster_colors,
    labels = function(x) str_wrap(str_replace_all(x, "_", " "), width = 25)
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = axis_line_width),
    
    axis.title = element_blank(),      # remove axis titles
    axis.text = element_text(size = axis_text_size),
    axis.line = element_line(linewidth = axis_line_width),
    axis.ticks = element_line(linewidth = tick_line_width),
    axis.ticks.length = unit(tick_length_cm, "cm"),
    
    legend.position = "none",          # remove legend
    
    plot.title = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(x = NULL, y = NULL, title = NULL)

DimPlot(
  integrated_labelv3,
  group.by = "cell_type",
  reduction = "umap",
  pt.size = point_size,
  label = FALSE
) +
  scale_color_manual(
    values = cluster_colors,
    labels = function(x) str_wrap(str_replace_all(x, "_", " "), width = 25)
  ) +
  theme_classic() +
  theme(
    panel.border = element_blank(),   # remove outer border
    axis.line = element_blank(),      # remove axis lines
    axis.ticks = element_blank(),     # remove ticks
    axis.text = element_blank(),      # remove tick labels
    axis.title = element_blank(),     # remove axis titles
    
    legend.position = "none",
    
    plot.title = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(x = NULL, y = NULL, title = NULL)

ggsave(
  "fig5_UMAP_no_border_08.tiff",
  plot = last_plot(),
  width = 9,    
  height = 8,      
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  compression = "lzw"
)

### ----------------------------------------------
### separate UMAP WT KO and DKO
### ----------------------------------------------
{
  {
    {
      legend_title_size <- 16
      legend_text_size  <- 24
      legend_key_lines  <- 1.2
      legend_dot_size   <- 5
      
      axis_title_size   <- 24
      axis_text_size    <- 24
      
      border_width      <- 1.5
      axis_line_width   <- 0.5
      tick_line_width   <- 0.5
      tick_length_cm    <- 0.25
      
      point_size        <- 0.2
    }
    
    integrated_labelv3 <- AFF3_integrated_final
    
    unique(integrated_labelv3$sample)
    
    obj_WT  <- subset(integrated_labelv3, subset = sample == "WT")
    obj_KO  <- subset(integrated_labelv3, subset = sample == "KO")
    obj_DKO <- subset(integrated_labelv3, subset = sample == "DKO")
  }
  
  plot_umap_by_group <- function(obj, plot_title = NULL) {
    DimPlot(
      obj,
      group.by = "cell_type",
      reduction = "umap",
      pt.size = point_size,
      label = FALSE
    ) +
      scale_color_manual(
        values = cluster_colors,
        labels = function(x) str_wrap(str_replace_all(x, "_", " "), width = 25)
      ) +
      theme_classic() +
      theme(
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(linewidth = axis_line_width),
        axis.ticks = element_line(linewidth = tick_line_width),
        axis.ticks.length = unit(tick_length_cm, "cm"),
        legend.position = "none",
        legend.text = element_text(size = legend_text_size),
        legend.title = element_text(size = legend_title_size, face = "bold"),
        legend.key.size = unit(legend_key_lines, "lines"),
        legend.spacing.y = unit(1, "lines"),
        plot.title = element_blank(),
        plot.margin = margin(15, 15, 15, 15)
      ) +
      labs(
        x = "UMAP 1",
        y = "UMAP 2",
        title = plot_title
      ) +
      guides(color = guide_legend(
        override.aes = list(size = legend_dot_size)
      ))
  }
  
  # Generate plots
  {
    p_WT  <- plot_umap_by_group(obj_WT,  "WT")
    p_KO  <- plot_umap_by_group(obj_KO,  "KO")
    p_DKO <- plot_umap_by_group(obj_DKO, "DKO")
  }
  
  p_WT | p_KO | p_DKO
}


# CHecking gene expression pattern using feature plot side by side
{
  gene <- "SFRP2"   # change gene here only
  
  plot_feature <- function(obj, genotype_label) {
    FeaturePlot(
      subset(obj, subset = sample == genotype_label),
      features  = gene,
      reduction = "umap",
      order     = TRUE,
      pt.size   = 0.25,
      keep.scale = "all"
    ) +
      scale_color_gradient(
        low  = "lightgrey",
        high = "#4B2EFF",
        limits = c(0, 5.5)
      ) +
      theme_classic() +
      theme(
        axis.title = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 11),
        axis.line  = element_line(linewidth = 1.2),
        axis.ticks = element_line(linewidth = 1.2)
      ) +
      labs(
        x = "UMAP 1",
        y = "UMAP 2",
        title = genotype_label
      )
  }
  
  p_WT  <- plot_feature(integrated_labelv3, "WT")
  p_KO  <- plot_feature(integrated_labelv3, "KO")
  p_DKO <- plot_feature(integrated_labelv3, "DKO")
  
  p_combined <- (p_WT | p_KO | p_DKO) +
    plot_annotation(title = paste(gene, "expression"))
  
  p_combined
}
{
  # Remove borders, axes, ticks for all three
  clean_theme <- theme(
    panel.border = element_blank(),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    axis.text    = element_blank(),
    axis.title   = element_blank(),
    legend.position = "none",
    plot.title      = element_blank(),
    
    panel.background = element_rect(fill = "#F7F7F7", color = NA),
    plot.background  = element_rect(fill = "#F7F7F7", color = NA),
    legend.background  = element_rect(fill = "#F7F7F7", color = NA)
  )
  
  p_WT  <- p_WT  + clean_theme
  p_KO  <- p_KO  + clean_theme
  p_DKO <- p_DKO + clean_theme
  
  p_combined <- (p_WT | p_KO | p_DKO)
}

ggsave(
  "SFRP2_WT_lab_grey.tiff",
  plot = p_WT,
  width = 3,
  height = 2.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)





