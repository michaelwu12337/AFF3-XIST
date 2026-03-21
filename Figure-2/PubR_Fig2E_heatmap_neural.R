{
  setwd("C:/Users/Michael Wu/OneDrive/HKU PhD/Publication ready/Figure2E_Heatmap")
  library(DESeq2)
  library(tidyverse)
  library(tidyr)
  library(tibble)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggnewscale)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(scales)
  library(grid)
  library(showtext)
  library(extrafont)
}

# Raw data processing -----------------------------------------------------
{
  F_KO1_counts <- read.table(file = "AFF3_F_KO_1_mm_counts_revised.txt", header = T)
  F_KO2_counts <- read.table(file = "AFF3_F_KO_2_mm_counts_revised.txt", header = T)
  F_WT1_counts <- read.table(file = "AFF3_F_WT_1_mm_counts_revised.txt", header = T)
  F_WT2_counts <- read.table(file = "AFF3_F_WT_2_mm_counts_revised.txt", header = T)
}
{
  counts_Data <- data.frame(AFF3_F_WT_1 = F_WT1_counts$counts)
  counts_Data$AFF3_F_WT_2 <- F_WT2_counts$counts
  counts_Data$AFF3_F_KO_1 <- F_KO1_counts$counts
  counts_Data$AFF3_F_KO_2 <- F_KO2_counts$counts
  rownames(counts_Data) <- F_WT1_counts$Geneid
}
{
  col_Data <- data.frame(condition = c("WT","WT","KO","KO"))
  rownames(col_Data) <- c("AFF3_F_WT_1", 
                          "AFF3_F_WT_2", 
                          "AFF3_F_KO_1", 
                          "AFF3_F_KO_2")
  condition <- factor(col_Data$condition)
}
{
  all(colnames(counts_Data) %in% rownames(col_Data))
  all(colnames(counts_Data) == rownames(col_Data))
}


# DESeq2 ------------------------------------------------------------------
{
  dds <- DESeqDataSetFromMatrix(countData = counts_Data,
                                colData = col_Data,
                                design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = "WT")
  dds <- DESeq(dds)
  dds_vst <- vst(dds, blind=FALSE)
  vst_counts <- assay(dds_vst)
}

# Subseting --------------------------------------------------------------

#New neural lineage genes old version
gene_list_neural_new <- c("PAX3", "PAX7", "OTX2", "PAX6", "EMX2",
                          "SRGAP3", "NAV1", "NRXN3", "EPHB1", "EPHB2")

#check if all genes are in the vst matrix
gene_list_neural_new[!gene_list_neural_new %in% rownames(vst_counts)]

#matrix creation and further subseting to keep all gene list = 10
{
  counts_neural <- vst_counts[gene_list_neural_new, ]
  z_matrix_neural <- t(scale(t(counts_neural)))
}


# Plotting ----------------------------------------------------------------
{
  df_long <- as.data.frame(z_matrix_neural) %>%   
    rownames_to_column("Gene") %>%
    pivot_longer(
      cols = -Gene,
      names_to = "Sample",
      values_to = "Zscore"
    )
  df_long$Sample <- factor(df_long$Sample, levels = c("AFF3_F_WT_1", 
                                                      "AFF3_F_WT_2", 
                                                      "AFF3_F_KO_1", 
                                                      "AFF3_F_KO_2"))
  df_long$Gene <- factor(df_long$Gene, levels = rev(unique(df_long$Gene)))
}
{
  lim <- max(abs(df_long$Zscore), na.rm = TRUE)
  q_lo <- quantile(abs(df_long$Zscore), 0.70, na.rm = TRUE) 
  q_hi <- quantile(abs(df_long$Zscore), 0.95, na.rm = TRUE) 
}
{
  stops <- c(-lim, -0.91, -0.792, 0, 0.792, 0.91, lim)
  vals  <- rescale(stops, to = c(0, 1))
  cols_new <- c("#2D004BFF",  
                "#542788FF", 
                "#8073ACFF", 
                "#B2ABD2FF", 
                "#D8DAEBFF",  
                "#F7F7F7FF", 
                "#FBE3E3",  
                "#F8B8B0", 
                "#EF7260",  
                "#CC2E2E",  
                "#7F0000")
  new_labels <- c("WT_1", "WT_2", "KO_1", "KO_2")
}



ggplot(df_long, aes(Sample, Gene, fill = Zscore)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = cols_new,
    values = vals, 
    limits = c(-lim, lim),
    oob = squish, 
    name = "Normalized expression",
    breaks = c(-1, 0, 1),
    guide = guide_colorbar(
      barwidth = 13,  
      barheight = 2.15, 
      ticks.colour = "black",
      direction = "horizontal",  
      title.position = "top",    
      title.hjust = 0
    )
  ) +
  scale_y_discrete(
    labels = rev(gene_list_neural_new), 
    position = "right",
    expand = c(0, 0)
  ) +
  scale_x_discrete(expand = c(0, 0),
                   labels = new_labels) +   
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 10, base_family = "Aptos") +
  theme(
    panel.grid = element_blank(),
    axis.text.y.right = element_text(size = 32,
                                     margin = margin(l = 2.5),
                                     hjust = 0,
                                     family = "Aptos", 
                                     face = "bold.italic"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 32),
    axis.text.x = element_text(
      size = 30,
      margin = margin(t = 2.5),
      angle = 90,   
      vjust = 0.5,
      hjust = 1,
      family = "Aptos", 
      face = "bold"
    ),
    legend.title = element_text(family = "Aptos", size = 28),
    legend.text = element_text(family = "Aptos", size = 28),
    legend.position = "bottom",   
    legend.box.margin = margin(t = 13),  
    legend.justification = "left",     
    legend.box.just = "left",
  ) +
  labs(x = NULL, y = NULL)


ggsave(
  "heatmap_neural.tiff",
  plot = last_plot(),
  width = 4.05,          
  height = 12.5,      
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  compression = "lzw"
)