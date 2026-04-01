{
  setwd("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Heat map/XIST WT & KO")
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
  library(extrafont)
}

loadfonts()

# Raw data processing -----------------------------------------------------
{
  F_KO1_counts <- read.table(file = "AFF3_F_KO_1_mm_counts_revised.txt", header = T)
  F_KO2_counts <- read.table(file = "AFF3_F_KO_2_mm_counts_revised.txt", header = T)
  F_WT1_counts <- read.table(file = "AFF3_F_WT_1_mm_counts_revised.txt", header = T)
  F_WT2_counts <- read.table(file = "AFF3_F_WT_2_mm_counts_revised.txt", header = T)
  F_DKO1_counts <- read.table(file = "AFF3_F_KO_Xt_KO_A3_1_mm_counts_revised.txt", header = T)
  F_DKO2_counts <- read.table(file = "AFF3_F_KO_Xt_KO_A3_2_mm_counts_revised.txt", header = T)
}
{
  counts_Data <- data.frame(AFF3_F_WT_1 = F_WT1_counts$counts)
  counts_Data$AFF3_F_WT_2 <- F_WT2_counts$counts
  counts_Data$AFF3_F_KO_1 <- F_KO1_counts$counts
  counts_Data$AFF3_F_KO_2 <- F_KO2_counts$counts
  counts_Data$AFF3_F_DKO_1 <- F_DKO1_counts$counts
  counts_Data$AFF3_F_DKO_2 <- F_DKO2_counts$counts
  rownames(counts_Data) <- F_WT1_counts$Geneid
}
head(counts_Data)
{
  col_Data <- data.frame(
    condition = factor(
      c("WT","WT","KO","KO","DKO","DKO"),
      levels = c("WT", "KO", "DKO")
    )
  )
  rownames(col_Data) <- c("AFF3_F_WT_1", "AFF3_F_WT_2", 
                          "AFF3_F_KO_1", "AFF3_F_KO_2",
                          "AFF3_F_DKO_1", "AFF3_F_DKO_2")
}
head(col_Data)
{
  all(colnames(counts_Data) %in% rownames(col_Data))
  all(colnames(counts_Data) == rownames(col_Data))
}

# DESeq2 ------------------------------------------------------------------
{
  dds <- DESeqDataSetFromMatrix(countData = counts_Data,
                                colData = col_Data,
                                design = ~ condition)
  dds <- DESeq(dds)
  counts <- assay(dds, normalized = F)
  norm_counts <- assay(dds, normalized = T)
  dds_vst <- vst(dds, blind=FALSE)
  vst_counts <- assay(dds_vst)
}

gene_list <- c("ALDH4A1", "ARNTL", "NAV1", "KLF6", "CAMK1G", "PLEKHA6", "EFS", "EPHB2",    
"GSN-AS1", "BTG2", "IFI16", "BCL2", "CDON", "PLEKHA7", "TLE3", "PAG1", "HHAT", "SDC4", "PKP1", 
"JAKMIP1", "TEAD1", "SLIT3", "PTPRT", "PRR9", "CD34", "HUNK", "PHC2", "STOM", "STRA6", "NRXN3",
"C5", "SORCS2", "KCNK17", "SKAP2", "TNFSF15", "PKNOX2", "PRELP", "EPHB3", "EPHB1", "TBC1D9",   
"KRT16", "KALRN", "TENM4", "DSCAML1", "KRT19", "KRT80", "CACNA1C", "KRT7", "EYA2", "KRT87P",
"FMOD", "TMPRSS4", "NCAM1", "RPH3A", "ALK", "ADAMTS18", "RIMBP2", "PLXNA2", "ALPK2", "SYBU",
"RASSF10", "DGKG", "C10orf71", "NCAM1-AS1")

{
  vst_counts_filtered <- vst_counts[gene_list, ]
  counts_filtered <- counts[gene_list, ]
  norm_counts_filtered <- norm_counts[gene_list, ]
  z_matrix_filtered <- t(scale(t(vst_counts_filtered)))
  z_matrix_filtered <- z_matrix_filtered[complete.cases(z_matrix_filtered), ]
  
  pos_WT <- z_matrix_filtered[
    apply(z_matrix_filtered[, c("AFF3_F_WT_1", "AFF3_F_WT_2")] > 0, 1, any), ]
  
  rows_to_remove <- c("ARNTL", "C10orf71")
  pos_WT <- pos_WT[ !rownames(pos_WT) %in% rows_to_remove, ]
  
  rows_to_add <- c("NRXN3", "CDON", "TLE3")
  pos_WT <- unique(rbind(
    pos_WT,
    z_matrix_filtered[rows_to_add, , drop = FALSE]
  ))
  pos_WT <- pos_WT[ rownames(z_matrix_filtered)[rownames(z_matrix_filtered) %in% rownames(pos_WT)], , drop = FALSE ]
  gene_list <- rownames(pos_WT)
  gene_list <- c("NAV1", "EFS", "EPHB2", "GSN-AS1", "CDON", "PLEKHA7", "TLE3", "HHAT", "PTPRT", 
                 "PHC2", "NRXN3", "C5", "SORCS2", "KCNK17", "PKNOX2", "TENM4", "RPH3A", "RASSF10")
}

{
  {
    df_long <- as.data.frame(pos_WT) %>%     #<-------------change dataset here
      rownames_to_column("Gene") %>%
      pivot_longer(
        cols = -Gene,
        names_to = "Sample",
        values_to = "Zscore"
      )
    df_long$Sample <- factor(df_long$Sample, levels = c("AFF3_F_DKO_1", "AFF3_F_DKO_2", 
                                                        "AFF3_F_KO_1", "AFF3_F_KO_2",
                                                        "AFF3_F_WT_1", "AFF3_F_WT_2"))
    df_long$Gene <- factor(df_long$Gene, levels = rev(unique(df_long$Gene)))
  }
  {
    new_labels <- c("DKO 2", "DKO 1", "KO 2", "KO 1", "WT 2", "WT 1")
    genes_to_label <- rownames(pos_WT) 
    gene_labels <- setNames(
      ifelse(levels(df_long$Gene) %in% genes_to_label, levels(df_long$Gene), ""),
      levels(df_long$Gene)
    )
    df_long$label_flag <- ifelse(df_long$Gene %in% genes_to_label, TRUE, FALSE)
    label_pos <- df_long %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarize(
        x_mid = median(as.numeric(Sample))
      ) %>%
      dplyr::filter(Gene %in% genes_to_label)
    label_pos$Gene <- factor(label_pos$Gene, levels = levels(df_long$Gene))
  }
  {
    lim <- max(abs(df_long$Zscore), na.rm = TRUE)
    q_lo <- quantile(abs(df_long$Zscore), 0.70, na.rm = TRUE) 
    q_hi <- quantile(abs(df_long$Zscore), 0.95, na.rm = TRUE) 
  }
  {
    stops <- c(-lim, -1.03, -0.64, -0.4, 0.14, 0.83, lim)
    vals  <- rescale(stops, to = c(0, 1))
    cols  <- c("#2C7BB6", 
               "#E6EEF6", 
               "white", 
               "#FDE0DD", 
               "#B2182B")
  }
}

ggplot(df_long, aes(Sample, Gene, fill = Zscore)) +
  geom_tile() +
  geom_segment(
    data = label_pos,
    inherit.aes = FALSE,
    aes(
      # start & end x just outside the last column of the heatmap
      x    = nlevels(df_long$Sample) + 0.5,
      xend = nlevels(df_long$Sample) + 0.65,
      y = Gene, yend = Gene
    ),
    color = "black", linewidth = 0.3, lineend = "round"
  ) +
  scale_fill_gradientn(
    colors = cols, 
    values = vals, 
    limits = c(-lim, lim),
    oob = squish, 
    name = "Normalized \nexpression",
    breaks = c(-1.3, 1.3),
    guide = guide_colorbar(
      barwidth = 2,
      barheight = 10, 
      # frame.colour = "black",
      ticks.colour = "black" 
    )
  ) +
  scale_y_discrete(
    labels = gene_labels,
    position = "right",
    expand = c(0, 0)
  ) +
  scale_x_discrete(labels = new_labels,
                   expand = c(0, 0)) +   
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x.top  = element_text(
      size = 25,
      angle = 45,        # <-- rotate labels
      hjust = 0,         # adjust horizontal alignment
      vjust = 0.02,       # adjust vertical position
      margin = margin(l = 1),
      face = "italic"
    ),
    # axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.y = element_blank(),   
    # axis.title.x = element_text(size = 12),
    # axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 25),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24)
  ) +
  labs(x = NULL, y = NULL)


ggsave(
  "heatmap_combined_short_WTtop_edit.tiff",
  plot = last_plot(),
  width = 11,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Fig4_heatmap.svg",
  plot     = last_plot(),
  device   = svglite::svglite,
  width    = 12,
  height   = 5,
  units    = "in"
) 
