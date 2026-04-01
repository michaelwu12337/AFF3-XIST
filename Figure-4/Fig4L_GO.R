{
  setwd("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Heat map/GO plot for XIST WT&KO&DKO")
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
}


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

counts_Data1 <- counts_Data[, !colnames(counts_Data) %in% c("AFF3_F_DKO_1", "AFF3_F_DKO_2")]
col_Data1 <- col_Data %>% filter(!rownames(col_Data) %in% c("AFF3_F_DKO_1", "AFF3_F_DKO_2"))

counts_Data2 <- counts_Data[, !colnames(counts_Data) %in% c("AFF3_F_KO_1", "AFF3_F_KO_2")]
col_Data2 <- col_Data %>% filter(!rownames(col_Data) %in% c("AFF3_F_KO_1", "AFF3_F_KO_2"))

# DESeq2 ------------------------------------------------------------------
{
  dds <- DESeqDataSetFromMatrix(countData = counts_Data,
                                colData = col_Data,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  res_KO_WT <- results(dds, contrast=c("condition","KO","WT"))
  results_KO_WT <- results(dds)
  res_DKO_KO <- results(dds, contrast=c("condition","DKO","WT"))
}
{
  filtered_KO_up <- res_KO_WT[res_KO_WT$log2FoldChange > 1 &
                                  res_KO_WT$padj < 0.05 &
                                  !is.na(res_KO_WT$log2FoldChange) &
                                  !is.na(res_KO_WT$padj), ]
  genes_KO_up <- rownames(filtered_KO_up)
  filtered_KO_down <- res_KO_WT[res_KO_WT$log2FoldChange < -1 &
                                res_KO_WT$padj < 0.05 &
                                !is.na(res_KO_WT$log2FoldChange) &
                                !is.na(res_KO_WT$padj), ]
  genes_KO_down <- rownames(filtered_KO_down)
  
  filtered_DKO_up <- res_DKO_KO[res_DKO_KO$log2FoldChange > 1 &
                                res_DKO_KO$padj < 0.05 &
                                !is.na(res_DKO_KO$log2FoldChange) &
                                !is.na(res_DKO_KO$padj), ]
  genes_DKO_up <- rownames(filtered_DKO_up)
  filtered_DKO_down <- res_DKO_KO[res_DKO_KO$log2FoldChange > 1 &
                                  res_DKO_KO$padj < 0.05 &
                                  !is.na(res_DKO_KO$log2FoldChange) &
                                  !is.na(res_DKO_KO$padj), ]
  genes_DKO_down <- rownames(filtered_DKO_down)
}



# GO Analysis -------------------------------------------------------------
{
  GO_KO_up <- enrichGO(gene = genes_KO_up, 
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "SYMBOL", 
                       ont = "ALL")
  df_KO_up <- as.data.frame(GO_KO_up)
  df_KO_up <- df_KO_up %>%
    mutate(
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue), 
      Description = str_wrap(Description, width = 40)
    )
  
  GO_KO_down <- enrichGO(gene = genes_KO_down, 
                         OrgDb = "org.Hs.eg.db", 
                         keyType = "SYMBOL", 
                         ont = "ALL")
  df_KO_down <- as.data.frame(GO_KO_down)
  df_KO_down <- df_KO_down %>%
    mutate(
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue), 
      Description = str_wrap(Description, width = 40)
    )
}
{
  gene_list <- c("NAV1", "EFS", "EPHB2", "GSN-AS1", "CDON", 
                 "PLEKHA7", "TLE3", "HHAT", "PTPRT", 
                 "PHC2", "NRXN3", "C5", "SORCS2", "KCNK17", 
                 "PKNOX2", "TENM4", "RPH3A", "RASSF10")
  gene_list_GO <- c("EPHB2", "CDON", "TENM4", "RASSF10", "PLEKHA7", "C5", 
                    "SORCS2", "RPH3A", "NRXN3", "PLEKHA7", "PTPRT")
  GO_selected <- enrichGO(gene = gene_list, 
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "SYMBOL", 
                          ont = "ALL",
                          pvalueCutoff = 1,
                          pAdjustMethod = "none",
                          qvalueCutoff = 1)
  df_selected <- as.data.frame(GO_selected)
  df_selected <- df_selected %>%
    mutate(
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue), 
      Description = str_wrap(Description, width = 40)
    )
}

df_filtered <- df_selected %>%
  filter(grepl("NAV1", geneID))

# Plotting ----------------------------------------------------------------

pl <- ggplot(df_selected %>% slice_max(order_by = neglog10p, n = 10), 
             aes(x = neglog10p, 
                 y = reorder(Description, neglog10p), 
                 fill = log2_OE)) +
  geom_col(color = "black") +
  labs(x = "-log10(p-value)", y = NULL, fill = "log2(O/E)") +
  scale_fill_gradient(low = "#E6E0DC", 
                      high = "salmon",   #F0AA8F
                      labels = scales::number_format(accuracy = 0.1),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

selected_terms <- c(
  "positive regulation of neurogenesis",
  "neuron spine",
  "synaptic membrane",
  "pore complex",
  "neuron to neuron synapse"
)

df_selected_filtered <- df_selected %>%
  filter(Description %in% selected_terms) %>%
  mutate(
    Description = str_to_title(Description),  # Capitalize the first letter of each word
    Description = str_wrap(Description, width = 22),  # Apply string wrapping
  )

{
  pl_selected <- ggplot(df_selected_filtered, 
                        aes(x = neglog10p, 
                            y = reorder(Description, neglog10p), 
                            fill = log2_OE)) +
    geom_col(color = "black") +
    labs(x = "-log10 (p-value)", y = NULL, fill = "log2(O/E)") +
    scale_fill_gradient(low = "#E8F1F2",
                        high = "#2A6F73",
                        labels = scales::number_format(accuracy = 0.1),
                        breaks = c(3.0, 5.0, 7.0),
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = "black")) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 28),
      axis.title.x = element_text(size = 28),
      axis.text.x = element_text(size = 24),
      legend.title = element_text(size = 28),
      legend.text = element_text(size = 24)
    )
  pl_selected
}

ggsave("GO_XIST.tiff",
       plot = pl,
       device = "tiff",
       dpi = 900,
       width = 10, height = 4, units = "in",
       compression = "lzw")

ggsave(
  filename = "Fig4_XIST_GO_edit.svg",
  plot     = pl_selected,
  device   = svglite::svglite,
  width    = 11,
  height   = 5.5,
  units    = "in"
)

