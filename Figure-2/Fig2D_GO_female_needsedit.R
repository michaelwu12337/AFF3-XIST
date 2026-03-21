{
  setwd("C:/Users/Michael Wu/OneDrive/HKU PhD/Publication ready/Figure2D_GO")
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
  library(extrafont)
}

# DESeq2 female -----------------------------------------------------------
# Read count files
{
  d8_WT1_mm_counts_E <- read.table(file = "AFF3_F_WT_1_mm_counts_revised.txt", header = TRUE)
  d8_WT2_mm_counts_E <- read.table(file = "AFF3_F_WT_2_mm_counts_revised.txt", header = TRUE)
  d8_KO1_mm_counts_E <- read.table(file = "AFF3_F_KO_1_mm_counts_revised.txt", header = TRUE)
  d8_KO2_mm_counts_E <- read.table(file = "AFF3_F_KO_2_mm_counts_revised.txt", header = TRUE)
}

# Generate countData
{
  countData <- data.frame(
    b1W = d8_WT1_mm_counts_E$counts,
    b2W = d8_WT2_mm_counts_E$counts,
    b1K = d8_KO1_mm_counts_E$counts,
    b2K = d8_KO2_mm_counts_E$counts
  )
  rownames(countData) <- d8_WT1_mm_counts_E$Geneid
}

# Generate colData
{
  colData <- data.frame(
    condition = factor(c("WT", "WT", "KO", "KO"), levels = c("WT", "KO"))
  )
  rownames(colData) <- c("b1W", "b2W", "b1K", "b2K")
}

{
  all(colnames(countData) %in% rownames(colData))
  all(colnames(countData) == rownames(colData))
}

{
  dds_F <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ condition)
  keep_F <- rowSums(counts(dds_F)) >= 10
  dds_F <- dds_F[keep_F,]
  
  dds_F$Condition <- relevel(dds_F$condition, ref = "WT")
  
  dds_F <- estimateSizeFactors(dds_F)
  sizeFactors(dds_F)
  dds_F <- estimateDispersions(dds_F)
  
  dds_F <- DESeq(dds_F)
  res_F <- results(dds_F)
  
  sigs_F <- na.omit(res_F)
  sigs_F <- sigs_F[sigs_F$padj < 0.05,]
}


# GO Analysis Female upregulated ----------------------------------------
{
  genes_F_up <- rownames(sigs_F[sigs_F$log2FoldChange > 1,])
  GO_F_up <- enrichGO(gene = genes_F_up, 
                      OrgDb = "org.Hs.eg.db", 
                      keyType = "SYMBOL", 
                      ont = "ALL")
  GO_F_up_df <- as.data.frame(GO_F_up)
  GO_F_up_df <- GO_F_up_df %>%
    mutate(
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue), 
      Description = str_wrap(Description, width = 40)
    )
}

# GO Analysis Female downregulated ----------------------------------------
{
  genes_F_down <- rownames(sigs_F[sigs_F$log2FoldChange < -1,])
  GO_F_down <- enrichGO(gene = genes_F_down, 
                        OrgDb = "org.Hs.eg.db", 
                        keyType = "SYMBOL", 
                        ont = "ALL")
  GO_F_down_df <- as.data.frame(GO_F_down)
  GO_F_down_df <- GO_F_down_df %>%
    mutate(
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue), 
      Description = str_wrap(Description, width = 40)
    )
}

# Diverging plot female ---------------------------------------------------

F_up_df <- subset(GO_F_up_df, 
                  select = c(Description,
                             log2_OE,
                             neglog10p
                  ))
F_down_df <- subset(GO_F_down_df, 
                    select = c(Description,
                               log2_OE,
                               neglog10p
                    ))
F_up_df <- F_up_df %>%
  mutate(
    group = "Upregulated",
    Description = str_wrap(Description, width = 30)
  ) %>%
  slice_max(order_by = neglog10p, n = 8)
F_down_df <- F_down_df %>%
  slice_max(order_by = abs(neglog10p), n = 8) %>%
  mutate(
    group = "Downregulated",
    Description = str_wrap(Description, width = 30), 
    neglog10p = -neglog10p
  ) 
F_combined <- bind_rows(F_up_df, F_down_df)

ord_down_F <- F_combined %>%
  filter(group == "Downregulated") %>%
  arrange(desc(abs(neglog10p))) %>%
  pull(Description)
ord_up_F <- F_combined %>%
  filter(group == "Upregulated") %>%
  arrange(abs(neglog10p)) %>%   
  pull(Description)
ord_all_F <- c(ord_up_F, ord_down_F)

F_combined <- F_combined %>%
  mutate(Description = factor(Description, levels = ord_all_F))
blue_order <- c(
  "homophilic cell adhesion via plasma\nmembrane adhesion molecules",
  "cell-cell adhesion via plasma-membrane\nadhesion molecules",
  "axon guidance",
  "neuron projection guidance",
  "synapse assembly",
  "regulation of synapse organization",
  "axonogenesis",
  "forebrain development"
)
F_down_list <- F_combined %>%
  filter(group == "Downregulated") %>%
  mutate(y_blue = factor(Description, levels = blue_order))
terms_F <- c(
  "embryonic organ development",
  "mesenchyme development",
  "regulation of synapse organization",
  "neuron projection guidance",
  "axonogenesis"
)
terms_F <- str_wrap(terms_F, width = 30)

F_selected <- F_combined %>%
  filter(Description %in% terms_F)
F_selected_up <- subset(F_selected, group == "Upregulated")
F_selected_down <- subset(F_selected, group == "Downregulated")

F_diverge <- ggplot() +
  geom_col(
    data = subset(F_selected_down, group == "Downregulated"),
    # data = F_down_list,
    aes(x = neglog10p, y = Description, fill = log2_OE), #y = y_blue OR Description OR reorder(Description, rev(log2_OE))
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#D4DEE0", high = "#3C8DBC",
    name = "",
    breaks = range(F_selected_down$log2_OE),
    labels = scales::number_format(accuracy = 0.1),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", 
                           direction = "vertical",
                           barheight = unit(50, "pt"), barwidth  = unit(12, "pt"),
                           title.theme = element_text(margin = margin(b = 2))
    )
  ) +
  ggnewscale::new_scale("fill") +
  geom_col(
    data = subset(F_selected_up, group == "Upregulated"),
    aes(x = neglog10p, y = Description, fill = log2_OE), #y = Description OR reorder(Description, log2_OE)
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#EBD6D0", high = "salmon",
    name = "log2(O/E)",
    breaks = range(F_selected_up$log2_OE),
    labels = scales::number_format(accuracy = 0.1),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", 
                           direction = "vertical",
                           barheight = unit(50, "pt"), barwidth  = unit(12, "pt"),
                           title.theme = element_text(margin = margin(b = 15))
    )
  ) +
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "-log10(p-value)", y = NULL) +
  scale_x_continuous(breaks = c(-20, 0, 10)) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    legend.title = element_text(size = 28), 
    legend.text = element_text(size = 28)
  )
F_diverge


# Plot saving -------------------------------------------------------------
ggsave("GO_Female_Selected_textedit.tiff",
       plot = F_diverge,
       device = "tiff",
       dpi = 900,
       width = 10, height = 4, units = "in",
       compression = "lzw")

