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

# DESeq2 male -------------------------------------------------------------
# Read count files
{
  d8_WT1_mm_counts_E <- read.table(file = "AFF3_M_WT_1_mm_counts_revised.txt", header = TRUE)
  d8_WT2_mm_counts_E <- read.table(file = "AFF3_M_WT_2_mm_counts_revised.txt", header = TRUE)
  d8_KO1_mm_counts_E <- read.table(file = "AFF3_M_KO_1_1_mm_counts_revised.txt", header = TRUE)
  d8_KO2_mm_counts_E <- read.table(file = "AFF3_M_KO_1_2_mm_counts_revised.txt", header = TRUE)
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
  dds_M <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ condition)
  keep_M <- rowSums(counts(dds_M)) >= 10
  dds_M <- dds_M[keep_M,]
  
  dds_M$Condition <- relevel(dds_M$condition, ref = "WT")
  
  dds_M <- estimateSizeFactors(dds_M)
  sizeFactors(dds_M)
  dds_M <- estimateDispersions(dds_M)
  
  dds_M <- DESeq(dds_M)
  res_M <- results(dds_M)
  
  sigs_M <- na.omit(res_M)
  sigs_M <- sigs_M[sigs_M$padj < 0.05,]
}

# GO Analysis Male upregulated--------------------------------------------------------
{
  genes_M_up <- rownames(sigs_M[sigs_M$log2FoldChange > 1,])
  GO_M_up <- enrichGO(gene = genes_M_up, 
                      OrgDb = "org.Hs.eg.db", 
                      keyType = "SYMBOL", 
                      ont = "ALL")
  GO_M_up_df <- as.data.frame(GO_M_up)
  GO_M_up_df <- GO_M_up_df %>%
    mutate(
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue), 
      Description = str_wrap(Description, width = 40)
    )
}

# GO Analysis Male downregulated ------------------------------------------
{
  genes_M_down <- rownames(sigs_M[sigs_M$log2FoldChange < -1,])
  GO_M_down <- enrichGO(gene = genes_M_down, 
                        OrgDb = "org.Hs.eg.db", 
                        keyType = "SYMBOL", 
                        ont = "ALL")
  GO_M_down_df <- as.data.frame(GO_M_down)
  GO_M_down_df <- GO_M_down_df %>%
    mutate(
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue), 
      Description = str_wrap(Description, width = 40)
    )
}

# Diverging plot male---------------------------------------------------------
M_up_df <- subset(GO_M_up_df, 
                  select = c(Description,
                             log2_OE,
                             neglog10p
                  ))
M_down_df <- subset(GO_M_down_df, 
                    select = c(Description,
                               log2_OE,
                               neglog10p
                    ))
M_up_df <- M_up_df %>%
  mutate(
    group = "Upregulated",
    Description = str_wrap(Description, width = 30)
  ) %>%
  slice_max(order_by = neglog10p, n = 8)
M_down_df <- M_down_df %>%
  slice_max(order_by = abs(neglog10p), n = 8) %>%
  mutate(
    group = "Downregulated",
    Description = str_wrap(Description, width = 30), 
    neglog10p = -neglog10p
  ) 
M_combined <- bind_rows(M_up_df, M_down_df)

ord_down <- M_combined %>%
  filter(group == "Downregulated") %>%
  arrange(desc(abs(neglog10p))) %>%
  pull(Description)
ord_up <- M_combined %>%
  filter(group == "Upregulated") %>%
  arrange(abs(neglog10p)) %>%   
  pull(Description)
ord_all <- c(ord_up, ord_down)

M_combined <- M_combined %>%
  mutate(Description = factor(Description, levels = ord_all))
terms_M <- c(
  "skin epidermis development",
  "retina layer formation",
  "somitogenesis",
  "hemoglobin complex",
  "oxygen transport"
)
M_selected <- M_combined %>%
  filter(Description %in% terms_M)
# M_selected$Description <- str_replace(M_selected$Description,
#                                       "^(\\w)", 
#                                       toupper)
M_selected_up <- subset(M_selected, group == "Upregulated")
M_selected_down <- subset(M_selected, group == "Downregulated")

M_diverge <- ggplot() +
  geom_col(
    data = subset(M_selected_down, group == "Downregulated"),
    aes(x = neglog10p, y = Description, fill = log2_OE), #y = Description OR reorder(Description, rev(log2_OE))
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#D4DEE0", high = "#3C8DBC",
    name = "",
    breaks = range(M_selected_down$log2_OE),
    labels = scales::number_format(accuracy = 0.1),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", 
                           direction = "vertical",
                           barheight = unit(50, "pt"), barwidth  = unit(12, "pt"),
                           title.theme = element_text(margin = margin(b = 2))
    )
  ) +
  ggnewscale::new_scale("fill") +
  geom_col(
    data = subset(M_selected_up, group == "Upregulated"),
    aes(x = neglog10p, y = Description, fill = log2_OE), #y = Description OR reorder(Description, log2_OE)
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#EBD6D0", high = "salmon",
    name = "log2(O/E)",
    breaks = range(M_selected_up$log2_OE),
    labels = scales::number_format(accuracy = 0.1),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", 
                           direction = "vertical",
                           barheight = unit(50, "pt"), barwidth  = unit(12, "pt"),
                           title.theme = element_text(margin = margin(b = 15))
    )
  ) +
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "-log10(p-value)", y = NULL) +
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
M_diverge


# Plot saving -------------------------------------------------------------
ggsave("GO_Male_Selected_textedit.tiff",
       plot = M_diverge,
       device = "tiff",
       dpi = 900,
       width = 10, height = 4, units = "in",
       compression = "lzw")


