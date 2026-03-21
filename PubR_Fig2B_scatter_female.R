{
  setwd("C:/Users/Michael Wu/OneDrive/HKU PhD/Publication ready/Figure2B_Scatterplot")
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(extrafont) 
}

# Read count files
{
  d8_WT1_mm_counts_E <- read.table(file = "AFF3_F_WT_1_mm_counts_revised.txt", header = TRUE)
  d8_WT2_mm_counts_E <- read.table(file = "AFF3_F_WT_2_mm_counts_revised.txt", header = TRUE)
  d8_KO1_mm_counts_E <- read.table(file = "AFF3_F_KO_1_mm_counts_revised.txt", header = TRUE)
  d8_KO2_mm_counts_E <- read.table(file = "AFF3_F_KO_2_mm_counts_revised.txt", header = TRUE)
}

# Generate countData
{
  countData_d8_E <- data.frame(
    b1W = d8_WT1_mm_counts_E$counts,
    b2W = d8_WT2_mm_counts_E$counts,
    b1K = d8_KO1_mm_counts_E$counts,
    b2K = d8_KO2_mm_counts_E$counts
  )
  rownames(countData_d8_E) <- d8_WT1_mm_counts_E$Geneid
}

# Generate colData
{
  colData <- data.frame(
    condition = factor(c("WT", "WT", "KO", "KO"), levels = c("WT", "KO"))
  )
  rownames(colData) <- c("b1W", "b2W", "b1K", "b2K")
}

# Create DESeq2 object and run analysis
{
  dds_d8_AFF3 <- DESeqDataSetFromMatrix(
    countData = countData_d8_E,
    colData = colData,
    design = ~ condition
  )
  dds_d8_AFF3 <- DESeq(dds_d8_AFF3)
}

# Get results
res_d8_AFF3 <- results(dds_d8_AFF3, contrast = c("condition", "KO", "WT"))
res_d8_AFF3 <- as.data.frame(res_d8_AFF3)

# Add gene names
res_d8_AFF3$gene <- rownames(res_d8_AFF3)

# Calculate mean log2 counts for WT and KO
norm_counts <- counts(dds_d8_AFF3, normalized = TRUE)

# Calculate mean counts (adding pseudocount to avoid log(0))
res_d8_AFF3$WT_mean <- rowMeans(norm_counts[, c("b1W", "b2W")]) + 1
res_d8_AFF3$KO_mean <- rowMeans(norm_counts[, c("b1K", "b2K")]) + 1

# Convert to log2 scale
res_d8_AFF3$WT_log2 <- log2(res_d8_AFF3$WT_mean)
res_d8_AFF3$KO_log2 <- log2(res_d8_AFF3$KO_mean)

# Define regulation status with specific colors
res_d8_AFF3$regulation <- "Not DEG"
res_d8_AFF3$regulation[res_d8_AFF3$padj < 0.05 & res_d8_AFF3$log2FoldChange > 1] <- "Up-regulated"
res_d8_AFF3$regulation[res_d8_AFF3$padj < 0.05 & res_d8_AFF3$log2FoldChange < -1] <- "Down-regulated"

# Create two groups for non-DEGs based on expression level
non_deg_data <- res_d8_AFF3[res_d8_AFF3$regulation == "Not DEG", ]
median_expr <- median(non_deg_data$WT_log2 + non_deg_data$KO_log2, na.rm = TRUE)

res_d8_AFF3$non_deg_group <- "Not DEG - Low"
res_d8_AFF3$non_deg_group[res_d8_AFF3$regulation == "Not DEG" & 
                            (res_d8_AFF3$WT_log2 + res_d8_AFF3$KO_log2) > median_expr] <- "Not DEG - High"
res_d8_AFF3$non_deg_group[res_d8_AFF3$regulation != "Not DEG"] <- res_d8_AFF3$regulation[res_d8_AFF3$regulation != "Not DEG"]

# Identify genes to label
genes_to_label <- c("SIX3", "PAX6", "PAX7", "XIST")
res_d8_AFF3$label <- ifelse(res_d8_AFF3$gene %in% genes_to_label, res_d8_AFF3$gene, "")

# Create the scatterplot
ggplot(res_d8_AFF3, aes(x = KO_log2, y = WT_log2)) +
  # Plot non-DEG groups with different shapes but no legend
  geom_point(data = subset(res_d8_AFF3, regulation == "Not DEG" & non_deg_group == "Not DEG - Low"),
             color = "gray70", alpha = 0.6, size = 1) +
  geom_point(data = subset(res_d8_AFF3, regulation == "Not DEG" & non_deg_group == "Not DEG - High"),
             color = "gray70", alpha = 0.6, size = 1) +
  
  # Plot DEGs with specified colors
  geom_point(data = subset(res_d8_AFF3, regulation == "Up-regulated"),
             color = "#85312B", alpha = 0.8, size = 1.5) +
  geom_point(data = subset(res_d8_AFF3, regulation == "Down-regulated"),
             color = "#373671", alpha = 0.8, size = 1.5) +
  
  # Add diagonal line (y = x)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue", alpha = 0.5) +
  
  # Label representative genes
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
    title = "RNA-seq (WT vs KO)",
    x = "KO (Log2 values)",
    y = "WT (Log2 values)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Aptos", size = 50),  
    axis.title = element_text(family = "Aptos", size = 50), 
    axis.text = element_text(family = "Aptos", size = 50) 
  ) +
  coord_fixed(ratio = 1)


ggsave("Scatterplot WT vs KO female.tiff",
       plot = last_plot(),
       device = "tiff",
       dpi = 900,
       width = 9, height = 9, units = "in",
       compression = "lzw")


