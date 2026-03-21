{
  setwd("C:/Users/Michael Wu/OneDrive/HKU PhD/Publication ready/Figure2C_Venn")
  library(VennDiagram)
  library(grid)
  library(gridSVG)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(extrafont)
}

###For female AFF3 WT & AFF3 KO
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

{
  # Get results
  res_d8_AFF3 <- results(dds_d8_AFF3, contrast = c("condition", "KO", "WT"))
  res_d8_AFF3 <- as.data.frame(res_d8_AFF3)
  
  #Calculate the up-regulated genes
  filtered_res_d8_AFF3_up <- res_d8_AFF3[res_d8_AFF3$log2FoldChange > 0.5 & res_d8_AFF3$padj < 0.05 & !is.na(res_d8_AFF3$log2FoldChange) & !is.na(res_d8_AFF3$padj), ]
  ordered_res_d8_AFF3_up <- filtered_res_d8_AFF3_up[order(filtered_res_d8_AFF3_up$log2FoldChange), ]
  gene_names_d8_AFF3_up <- rownames(ordered_res_d8_AFF3_up)
  
  #Calculate the down-regulated genes
  filtered_res_d8_AFF3_down <- res_d8_AFF3[res_d8_AFF3$log2FoldChange < -0.5 & res_d8_AFF3$padj < 0.05 & !is.na(res_d8_AFF3$log2FoldChange) & !is.na(res_d8_AFF3$padj), ]
  ordered_res_d8_AFF3_down <- filtered_res_d8_AFF3_down[order(filtered_res_d8_AFF3_down$log2FoldChange), ]
  gene_names_d8_AFF3_down <- rownames(ordered_res_d8_AFF3_down)
}



###For female AFF3 WT & AFF3 KO
# Read count files
{
  d8_WT1_mm_counts_M <- read.table(file = "AFF3_M_WT_1_mm_counts_revised.txt", header = TRUE)
  d8_WT2_mm_counts_M <- read.table(file = "AFF3_M_WT_2_mm_counts_revised.txt", header = TRUE)
  d8_KO1_mm_counts_M <- read.table(file = "AFF3_M_KO_1_1_mm_counts_revised.txt", header = TRUE)
  d8_KO2_mm_counts_M <- read.table(file = "AFF3_M_KO_1_2_mm_counts_revised.txt", header = TRUE)
}

# Generate countData
{
  countData_d8_M <- data.frame(
    b1W = d8_WT1_mm_counts_M$counts,
    b2W = d8_WT2_mm_counts_M$counts,
    b1K = d8_KO1_mm_counts_M$counts,
    b2K = d8_KO2_mm_counts_M$counts
  )
  rownames(countData_d8_M) <- d8_WT1_mm_counts_M$Geneid
}

# Generate colData
{
  colData_M <- data.frame(
    condition = factor(c("WT", "WT", "KO", "KO"), levels = c("WT", "KO"))
  )
  rownames(colData_M) <- c("b1W", "b2W", "b1K", "b2K")
}

# Create DESeq2 object and run analysis
{
  dds_d8_MA <- DESeqDataSetFromMatrix(
    countData = countData_d8_M,
    colData = colData_M,
    design = ~ condition
  )
  dds_d8_MA <- DESeq(dds_d8_MA)
}

{
  # Get results
  res_d8_MA <- results(dds_d8_MA, contrast = c("condition", "KO", "WT"))
  res_d8_MA <- as.data.frame(res_d8_MA)
  
  #Calculate the up-regulated genes
  filtered_res_d8_MA_up <- res_d8_MA[res_d8_MA$log2FoldChange > 0.5 & res_d8_MA$padj < 0.05 & !is.na(res_d8_MA$log2FoldChange) & !is.na(res_d8_MA$padj), ]
  ordered_res_d8_MA_up <- filtered_res_d8_MA_up[order(filtered_res_d8_MA_up$log2FoldChange), ]
  gene_names_d8_MA_up <- rownames(ordered_res_d8_MA_up)
  
  #Calculate the down-regulated genes
  filtered_res_d8_MA_down <- res_d8_MA[res_d8_MA$log2FoldChange < -0.5 & res_d8_MA$padj < 0.05 & !is.na(res_d8_MA$log2FoldChange) & !is.na(res_d8_MA$padj), ]
  ordered_res_d8_MA_down <- filtered_res_d8_MA_down[order(filtered_res_d8_MA_down$log2FoldChange), ]
  gene_names_d8_MA_down <- rownames(ordered_res_d8_MA_down)
}


# Create Venn diagram for downregulated genes
venn_plot_down <- venn.diagram(
  x = list(
    AFF3 = gene_names_d8_AFF3_down,
    MA = gene_names_d8_MA_down
  ),
  category.names = c("Female Down", "Male Down"),
  filename = NULL,
  output = TRUE,
  
  # Colors and transparency
  fill = c("#BE9EB0", "#A0B5D1"),
  alpha = 0.5,
  
  # Circle outlines
  lwd = 2,
  lty = 'solid',
  
  # Label formatting
  cex = 3,
  fontface = "bold",
  fontfamily = "Aptos",
  
  # Category label positioning
  cat.cex = 3.6,
  cat.fontface = "bold",
  cat.fontfamily = "Aptos",
  cat.pos = c(-25, 10),
  cat.dist = c(0.05, 0.06),
  # cat.default.pos = "outer",
  
  # Display only counts
  print.mode = c("raw"),
  sigdigs = 3,
  
  # Make circles oval
  rotation.degree = 0,
  compression = "euler",
  
  # Additional parameters
  ext.text = FALSE,
  ext.line.lwd = 0,
  ext.dist = 0,
  ext.length = 0,
  ext.pos = 0,
  
  # Margin adjustment
  margin = 0.04
)

# Save downregulated Venn diagram as SVG with transparent background
{
  svg("venn_downregulated.svg", bg = "transparent", width = 10.5, height = 6)
  grid.newpage()
  grid.draw(venn_plot_down)
  dev.off()
}
