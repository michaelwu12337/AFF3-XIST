#Dotplot or barplot for Bulk-RNA-seq dataset
{
  setwd("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/new plot")
  library(DESeq2)
  library(ggplot2)
  library(tidyverse)
  library(RColorBrewer)
  library(extrafont)
  library(dplyr)
  library(grid)      
  library(patchwork) 
}


# 1. Define gene sets as NAMED list

gene_sets <- list(
  Advanced_Mesoderm = c("VCAN", "COL6A1", "BMP5", "COL6A2", "LIX1", "PITX1", "FN1", "COL6A3", "RGS4", "TENM4", 
                        "SDK1", "PRRX1", "CYB5A", "MAP1B", "BAMBI", "HAND2", "FEZ1", "FGFR2", "PGM5", "COL3A1"),
  
  Amnion = c("KRT19", "WFDC2", "TINAGL1", "ARL4C", "ID1", "GAPDH", "GABRP", "TBCA", 
             "CLDN6", "SLC2A3", "VTCN1", "SPINT2", "KCNK12", "CLDN4", "PRSS8", "TRIML2", 
             "TFAP2A", "SINHCAF", "PRR15", "EPCAM"),
  
  Blood_Progenitors = c("GPX1", "RGS10", "PLEK", "ITGA2B", "PFN1", "PF4", "TESC", "TAGLN2", "FERMT3", 
                        "PRKAR2B", "CCND3", "SH3BGRL3", "CMTM5", "IFITM2", "NFE2", "FCER1G", "GMPR", "S100A4", "TALDO1"),
  
  Hemogenic_Endothelium = c("GNG11", "ECSCR", "VAMP5", "CALCRL", "RAMP2", "IGFBP4", "PRCP", "S100A16", "CDH5", "KDR", 
                            "MEF2C", "PECAM1", "TIE1", "LIMCH1", "CCDC85B", "IGFBP2", "ICAM2", "MARCKS", "PLVAP", "ESAM"),
  
  Yolk_Sac_Endoderm = c("APOA1", "APOA2", "TTR", "RBP4", "LGALS3", "AFP", "SERPINA1", "APOC1", "FGB", "FABP1", 
                        "GSTA2", "APOC3", "AHSG", "VTN", "S100A14", "ACAA1", "APOE", "GSTA1", "FGG", "DUSP9"),
  
  Yolk_Sac_Mesoderm = c("LUM", "H19", "IGF2", "LUM", "SPARC", "FRZB", "COL1A1", "TGFBI", "COL1A2", "HAND1", 
                        "EPHA7", "COL3A1", "MEST", "PTN", "DSG2", "IGFBP3", "IGDCC3", "POSTN", "KRT18", "CD248", "PARM1"),
  
  Retinal_Progenitors = c("SIX6", "FEZF2", "PAX6", "SIX3", "RAX", "FEZF1", "EMX2", "ZNF593OS", 
                          "FEZF1-AS1", "CLEC19A", "GRP", "LHX5-AS1", "HES5", "BCO1", "LHX2", 
                          "CLDN1", "RMST", "NOS2", "EMX2OS"),
  
  Neural_Progenitors = c("LRRC53", "NEUROD1", "LINC02217", "PPP1R17", "POU4F1", "NHLH1", "KLHL40", "NEUROD4", "ONECUT1", "INSM1", 
                         "TLX3", "ENSG00000259203", "B3GALT2", "ELAVL3", "SVOP", "THSD7B", "ATP2B2", "P2RX3", "NEUROG1", "SIX1")
)

gene_cluster <- read.csv("celltype_v2_top50_markers.csv")
{
  genes_amnion <- gene_cluster %>%
    filter(celltype == "Amnion") %>%
    pull(gene)
  genes_blood <- gene_cluster %>%
    filter(celltype == "Blood Progenitors") %>%
    pull(gene)
  genes_cardiac <- gene_cluster %>%
    filter(celltype == "Cardiac Mesoderm") %>%
    pull(gene)
  genes_YSE <- gene_cluster %>%
    filter(celltype == "YSE") %>%
    pull(gene)
  genes_neural_new <- gene_cluster %>%
    filter(celltype == "Neural Progenitors") %>%
    pull(gene)
  genes_neuron <- gene_cluster %>%
    filter(celltype == "Neurons") %>%
    pull(gene)
  genes_PGC <- gene_cluster %>%
    filter(celltype == "PGC") %>%
    pull(gene)
  genes_PSC_EPI <- gene_cluster %>%
    filter(celltype == "PSC/Epi") %>%
    pull(gene)
  genes_TE <- gene_cluster %>%
    filter(celltype == "TE") %>%
    pull(gene)
  genes_VEC <- gene_cluster %>%
    filter(celltype == "VEC") %>%
    pull(gene)
}
cat(paste(shQuote(genes_amnion), collapse = ", "))

gene_sets_new <- list(
  Amnion = genes_amnion,
  Blood_Progenitors = genes_blood,
  Cardiac_Mesoderm = genes_cardiac,
  YSE = genes_YSE,
  Neural_Progenitors = genes_neural_new,
  Neurons = genes_neuron,
  PGC = genes_PGC,
  Epi = genes_PSC_EPI,
  TE = genes_TE,
  VEC = genes_VEC
)

# 1a. Data import and DESeq2
{
  d8_ko_1_1_mm_counts_F <- read.table(file = "d8_ko_1_1_mm_counts_revised.txt", header = T)
  d8_ko_1_2_mm_counts_F <- read.table(file = "d8_ko_1_2_mm_counts_revised.txt", header = T)
  d8_ko_2_1_mm_counts_F <- read.table(file = "d8_ko_2_1_mm_counts_revised.txt", header = T)
  d8_ko_2_2_mm_counts_F <- read.table(file = "d8_ko_2_2_mm_counts_revised.txt", header = T)
  d8_wt_1_mm_counts_F <- read.table(file = "d8_wt_1_mm_counts_revised.txt", header = T)
  d8_wt_2_mm_counts_F <- read.table(file = "d8_wt_2_mm_counts_revised.txt", header = T)
}
{
  countData <- data.frame(ko_1_1 = d8_ko_1_1_mm_counts_F$counts)
  countData$ko_1_2 <- d8_ko_1_2_mm_counts_F$counts
  countData$ko_2_1 <- d8_ko_2_1_mm_counts_F$counts
  countData$ko_2_2 <- d8_ko_2_2_mm_counts_F$counts
  countData$wt_1 <- d8_wt_1_mm_counts_F$counts
  countData$wt_2 <- d8_wt_2_mm_counts_F$counts
  rownames(countData) <- d8_wt_2_mm_counts_F$Geneid
}

{
  colData <- data.frame(
    condition = factor(
      c("KO","KO","KO","KO","WT","WT"),
      levels = c("KO", "WT")
    )
  )
  rownames(colData) <- c("ko_1_1", "ko_1_2",
                         "ko_2_1", "ko_2_2",
                         "wt_1", "wt_2")
}

{
  all(colnames(countData) %in% rownames(colData))
  all(colnames(countData) == rownames(colData))
}
{
  dds_d8_F <- DESeqDataSetFromMatrix(countData = countData, 
                                     colData = colData, 
                                     design = ~ condition)
  dds_d8_F <- DESeq(dds_d8_F)
  res_d8_F <- results(dds_d8_F, contrast=c("condition","KO","WT"))
}


# 2. Extract normalized counts (use consistent DESeqDataSet name)
norm_counts <- counts(dds_d8_F, normalized = TRUE) %>%  # Changed dds_d8_F to dds
  as.data.frame() %>%
  rownames_to_column("gene")

# 3. Prepare metadata - extract condition from sample names
metadata <- as.data.frame(colData(dds_d8_F)) %>%
  rownames_to_column("sample") %>%
  mutate(condition = case_when(
    grepl("w", sample) ~ "WT",  # Samples ending with "W"
    grepl("k", sample) ~ "KO"    # Samples containing "K"
  )) %>%
  dplyr::select(sample, condition)

# 4. Calculate expression metrics (simplified and robust)
results_df <- map_dfr(names(gene_sets_new), function(set_name) {
  genes <- gene_sets_new[[set_name]]
  pseudocount <- 0.1  # Prevents log(0) errors
  
  # Get expression data for the gene set
  expr_data <- norm_counts %>%
    dplyr::filter(gene %in% genes) %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    left_join(metadata, by = "sample")
  
  # Calculate average expression per condition
  set_avgs <- expr_data %>%
    group_by(condition) %>%
    summarise(avg_expression = mean(expression, na.rm = TRUE), .groups = "drop")
  
  # Extract values
  wt_avg <- set_avgs$avg_expression[set_avgs$condition == "WT"]
  ko_avg <- set_avgs$avg_expression[set_avgs$condition == "KO"]
  
  # Create result tibble
  tibble(
    GeneSet = set_name,
    condition = c("WT", "KO"),
    avg_expression = c(wt_avg, ko_avg),
    log_ratio = c(
      log2((wt_avg + pseudocount) / (ko_avg + pseudocount)),
      log2((ko_avg + pseudocount) / (wt_avg + pseudocount))
    )
  )
})

results_df <- results_df %>%
  mutate(
    GeneSet = factor(GeneSet, levels = c(
      "Amnion",
      "Blood_Progenitors",
      "Cardiac_Mesoderm",
      "YSE",
      "Neural_Progenitors",
      "Neurons",
      "PGC",
      "Epi",
      "TE",
      "VEC"
    )),
    condition = factor(condition, levels = c("WT", "KO"))
  )

max_abs_log <- max(abs(results_df$log_ratio), na.rm = TRUE)

p_bubble <- ggplot(results_df, aes(x = condition, y = GeneSet)) +
  geom_point(
    aes(size = avg_expression, fill = log_ratio),
    colour = "black",
    stroke  = 0.9,       # thicker border around circles
    shape   = 21,
    alpha   = 0.95
  ) +
  scale_size_continuous(
    name  = "Average \nExpression",
    range = c(6, 14),
    breaks = pretty(results_df$avg_expression, n = 3)
  ) +
  scale_fill_distiller(
    name      = "Log2 Ratio\n(KO/WT)",
    type      = "div",
    palette   = "PuOr",
    direction = -1,
    limits    = c(-max_abs_log, max_abs_log)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    # panel.background = element_rect(fill = "grey98", colour = "grey70", linewidth = 0.5),
    panel.background = element_rect(fill = "grey98", colour = NA),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, colour = "black", face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(face = "bold", size = 20),
    legend.text  = element_text(size = 20),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = margin(15, 15, 15, 15),
    # plot.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.5),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  guides(
    size = guide_legend(override.aes = list(
        size = c(4, 7, 10),  # smaller sizes for legend bubbles
        fill = "grey70", 
        colour = "black"),
      keyheight = unit(0.6, "lines"),     # reduce vertical spacing
      keywidth  = unit(0.6, "lines")
    ),
    fill = guide_colorbar(
      barwidth = 6,
      barheight = 1,
      ticks = FALSE,
      label = FALSE,
      frame.colour = NA      # <-- remove outline around color bar
    )
  )


# Create the diverging plot
# First, calculate the percentages and create deg_df
padj_threshold <- 0.05

# Calculate DEG percentages
deg_df <- map_dfr(names(gene_sets_new), function(set_name) {
  genes <- gene_sets_new[[set_name]]
  
  # Get results for genes in set
  set_res <- as.data.frame(res_d8_F) %>%
    rownames_to_column("gene") %>%
    filter(gene %in% genes) %>%
    mutate(
      direction = case_when(
        padj < padj_threshold & log2FoldChange > 0 ~ "Up",
        padj < padj_threshold & log2FoldChange < 0 ~ "Down",
        TRUE ~ "Non-sig"
      )
    )
  
  # Count DEGs
  deg_counts <- table(set_res$direction)
  
  # Calculate percentages
  tibble(
    GeneSet = set_name,
    Direction = c("Down", "Up"),
    Percentage = c(
      ifelse(is.na(deg_counts["Down"]), 0, (deg_counts["Down"] / length(genes)) * 100),
      ifelse(is.na(deg_counts["Up"]), 0, (deg_counts["Up"] / length(genes)) * 100)
    )
  )
})

# Now add x_position column to deg_df
deg_df <- deg_df %>%
  mutate(x_position = ifelse(Direction == "Down", -Percentage, Percentage))

# Create plot with dynamic limits
max_x <- max(abs(deg_df$x_position)) * 1.1

# Calculate position for labels
n_sets <- length(unique(deg_df$GeneSet))

plot_order <- c(
  "Amnion",
  "Blood_Progenitors",
  "Cardiac_Mesoderm",
  "YSE",
  "Neural_Progenitors",
  "Neurons",
  "PGC",
  "PSC_Epi",
  "TE",
  "VEC"
)

deg_df$GeneSet <- factor(deg_df$GeneSet, levels = plot_order)
label_y <- length(unique(deg_df$GeneSet)) + 2.5   # push OUTSIDE the panel

p_diverge <- ggplot(deg_df, aes(y = GeneSet, x = x_position, fill = Direction)) +
  geom_col(width = 0.8) +
  geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
  scale_fill_manual(values = c(
    Down = "#3C8DBC",
    Up   = "#E36B6B"
  )) +
  scale_x_continuous(
    name   = "% Differentially Expressed Genes",
    # labels = function(x) abs(x),
    labels = c("80", "40", "0", "40", "80"),
    # breaks = seq(-100, 100, 20),
    breaks = c(-80, -40, 0, 40, 80),
    limits = c(-max_x, max_x)
  ) +
  # normal expansion; place UP/DOWN using y = Inf 
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  labs(
    y = NULL,   # remove y-axis label
    title = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.3),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
    axis.text.x = element_text(size = 20, colour = "black"),
    # axis.text.y = element_text(family = "Aptos", face = "bold", size = 16, colour = "black"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 20, 
                              face = "bold", 
                              hjust = 0.5,
                              margin = margin(t = 5, b = 30)),
    legend.position = "none",
    
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),

    # plot.margin = margin(10, 10, 10, 0),
    
    plot.margin = margin(t = 30, r = 15, b = 15, l = 15)
  ) +
  coord_cartesian(clip = "off") +
  # DOWN label truly outside (above) the panel
  annotate(
    "text",
    x = -max_x * 0.45,
    y = Inf,
    label = "DOWN",
    vjust = -1,   
    size  = 7,
    colour = "#3C8DBC",
    fontface = "bold"
  ) +
  # UP label outside the panel
  annotate(
    "text",
    x =  max_x * 0.45,
    y =  Inf,
    label = "UP",
    vjust = -1,
    size  = 7,
    colour = "#E36B6B",
    fontface = "bold"
  )

p_combined <- p_bubble + p_diverge +
  plot_layout(widths = c(1, 2))  
p_combined

ggsave(
  "Bubble_and_Diverging_edit_03.tiff",
  plot = p_combined,
  width = 9.8,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


