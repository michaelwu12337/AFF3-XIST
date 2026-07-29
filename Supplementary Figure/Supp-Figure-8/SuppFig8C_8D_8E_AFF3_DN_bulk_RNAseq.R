{
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  library(pheatmap)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  table_dir <- "generated_tables"
  fig8c_dir <- file.path("generated_figures", "SuppFig8C")
  fig8d_dir <- file.path("generated_figures", "SuppFig8D")
  fig8e_dir <- file.path("generated_figures", "SuppFig8E")

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig8c_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig8d_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig8e_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load WT and AFF3 dominant-negative count tables ------------------------
{
  WT_1_1_counts <- read.table(
    file.path(source_dir, "D1_1_mm_counts_revised.txt"),
    header = TRUE
  )
  WT_1_2_counts <- read.table(
    file.path(source_dir, "D1_2_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_1_1_counts <- read.table(
    file.path(source_dir, "D2_1_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_1_2_counts <- read.table(
    file.path(source_dir, "D2_2_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_2_1_counts <- read.table(
    file.path(source_dir, "D3_1_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_2_2_counts <- read.table(
    file.path(source_dir, "D3_2_mm_counts_revised.txt"),
    header = TRUE
  )

  count_tables <- list(
    WT_1_1_counts,
    WT_1_2_counts,
    DN_1_1_counts,
    DN_1_2_counts,
    DN_2_1_counts,
    DN_2_2_counts
  )
  stopifnot(all(vapply(
    count_tables[-1],
    function(x) identical(x$Geneid, count_tables[[1]]$Geneid),
    logical(1)
  )))
}

# Supplementary Figure 8C–D: DESeq2, scatterplot, and GO analysis --------
{
  counts_Data <- data.frame(
    AFF3_WT_1_1 = WT_1_1_counts$counts,
    AFF3_WT_1_2 = WT_1_2_counts$counts,
    AFF3_DN_1_1 = DN_1_1_counts$counts,
    AFF3_DN_1_2 = DN_1_2_counts$counts,
    AFF3_DN_2_1 = DN_2_1_counts$counts,
    AFF3_DN_2_2 = DN_2_2_counts$counts
  )
  rownames(counts_Data) <- WT_1_1_counts$Geneid

  col_Data <- data.frame(
    condition = c("WT", "WT", "DN", "DN", "DN", "DN")
  )
  rownames(col_Data) <- colnames(counts_Data)

  dds <- DESeqDataSetFromMatrix(
    countData = counts_Data,
    colData = col_Data,
    design = ~ condition
  )
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", "DN", "WT"))
  res <- as.data.frame(res)
  res$gene <- rownames(res)

  norm_counts <- counts(dds, normalized = TRUE)
  res$WT_mean <- rowMeans(
    norm_counts[, c("AFF3_WT_1_1", "AFF3_WT_1_2")]
  ) + 1
  res$DN_mean <- rowMeans(
    norm_counts[, c(
      "AFF3_DN_1_1",
      "AFF3_DN_1_2",
      "AFF3_DN_2_1",
      "AFF3_DN_2_2"
    )]
  ) + 1
  res$WT_log2 <- log2(res$WT_mean)
  res$DN_log2 <- log2(res$DN_mean)

  res$regulation <- "Not DEG"
  res$regulation[
    res$padj < 0.05 & res$log2FoldChange > 0.5
  ] <- "Up-regulated"
  res$regulation[
    res$padj < 0.05 & res$log2FoldChange < -0.5
  ] <- "Down-regulated"

  non_deg_data <- res[res$regulation == "Not DEG", ]
  median_expr <- median(
    non_deg_data$WT_log2 + non_deg_data$DN_log2,
    na.rm = TRUE
  )

  res$non_deg_group <- "Not DEG - Low"
  res$non_deg_group[
    res$regulation == "Not DEG" &
      (res$WT_log2 + res$DN_log2) > median_expr
  ] <- "Not DEG - High"
  res$non_deg_group[
    res$regulation != "Not DEG"
  ] <- res$regulation[res$regulation != "Not DEG"]

  DEG_list <- res %>%
    filter(regulation %in% c("Up-regulated", "Down-regulated")) %>%
    arrange(padj)

  write.csv(
    res,
    file.path(table_dir, "SuppFig8C_DESeq2_all_results.csv"),
    row.names = TRUE
  )
  write.csv(
    DEG_list %>% filter(regulation == "Up-regulated"),
    file.path(table_dir, "SuppFig8C_upregulated_genes.csv"),
    row.names = TRUE
  )
  write.csv(
    DEG_list %>% filter(regulation == "Down-regulated"),
    file.path(table_dir, "SuppFig8C_downregulated_genes.csv"),
    row.names = TRUE
  )
}

# Supplementary Figure 8C: WT-versus-DN scatterplot ----------------------
{
  fig8c_plot <- ggplot(res, aes(x = DN_log2, y = WT_log2)) +
    geom_point(
      data = subset(
        res,
        regulation == "Not DEG" & non_deg_group == "Not DEG - Low"
      ),
      color = "gray70",
      alpha = 0.6,
      size = 1
    ) +
    geom_point(
      data = subset(
        res,
        regulation == "Not DEG" & non_deg_group == "Not DEG - High"
      ),
      color = "gray70",
      alpha = 0.6,
      size = 1
    ) +
    geom_point(
      data = subset(res, regulation == "Up-regulated"),
      color = "#85312B",
      alpha = 0.8,
      size = 1.5
    ) +
    geom_point(
      data = subset(res, regulation == "Down-regulated"),
      color = "#373671",
      alpha = 0.8,
      size = 1.5
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "blue",
      alpha = 0.5
    ) +
    theme_minimal() +
    labs(
      title = "RNA-seq (WT vs DN)",
      x = "DN (Log2 values)",
      y = "WT (Log2 values)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20)
    ) +
    coord_fixed(ratio = 1)

  ggsave(
    file.path(fig8c_dir, "SuppFig8C_AFF3_WT_vs_DN_scatterplot.tiff"),
    plot = fig8c_plot,
    device = "tiff",
    dpi = 900,
    width = 9,
    height = 9,
    units = "in",
    compression = "lzw"
  )
}

# Supplementary Figure 8D: GO analysis of up- and downregulated genes ----
{
  DEG_up <- DEG_list %>%
    filter(non_deg_group == "Up-regulated")
  DEG_down <- DEG_list %>%
    filter(non_deg_group == "Down-regulated")

  up_entrez <- bitr(
    DEG_up$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  down_entrez <- bitr(
    DEG_down$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  go_up <- enrichGO(
    up_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  go_down <- enrichGO(
    down_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  go_up_df <- as.data.frame(go_up)
  go_down_df <- as.data.frame(go_down)

  write.csv(
    go_up_df,
    file.path(table_dir, "SuppFig8D_GO_up_all_terms.csv"),
    row.names = FALSE
  )
  write.csv(
    go_down_df,
    file.path(table_dir, "SuppFig8D_GO_down_all_terms.csv"),
    row.names = FALSE
  )

  go_up_bp <- go_up_df %>%
    filter(ONTOLOGY == "BP") %>%
    mutate(
      neglog10_padj = -log10(p.adjust),
      Direction = "Up in DN"
    )

  go_down_bp <- go_down_df %>%
    filter(ONTOLOGY == "BP") %>%
    mutate(
      neglog10_padj = -log10(p.adjust),
      Direction = "Down in DN"
    )

  fig8d_up_plot <- go_up_bp %>%
    slice_min(order_by = p.adjust, n = 15) %>%
    ggplot(aes(
      x = neglog10_padj,
      y = reorder(Description, neglog10_padj),
      fill = FoldEnrichment
    )) +
    geom_col(color = "black") +
    labs(
      x = "-log10 (p.adj)",
      y = NULL,
      fill = "Fold enrichment",
      title = "GO BP enriched in upregulated genes"
    ) +
    scale_fill_gradient(
      low = "#EFE6D8",
      high = "#85312B"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )

  fig8d_down_plot <- go_down_bp %>%
    slice_min(order_by = p.adjust, n = 15) %>%
    ggplot(aes(
      x = neglog10_padj,
      y = reorder(
        stringr::str_wrap(Description, width = 40),
        neglog10_padj
      ),
      fill = FoldEnrichment
    )) +
    geom_col(color = "black") +
    labs(
      x = "-log10 (p.adj)",
      y = NULL,
      fill = "Fold enrichment",
      title = "GO BP enriched in downregulated genes"
    ) +
    scale_fill_gradient(
      low = "#EFE6D8",
      high = "#373671"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )

  ggsave(
    file.path(fig8d_dir, "SuppFig8D_GO_up.tiff"),
    plot = fig8d_up_plot,
    device = "tiff",
    dpi = 900,
    width = 10,
    height = 5,
    units = "in",
    compression = "lzw"
  )
  ggsave(
    file.path(fig8d_dir, "SuppFig8D_GO_down.tiff"),
    plot = fig8d_down_plot,
    device = "tiff",
    dpi = 900,
    width = 10,
    height = 5,
    units = "in",
    compression = "lzw"
  )
}

# Supplementary Figure 8E: DN mutant 1 versus WT VST matrix --------------
{
  counts_Data_DN1 <- data.frame(
    AFF3_WT_1_1 = WT_1_1_counts$counts,
    AFF3_WT_1_2 = WT_1_2_counts$counts,
    AFF3_DN_1_1 = DN_1_1_counts$counts,
    AFF3_DN_1_2 = DN_1_2_counts$counts
  )
  rownames(counts_Data_DN1) <- WT_1_1_counts$Geneid

  col_Data_DN1 <- data.frame(
    condition = c("WT", "WT", "DN", "DN")
  )
  rownames(col_Data_DN1) <- colnames(counts_Data_DN1)

  dds_DN1 <- DESeqDataSetFromMatrix(
    countData = counts_Data_DN1,
    colData = col_Data_DN1,
    design = ~ condition
  )
  dds_DN1$condition <- relevel(dds_DN1$condition, ref = "WT")
  dds_DN1 <- DESeq(dds_DN1)

  dds_vst_DN1 <- vst(dds_DN1, blind = FALSE)
  vst_counts_DN1 <- assay(dds_vst_DN1)
}

# Marker lists used in Supplementary Figure 8E --------------------------
{
  gene_list_neural_progenitor <- gsub(" ", "", c(
    "SOX1", "SOX2", "PAX6", "NES",
    "HES1", "HES5",
    "VIM",
    "FABP7",
    "ID1", "ID3"
  ))

  gene_list_mature_neuron_core <- gsub(" ", "", c(
    "RBFOX3",
    "MAP2",
    "NEFL", "NEFM", "NEFH",
    "MAPT",
    "STMN2",
    "SYP",
    "SYN1", "SYN2",
    "SNAP25",
    "SYT1",
    "VAMP2",
    "RAB3A"
  ))

  gene_list_excitatory_neuron <- gsub(" ", "", c(
    "SLC17A7",
    "GLS",
    "GRIN1",
    "GRIA2",
    "CAMK2A",
    "SATB2",
    "BCL11B"
  ))

  gene_list_inhibitory_final <- c(
    "GAD1",
    "GAD2",
    "CALB1",
    "CNR1",
    "CXCL14",
    "TAC3",
    "PVALB",
    "VIP",
    "LAMP5",
    "NDNF"
  )

  stopifnot(all(c(
    gene_list_neural_progenitor,
    gene_list_mature_neuron_core,
    gene_list_excitatory_neuron,
    gene_list_inhibitory_final
  ) %in% rownames(vst_counts_DN1)))
}

# Row-z-scored expression matrices ---------------------------------------
{
  counts_neural_progenitor <- vst_counts_DN1[
    gene_list_neural_progenitor,
  ]
  z_matrix_neural_progenitor <- t(scale(t(counts_neural_progenitor)))
  z_matrix_neural_progenitor_cleaned <-
    z_matrix_neural_progenitor[
      !apply(z_matrix_neural_progenitor, 1, function(x) {
        any(is.na(x) | is.infinite(x))
      }),
    ]
  z_matrix_neural_progenitor_cleaned <-
    z_matrix_neural_progenitor_cleaned[
      ,
      !apply(z_matrix_neural_progenitor_cleaned, 2, function(x) {
        any(is.na(x) | is.infinite(x))
      })
    ]

  counts_mature_neuron_core <- vst_counts_DN1[
    gene_list_mature_neuron_core,
  ]
  z_matrix_mature_neuron_core <- t(scale(t(counts_mature_neuron_core)))
  z_matrix_mature_neuron_core_cleaned <-
    z_matrix_mature_neuron_core[
      !apply(z_matrix_mature_neuron_core, 1, function(x) {
        any(is.na(x) | is.infinite(x))
      }),
    ]
  z_matrix_mature_neuron_core_cleaned <-
    z_matrix_mature_neuron_core_cleaned[
      ,
      !apply(z_matrix_mature_neuron_core_cleaned, 2, function(x) {
        any(is.na(x) | is.infinite(x))
      })
    ]

  counts_excitatory_neuron <- vst_counts_DN1[
    gene_list_excitatory_neuron,
  ]
  z_matrix_excitatory_neuron <- t(scale(t(counts_excitatory_neuron)))

  counts_inhibitory_final <- vst_counts_DN1[
    gene_list_inhibitory_final,
  ]
  z_matrix_inhibitory_final <- t(scale(t(counts_inhibitory_final)))
  z_matrix_inhibitory_final_cleaned <-
    z_matrix_inhibitory_final[
      !apply(z_matrix_inhibitory_final, 1, function(x) {
        any(is.na(x) | is.infinite(x))
      }),
    ]
  z_matrix_inhibitory_final_cleaned <-
    z_matrix_inhibitory_final_cleaned[
      ,
      !apply(z_matrix_inhibitory_final_cleaned, 2, function(x) {
        any(is.na(x) | is.infinite(x))
      })
    ]

  write.csv(
    z_matrix_neural_progenitor_cleaned,
    file.path(table_dir, "SuppFig8E_neural_progenitor_z_scores.csv"),
    row.names = TRUE
  )
  write.csv(
    z_matrix_mature_neuron_core_cleaned,
    file.path(table_dir, "SuppFig8E_mature_neuron_z_scores.csv"),
    row.names = TRUE
  )
  write.csv(
    z_matrix_excitatory_neuron,
    file.path(table_dir, "SuppFig8E_excitatory_neuron_z_scores.csv"),
    row.names = TRUE
  )
  write.csv(
    z_matrix_inhibitory_final_cleaned,
    file.path(table_dir, "SuppFig8E_inhibitory_neuron_z_scores.csv"),
    row.names = TRUE
  )
}

# Supplementary Figure 8E heatmaps ---------------------------------------
{
  cols_new <- rev(RColorBrewer::brewer.pal(11, "PuOr"))

  plot_lineage_heatmap <- function(z_matrix, file_name) {
    tiff(
      file.path(fig8e_dir, file_name),
      width = 7,
      height = 14,
      units = "in",
      res = 300
    )

    pheatmap(
      z_matrix,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      color = cols_new,
      main = "",
      fontsize = 12,
      border_color = "grey",
      scale = "none",
      cellwidth = 20,
      cellheight = 20,
      legend = TRUE,
      legend.title = "Normalized expression",
      legend.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    )

    dev.off()
  }

  plot_lineage_heatmap(
    z_matrix_neural_progenitor_cleaned,
    "SuppFig8E_neural_progenitor.tiff"
  )
  plot_lineage_heatmap(
    z_matrix_mature_neuron_core_cleaned,
    "SuppFig8E_mature_neuron.tiff"
  )
  plot_lineage_heatmap(
    z_matrix_excitatory_neuron,
    "SuppFig8E_excitatory_neuron.tiff"
  )
  plot_lineage_heatmap(
    z_matrix_inhibitory_final_cleaned,
    "SuppFig8E_inhibitory_neuron.tiff"
  )
}
