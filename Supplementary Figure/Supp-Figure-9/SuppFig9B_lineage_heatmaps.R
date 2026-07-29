{
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(scales)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "SuppFig9B")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Read the old and new WT/KO count tables --------------------------------
{
  old_WT1_counts <- read.table(
    file.path(source_dir, "AFF3_F_WT_1_mm_counts_revised.txt"),
    header = TRUE
  )
  old_WT2_counts <- read.table(
    file.path(source_dir, "AFF3_F_WT_2_mm_counts_revised.txt"),
    header = TRUE
  )
  old_KO1_counts <- read.table(
    file.path(source_dir, "AFF3_F_KO_1_mm_counts_revised.txt"),
    header = TRUE
  )
  old_KO2_counts <- read.table(
    file.path(source_dir, "AFF3_F_KO_2_mm_counts_revised.txt"),
    header = TRUE
  )

  new_WT1_counts <- read.table(
    file.path(source_dir, "WT_1_mm_counts_revised.txt"),
    header = TRUE
  )
  new_WT2_counts <- read.table(
    file.path(source_dir, "WT_2_mm_counts_revised.txt"),
    header = TRUE
  )
  new_KO1_counts <- read.table(
    file.path(source_dir, "KO_1_mm_counts_revised.txt"),
    header = TRUE
  )
  new_KO2_counts <- read.table(
    file.path(source_dir, "KO_2_mm_counts_revised.txt"),
    header = TRUE
  )
}

# DESeq2 VST transformation ----------------------------------------------
{
  make_WT_KO_vst <- function(WT1_counts, WT2_counts, KO1_counts, KO2_counts) {
    stopifnot(
      identical(WT1_counts$Geneid, WT2_counts$Geneid),
      identical(WT1_counts$Geneid, KO1_counts$Geneid),
      identical(WT1_counts$Geneid, KO2_counts$Geneid)
    )

    counts_Data <- data.frame(
      AFF3_F_WT_1 = WT1_counts$counts,
      AFF3_F_WT_2 = WT2_counts$counts,
      AFF3_F_KO_1 = KO1_counts$counts,
      AFF3_F_KO_2 = KO2_counts$counts
    )
    rownames(counts_Data) <- WT1_counts$Geneid

    col_Data <- data.frame(
      condition = factor(
        c("WT", "WT", "KO", "KO"),
        levels = c("WT", "KO")
      )
    )
    rownames(col_Data) <- colnames(counts_Data)

    dds <- DESeqDataSetFromMatrix(
      countData = counts_Data,
      colData = col_Data,
      design = ~ condition
    )
    dds <- DESeq(dds)
    dds_vst <- vst(dds, blind = FALSE)
    assay(dds_vst)
  }

  old_vst_counts <- make_WT_KO_vst(
    old_WT1_counts,
    old_WT2_counts,
    old_KO1_counts,
    old_KO2_counts
  )
  new_vst_counts <- make_WT_KO_vst(
    new_WT1_counts,
    new_WT2_counts,
    new_KO1_counts,
    new_KO2_counts
  )
}

# Lineage-marker lists used in the v14 figure ----------------------------
{
  neural_genes <- c(
    "PAX3", "PAX7", "OTX2", "PAX6", "EMX2",
    "SRGAP3", "NAV1", "NRXN3", "EPHB1", "EPHB2"
  )
  hematopoietic_genes <- c(
    "RUNX1", "ERG", "HOXA10", "KDR", "IGFBP4",
    "ETV6", "CLEC11A", "GYPA", "LMO2", "HBZ"
  )
  cardiac_genes <- c(
    "NKX2-5", "HAND1", "HAND2", "TBX1", "ISL1",
    "NFATC1", "WT1", "THY1", "VIM", "TNNT2"
  )
  germ_layer_genes <- c(
    "SOX17", "FOXA2", "TFAP2A", "GATA3", "MESP1",
    "EOMES", "CDH11", "MDK", "COL3A1", "LUM"
  )
}

# Lineage heatmap function ------------------------------------------------
{
  cols_new <- c(
    "#2D004BFF",
    "#542788FF",
    "#8073ACFF",
    "#B2ABD2FF",
    "#D8DAEBFF",
    "#F7F7F7FF",
    "#FBE3E3",
    "#F8B8B0",
    "#EF7260",
    "#CC2E2E",
    "#7F0000"
  )

  plot_lineage_heatmap <- function(
      vst_counts,
      gene_list,
      lower_inner_stop,
      upper_inner_stop,
      output_name) {
    missing_genes <- setdiff(gene_list, rownames(vst_counts))
    if (length(missing_genes) > 0) {
      stop(
        paste(
          "Missing lineage-marker genes:",
          paste(missing_genes, collapse = ", ")
        )
      )
    }

    counts_filtered <- vst_counts[gene_list, ]
    z_matrix_filtered <- t(scale(t(counts_filtered)))

    df_long <- as.data.frame(z_matrix_filtered) %>%
      rownames_to_column("Gene") %>%
      pivot_longer(
        cols = -Gene,
        names_to = "Sample",
        values_to = "Zscore"
      )
    df_long$Sample <- factor(
      df_long$Sample,
      levels = c(
        "AFF3_F_WT_1",
        "AFF3_F_WT_2",
        "AFF3_F_KO_1",
        "AFF3_F_KO_2"
      )
    )
    df_long$Gene <- factor(
      df_long$Gene,
      levels = rev(unique(df_long$Gene))
    )

    lim <- max(abs(df_long$Zscore), na.rm = TRUE)
    stops <- c(
      -lim,
      lower_inner_stop[1],
      lower_inner_stop[2],
      0,
      upper_inner_stop[1],
      upper_inner_stop[2],
      lim
    )
    vals <- rescale(stops, to = c(0, 1))
    new_labels <- c("WT_1", "WT_2", "KO_1", "KO_2")

    heatmap_plot <- ggplot(
      df_long,
      aes(Sample, Gene, fill = Zscore)
    ) +
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
        labels = rev(gene_list),
        position = "right",
        expand = c(0, 0)
      ) +
      scale_x_discrete(
        expand = c(0, 0),
        labels = new_labels
      ) +
      coord_fixed(ratio = 1) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.text.y.right = element_text(
          size = 32,
          margin = margin(l = 2.5),
          hjust = 0,
          face = "bold.italic"
        ),
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
          face = "bold"
        ),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        legend.position = "bottom",
        legend.box.margin = margin(t = 13),
        legend.justification = "left",
        legend.box.just = "left"
      ) +
      labs(x = NULL, y = NULL)

    ggsave(
      file.path(figure_dir, paste0(output_name, ".tiff")),
      plot = heatmap_plot,
      width = 4.05,
      height = 12.5,
      units = "in",
      dpi = 600,
      limitsize = FALSE,
      compression = "lzw"
    )
  }
}

# Supplementary Figure 9B: old experiment heatmaps -----------------------
{
  plot_lineage_heatmap(
    old_vst_counts,
    neural_genes,
    lower_inner_stop = c(-0.91, -0.792),
    upper_inner_stop = c(0.792, 0.91),
    output_name = "SuppFig9B_old_neural"
  )
  plot_lineage_heatmap(
    old_vst_counts,
    hematopoietic_genes,
    lower_inner_stop = c(-1.19, -0.885),
    upper_inner_stop = c(0.885, 1.19),
    output_name = "SuppFig9B_old_hematopoietic"
  )
  plot_lineage_heatmap(
    old_vst_counts,
    cardiac_genes,
    lower_inner_stop = c(-1.06, -0.837),
    upper_inner_stop = c(0.837, 1.06),
    output_name = "SuppFig9B_old_cardiac"
  )
  plot_lineage_heatmap(
    old_vst_counts,
    germ_layer_genes,
    lower_inner_stop = c(-1.19, -0.812),
    upper_inner_stop = c(0.812, 1.19),
    output_name = "SuppFig9B_old_germ_layer"
  )
}

# Supplementary Figure 9B: new experiment heatmaps -----------------------
{
  plot_lineage_heatmap(
    new_vst_counts,
    neural_genes,
    lower_inner_stop = c(-0.91, -0.792),
    upper_inner_stop = c(0.792, 0.91),
    output_name = "SuppFig9B_new_neural"
  )
  plot_lineage_heatmap(
    new_vst_counts,
    hematopoietic_genes,
    lower_inner_stop = c(-1.19, -0.885),
    upper_inner_stop = c(0.885, 1.19),
    output_name = "SuppFig9B_new_hematopoietic"
  )
  plot_lineage_heatmap(
    new_vst_counts,
    cardiac_genes,
    lower_inner_stop = c(-1.06, -0.837),
    upper_inner_stop = c(0.837, 1.06),
    output_name = "SuppFig9B_new_cardiac"
  )
  plot_lineage_heatmap(
    new_vst_counts,
    germ_layer_genes,
    lower_inner_stop = c(-1.19, -0.812),
    upper_inner_stop = c(0.812, 1.19),
    output_name = "SuppFig9B_new_germ_layer"
  )
}
