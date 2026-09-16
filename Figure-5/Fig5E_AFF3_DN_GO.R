{
  library(DESeq2)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "Fig5E")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the WT and dominant-negative count tables -------------------------
{
  WT_1_counts <- read.table(
    file.path(source_dir, "Fig5D_WT_1_mm_counts_revised.txt"),
    header = TRUE
  )
  WT_2_counts <- read.table(
    file.path(source_dir, "Fig5D_WT_2_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_1_1_counts <- read.table(
    file.path(source_dir, "Fig5D_DN_1_1_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_1_2_counts <- read.table(
    file.path(source_dir, "Fig5D_DN_1_2_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_2_1_counts <- read.table(
    file.path(source_dir, "Fig5D_DN_2_1_mm_counts_revised.txt"),
    header = TRUE
  )
  DN_2_2_counts <- read.table(
    file.path(source_dir, "Fig5D_DN_2_2_mm_counts_revised.txt"),
    header = TRUE
  )

  count_tables <- list(
    WT_1_counts,
    WT_2_counts,
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

# DESeq2 analysis used for the WT-versus-DN comparison ------------------
{
  counts_data <- data.frame(
    AFF3_WT_1 = WT_1_counts$counts,
    AFF3_WT_2 = WT_2_counts$counts,
    AFF3_DN_1_1 = DN_1_1_counts$counts,
    AFF3_DN_1_2 = DN_1_2_counts$counts,
    AFF3_DN_2_1 = DN_2_1_counts$counts,
    AFF3_DN_2_2 = DN_2_2_counts$counts
  )
  rownames(counts_data) <- WT_1_counts$Geneid

  col_data <- data.frame(
    condition = c("WT", "WT", "DN", "DN", "DN", "DN")
  )
  rownames(col_data) <- colnames(counts_data)

  dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", "DN", "WT"))
  res <- as.data.frame(res)
  res$gene <- rownames(res)

  res$regulation <- "Not DEG"
  res$regulation[
    !is.na(res$padj) &
      res$padj < 0.05 &
      res$log2FoldChange > 0.5
  ] <- "Up-regulated"
  res$regulation[
    !is.na(res$padj) &
      res$padj < 0.05 &
      res$log2FoldChange < -0.5
  ] <- "Down-regulated"
}

# GO analysis follows the Supplementary Figure 8 workflow ---------------
{
  DEG_up <- res %>%
    filter(regulation == "Up-regulated")
  DEG_down <- res %>%
    filter(regulation == "Down-regulated")

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

  go_up_bp <- as.data.frame(go_up) %>%
    filter(ONTOLOGY == "BP") %>%
    mutate(Direction = "Upregulated in DN")
  go_down_bp <- as.data.frame(go_down) %>%
    filter(ONTOLOGY == "BP") %>%
    mutate(Direction = "Downregulated in DN")
}

# Five manually selected terms from each direction ----------------------
{
  selected_up_terms <- c(
    "regulation of neural precursor cell proliferation",
    "forebrain regionalization",
    "axonogenesis",
    "visual system development",
    "sensory system development",
    "eye development",
    "telencephalon regionalization",
    "regulation of neurogenesis",
    "regulation of neuroblast proliferation"
  )
  selected_down_terms <- c(
    "response to BMP",
    "cellular response to BMP stimulus",
    "endoderm development",
    "cell fate commitment",
    "formation of primary germ layer",
    "gastrulation",
    "mesoderm development",
    "mesoderm morphogenesis",
    "anterior/posterior axis specification"
  )

  selected_up <- go_up_bp %>%
    filter(Description %in% selected_up_terms)
  selected_down <- go_down_bp %>%
    filter(Description %in% selected_down_terms)

  stopifnot(setequal(selected_up$Description, selected_up_terms))
  stopifnot(setequal(selected_down$Description, selected_down_terms))

  selected_go <- bind_rows(selected_up, selected_down) %>%
    mutate(
      signed_neglog10_padj = if_else(
        Direction == "Upregulated in DN",
        -log10(p.adjust),
        log10(p.adjust)
      )
    )

  negative_axis_limit <- ceiling(
    max(-selected_go$signed_neglog10_padj) / 5
  ) * 5
  positive_axis_limit <- 6

  selected_go <- selected_go %>%
    mutate(
      display_score = if_else(
        signed_neglog10_padj < 0,
        signed_neglog10_padj / negative_axis_limit,
        signed_neglog10_padj / positive_axis_limit
      )
    )

  term_order <- selected_go %>%
    arrange(signed_neglog10_padj) %>%
    pull(Description)
  selected_go$Description <- factor(
    selected_go$Description,
    levels = term_order
  )
}

# Figure 5D: compact mirrored GO plot -----------------------------------
{
  fig6d_plot <- ggplot(
    selected_go,
    aes(
      x = display_score,
      y = Description
    )
  ) +
    geom_col(
      data = filter(selected_go, Direction == "Upregulated in DN"),
      aes(fill = FoldEnrichment),
      width = 0.80,
      color = "black",
      linewidth = 0.25
    ) +
    scale_fill_gradient(
      name = "Fold enrichment (up)",
      low = "#EFE6D8",
      high = "#85312B",
      guide = guide_colorbar(order = 1)
    ) +
    new_scale_fill() +
    geom_col(
      data = filter(selected_go, Direction == "Downregulated in DN"),
      aes(fill = FoldEnrichment),
      width = 0.80,
      color = "black",
      linewidth = 0.25
    ) +
    scale_fill_gradient(
      name = "Fold enrichment (down)",
      low = "#EFE6D8",
      high = "#373671",
      guide = guide_colorbar(order = 2)
    ) +
    geom_vline(
      xintercept = 0,
      color = "black",
      linewidth = 0.35
    ) +
    geom_text(
      data = filter(selected_go, Direction == "Upregulated in DN"),
      aes(x = -0.015, y = Description, label = Description),
      inherit.aes = FALSE,
      hjust = 1,
      size = 5,
      fontface = "bold",
      show.legend = FALSE
    ) +
    geom_text(
      data = filter(selected_go, Direction == "Downregulated in DN"),
      aes(x = 0.015, y = Description, label = Description),
      inherit.aes = FALSE,
      hjust = 0,
      size = 5,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_x_continuous(
      name = "Signed -log10(adjusted P)",
      limits = c(-1, 1),
      breaks = c(-1, -2 / 3, -1 / 3, 0, 1 / 3, 2 / 3, 1),
      labels = c(
        -negative_axis_limit,
        -2 * negative_axis_limit / 3,
        -negative_axis_limit / 3,
        0,
        positive_axis_limit / 3,
        2 * positive_axis_limit / 3,
        positive_axis_limit
      ),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(y = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 8) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = grid::unit(22, "pt"),
      plot.margin = margin(4, 8, 4, 8)
    )

  ggsave(
    file.path(figure_dir, "Fig5E_AFF3_DN_GO.tiff"),
    plot = fig6d_plot,
    device = "tiff",
    dpi = 600,
    width = 10.8,
    height = 4.8,
    units = "in",
    compression = "lzw"
  )
}
