{
  library(VennDiagram)
  library(grid)
  library(svglite)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  table_dir <- "generated_tables"
  fig5e_dir <- file.path("generated_figures", "Fig5E")
  fig5f_dir <- file.path("generated_figures", "Fig5F")

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig5e_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig5f_dir, recursive = TRUE, showWarnings = FALSE)
}

# Gene sets used in the v46 Venn diagram ---------------------------------
{
  set_xist <- read.csv(
    file.path(source_dir, "XIST_binding_933_genes.csv"),
    check.names = FALSE
  )$Gene

  set_dko <- read.csv(
    file.path(source_dir, "BO_DKO_rescued_genes_1599.csv"),
    check.names = FALSE
  )$x

  set_xist <- set_xist |>
    as.character() |>
    trimws() |>
    toupper() |>
    unique()

  set_dko <- set_dko |>
    as.character() |>
    trimws() |>
    toupper() |>
    unique()
}

# Figure 5E: DKO-rescued genes and autosomal XIST targets ----------------
{
  venn_plot <- venn.diagram(
    x = list(
      XIST_binding = set_xist,
      genes_binding = set_dko
    ),
    category.names = c("XIST binding", "Genes binding"),
    filename = NULL,
    output = TRUE,
    fill = c("#FFCDBA", "#FFC566"),
    alpha = 0.5,
    lwd = 2,
    lty = "solid",
    cex = 1.8,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 2.4,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    cat.pos = c(-40, 40),
    cat.dist = c(0.08, 0.08),
    print.mode = "raw",
    margin = 0.1
  )

  svglite(
    file.path(fig5e_dir, "Fig5E_DKO_rescued_XIST_target_Venn.svg"),
    bg = "transparent",
    width = 8,
    height = 6
  )
  grid.newpage()
  pushViewport(
    viewport(
      width = unit(1, "npc"),
      height = unit(0.7, "npc")
    )
  )
  grid.draw(venn_plot)
  popViewport()
  dev.off()

  XIST_only <- setdiff(set_xist, set_dko)
  DKO_only <- setdiff(set_dko, set_xist)
  overlap <- intersect(set_xist, set_dko)

  max_len <- max(length(XIST_only), length(DKO_only), length(overlap))
  venn_genes <- data.frame(
    XIST_only = c(XIST_only, rep(NA, max_len - length(XIST_only))),
    DKO_only = c(DKO_only, rep(NA, max_len - length(DKO_only))),
    overlap = c(overlap, rep(NA, max_len - length(overlap)))
  )

  venn_summary <- data.frame(
    gene_set = c(
      "Autosomal XIST targets",
      "DKO rescued genes",
      "Shared genes",
      "XIST targets only",
      "DKO rescued only"
    ),
    gene_count = c(
      length(set_xist),
      length(set_dko),
      length(overlap),
      length(XIST_only),
      length(DKO_only)
    )
  )

  write.csv(
    venn_genes,
    file.path(table_dir, "Fig5E_Venn_gene_sets.csv"),
    row.names = FALSE
  )
  write.csv(
    venn_summary,
    file.path(table_dir, "Fig5E_Venn_summary.csv"),
    row.names = FALSE
  )
}

# Figure 5F: GO analysis of the 44 overlapping genes ---------------------
{
  GO_overlap <- enrichGO(
    gene = overlap,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pAdjustMethod = "BH"
  )

  GO_overlap_df <- as.data.frame(GO_overlap) |>
    mutate(
      Description_raw = Description,
      FoldEnrichment = as.numeric(FoldEnrichment),
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue)
    )

  patterns_keep <- c(
    "cell fate specification",
    "motor neuron apoptotic",
    "skeletal muscle organ development",
    "synaptic membrane",
    "pore complex"
  )

  GO_overlap_df_sel <- GO_overlap_df |>
    filter(
      stringr::str_detect(
        Description_raw,
        paste(patterns_keep, collapse = "|")
      )
    ) |>
    mutate(
      Description = stringr::str_wrap(Description_raw, width = 40)
    )

  fig5f_plot <- ggplot(
    GO_overlap_df_sel |>
      slice_max(order_by = neglog10p, n = 5),
    aes(
      x = neglog10p,
      y = reorder(Description, neglog10p),
      fill = log2_OE
    )
  ) +
    geom_col(
      color = "black",
      width = 0.8,
      linewidth = 1.4
    ) +
    labs(
      x = "-log10 (p-value)",
      y = NULL,
      fill = "log2 (O/E)"
    ) +
    scale_fill_gradient(
      low = "#E6EDF2",
      high = "#4C6A92",
      labels = scales::number_format(accuracy = 0.1),
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    theme_classic() +
    theme(
      axis.line = element_line(linewidth = 1.4, colour = "black"),
      axis.text.y = element_text(size = 30, face = "bold"),
      axis.title.x = element_text(size = 30, face = "bold"),
      axis.text.x = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 30),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(0.8, "cm")
    )

  ggsave(
    file.path(fig5f_dir, "Fig5F_DKO_rescued_XIST_target_GO.tiff"),
    plot = fig5f_plot,
    width = 15,
    height = 6,
    units = "in",
    dpi = 600,
    limitsize = FALSE,
    compression = "lzw"
  )

  write.csv(
    GO_overlap_df,
    file.path(table_dir, "Fig5F_GO_all_terms.csv"),
    row.names = FALSE
  )
  write.csv(
    GO_overlap_df_sel,
    file.path(table_dir, "Fig5F_GO_selected_terms.csv"),
    row.names = FALSE
  )
}
