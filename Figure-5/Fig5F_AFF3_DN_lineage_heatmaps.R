{
  library(DESeq2)
  library(pheatmap)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "Fig5F")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the WT and dominant-negative mutant 1 count tables ----------------
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

  count_tables <- list(
    WT_1_counts,
    WT_2_counts,
    DN_1_1_counts,
    DN_1_2_counts
  )
  stopifnot(all(vapply(
    count_tables[-1],
    function(x) identical(x$Geneid, count_tables[[1]]$Geneid),
    logical(1)
  )))
}

# VST expression matrix used in the original lineage heatmaps ------------
{
  counts_data <- data.frame(
    AFF3_WT_1_1 = WT_1_counts$counts,
    AFF3_WT_1_2 = WT_2_counts$counts,
    AFF3_DN_1_1 = DN_1_1_counts$counts,
    AFF3_DN_1_2 = DN_1_2_counts$counts
  )
  rownames(counts_data) <- WT_1_counts$Geneid

  col_data <- data.frame(
    condition = factor(c("WT", "WT", "DN", "DN"), levels = c("WT", "DN"))
  )
  rownames(col_data) <- colnames(counts_data)

  dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq(dds)
  vst_counts <- assay(vst(dds, blind = FALSE))
}

# Preserve the existing clustered row order after removing selected genes -
{
  neural_progenitor_order <- c(
    "ID3", "PAX6", "HES5", "VIM", "HES1", "FABP7", "SOX1", "SOX2"
  )
  neural_progenitor_expanded_order <- c(
    neural_progenitor_order,
    "SOX3", "GLI3", "LHX2", "SOX9", "NOTCH1", "MSI1", "FGFBP3",
    "OTX2", "ZIC2"
  )
  excitatory_neuron_order <- c(
    "SATB2", "GRIN1", "GRIA2", "SLC17A7", "CAMK2A"
  )
  excitatory_neuron_expanded_order <- c(
    "SATB2",
    "GRIN1", "GRIA2", "GRIA1", "GRIN2B",
    "SHANK3",
    "SLC17A7", "CAMK2A"
  )
  inhibitory_neuron_order <- c(
    "GAD2", "GAD1", "TAC3", "CALB1", "CNR1", "PVALB"
  )
  inhibitory_neuron_expanded_order <- c(
    inhibitory_neuron_order,
    "ERBB4"
  )

  requested_genes <- unique(c(
    neural_progenitor_order,
    neural_progenitor_expanded_order,
    excitatory_neuron_order,
    excitatory_neuron_expanded_order,
    inhibitory_neuron_order,
    inhibitory_neuron_expanded_order
  ))
  stopifnot(all(requested_genes %in% rownames(vst_counts)))

  make_z_matrix <- function(gene_order) {
    expression_matrix <- vst_counts[gene_order, , drop = FALSE]
    z_matrix <- t(scale(t(expression_matrix)))
    stopifnot(!anyNA(z_matrix), all(is.finite(z_matrix)))
    z_matrix
  }

  z_neural_progenitor <- make_z_matrix(neural_progenitor_order)
  z_neural_progenitor_expanded <- make_z_matrix(
    neural_progenitor_expanded_order
  )
  z_excitatory_neuron <- make_z_matrix(excitatory_neuron_order)
  z_excitatory_neuron_expanded <- make_z_matrix(
    excitatory_neuron_expanded_order
  )
  z_inhibitory_neuron <- make_z_matrix(inhibitory_neuron_order)
  z_inhibitory_neuron_expanded <- make_z_matrix(
    inhibitory_neuron_expanded_order
  )
}

# Figure 5E: three separate lineage heatmaps ------------------------------
{
  heatmap_colors <- rev(RColorBrewer::brewer.pal(11, "PuOr"))

  plot_lineage_heatmap <- function(
    z_matrix,
    file_name,
    figure_height = 4.3
  ) {
    tiff(
      file.path(figure_dir, file_name),
      width = 3.2,
      height = figure_height,
      units = "in",
      res = 600,
      compression = "lzw"
    )

    pheatmap(
      z_matrix,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      color = heatmap_colors,
      main = "",
      fontsize = 12,
      border_color = "grey",
      scale = "none",
      cellwidth = 20,
      cellheight = 20,
      legend = TRUE
    )

    dev.off()
  }

  plot_lineage_heatmap(
    z_neural_progenitor,
    "Fig5F_neural_progenitors.tiff"
  )
  plot_lineage_heatmap(
    z_neural_progenitor_expanded,
    "Fig5F_neural_progenitors_expanded.tiff",
    figure_height = 6.6
  )
  plot_lineage_heatmap(
    z_excitatory_neuron,
    "Fig5F_excitatory_neurons.tiff"
  )
  plot_lineage_heatmap(
    z_excitatory_neuron_expanded,
    "Fig5F_excitatory_neurons_expanded.tiff"
  )
  plot_lineage_heatmap(
    z_inhibitory_neuron,
    "Fig5F_inhibitory_neurons.tiff"
  )
  plot_lineage_heatmap(
    z_inhibitory_neuron_expanded,
    "Fig5F_inhibitory_neurons_expanded.tiff"
  )
}
