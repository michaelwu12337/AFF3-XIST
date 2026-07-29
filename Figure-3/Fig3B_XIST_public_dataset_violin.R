{
  library(Seurat)
  library(anndataR)
  library(Matrix)
  library(ggplot2)
}
{
  ad <- read_h5ad("source_files/adata_no_W51.h5ad")
}
{
  as_dgC <- function(m) {
    if (inherits(m, "dgCMatrix")) return(m)
    if (is.matrix(m)) return(Matrix(m, sparse = TRUE))
    Matrix(as.matrix(m), sparse = TRUE)
  }
  
  X <- ad$X                      
  counts <- NULL
  if (!is.null(ad$layers) && "counts" %in% names(ad$layers)) {
    counts <- ad$layers[["counts"]]  
  }
  
  genes <- rownames(ad$var)
  cells <- rownames(ad$obs)
  
  data_mat <- t(as_dgC(X))
  colnames(data_mat) <- cells
  rownames(data_mat) <- genes
  
  if (!is.null(counts)) {
    counts_mat <- t(as_dgC(counts))
    colnames(counts_mat) <- cells
    rownames(counts_mat) <- genes
    
    obj <- CreateSeuratObject(counts = counts_mat, meta.data = as.data.frame(ad$obs))
    obj <- SetAssayData(obj, assay = "RNA", slot = "data", new.data = data_mat)
  } else {
    obj <- CreateSeuratObject(counts = data_mat, meta.data = as.data.frame(ad$obs))
  }
  
  obj
}
{
  p <- VlnPlot(
    obj,
    features = "XIST",
    group.by = "week_stage",
    slot = "counts",
    pt.size = 0,
    adjust = 0.55,     # thinner violin
    raster = TRUE
  ) +
    NoLegend() +
    theme_classic(base_size = 14) +
    theme(
      plot.title   = element_blank(),
      axis.title.x = element_text(size = 17),
      axis.title.y = element_text(size = 17),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      axis.line    = element_line(linewidth = 0.6),
      axis.ticks   = element_line(linewidth = 0.6),
      axis.ticks.length = unit(5, "pt"),
      panel.grid   = element_blank()
    ) +
    labs(title = "XIST", 
         x = "Week Stage", 
         y = "XIST Raw Counts")+
    guides(fill = "none", color = "none")
}

ggsave(
  filename = "Fig3B_XIST_public_dataset_violin.svg",
  plot     = p,
  device   = svglite::svglite,
  width    = 3,
  height   = 4.5,
  units    = "in"
)
