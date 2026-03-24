{
  setwd("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/scRNA/scRNA/notebook")
  library(SeuratDisk)
  library(Seurat)
  library(zellkonverter)
  library(anndataR)
  library(SingleCellExperiment)
  library(Matrix)
  library(ggplot2)
}

## 1) Path to h5ad
h5ad <- "~/Library/CloudStorage/OneDrive-Personal/HKU PhD/scRNA/scRNA/data/adata_no_W51.h5ad"

## 2) Read h5ad -> AnnData (R-native)
ad <- read_h5ad(h5ad)

## 3) Inspect structure (safe checks)
ad
ad$X
names(ad$layers)
head(ad$obs)
head(ad$var)
names(ad$obsm)

# helper: convert dense -> sparse if needed
{
  as_dgC <- function(m) {
    if (inherits(m, "dgCMatrix")) return(m)
    if (is.matrix(m)) return(Matrix(m, sparse = TRUE))
    Matrix(as.matrix(m), sparse = TRUE)
  }
  
  # 1) get cell x gene matrices from AnnData
  X <- ad$X                      # typically cells x genes
  counts <- NULL
  if (!is.null(ad$layers) && "counts" %in% names(ad$layers)) {
    counts <- ad$layers[["counts"]]   # typically cells x genes
  }
  
  # 2) get gene + cell names
  genes <- rownames(ad$var)
  cells <- rownames(ad$obs)
  
  # 3) coerce to sparse and transpose to genes x cells for Seurat
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
    # fall back: treat X as counts if no counts layer exists
    obj <- CreateSeuratObject(counts = data_mat, meta.data = as.data.frame(ad$obs))
  }
  
  obj
}


VlnPlot(
  obj,
  features = "nFeature_RNA",
  group.by = "week_stage",
  pt.size = 0
)

## Same as jupyter plot
VlnPlot(
  obj,
  features = "XIST",
  group.by = "week_stage",
  slot = "counts",
  pt.size = 0
)

## Modified
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
  
  p
}


ggsave(
  filename = "~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Figure 4 new/violin/Fig4A_violin.svg",
  plot     = p,
  device   = svglite::svglite,
  width    = 3,
  height   = 4.5,
  units    = "in"
)

