# Dotplot and Diverging Barplot for Bulk RNA-seq Dataset
source_dir <- "source_files"
output_dir <- file.path("generated_figures", "Fig2F")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load required packages
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(extrafont)
library(grid)      
library(patchwork) 


# 1. Define gene sets as NAMED list ----------------------------------------

gene_sets <- list(
  Skeletal_Muscle    = c("PAX3", "PAX7", "RBM24", "GAS1", "BOC", "CDON", "EPHB1", "ACTN3", "TNNT1", "SGCB", "DMD"),
  Amnion             = c("EZR", "KRT8", "SAT1", "KRT19", "GSTO1", "HPGD", "CDH1", "EPCAM", "CLDN4", "TMEM54", "GATA3", "KRT18", "BMP4", "DAPK1", "SMAGP", "TFAP2A", "PATJ", "CLDN10", "AHNAK", "STOM", "FOXP1", "RARRES2", "TPM1", "GRHL2", "SHANK2", "MIR205HG", "HAPLN1", "PLCE1", "BAMBI", "NAV2", "DSP", "UTRN", "ISL1", "ARHGAP18", "PKP2", "SLC7A3", "FTL", "DSC2", "HMGA2", "FMR1", "ATP6V0D2", "ERRFI1", "TDRP", "ID2", "ERBB3", "AMER2", "CLDN7", "SLC7A2", "CD55", "YTHDC1"),
  Blood_Progenitors  = c("LYL1", "IFITM2", "ARHGDIB", "EMP3", "TAGLN2", "GYPC", "MYL4", "USP15", "RUNX1", "GYPB", "GATA1", "B2M", "S100A6", "ANP32B", "TESC", "S100A4", "RGS10", "NFE2", "CCND3", "CLIC1", "BMP2K", "GMFG", "MIR4435-2HG", "SPN", "GYPE", "LMO2", "TLN1", "CYTOR", "CHCHD2", "TIMP3", "CST3", "GPX1", "HBG2", "LAPTM5", "DUT", "RHAG", "MBNL1", "UROD", "TIMP1", "HBG1", "RAC2", "PFN1", "GMPR", "ELF1", "FERMT3", "HBE1", "PRKAR2B", "PSTPIP2", "ZFPM2", "KLF1"),
  Cardiac_Mesoderm   = c("S100A11", "RHOBTB3", "HAS2", "FSCN1", "HAND1", "B2M", "SPARC", "BAMBI", "COL5A2", "NTS", "GYPC", "FN1", "S100A10", "BMPER", "HAPLN1", "PPIB", "COL1A2", "PDGFRA", "SERPINE2", "ANXA6", "FAM89A", "TMSB10", "IL6ST", "ALDH2", "GATA6", "ENSG00000283674", "BMP4", "TIMP1", "SIPA1L2", "DOK4", "EPB41L2", "IFI16", "SERPINH1", "GSTO1", "CDH11", "TCEAL9", "PPFIBP1", "WLS", "IFITM3", "ITGB1", "LRFN5", "MMP2", "UNC5C", "H19", "FREM1", "LIFR", "CKAP4", "HSP90B1", "DNAJC15", "KDR"),
  Epi                = c("TERF1", "ESRG", "GRID2", "ENSG00000254277", "CD24", "ENSG00000288749", "POU5F1", "L1TD1", "DNMT3B", "ENSG00000274090", "PLAAT3", "JARID2", "PODXL", "ENSG00000288895", "GABRB3", "PHC1", "DPPA4", "UGP2", "PTMA", "TCEA1", "ENSG00000289894", "ENSG00000289479", "LNCPRESS1", "PRDX1", "TMSB4X", "TDGF1", "ENSG00000203279", "IFITM1", "ENSG00000248605", "LINC00698", "TPI1", "ENSG00000234261", "TRIM24", "ENSG00000284294", "PMAIP1", "MAD2L2", "USP9X", "TARS1", "FOXH1", "SCGB3A2", "EMSLR", "NMRK2", "SNRPN", "ENSG00000291111", "SEPHS1", "CARHSP1", "USP28", "AASS", "EPCAM", "HSPD1"),
  Neural_Progenitors = c("DLK1", "SIX3", "GPC3", "TUBA1A", "ENSG00000228999", "MAPK10", "DACH1", "PAMR1", "SOX5", "GREB1L", "ADGRL3", "RAX", "EFNA5", "NRG1", "CSRP2", "CDH2", "CTNNA2", "MAP2", "TPBG", "FRZB", "SHROOM3", "FBN2", "EPHA4", "TRPM3", "SYT1", "LIX1", "LHX2", "LMAN1", "FEZF1", "FZD3", "TMSB15A", "SSBP2", "FZD2", "NAV3", "MIR325HG", "LSAMP", "ZFHX4", "CNTNAP2", "NLGN1", "PRTG", "CDH6", "RPL6", "FHOD3", "RMST", "GPC6", "FEZF1-AS1", "SFRP1", "NNAT", "CDK14", "ST6GALNAC5"),
  Neurons            = c("SFRP1", "FEZ1", "FERMT2", "PTMS", "ZIC2", "JAKMIP2-AS1", "MAPK10", "NR6A1", "NNAT", "LHX5-AS1", "LMO1", "POU3F1", "FGFBP3", "DPYSL3", "SOX2", "CNTNAP2", "EPHA4", "GPC3", "FJX1", "PAMR1", "PTPN13", "DLGAP1", "GREB1L", "ASPH", "SHISA2", "CDH2", "LGR4", "FGFR1", "TNNT1", "CNKSR2", "NEBL", "RBFOX1", "HESX1", "OTX2", "PSIP1", "MAP2", "FHOD3", "LINC02751", "LRIG1", "DEK", "TUBA1A", "SOX11", "CRABP1", "XACT", "CHN1", "CYP26A1", "PTPRD", "LHX5", "NRP2", "NME4"),
  PGC                = c("GMPR", "COL23A1", "PCSK1N", "CERT1", "PLPP1", "POU5F1", "GABRA3", "TRPC5", "FRMD6", "PDPN", "EPB41L2", "APBB2", "SLC25A16", "LINC00937", "CACNA2D2", "AKAP12", "IFI16", "PHLDA3", "ASRGL1", "GPX2", "PRDM1", "PCAT14", "MKRN1", "NANOG", "DMD", "SLC16A10", "WASL", "NANOS3", "LAMA4", "S100A10", "TOMM7", "ZBTB44", "RTL4", "CD81", "NRK", "SOCS4", "TNRC6C", "RNF214", "NFE2L3", "IL6ST", "SLC4A8", "CEBPZ", "SMAGP", "PLBD1", "TXNDC12", "SEC63", "SUGCT", "CACNA2D3", "TMEM132D", "LINGO2"),
  TE                 = c("S100A11", "TPM1", "EPAS1", "ANXA2", "GABRP", "DSP", "KRT8", "KRT19", "TMSB10", "IGFBP7", "CD99", "MPZL1", "TBCA", "IGFBP3", "ANXA3", "RAB31", "SLC38A1", "AMOT", "KIAA1217", "ANKS1A", "ARHGAP29", "S100A10", "KRT18", "EPB41L3", "ATP2B1", "SESN3", "AHNAK", "GLB1", "CCBE1", "CALD1", "PKP2", "LCP1", "MYL12B", "LINC01924", "HAND1", "PTN", "TANC2", "P2RY6", "WLS", "SDK1", "CHN2", "SH3KBP1", "MID1", "RHOU", "LGALS3", "GATA3", "SUCO", "HAPLN1", "CTHRC1", "MEIS1"),
  VEC                = c("FLT1", "KDR", "TFPI", "GNG11", "EGFL7", "IGFBP4", "VAMP5", "IFI16", "IFITM2", "CALCRL", "FLI1", "VIM", "HLA-E", "SPTBN1", "PECAM1", "HAPLN1", "CLIC1", "NRP1", "PLK2", "ERG", "PLVAP", "ECSCR", "IFITM3", "CD99", "TMSB10", "B2M", "ARHGDIB", "VAV3", "ITM2B", "TIMP1", "ADGRL4", "PPFIBP1", "RNASE1", "LDB2", "RAP1B", "MEF2C", "TIMP3", "JAK1", "S100A10", "CDH5", "GYPC", "RASGRP3", "ARPC1B", "ETS1", "ANXA2", "GPX1", "S100A11", "DYSF", "ELK3", "ITGB1"),
  YSE                = c("S100A10", "FN1", "S100A11", "SERPINE2", "PPIB", "FSCN1", "RHOBTB3", "SPARC", "HAS2", "CD99", "BAMBI", "TMEM141", "MGST2", "TMSB10", "RGS5", "COL4A1", "APOC1", "SERPINH1", "ANXA6", "CALD1", "ALDH2", "APOA1", "HSP90B1", "RRBP1", "TBX3", "ANXA2", "CST3", "IL6ST", "COL4A2", "IFI16", "LIFR", "CAMK2D", "CD63", "LPGAT1", "FAM114A1", "ARSL", "NTS", "TTC3", "KRT8", "GATA6", "FGD6", "KRT19", "GSTO1", "TMED9", "FNDC3B", "DNAJC15", "PCDH7", "BEX1", "MACROH2A1", "RHOC")
)

# Preserve order of gene sets dynamically
set_order <- names(gene_sets)


# 2. Data Import and DESeq2 Pipeline ---------------------------------------

d8_ko_1_1 <- read.table(file = file.path(source_dir, "d8_ko_1_1_mm_counts_revised.txt"), header = TRUE)
d8_ko_1_2 <- read.table(file = file.path(source_dir, "d8_ko_1_2_mm_counts_revised.txt"), header = TRUE)
d8_ko_2_1 <- read.table(file = file.path(source_dir, "d8_ko_2_1_mm_counts_revised.txt"), header = TRUE)
d8_ko_2_2 <- read.table(file = file.path(source_dir, "d8_ko_2_2_mm_counts_revised.txt"), header = TRUE)
d8_wt_1   <- read.table(file = file.path(source_dir, "d8_wt_1_mm_counts_revised.txt"),   header = TRUE)
d8_wt_2   <- read.table(file = file.path(source_dir, "d8_wt_2_mm_counts_revised.txt"),   header = TRUE)

countData <- data.frame(
  ko_1_1 = d8_ko_1_1$counts,
  ko_1_2 = d8_ko_1_2$counts,
  ko_2_1 = d8_ko_2_1$counts,
  ko_2_2 = d8_ko_2_2$counts,
  wt_1   = d8_wt_1$counts,
  wt_2   = d8_wt_2$counts,
  row.names = d8_wt_2$Geneid
)

colData <- data.frame(
  condition = factor(c("KO", "KO", "KO", "KO", "WT", "WT"), levels = c("KO", "WT")),
  row.names = colnames(countData)
)

stopifnot(all(colnames(countData) == rownames(colData)))

dds_d8_F <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds_d8_F <- DESeq(dds_d8_F)
res_d8_F <- results(dds_d8_F, contrast = c("condition", "KO", "WT"))


# 3. Data Preparation for Bubble Plot -------------------------------------

norm_counts <- counts(dds_d8_F, normalized = TRUE) %>% 
  as.data.frame() %>%
  rownames_to_column("gene")

metadata <- as.data.frame(colData(dds_d8_F)) %>%
  rownames_to_column("sample")

# Calculate average expression and log2 ratio for each gene set
results_df <- map_dfr(names(gene_sets), function(set_name) {
  genes <- gene_sets[[set_name]]
  pseudocount <- 0.1
  
  expr_data <- norm_counts %>%
    dplyr::filter(gene %in% genes) %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    left_join(metadata, by = "sample")
  
  set_avgs <- expr_data %>%
    group_by(condition) %>%
    summarise(avg_expression = mean(expression, na.rm = TRUE), .groups = "drop")
  
  wt_avg <- set_avgs$avg_expression[set_avgs$condition == "WT"]
  ko_avg <- set_avgs$avg_expression[set_avgs$condition == "KO"]
  
  if (length(wt_avg) == 0) wt_avg <- 0
  if (length(ko_avg) == 0) ko_avg <- 0
  
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
    GeneSet   = factor(GeneSet, levels = set_order),
    condition = factor(condition, levels = c("WT", "KO"))
  )

max_abs_log <- max(abs(results_df$log_ratio), na.rm = TRUE)


# 4. Create Bubble Plot ----------------------------------------------------

p_bubble <- ggplot(results_df, aes(x = condition, y = GeneSet)) +
  geom_point(
    aes(size = avg_expression, fill = log_ratio),
    colour = "black",
    stroke = 0.7,
    shape  = 21,
    alpha  = 0.95
  ) +
  scale_size_continuous(
    name   = "Average \nExpression",
    range  = c(3.5, 8),
    breaks = pretty(results_df$avg_expression, n = 3)
  ) +
  scale_fill_distiller(
    name      = "Log2 Ratio\n(KO/WT)",
    type      = "div",
    palette   = "PuOr",
    direction = -1,
    limits    = c(-max_abs_log, max_abs_log)
  ) +
  scale_y_discrete(
    labels = function(label) {
      label %>%
        stringr::str_replace_all("^Skeletal_Muscle$", "Musculo\nskeletal") %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_wrap(width = 12)
    }
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "grey98", colour = NA),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 12, colour = "black", face = "plain"),
    axis.text.y      = element_text(
      size = 9.5, colour = "black", face = "plain", lineheight = 0.78
    ),
    axis.title       = element_blank(),
    legend.title     = element_text(face = "plain", size = 11),
    legend.text      = element_text(size = 10),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    aspect.ratio     = 2.05,
    plot.margin      = margin(8, 0, 8, 8),
    plot.background  = element_rect(fill = "white", colour = NA)
  ) +
  guides(
    size = guide_legend(
      override.aes = list(
        fill   = "grey70",
        colour = "black"
      ),
      keyheight = unit(0.6, "lines"),
      keywidth  = unit(0.6, "lines")
    ),
    fill = guide_colorbar(
      barwidth     = 4.5,
      barheight    = 0.7,
      ticks        = FALSE,
      label        = FALSE,
      frame.colour = NA
    )
  )


# 5. Prepare Data & Plot Diverging Bar Chart -------------------------------

padj_threshold <- 0.05

deg_df <- map_dfr(names(gene_sets), function(set_name) {
  genes <- gene_sets[[set_name]]
  
  set_res <- as.data.frame(res_d8_F) %>%
    rownames_to_column("gene") %>%
    filter(gene %in% genes) %>%
    mutate(
      direction = case_when(
        padj < padj_threshold & log2FoldChange > 0 ~ "Up",
        padj < padj_threshold & log2FoldChange < 0 ~ "Down",
        TRUE ~ "Non-sig"
      ),
      direction = factor(direction, levels = c("Down", "Up", "Non-sig"))
    )
  
  deg_counts <- table(set_res$direction)
  total_genes <- max(length(genes), 1)
  
  tibble(
    GeneSet    = set_name,
    Direction  = c("Down", "Up"),
    Percentage = c(
      (deg_counts["Down"] / total_genes) * 100,
      (deg_counts["Up"] / total_genes) * 100
    )
  )
})

deg_df <- deg_df %>%
  mutate(
    x_position = ifelse(Direction == "Down", -Percentage, Percentage),
    GeneSet    = factor(GeneSet, levels = set_order)
  )

max_x <- max(abs(deg_df$x_position), na.rm = TRUE) * 1.1
if (max_x == 0) max_x <- 10

p_diverge <- ggplot(deg_df, aes(y = GeneSet, x = x_position, fill = Direction)) +
  geom_col(width = 0.8) +
  geom_vline(xintercept = 0, color = "black", linewidth = 1.2) +
  scale_fill_manual(values = c(Down = "#3C8DBC", Up = "#E36B6B")) +
  guides(fill = "none") +
  scale_x_continuous(
    name   = "% Differentially\nExpressed Genes",
    labels = abs,
    breaks = pretty(c(-max_x, max_x), n = 3),
    limits = c(-max_x, max_x)
  ) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  labs(y = NULL, title = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.border       = element_rect(color = "black", linewidth = 0.3),
    panel.grid.major   = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_line(colour = "grey92", linewidth = 0.15),
    axis.text.x        = element_text(size = 11, colour = "black"),
    axis.title.x       = element_text(size = 11, face = "plain"),
    axis.title.y       = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    legend.position    = "none",
    plot.margin        = margin(t = 18, r = 8, b = 8, l = 0)
  ) +
  coord_cartesian(clip = "off") +
  annotate(
    "text",
    x = -max_x * 0.45,
    y = Inf,
    label = "DOWN",
    vjust = -1,
    size  = 4,
    colour = "#3C8DBC",
    fontface = "plain"
  ) +
  annotate(
    "text",
    x = max_x * 0.45,
    y = Inf,
    label = "UP",
    vjust = -1,
    size  = 4,
    colour = "#E36B6B",
    fontface = "plain"
  )


# 6. Combine and Export Plots ---------------------------------------------

p_combined <- (
  p_bubble + p_diverge +
    plot_layout(widths = c(0.78, 0.8), guides = "collect")
) & theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "Fig2F_bubble_and_diverging_compact.tiff"),
  plot        = p_combined,
  width       = 5.8,
  height      = 4.8,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)
