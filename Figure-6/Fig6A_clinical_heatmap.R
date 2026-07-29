{
  library(tidyverse)
  library(pheatmap)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "Fig6A")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Clinical classifications used in the v46 heatmap ----------------------
{
  clinical_df <- read.csv(
    file.path(source_dir, "Fig6A_clinical_heatmap_matrix.csv"),
    na.strings = c("", "NA"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  phenotype_cols <- setdiff(
    colnames(clinical_df),
    c("Patient", "Sex", "Zygosity")
  )

  clinical_df <- clinical_df %>%
    mutate(
      Sex = factor(Sex, levels = c("F", "M")),
      Zygosity = factor(
        Zygosity,
        levels = c(
          "Het LoF/+",
          "Compound LoF/missense",
          "Hom LoF/LoF"
        )
      )
    )

  mat <- clinical_df %>%
    dplyr::select(all_of(phenotype_cols)) %>%
    as.matrix()

  rownames(mat) <- clinical_df$Patient

  ann_row <- clinical_df %>%
    dplyr::select(Sex, Zygosity) %>%
    as.data.frame()

  rownames(ann_row) <- clinical_df$Patient

  ann_colors <- list(
    Sex = c("F" = "#C44E52", "M" = "#4C72B0"),
    Zygosity = c(
      "Het LoF/+" = "#7F7F7F",
      "Compound LoF/missense" = "#DD8452",
      "Hom LoF/LoF" = "#9467BD"
    )
  )

  # Preserve the v46 treatment of missing values as "Not applicable".
  mat_plot <- mat
  mat_plot[is.na(mat_plot)] <- 2
}

# Figure 6A --------------------------------------------------------------
{
  tiff(
    file.path(figure_dir, "Fig6A_clinical_heatmap.tiff"),
    width = 5.8,
    height = 4.5,
    units = "in",
    res = 600
  )

  pheatmap(
    mat_plot,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    annotation_names_row = FALSE,
    color = c("white", "#7A1F1F", "grey85"),
    breaks = c(-0.5, 0.5, 1.5, 2.5),
    border_color = "grey70",
    fontsize = 8,
    fontsize_row = 9,
    fontsize_col = 8,
    angle_col = 45,
    legend_breaks = c(0, 1, 2),
    legend_labels = c("Absent", "Present", "Not applicable"),
    main = ""
  )

  dev.off()
}
