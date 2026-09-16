{
  library(tidyverse)
  library(pheatmap)
  library(grid)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "Fig5A")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Clinical classifications for the LoF and DN comparison ----------------
{
  clinical_df <- read.csv(
    file.path(source_dir, "Fig5A_clinical_heatmap_matrix.csv"),
    na.strings = c("", "NA"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  phenotype_cols <- c(
    "DD_ID_severe",
    "DD_ID_moderate",
    "DD_ID_mild",
    "Seizure_abnormal_EEG",
    "Brain_MRI_abnormality",
    "Muscle_tone_abnormality",
    "Skeletal_abnormality"
  )

  clinical_df <- clinical_df %>%
    mutate(
      Variant_group = factor(Variant_group, levels = c("DN", "LoF")),
      Sex = factor(Sex, levels = c("M", "F"))
    ) %>%
    arrange(Variant_group, Sex, Plot_order)

  mat <- clinical_df %>%
    dplyr::select(all_of(phenotype_cols)) %>%
    as.matrix()

  rownames(mat) <- clinical_df$Patient

  ann_row <- clinical_df %>%
    dplyr::select(Variant_group) %>%
    as.data.frame()

  rownames(ann_row) <- clinical_df$Patient
  colnames(ann_row) <- "Variant group"

  ann_col_sex <- clinical_df %>%
    dplyr::select(Sex) %>%
    as.data.frame()

  rownames(ann_col_sex) <- clinical_df$Patient

  ann_colors <- list(
    "Variant group" = c("LoF" = "#C44E52", "DN" = "#4C72B0")
  )

  sex_colors <- list(
    Sex = c("M" = "#4C72B0", "F" = "#C44E52")
  )

  # Grey is kept for values that were not reported or could not be classified.
  mat_plot <- mat
  mat_plot[is.na(mat_plot)] <- 2
}

# Figure 5A --------------------------------------------------------------
{
  draw_fig6a <- function(file_name, canvas_width, cell_width) {
    tiff(
      file.path(figure_dir, file_name),
      width = canvas_width,
      height = 7.2,
      units = "in",
      res = 600
    )

    heatmap_plot <- pheatmap(
      mat_plot,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = ann_row,
      annotation_colors = ann_colors,
      annotation_names_row = FALSE,
      gaps_row = 16,
      color = c("white", "#7A1F1F", "grey85"),
      breaks = c(-0.5, 0.5, 1.5, 2.5),
      border_color = "grey70",
      fontsize = 8,
      fontsize_row = 9,
      fontsize_col = 8,
      cellwidth = cell_width,
      cellheight = 14,
      angle_col = 45,
      labels_col = c(
        "Severe DD/ID",
        "Moderate DD/ID",
        "Mild DD/ID",
        "Seizure/Epilepsy",
        "Abnormal brain MRI",
        "Abnormal muscle tone",
        "Skeletal abnormalities"
      ),
      legend_breaks = c(0, 1, 2),
      legend_labels = c("Absent", "Present", "Not reported/unknown"),
      main = "",
      silent = TRUE
    )

    # Keep the group strip the same width as the heatmap cells.
    annotation_layout <- heatmap_plot$gtable$layout[
      heatmap_plot$gtable$layout$name == "row_annotation",
      ,
      drop = FALSE
    ]
    if (nrow(annotation_layout) == 1) {
      heatmap_plot$gtable$widths[annotation_layout$l] <- unit(
        cell_width,
        "bigpts"
      )
    }

    grid.newpage()
    grid.draw(heatmap_plot$gtable)

    dev.off()
  }

  draw_fig6a_transposed <- function() {
    make_status_legend <- function() {
      labels <- c("Present", "Absent", "Not reported/\nunknown")
      fills <- c("#7A1F1F", "white", "grey85")

      legend_table <- gtable::gtable(
        widths = unit.c(
          unit(1, "in"),
          rep(unit(14, "bigpts"), 3),
          unit(1, "in")
        ),
        heights = unit.c(
          unit(0.28, "in"),
          unit(14, "bigpts"),
          unit(0.85, "in")
        )
      )

      legend_table <- gtable::gtable_add_grob(
        legend_table,
        textGrob(
          "Phenotype status",
          x = unit(1, "in"),
          just = "left",
          gp = gpar(fontsize = 9, fontface = "bold")
        ),
        t = 1,
        l = 1,
        r = 5
      )

      for (i in seq_along(labels)) {
        legend_table <- gtable::gtable_add_grob(
          legend_table,
          rectGrob(
            width = unit(1, "npc"),
            height = unit(1, "npc"),
            gp = gpar(fill = fills[i], col = "grey50", lwd = 1)
          ),
          t = 2,
          l = i + 1
        )
        legend_table <- gtable::gtable_add_grob(
          legend_table,
          textGrob(
            labels[i],
            x = unit(0.5, "npc"),
            y = unit(1, "npc"),
            just = c("right", "top"),
            rot = 45,
            gp = gpar(fontsize = 8)
          ),
          t = 3,
          l = i + 1,
          clip = "off"
        )
      }

      legend_table
    }

    make_sex_legend <- function() {
      labels <- c("Male", "Female")
      fills <- c("#4C72B0", "#C44E52")

      legend_table <- gtable::gtable(
        widths = unit.c(
          unit(0.75, "in"),
          rep(unit(14, "bigpts"), 2),
          unit(0.75, "in")
        ),
        heights = unit.c(
          unit(0.28, "in"),
          unit(14, "bigpts"),
          unit(0.85, "in")
        )
      )

      legend_table <- gtable::gtable_add_grob(
        legend_table,
        textGrob(
          "Gender",
          x = unit(0.75, "in"),
          just = "left",
          gp = gpar(fontsize = 9, fontface = "bold")
        ),
        t = 1,
        l = 1,
        r = 4
      )

      for (i in seq_along(labels)) {
        legend_table <- gtable::gtable_add_grob(
          legend_table,
          rectGrob(
            width = unit(1, "npc"),
            height = unit(1, "npc"),
            gp = gpar(fill = fills[i], col = "grey50", lwd = 1)
          ),
          t = 2,
          l = i + 1
        )
        legend_table <- gtable::gtable_add_grob(
          legend_table,
          textGrob(
            labels[i],
            x = unit(0.5, "npc"),
            y = unit(1, "npc"),
            just = c("right", "top"),
            rot = 45,
            gp = gpar(fontsize = 8)
          ),
          t = 3,
          l = i + 1,
          clip = "off"
        )
      }

      legend_table
    }

    tiff(
      file.path(figure_dir, "Fig5A_clinical_heatmap_transposed.tiff"),
      width = 11,
      height = 4.8,
      units = "in",
      res = 600
    )

    heatmap_plot <- pheatmap(
      t(mat_plot),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_col = ann_col_sex,
      annotation_colors = sex_colors,
      annotation_names_col = FALSE,
      gaps_col = 16,
      color = c("white", "#7A1F1F", "grey85"),
      breaks = c(-0.5, 0.5, 1.5, 2.5),
      border_color = "grey70",
      fontsize = 8,
      fontsize_row = 9,
      fontsize_col = 8,
      cellwidth = 14,
      cellheight = 14,
      angle_col = 45,
      labels_row = c(
        "Severe DD/ID",
        "Moderate DD/ID",
        "Mild DD/ID",
        "Seizure/Epilepsy",
        "Abnormal brain MRI",
        "Abnormal muscle tone",
        "Skeletal abnormalities"
      ),
      legend = FALSE,
      annotation_legend = FALSE,
      main = "",
      silent = TRUE
    )

    # Keep the sex strip the same height as the square heatmap cells.
    annotation_layout <- heatmap_plot$gtable$layout[
      heatmap_plot$gtable$layout$name == "col_annotation",
      ,
      drop = FALSE
    ]
    if (nrow(annotation_layout) == 1) {
      heatmap_plot$gtable$heights[annotation_layout$t] <- unit(
        14,
        "bigpts"
      )
    }

    status_legend <- make_status_legend()
    sex_legend <- make_sex_legend()

    bottom_legends <- gtable::gtable(
      widths = unit.c(
        unit(1, "null"),
        unit(3.2, "in"),
        unit(0.15, "in"),
        unit(2.3, "in"),
        unit(1, "null")
      ),
      heights = unit(1.6, "in")
    )
    bottom_legends <- gtable::gtable_add_grob(
      bottom_legends,
      status_legend,
      t = 1,
      l = 2
    )
    bottom_legends <- gtable::gtable_add_grob(
      bottom_legends,
      sex_legend,
      t = 1,
      l = 4
    )

    heatmap_gtable <- gtable::gtable_add_rows(
      heatmap_plot$gtable,
      heights = unit(1.6, "in")
    )
    heatmap_gtable <- gtable::gtable_add_grob(
      heatmap_gtable,
      bottom_legends,
      t = length(heatmap_gtable$heights),
      l = 1,
      r = length(heatmap_gtable$widths)
    )

    grid.newpage()
    grid.draw(heatmap_gtable)

    dev.off()
  }

  draw_fig6a(
    file_name = "Fig5A_clinical_heatmap.tiff",
    canvas_width = 4.8,
    cell_width = 14
  )

  draw_fig6a(
    file_name = "Fig5A_clinical_heatmap_wide.tiff",
    canvas_width = 6.5,
    cell_width = 26
  )

  draw_fig6a_transposed()
}
