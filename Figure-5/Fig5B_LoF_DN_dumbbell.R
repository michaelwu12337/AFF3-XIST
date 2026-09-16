{
  library(tidyverse)
  library(ggplot2)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "Fig5B")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Percentages for the LoF and DN comparison -----------------------------
{
  kinssh_bp <- read.csv(
    file.path(source_dir, "Fig5B_LoF_DN_dumbbell_summary.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  phenotype_order <- c(
    "Severe DD/ID",
    "Seizure/EEG",
    "Abnormal brain MRI",
    "Abnormal muscle tone"
  )

  kinssh_bp <- kinssh_bp %>%
    filter(Label %in% phenotype_order) %>%
    mutate(
      Label = factor(Label, levels = rev(phenotype_order))
    )
}

# Figure 5B --------------------------------------------------------------
{
  fig6b_plot <- ggplot(
    kinssh_bp,
    aes(y = Label)
  ) +
    geom_segment(
      aes(
        x = LoF_pct,
        xend = DN_pct,
        yend = Label
      ),
      color = "grey70",
      linewidth = 1.1
    ) +
    geom_point(
      aes(x = LoF_pct, colour = "LoF"),
      size = 7
    ) +
    geom_point(
      aes(x = DN_pct, colour = "DN"),
      size = 7
    ) +
    scale_colour_manual(
      name = NULL,
      breaks = c("LoF", "DN"),
      values = c(
        "LoF" = "#C44E52",
        "DN" = "#4C72B0"
      )
    ) +
    scale_x_continuous(
      name = "Patients with phenotype (%)",
      limits = c(-5, 105),
      breaks = seq(0, 100, 20),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(
      labels = function(label) stringr::str_wrap(label, width = 16)
    ) +
    labs(
      y = NULL,
      title = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      plot.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20),
      legend.position = "top",
      legend.justification = "left"
    )

  ggsave(
    file.path(figure_dir, "Fig5B_LoF_DN_dumbbell.tiff"),
    plot = fig6b_plot,
    width = 9,
    height = 8,
    units = "in",
    dpi = 600,
    limitsize = FALSE,
    compression = "lzw"
  )
}
