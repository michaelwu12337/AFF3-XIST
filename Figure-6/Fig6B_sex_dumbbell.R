{
  library(tidyverse)
  library(ggplot2)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "Fig6B")

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Sex-stratified percentages used in the v46 dumbbell plot --------------
{
  kinssh_bp <- read.csv(
    file.path(source_dir, "Fig6B_sex_dumbbell_summary.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Figure 6B --------------------------------------------------------------
{
  fig6b_plot <- ggplot(
    kinssh_bp,
    aes(y = fct_reorder(Variable, Female_pct - Male_pct))
  ) +
    geom_segment(
      aes(
        x = Male_pct,
        xend = Female_pct,
        y = Variable,
        yend = Variable
      ),
      color = "grey70",
      linewidth = 1.1
    ) +
    geom_point(
      aes(x = Male_pct, y = Variable, colour = "Male"),
      size = 7
    ) +
    geom_point(
      aes(x = Female_pct, y = Variable, colour = "Female"),
      size = 7
    ) +
    scale_colour_manual(
      name = NULL,
      values = c(
        "Female" = "#E36B6B",
        "Male" = "#3C8DBC"
      )
    ) +
    scale_x_continuous(
      name = "Patients with symptom (%)",
      limits = c(-5, 60),
      breaks = seq(0, 60, 20),
      expand = expansion(mult = c(0, 0.05))
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
    file.path(figure_dir, "Fig6B_sex_dumbbell.tiff"),
    plot = fig6b_plot,
    width = 12,
    height = 6,
    units = "in",
    dpi = 600,
    limitsize = FALSE,
    compression = "lzw"
  )
}
