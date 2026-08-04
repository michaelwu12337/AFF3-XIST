{
  library(ggplot2)
  library(dplyr)
  library(tidyr)
}

# Input and output directories -------------------------------------------
{
  source_dir <- "source_files"
  figure_dir <- file.path("generated_figures", "Fig4E")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

# Morphology counts for Figure 4E from the quantification slides ---------
{
  morphology_counts <- read.csv(
    file.path(source_dir, "Fig4E_gastruloid_morphology_counts.csv"),
    check.names = FALSE
  )

  stopifnot(
    all(
      morphology_counts$Polarization +
        morphology_counts[["No polarization"]] ==
        morphology_counts$Total
    )
  )
}

# Match the Figure 3F percentage-bar aesthetics --------------------------
{
  fill_cols <- c(
    "Polarization" = "#6F5A8A",
    "No polarization" = "#DDE6F6"
  )

  make_percentage_plot <- function(plot_data, group_levels) {
    plot_long <- plot_data %>%
      pivot_longer(
        cols = c("Polarization", "No polarization"),
        names_to = "Category",
        values_to = "Count"
      ) %>%
      group_by(Group) %>%
      mutate(Percent = Count / sum(Count) * 100) %>%
      ungroup() %>%
      mutate(
        Group = factor(Group, levels = group_levels),
        Category = factor(
          Category,
          levels = c("No polarization", "Polarization")
        )
      )

    ggplot(plot_long, aes(x = Group, y = Percent, fill = Category)) +
      geom_col(
        width = 0.6,
        color = "black",
        linewidth = 0.8
      ) +
      scale_y_continuous(
        breaks = c(0, 50, 100),
        expand = c(0, 0)
      ) +
      coord_cartesian(ylim = c(0, 100), expand = FALSE) +
      scale_fill_manual(values = fill_cols) +
      guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
      labs(
        x = NULL,
        y = "Gastruloids (%)",
        fill = NULL
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(color = "black", size = 18),
        axis.line = element_line(linewidth = 1, color = "black"),
        axis.ticks = element_line(linewidth = 1, color = "black"),
        axis.ticks.length = unit(6, "pt"),
        legend.position = "top",
        legend.text = element_text(size = 18),
        legend.key.height = unit(1.0, "lines"),
        legend.key.width = unit(1, "lines")
      )
  }
}

# WT, KO, and DKO ---------------------------------------------------------
{
  dko_plot <- make_percentage_plot(
    morphology_counts %>% filter(Experiment == "DKO"),
    group_levels = c("WT", "KO", "DKO")
  )

  ggsave(
    filename = file.path(
      figure_dir,
      "Fig4E_DKO_gastruloid_morphology_percentages.svg"
    ),
    plot = dko_plot,
    device = svglite::svglite,
    width = 3,
    height = 3.7,
    units = "in"
  )
}

# WT, KO, and KO+X1 -------------------------------------------------------
{
  kox1_plot <- make_percentage_plot(
    morphology_counts %>% filter(Experiment == "KO+X1"),
    group_levels = c("WT", "KO", "KO+X1")
  )

  ggsave(
    filename = file.path(
      figure_dir,
      "Fig4E_KOX1_gastruloid_morphology_percentages.svg"
    ),
    plot = kox1_plot,
    device = svglite::svglite,
    width = 3,
    height = 3.7,
    units = "in"
  )
}
