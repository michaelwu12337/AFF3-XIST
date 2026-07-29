{
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(dplyr)
  library(tidyr)
}

## ---- cell-percentage data ----
df <- read.csv(
  "source_files/Fig3F_RNA_FISH_cell_percentages.csv",
  check.names = FALSE
)

## ---- long format for ggplot ----
df_long <- df %>%
  pivot_longer(cols = c(`XIST+`, `No signal`),
               names_to = "Category",
               values_to = "Percent") %>%
  mutate(
    Group = factor(Group, levels = c("WT", "KO")),
    Category = factor(Category, levels = c("XIST+", "No signal"))
  )
df_long$Category <- factor(
  df_long$Category,
  levels = c("No signal", "XIST+")
)

## ---- colors (GraphPad-ish) ----
fill_cols <- c(
  "XIST+" = "#6F5A8A",  # muted purple
  "No signal"   = "#DDE6F6"   # light blue-grey
)

## ---- plot ----
{
  p <- ggplot(df_long, aes(x = Group, y = Percent, fill = Category)) +
    geom_col(
      width = 0.6,
      color = "black",
      linewidth = 0.8
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 50, 100),
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = fill_cols) +
    labs(
      title = "XIST RNA FISH",
      x = NULL,
      y = "Cells (%)",
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
      legend.key.width  = unit(1, "lines")
    )
}

## ---- save (optional) ----
ggsave(
  filename = "Fig3F_RNA_FISH_cell_percentages.svg",
  plot     = p,
  device   = svglite::svglite,
  width    = 3,
  height   = 3.7,
  units    = "in"
)


## ===============================
## RNA FISH intensity
## ===============================
{
  ## Intensity data ----
  df1 <- read.csv(
    "source_files/Fig3F_RNA_FISH_intensity.csv"
  ) %>%
    mutate(Group = factor(Group, levels = c("WT", "KO")),
           Value_div_1000 = Value / 1000)  # Create a new column with values divided by 1000
  
  ## ===============================
  ## Colors (filled bars; GraphPad-like) ----
  {
    wt_fill <- "#F2B5AE"   # light red fill
    wt_pt   <- "#E41A1C"   # red points
    
    ko_fill <- "#E6D6C3"   # light brown fill
    ko_pt   <- "#A67C52"   # brown points
  }
  
  ## ===============================
  ## Summary for bars (mean ± SEM) ----
  sum_df <- df1 %>%
    group_by(Group) %>%
    summarise(
      mean = mean(Value_div_1000, na.rm = TRUE),  # Use new column
      sem  = sd(Value_div_1000, na.rm = TRUE) / sqrt(sum(!is.na(Value_div_1000))),
      .groups = "drop"
    )
  
  ## ===============================
  ## Stats (rstatix -> stat_pvalue_manual) ----
  {
    stat_test <- df1 %>%
      t_test(Value_div_1000 ~ Group) %>%  # Use new column
      add_significance("p") %>%
      add_xy_position(x = "Group")   # creates xmin/xmax
    
    # set annotation height nicely above the highest point
    y_max <- max(df1$Value_div_1000, na.rm = TRUE)  # Use new column
    stat_test$y.position <- y_max * 1.12
  }
}

## ===============================
## Plot (match monoallelic aesthetics) ----
{
  p <- ggplot() +
    ## bars (mean) - black border like monoallelic plot
    geom_col(
      data = sum_df,
      aes(x = Group, y = mean, fill = Group),
      width = 0.6,
      color = "black",
      linewidth = 0.8
    ) +
    ## error bars (SEM)
    geom_errorbar(
      data = sum_df,
      aes(x = Group, ymin = mean - sem, ymax = mean + sem),
      width = 0.18,
      linewidth = 0.8,
      color = "black"
    ) +
    ## points: WT triangles, KO diamonds (all data points shown)
    geom_jitter(
      data = df1 %>% filter(Group == "WT"),
      aes(x = Group, y = Value_div_1000),  # Use new column
      shape = 17, size = 3,
      color = wt_pt,
      width = 0.08, height = 0
    ) +
    geom_jitter(
      data = df1 %>% filter(Group == "KO"),
      aes(x = Group, y = Value_div_1000),  # Use new column
      shape = 18, size = 3,
      color = ko_pt,
      width = 0.08, height = 0
    ) +
    ## fills
    scale_fill_manual(values = c(WT = wt_fill, KO = ko_fill), guide = "none") +
    ## axes / labels (GraphPad-ish, like monoallelic plot)
    labs(
      title = "XIST Intensity",
      x = NULL,
      y = "Mean XIST Intensity \n(in ROI)"
    ) +
    scale_y_continuous(
      limits = c(0, y_max * 1.20),
      breaks = c(0, 5, 10, 15, 20),  # Adjust as needed
      expand = c(0, 0)
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(color = "black", size = 18),
      axis.text.y = element_text(color = "black", size = 13),
      axis.line = element_line(linewidth = 1, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(6, "pt")
    ) +
    ## significance (stat_pvalue_manual, no tips)
    stat_pvalue_manual(
      stat_test,
      label = "p.signif",
      tip.length = 0,
      bracket.size = 0.8,
      label.size = 7
    )
}


ggsave(
  filename = "Fig3F_RNA_FISH_intensity.svg",
  plot     = p,
  device   = svglite::svglite,
  width    = 3,
  height   = 3.7,
  units    = "in"
)
