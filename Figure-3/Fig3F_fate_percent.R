{
  setwd("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Figure 3 new/Barplot lineage fate percent")
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
}

# PSC-origin contribution within each cluster (PSC_origin_fraction_%)
pscfrac_df <- tibble(
  Cluster = c(
    "Amnion",
    "Blood Progenitors",
    "Cardiac Mesoderm",
    "Neural Progenitors",
    "Neurons",
    "Epi",
    "TB"
  ),
  WT_percent = c(
    0.00,
    0.00,
    0.00,
    43.24,
    32.43,
    24.32,
    0.00
  ),
  KO_percent = c(
    2.5,
    42.5,
    12.5,
    0.00,
    0.00,
    5.0,
    37.5
  )
)

psc_long <- pscfrac_df %>%
  pivot_longer(
    cols = c(WT_percent, KO_percent),
    names_to = "Condition",
    values_to = "Percent"
  ) %>%
  mutate(
    Condition = recode(Condition, WT_percent = "WT", KO_percent = "KO"),
    Cluster = factor(Cluster, levels = c(
      "Amnion",
      "Blood Progenitors",
      "Cardiac Mesoderm",
      "Neural Progenitors",
      "Neurons",
      "Epi",
      "TB"
    ))
  )

# no zeros, two separate plots
{
  # ---- axis settings ----
  XMAX_WT <- 50
  XBREAKS_WT <- c(0, 15, 30, 45)
  
  XMAX_KO <- 45
  XBREAKS_KO <- c(0, 15, 30, 45)
  
  # ---- BAR COLOR CONTROL ----
  KO_FILL  <- "#F4B183"   
  WT_FILL  <- "#9ED3E6" 
  
  WT_BORDER <- NA     
  KO_BORDER <- NA
  
  BAR_LW  <- 0.5           
  ALPHA   <- 0.7
  
  # ---- AXIS CONTROL ----
  AXIS_LW <- 0.5
  TICK_LW <- 0.5
  
  # ---- TEXT CONTROL ----
  BASE_SIZE        <- 20
  AXIS_TEXT_SIZE_X <- 20
  AXIS_TEXT_SIZE_Y <- 20
  VALUE_TEXT_SIZE  <- 6.5
  
  gp_theme_thin <- theme_classic(base_size = BASE_SIZE) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = AXIS_TEXT_SIZE_X, face = "bold", color = "black"),
      axis.text.y  = element_text(size = AXIS_TEXT_SIZE_Y, face = "bold", color = "black"),
      axis.line    = element_line(linewidth = AXIS_LW, color = "black"),
      axis.ticks   = element_line(linewidth = TICK_LW, color = "black"),
      axis.ticks.length = unit(5, "pt"),
      legend.position = "none",
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
  
  # ---- WT plot ----
  p_WT <- psc_long %>%
    filter(Condition == "WT", Percent > 0) %>%
    mutate(Cluster = reorder(Cluster, Percent)) %>%   # reorder by value
    ggplot(aes(x = Cluster, y = Percent)) +
    geom_col(
      width = 0.75,
      fill = WT_FILL,
      color = WT_BORDER,
      linewidth = BAR_LW,
      alpha = ALPHA
    ) +
    geom_text(
      aes(label = sprintf("%.2f", Percent)),
      hjust = -0.15,
      size = VALUE_TEXT_SIZE,
      fontface = "bold"
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, XMAX_WT),
      breaks = XBREAKS_WT,
      expand = expansion(mult = c(0, 0.12))
    ) +
    gp_theme_thin
  
  # ---- KO plot ----
  p_KO <- psc_long %>%
    filter(Condition == "KO", Percent > 0) %>%
    mutate(Cluster = reorder(Cluster, Percent)) %>%   # reorder by value
    ggplot(aes(x = Cluster, y = Percent)) +
    geom_col(
      width = 0.75,
      fill = KO_FILL,
      color = KO_BORDER,
      linewidth = BAR_LW,
      alpha = ALPHA
    ) +
    geom_text(
      aes(label = sprintf("%.2f", Percent)),
      hjust = -0.15,
      size = VALUE_TEXT_SIZE,
      fontface = "bold"
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, XMAX_KO),
      breaks = XBREAKS_KO,
      expand = expansion(mult = c(0, 0.12))
    ) +
    gp_theme_thin
  
  p_WT + p_KO
}

ggsave(
  "fig3_fate_percent_WT_02.tiff",
  plot = p_WT,
  width = 6,    
  height = 3,      
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  compression = "lzw"
)

ggsave(
  "fig3_fate_percent_KO_02.tiff",
  plot = p_KO,
  width = 7.5,    
  height = 5,      
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  compression = "lzw"
)