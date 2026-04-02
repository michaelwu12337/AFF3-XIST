{
  setwd("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Figure 5 new/Barplot")
  library(ggplot2)
}


# For total neural lineage percentage -------------------------------------
df <- data.frame(
  group = factor(c("WT", "KO", "DKO"), levels = c("WT", "KO", "DKO")),
  pct   = c(40, 21, 37.5)
)

bar_cols <- c(WT = "#ff1a1a", KO = "#b07d3c", DKO = "#ff8c00")
shape_map <- c(WT = 15, KO = 17, DKO = 18)

ggplot(df, aes(x = group, y = pct, fill = group)) +
  geom_col(
    width = 0.55,
    color = "black",     # outline color
    linewidth = 1.2      # outline thickness
  ) +
  geom_point(aes(shape = group), size = 4.2, color = "black") +
  scale_fill_manual(values = bar_cols, guide = "none") +
  scale_shape_manual(values = shape_map, guide = "none") +
  scale_y_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, 10),
    expand = c(0, 0)
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "% of Neuronal Lineage Cells"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 13, face = "bold", color = "black"),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(6, "pt")
  )

ggsave(
  "fig5_bar_02.tiff",
  plot = last_plot(),
  width = 3.5,    
  height = 4.5,      
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  compression = "lzw"
)
