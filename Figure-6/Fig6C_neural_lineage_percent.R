dir.create("generated_figures/Fig6C", recursive = TRUE, showWarnings = FALSE)

{
  library(ggplot2)
}


# For total neural lineage percentage -------------------------------------
df <- data.frame(
  group = factor(c("WT", "KO", "DKO"), levels = c("DKO", "KO", "WT")),
  pct   = c(40, 21, 37.5)
)

bar_cols <- c(
  WT  = "#E45757",   # softer red
  KO  = "#BFA06A",   # lighter, less muddy brown
  DKO = "#F4A261"    # muted orange (more pastel)
)
shape_map <- c(WT = 15, KO = 17, DKO = 19)

ggplot(df, aes(x = group, y = pct, fill = group)) +
  geom_col(
    width = 0.55,
    color = "black",
    linewidth = 1        
  ) +
  geom_point(aes(shape = group), size = 4.2, color = "black") +
  geom_text(
    aes(label = pct),
    hjust = -0.6,
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_manual(values = bar_cols, guide = "none") +
  scale_shape_manual(values = shape_map, guide = "none") +
  scale_y_continuous(
    limits = c(0, 45),
    breaks = seq(0, 45, 10),
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
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.line = element_line(size = 1, color = "black"),
    axis.ticks = element_line(size = 1, color = "black"),
    axis.ticks.length = unit(6, "pt")
  ) + 
  coord_flip()

ggsave(
  "generated_figures/Fig6C/Fig6C_neural_lineage_percent.tiff",
  plot = last_plot(),
  width = 6.5,    
  height = 2,      
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  compression = "lzw"
)
