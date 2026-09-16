# Figure 3F: combined XIST/XACT categories from the supplied approximate counts

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})


# 1. Paths and figure settings -------------------------------------------

source_file <- file.path("source_files", "Fig3F_RNA_FISH_XIST_XACT_counts.csv")
output_dir <- file.path("generated_figures", "Fig3F")
output_file <- file.path(
  output_dir,
  "Fig3F_RNA_FISH_XIST_XACT_allelic_status_requantified.tiff"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_order <- c("WT", "KO")
category_order <- c(
  "NO_XIST / Mono_XACT", "NO_XIST / Bi_XACT",
  "Mono_XIST / Mono_XACT", "Bi_XIST / Mono_XACT"
)
category_colors <- setNames(
  c("#E9DDA8", "#A9BDA0", "#D99D58", "#B85D56"),
  category_order
)

figure_width <- 5.4
figure_height <- 3.7


# 2. Load counts and calculate within-group percentages ------------------

count_data <- read.csv(source_file, check.names = FALSE)

stopifnot(
  identical(
    names(count_data),
    c("Group", "XIST_status", "XACT_status", "Approximate_count")
  ),
  all(count_data$Group %in% group_order),
  all(count_data$XIST_status %in% c("NO_XIST", "Mono_XIST", "Bi_XIST")),
  all(count_data$XACT_status %in% c("Mono_XACT", "Bi_XACT")),
  all(count_data$Approximate_count > 0)
)

group_totals <- count_data %>%
  group_by(Group) %>%
  summarise(Total = sum(Approximate_count), .groups = "drop")

total_by_group <- setNames(group_totals$Total, group_totals$Group)

plot_data <- count_data %>%
  mutate(Category = paste(XIST_status, XACT_status, sep = " / ")) %>%
  group_by(Group, Category) %>%
  summarise(Count = sum(Approximate_count), .groups = "drop") %>%
  complete(
    Group = group_order,
    Category = category_order,
    fill = list(Count = 0)
  ) %>%
  left_join(group_totals, by = "Group") %>%
  mutate(
    Percent = 100 * Count / Total,
    Group = factor(Group, levels = group_order),
    Category = factor(Category, levels = category_order)
  )

stopifnot(
  all(total_by_group[group_order] == c(WT = 6, KO = 10)),
  all(abs(tapply(plot_data$Percent, plot_data$Group, sum) - 100) < 1e-8)
)


# 3. Plot and export ------------------------------------------------------

xist_status_plot <- ggplot(
  plot_data,
  aes(x = Group, y = Percent, fill = Category)
) +
  geom_col(
    width = 0.40,
    position = position_stack(reverse = TRUE),
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = category_colors,
    breaks = category_order,
    labels = c(
      "No XIST / Mono XACT", "No XIST / Bi XACT",
      "Mono XIST / Mono XACT", "Bi XIST / Mono XACT"
    ),
    name = NULL,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 50, 100),
    expand = c(0, 0)
  ) +
  scale_x_discrete(expand = expansion(add = 0.38)) +
  labs(x = NULL, y = "Cells (%)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11, face = "plain"),
    legend.key.size = grid::unit(0.45, "cm"),
    axis.title.y = element_text(size = 18, face = "plain"),
    axis.text.x = element_text(size = 18, color = "black", face = "plain"),
    axis.text.y = element_text(size = 15, color = "black", face = "plain"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = grid::unit(6, "pt"),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  filename = output_file,
  plot = xist_status_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

print(plot_data %>% arrange(Group, Category))
message("Saved: ", output_file)
