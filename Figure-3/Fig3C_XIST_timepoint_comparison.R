{
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(dplyr)
  library(DESeq2)
  library(ggbreak)
}

# DESeq2 ------------------------------------------------------------------
{
  path <- "source_files"
  files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)
  
  for (f in files) {
    name <- tools::file_path_sans_ext(basename(f))
    assign(name, read.table(f, header = TRUE))
  }
  
  old_names <- ls(pattern = "^(AFF3_|d8_)")
  
  new_names <- old_names %>%
    gsub("^AFF3_", "d8_", .) %>%   # AFF3 → d8
    gsub("_f_", "_F_", .) %>%      # f → F
    gsub("_m_", "_M_", .) %>%      # m → M
    gsub("_KO_", "_ko_", .) %>%    # KO → ko
    gsub("_WT_", "_wt_", .)        # WT → wt
  
  for (i in seq_along(old_names)) {
    if (old_names[i] != new_names[i]) {
      assign(new_names[i], get(old_names[i]))
      rm(list = old_names[i])
    }
  }
}
{
  df_list <- mget(ls(pattern = "_counts_revised$"))
  
  counts_Data <- do.call(
    cbind,
    lapply(df_list, function(df) df$counts)
  )
  
  colnames(counts_Data) <- names(df_list)
  rownames(counts_Data) <- df_list[[1]]$Geneid
}
{
  sample_names <- colnames(counts_Data)
  
  col_Data <- data.frame(
    Time_point = as.numeric(sub("^d([0-9]+)_.*", "\\1", sample_names)),
    Condition  = toupper(sub("^d[0-9]+_[FM]_(wt|ko).*", "\\1", sample_names)),
    row.names  = sample_names,
    stringsAsFactors = FALSE
  )
  
  col_Data$Condition <- factor(col_Data$Condition, levels = c("WT", "KO"))
}
{
  all(colnames(counts_Data) %in% rownames(col_Data))
  all(colnames(counts_Data) == rownames(col_Data))
}
{
  dds <- DESeqDataSetFromMatrix(countData = counts_Data,
                                colData = col_Data,
                                design = ~ Condition)
  dds$Condition <- relevel(dds$Condition, ref = "WT")
  dds <- DESeq(dds)
}
{
  vst_dds <- vst(dds, blind=FALSE)
  XIST_vst <- assay(vst_dds)["XIST", , drop = FALSE]
  XIST_vst <- as.data.frame(t(XIST_vst))
  colnames(XIST_vst) <- "XIST_vst"
  
  rlog_dds <- rlog(dds, blind=FALSE)
  XIST_rlog <- assay(rlog_dds)["XIST", , drop = FALSE]
  XIST_rlog <- as.data.frame(t(XIST_rlog))
  colnames(XIST_rlog) <- "XIST_rlog"
  
  raw_counts <- counts(dds, normalized = FALSE)
  XIST_raw <- raw_counts["XIST", , drop = FALSE]
  XIST_raw <- as.data.frame(t(XIST_raw))
  colnames(XIST_raw) <- "XIST_raw"
  
  norm_counts <- counts(dds, normalized = TRUE)
  XIST_norm <- norm_counts["XIST", , drop = FALSE]
  XIST_norm <- as.data.frame(t(XIST_norm))
  colnames(XIST_norm) <- "XIST_norm"
  
  XIST_raw_log2  <- log2(XIST_raw  + 2)
  colnames(XIST_raw_log2) <- "XIST_raw_log2"
  
  XIST_norm_log2 <- log2(XIST_norm + 2)
  colnames(XIST_norm_log2) <- "XIST_norm_log2"
  
  XIST_norm_log10 <- log10(XIST_norm + 2)
  colnames(XIST_norm_log10) <- "XIST_norm_log10"
  
  all(rownames(XIST_raw) == rownames(XIST_norm))
  all(rownames(XIST_raw) == rownames(XIST_vst))
  all(rownames(XIST_raw) == rownames(XIST_rlog))
}
{
  XIST_all <- cbind(
    XIST_raw,
    XIST_raw_log2,   # log2(raw + 1)
    XIST_norm,
    XIST_norm_log2,  # log2(norm + 1)
    XIST_norm_log10, 
    XIST_vst,
    XIST_rlog
  )
  XIST_all_export <- data.frame(
    Sample = rownames(XIST_all),
    Gene = "XIST",
    XIST_all,
    row.names = NULL
  )
}
{
  xist_data <- XIST_all_export %>%
    dplyr::select(-Gene) %>%
    mutate(
      Time_point = case_when(
        grepl("^d0_", Sample)      ~ 0,
        grepl("^d4_", Sample)      ~ 4,
        grepl("^d8_M_", Sample)    ~ 8,
        grepl("^d8_F_", Sample)    ~ 12
      ),
      Condition = case_when(
        grepl("_ko_", Sample, ignore.case = TRUE) ~ "KO",
        grepl("_wt_", Sample, ignore.case = TRUE) ~ "WT"
      )
    ) %>%
    dplyr::select(Sample, Time_point, Condition, everything()) %>%
    arrange(Time_point, Condition)
  
  xist_data <- xist_data %>%
    mutate(
      Time_point = factor(as.character(Time_point),
                          levels = c("0", "4", "8", "12")),
      Condition  = factor(Condition, levels = c("WT", "KO")),
      GroupID    = paste0(Condition, "_", Time_point),
      
      # WT day 8 has the blue fill used in the figure
      FillGroup = ifelse(GroupID == "WT_8", "WT_8_special", as.character(Condition)),
      
      # four legend groups, kept seperate for fill and outline
      LegendGroup = case_when(
        FillGroup == "WT" ~ "WT_filled",                                 # orange filled
        FillGroup == "WT_8_special" ~ "WT_special_filled",               # blue filled
        FillGroup == "KO" & grepl("_8$", GroupID) ~ "KO_blue_outline",   # blue outline
        FillGroup == "KO" ~ "KO_orange_outline"                          # orange outline
      ),
      
      # this order also controls the dodge
      LegendGroup = factor(
        LegendGroup,
        levels = c("WT_filled",                 # F WT
                   "KO_orange_outline",         # F KO
                   "WT_special_filled",         # M WT
                   "KO_blue_outline")           # M KO
      )
    )
}

# Bar plot ----------------------------------------------------------------
{
  fill_vals <- c(
    WT_filled          = "#D8B365",  # F WT (BrBG brown)
    WT_special_filled  = "#5AB4AC",  # M WT (BrBG teal)
    KO_orange_outline  = NA,         # KO outline only
    KO_blue_outline    = NA
  )
  
  col_vals <- c(
    WT_filled          = "#A67C52",  # darker brown outline
    WT_special_filled  = "#3C8F8B",  # darker teal outline
    KO_orange_outline  = "#A67C52",  # F KO outline
    KO_blue_outline    = "#3C8F8B"   # M KO outline
  )
  
  legend_labels <- c(
    WT_filled          = "F WT",
    KO_orange_outline  = "F KO",
    WT_special_filled  = "M WT",
    KO_blue_outline    = "M KO"
  )
}

## mean and SE for the bar heights
{
  bar_df <- xist_data %>%
    group_by(Time_point, LegendGroup) %>%
    summarise(
      mean = mean(XIST_norm, na.rm = TRUE),
      se   = sd(XIST_norm, na.rm = TRUE) / sqrt(sum(!is.na(XIST_norm))),
      .groups = "drop"
    )
  
  ## bar version used in the final figure
  m <- ggplot(
    bar_df,
    aes(x = Time_point, y = mean, fill = LegendGroup, color = LegendGroup)
  ) +
    geom_col(
      width = 0.4,
      position = position_dodge(width = 0.8),
      alpha = 0.6
    ) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.12,
      position = position_dodge(width = 0.8),
      color = "black"
    ) +
    scale_fill_manual(
      values = fill_vals,
      labels = legend_labels,
      name   = "Condition"
    ) +
    scale_color_manual(
      values = col_vals,
      labels = legend_labels,
      name   = "Condition"
    ) +
    scale_x_discrete(
      labels = c("0" = "EPSC (0)",
                 "4" = "4",
                 "8" = "8",
                 "12" = "8")
    ) +
    labs(
      title = "",
      x = "Time point (Days)",
      y = "XIST Expression \n(WT vs KO)"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 16),
      legend.text = element_text(size = 16, colour = "black"),
      legend.title= element_text(size = 16, colour = "black")
    )
  
  ## stats still use the individual values in xist_data
  {
    stat_test_within <- xist_data %>%
      group_by(Time_point) %>%
      t_test(XIST_norm ~ Condition) %>%
      add_significance("p") %>%
      add_xy_position(x = "Time_point", dodge = 0.8)
    
    stat_test_within$y.position <- 680
    
    m <- m +
      stat_pvalue_manual(
        stat_test_within,
        label = "p.signif",
        tip.length = 0,
        label.vjust = -1
      )
  }
  
  ## keep the same gapped y axis
  m_broken <- m +
    scale_y_break(c(20, 160), scales = 0.7, space = 0.3) +
    scale_y_break(c(220, 550), scales = 0.8, space = 0.3) +
    scale_y_continuous(
      limits = c(0, 700),          # zero stays as the common baseline
      breaks = c(
        seq(0, 20, by = 5),
        seq(160, 220, by = 20),
        seq(550, 700, by = 50)
      ),
      expand = c(0, 0),
      oob = scales::squish
    ) +
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right  = element_blank(),
      axis.text.x.top    = element_blank(),
      axis.ticks.x.top   = element_blank(),
      axis.line.x.top    = element_blank()
    ) +
    stat_pvalue_manual(
      stat_test_within,
      label = "p.signif",
      tip.length = 0,
      label.vjust = -1
    )
  
  m_broken <- m_broken +
    geom_point(
      data = xist_data,
      aes(
        x = Time_point,
        y = XIST_norm,
        color = LegendGroup,   # same outline colors as the bars
        shape = Condition,     # WT vs KO shape
        group = LegendGroup    # dodge follows LegendGroup here
      ),
      position = position_jitterdodge(
        dodge.width  = 0.8,    # same width as the bars or points shift
        jitter.width = 0.10,
        jitter.height = 0
      ),
      size = 2.8,
      alpha = 0.85,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_shape_manual(values = c(WT = 17, KO = 18))
}

ggsave(
  filename = "Fig3C_XIST_timepoint_comparison.svg",
  plot     = m_broken,
  device   = svglite::svglite,
  width    = 7,
  height   = 4.5,
  units    = "in"
)
