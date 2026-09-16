{
  library(DESeq2)
  library(tidyverse)
  library(tidyr)
  library(tibble)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggnewscale)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
}

gene_list <- c(
  "SORCS2", "PAX7", "TTYH2", "RBFOX1", "EPHB1",
  "CDON", "RNF165", "NAV1", "PKNOX2"
)

GO_selected <- enrichGO(gene = gene_list,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pvalueCutoff = 1,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 1)
df_selected <- as.data.frame(GO_selected)
df_selected <- df_selected %>%
  mutate(
    log2_OE = log2(FoldEnrichment),
    neglog10p = -log10(p.adjust),
    Description = str_wrap(Description, width = 40),
    analysis_name = "Old/new shared rescued XIST-binding 9-gene GO analysis",
    analysis_description = "GO enrichment analysis for the 9 genes that are XIST-binding and rescued in both old and new experiments.",
    selection_criteria = "Genes are XIST-binding, shared between old and new rescued gene lists, and overlap the full XIST-binding list.",
    thresholds = "Shared rescue gene selection used p.adj < 0.05 and absolute log2FC > 0.5; GO uses clusterProfiler::enrichGO with ont = ALL and BH correction.",
    source_script = "Fig4I_shared_rescue_gene_GO.R"
  )

# Plotting ----------------------------------------------------------------

pl <- ggplot(df_selected %>% slice_max(order_by = neglog10p, n = 10, with_ties = FALSE),
             aes(x = neglog10p,
                 y = reorder(Description, neglog10p),
                 fill = log2_OE)) +
  geom_col(color = "black") +
  labs(x = "-log10(p.adj)", y = NULL, fill = "log2(O/E)") +
  scale_fill_gradient(low = "#E6EEF6",
                      high = "#2C7BB6",
                      labels = scales::number_format(accuracy = 0.1),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

dir.create("generated_figures/Fig4I", recursive = TRUE, showWarnings = FALSE)

ggsave("generated_figures/Fig4I/shared_rescued_XIST_9_gene_GO_top10.tiff",
       plot = pl,
       device = "tiff",
       dpi = 900,
       width = 10, height = 5, units = "in",
       compression = "lzw")

pl_top20 <- ggplot(df_selected %>% slice_max(order_by = neglog10p, n = 20, with_ties = FALSE),
                   aes(x = neglog10p,
                       y = reorder(Description, neglog10p),
                       fill = log2_OE)) +
  geom_col(color = "black") +
  labs(x = "-log10(p.adj)", y = NULL, fill = "log2(O/E)") +
  scale_fill_gradient(low = "#E6EEF6",
                      high = "#2C7BB6",
                      labels = scales::number_format(accuracy = 0.1),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

ggsave("generated_figures/Fig4I/shared_rescued_XIST_9_gene_GO_top20.tiff",
       plot = pl_top20,
       device = "tiff",
       dpi = 900,
       width = 10, height = 8, units = "in",
       compression = "lzw")

# Manual selected GO terms ------------------------------------------------

manual_terms <- c(
  "skeletal muscle cell differentiation",
  "skeletal muscle tissue development",
  "central nervous system neuron differentiation",
  "neural precursor cell proliferation",
  "positive regulation of nervous system development",
  "cell fate commitment"
)

df_manual <- df_selected %>%
  mutate(
    Description_clean = str_replace_all(Description, "\n", " ")
  ) %>%
  dplyr::filter(Description_clean %in% manual_terms) %>%
  mutate(
    manual_order = match(Description_clean, manual_terms),
    Description_manual_text = str_wrap(Description_clean, width = 30),
    Description_manual = factor(Description_manual_text, levels = rev(str_wrap(manual_terms, width = 30)))
  ) %>%
  arrange(manual_order)

missing_manual_terms <- setdiff(manual_terms, df_manual$Description_clean)
if (length(missing_manual_terms) > 0) {
  warning(paste("Missing manual GO terms:", paste(missing_manual_terms, collapse = ", ")))
}

pl_manual <- ggplot(df_manual,
                    aes(x = neglog10p,
                        y = Description_manual,
                        fill = log2_OE)) +
  geom_col(color = "black") +
  labs(x = "-log10(p.adj)", y = NULL, fill = "log2(O/E)") +
  scale_fill_gradient(low = "#E6EEF6",
                      high = "#2C7BB6",
                      labels = scales::number_format(accuracy = 0.1),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20)
  )

ggsave("generated_figures/Fig4I/shared_rescued_XIST_9_gene_GO_manual_selected_terms.tiff",
       plot = pl_manual,
       device = "tiff",
       dpi = 900,
       width = 11, height = 6.5, units = "in",
       compression = "lzw")

# Wider, flatter layout for the main-figure panel. The GO terms, enrichment
# values, ordering, and colour scale are identical to the original version.
wide_wrap_width <- 48

df_manual_wide <- df_manual %>%
  mutate(
    Description_manual_wide_text = str_wrap(Description_clean, width = wide_wrap_width),
    Description_manual_wide = factor(
      Description_manual_wide_text,
      levels = rev(str_wrap(manual_terms, width = wide_wrap_width))
    )
  )

pl_manual_wide <- ggplot(df_manual_wide,
                         aes(x = neglog10p,
                             y = Description_manual_wide,
                             fill = log2_OE)) +
  geom_col(color = "black") +
  labs(x = "-log10(p.adj)", y = NULL, fill = "log2(O/E)") +
  scale_fill_gradient(low = "#E6EEF6",
                      high = "#2C7BB6",
                      labels = scales::number_format(accuracy = 0.1),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20)
  )

ggsave("generated_figures/Fig4I/shared_rescued_XIST_9_gene_GO_manual_selected_terms_wide_flat.tiff",
       plot = pl_manual_wide,
       device = "tiff",
       dpi = 900,
       width = 14, height = 4.5, units = "in",
       compression = "lzw")
