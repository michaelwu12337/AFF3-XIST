# AFF3 safeguards female neural lineage by restraining XIST

The earliest neural lineage develops under sex-specific dosage-compensation machinery, yet how this shapes neural fate remains poorly defined. Mutations in the super-elongation complex factor AFF3 cause neurodevelopmental disease, but only loss-of-function variants show a female bias, whereas dominant-negative variants affect both sexes. Using human stem cell-based embryo models to dissect this vulnerability, we find that AFF3 is enriched in neuroectoderm and is required for neural lineage specification in females, but not in males. Mechanistically, AFF3 prevents premature upregulation of the X-inactivation mediator XIST during neural specification in both gastruloids and cerebral organoids. Loss of AFF3 leads to excessive XIST expression that represses autosomal neural programs. Pharmacogenetic attenuation of XIST restores neural specification in AFF3 loss-of-function models. These findings establish AFF3 as a safeguard of early female neural development, highlighting its role in temporally restraining XIST to permit neural lineage specification, with AFF3 loss driving sex-biased neurodevelopmental disease.


### Download link for processed files

```text
goes here
```

### Analysis script structure

Main-figure file names follow v66. Run scripts and notebooks with their Figure-N folder as the working directory.

Prism projects retain the final graph formatting and statistical analyses for Figures 2B, 3B–C, 3F, 4E and 4G. The corresponding Prism input tables are also provided as CSV files in source_files, with stored values and row order unchanged. Prism exclusion flags and graph settings remain in the original projects. Figure 3 R/Python workflows provide the upstream analysis or earlier plot layout; use the Prism projects for the final panel layout.

Figure 7E includes the original Prism data file and extracted single-value tables; separate plotting code has not been located. The new UMAP integration and annotation workflow is pending; the existing Figure 1 script is retained and is not a verified match to v66. Existing Figure 2, 3 and 6 UMAP/lineage workflows are retained without regenerating the missing integration or annotations.

Large processed single-cell objects are not included in Git. Place them in the source_files folders listed below once available. Generated figures and intermediate tables are not included.

## Figure-1

```text
Figure-1/
└── Fig1C_1D_scRNAseq_integration_highlight.py
```

---

## Figure-2

```text
Figure-2/
├── Fig2B_gastruloid_morphology.prism
├── Fig2C_RNAseq_scatter.ipynb
├── Fig2D_neural_and_skeletal_muscle_heatmaps.ipynb
├── Fig2F_lineage_bubble_and_DEG_proportions.R
└── Fig2H_2I_scRNA_lineage_tracking.ipynb
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-2/source_files/
├── AFF3_F_KO_1_mm_counts_revised.txt
├── AFF3_F_KO_2_mm_counts_revised.txt
├── AFF3_F_WT_1_mm_counts_revised.txt
├── AFF3_F_WT_2_mm_counts_revised.txt
├── AFF3_M_KO_1_1_mm_counts_revised.txt
├── AFF3_M_KO_1_2_mm_counts_revised.txt
├── AFF3_M_WT_1_mm_counts_revised.txt
├── AFF3_M_WT_2_mm_counts_revised.txt
├── Fig2B_Female_Prism.csv
├── Fig2B_Male_Prism.csv
├── adata_combined.h5ad
├── d8_ko_1_1_mm_counts_revised.txt
├── d8_ko_1_2_mm_counts_revised.txt
├── d8_ko_2_1_mm_counts_revised.txt
├── d8_ko_2_2_mm_counts_revised.txt
├── d8_wt_1_mm_counts_revised.txt
└── d8_wt_2_mm_counts_revised.txt
```

</details>

---

## Figure-3

```text
Figure-3/
├── Fig3B_3C_XIST_expression.prism
├── Fig3B_XIST_public_dataset_bar.ipynb
├── Fig3C_XIST_timepoint_comparison.R
├── Fig3D_XIST_UMAP.ipynb
├── Fig3F_RNA_FISH_allelic_status.R
└── Fig3F_RNA_FISH_allelic_status.prism
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-3/source_files/
├── AFF3_F_KO_1_mm_counts_revised.txt
├── AFF3_F_KO_2_mm_counts_revised.txt
├── AFF3_F_WT_1_mm_counts_revised.txt
├── AFF3_F_WT_2_mm_counts_revised.txt
├── AFF3_M_KO_1_1_mm_counts_revised.txt
├── AFF3_M_KO_1_2_mm_counts_revised.txt
├── AFF3_M_WT_1_mm_counts_revised.txt
├── AFF3_M_WT_2_mm_counts_revised.txt
├── Fig3B_XIST_expression_Prism.csv
├── Fig3C_XIST_timepoint_expression_Prism.csv
├── Fig3B_public_XIST_normalized_counts.csv
├── Fig3F_RNA_FISH_allelic_percentages_Prism.csv
├── Fig3F_RNA_FISH_XIST_XACT_counts.csv
├── Fig3F_RNA_FISH_cell_percentages.csv
├── Fig3F_RNA_FISH_intensity.csv
├── adata_combined.h5ad
├── d0_F_ko_1_1_mm_counts_revised.txt
├── d0_F_ko_1_2_mm_counts_revised.txt
├── d0_F_wt_1_mm_counts_revised.txt
├── d0_F_wt_2_mm_counts_revised.txt
├── d4_F_ko_1_1_mm_counts_revised.txt
├── d4_F_ko_1_2_mm_counts_revised.txt
├── d4_F_wt_1_mm_counts_revised.txt
└── d4_F_wt_2_mm_counts_revised.txt
```

</details>

---

## Figure-4

```text
Figure-4/
├── Fig4E_gastruloid_morphology.prism
├── Fig4G_NAV1_EPHB1_expression.prism
├── Fig4H_shared_rescue_gene_heatmap.ipynb
└── Fig4I_shared_rescue_gene_GO.R
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-4/source_files/
├── AFF3_F_KO_1_mm_counts_revised.txt
├── AFF3_F_KO_2_mm_counts_revised.txt
├── AFF3_F_KO_Xt_KO_A3_1_mm_counts_revised.txt
├── AFF3_F_KO_Xt_KO_A3_2_mm_counts_revised.txt
├── AFF3_F_WT_1_mm_counts_revised.txt
├── AFF3_F_WT_2_mm_counts_revised.txt
├── Fig4E_DKO_Prism.csv
├── Fig4E_X1_Prism.csv
├── Fig4G_EPHB1_DKO_Prism.csv
├── Fig4G_EPHB1_X1_Prism.csv
├── Fig4G_NAV1_DKO_Prism.csv
├── Fig4G_NAV1_EPHB1_normalized_counts.csv
├── Fig4G_NAV1_X1_Prism.csv
├── KOX1_1_mm_counts_revised.txt
├── KOX1_2_mm_counts_revised.txt
├── KO_1_mm_counts_revised.txt
├── KO_2_mm_counts_revised.txt
├── WT_1_mm_counts_revised.txt
├── WT_2_mm_counts_revised.txt
└── XIST_binding_933_genes.csv
```

</details>

---

## Figure-5

```text
Figure-5/
├── Fig5A_clinical_heatmap.R
├── Fig5B_LoF_DN_dumbbell.R
├── Fig5D_AFF3_WT_vs_DN_scatterplot.R
├── Fig5E_AFF3_DN_GO.R
└── Fig5F_AFF3_DN_lineage_heatmaps.R
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-5/source_files/
├── Fig5A_clinical_heatmap_matrix.csv
├── Fig5B_LoF_DN_dumbbell_summary.csv
├── Fig5D_DN_1_1_mm_counts_revised.txt
├── Fig5D_DN_1_2_mm_counts_revised.txt
├── Fig5D_DN_2_1_mm_counts_revised.txt
├── Fig5D_DN_2_2_mm_counts_revised.txt
├── Fig5D_WT_1_mm_counts_revised.txt
└── Fig5D_WT_2_mm_counts_revised.txt
```

</details>

---

## Figure-6

```text
Figure-6/
├── Fig6C_6D_brain_organoid_UMAP.R
├── Fig6C_neural_lineage_percent.R
└── Fig6E_6F_Venn_GO.R
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-6/source_files/
├── AFF3_integrated_final_v3.rds
├── BO_DKO_rescued_genes_1599.csv
└── XIST_binding_933_genes.csv
```

</details>

---

## Figure-7

```text
Figure-7/
└── Fig7E_transplantation_marker_quantification.pzfx
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-7/source_files/
├── Fig7E_EOMES_single_values.csv
├── Fig7E_SOX2_single_values.csv
└── Fig7E_TUJ1_single_values.csv
```

</details>

---

## Supplementary Figure

```text
Supplementary Figure/
├── Supp-Figure-1/
│   └── SuppFig1_dotplots_heatmaps.py
├── Supp-Figure-3/
│   └── SuppFig3A_3B_3C_TemporalVAE_query_projection.ipynb
├── Supp-Figure-4/
│   └── SuppFig4B_X_linked_heatmaps.R
├── Supp-Figure-5/
│   └── SuppFig5_individual_gene_UMAPs.R
├── Supp-Figure-6/
│   └── SuppFig6_individual_gene_UMAPs.R
├── Supp-Figure-7/
│   └── SuppFig7_imprinted_gene_UMAPs.R
├── Supp-Figure-8/
│   └── SuppFig8C_8D_8E_AFF3_DN_bulk_RNAseq.R
└── Supp-Figure-9/
    ├── SuppFig9A_RNAseq_scatterplots.R
    ├── SuppFig9B_lineage_heatmaps.R
    └── SuppFig9C_9D_rescue_analysis.R
```

<details>
<summary><strong>source_files</strong></summary>

```text
Supplementary Figure/
├── Supp-Figure-4/source_files/
│   ├── AFF3_F_KO_1_mm_counts_revised.txt
│   ├── AFF3_F_KO_2_mm_counts_revised.txt
│   ├── AFF3_F_WT_1_mm_counts_revised.txt
│   ├── AFF3_F_WT_2_mm_counts_revised.txt
│   ├── AFF3_M_KO_1_1_mm_counts_revised.txt
│   ├── AFF3_M_KO_1_2_mm_counts_revised.txt
│   ├── AFF3_M_WT_1_mm_counts_revised.txt
│   ├── AFF3_M_WT_2_mm_counts_revised.txt
│   └── X_linked_HGNC_genes_Ensembl115.csv
├── Supp-Figure-8/source_files/
│   ├── D1_1_mm_counts_revised.txt
│   ├── D1_2_mm_counts_revised.txt
│   ├── D2_1_mm_counts_revised.txt
│   ├── D2_2_mm_counts_revised.txt
│   ├── D3_1_mm_counts_revised.txt
│   └── D3_2_mm_counts_revised.txt
└── Supp-Figure-9/source_files/
    ├── AFF3_F_KO_1_mm_counts_revised.txt
    ├── AFF3_F_KO_2_mm_counts_revised.txt
    ├── AFF3_F_KO_Xt_KO_A3_1_mm_counts_revised.txt
    ├── AFF3_F_KO_Xt_KO_A3_2_mm_counts_revised.txt
    ├── AFF3_F_WT_1_mm_counts_revised.txt
    ├── AFF3_F_WT_2_mm_counts_revised.txt
    ├── KOX1_1_mm_counts_revised.txt
    ├── KOX1_2_mm_counts_revised.txt
    ├── KO_1_mm_counts_revised.txt
    ├── KO_2_mm_counts_revised.txt
    ├── WT+X1_1_mm_counts_revised.txt
    ├── WT+X1_2_mm_counts_revised.txt
    ├── WT_1_mm_counts_revised.txt
    ├── WT_2_mm_counts_revised.txt
    └── XIST_binding_933_genes.csv
```

</details>

---
