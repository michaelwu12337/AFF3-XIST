# AFF3 safeguards female neural lineage specification by restraining XIST

The earliest neural lineage develops under sex-specific dosage-compensation machinery, yet how this shapes neural fate remains poorly defined. Mutations in the super-elongation complex factor AFF3 cause neurodevelopmental disease, but only loss-of-function variants show a female bias, whereas dominant-negative variants affect both sexes. Using human stem cell-based embryo models to dissect this vulnerability, we find that AFF3 is enriched in neuroectoderm and is required for neural lineage specification in females, but not in males. Mechanistically, AFF3 prevents premature upregulation of the X-inactivation mediator XIST during neural specification in both gastruloids and cerebral organoids. Loss of AFF3 leads to excessive XIST expression that represses autosomal neural programs. Pharmacogenetic attenuation of XIST restores neural specification in AFF3 loss-of-function models. These findings establish AFF3 as a safeguard of early female neural development, highlighting its role in temporally restraining XIST to permit neural lineage specification, with AFF3 loss driving sex-biased neurodevelopmental disease.


### Download link for processed files

```text
goes here
```

### Analysis script structure

## Figure-1

```text
Figure-1/
└── Fig1C_1D_scRNAseq_integration_highlight.py
```

---

## Figure-2

```text
Figure-2/
├── Fig2C_RNAseq_scatter.R
├── Fig2D_neural_heatmap.R
└── Fig2F_2G_scRNA_lineage_tracking.ipynb
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
└── adata_combined.h5ad
```

</details>

---

## Figure-3

```text
Figure-3/
├── Fig3B_XIST_public_dataset_violin.R
├── Fig3C_XIST_timepoint_comparison.R
├── Fig3D_XIST_UMAP.ipynb
└── Fig3F_RNA_FISH_quantification.R
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
├── Fig3F_RNA_FISH_cell_percentages.csv
├── Fig3F_RNA_FISH_intensity.csv
├── adata_combined.h5ad
├── adata_no_W51.h5ad
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
├── Fig4G_shared_rescue_gene_heatmaps.R
└── Fig4H_shared_rescue_gene_GO.R
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
├── Fig5C_5D_brain_organoid_UMAP.R
└── Fig5E_5F_Venn_GO.R
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-5/source_files/
├── AFF3_integrated_final_v3.rds
├── BO_DKO_rescued_genes_1599.csv
└── XIST_binding_933_genes.csv
```

</details>

---

## Figure-6

```text
Figure-6/
├── Fig6A_clinical_heatmap.R
└── Fig6B_sex_dumbbell.R
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-6/source_files/
├── Fig6A_clinical_heatmap_matrix.csv
└── Fig6B_sex_dumbbell_summary.csv
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
