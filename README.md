# AFF3 safeguards female neural lineage specification by restraining XIST

Brain development spans the human lifespan, but the molecular control of the
earliest neural lineage remains poorly defined. Mutations in the super
elongation complex factor AFF3 cause KINSSHIP syndrome, a predominantly
female neurodevelopmental disorder, implicating sex-specific regulation of
early neural fate. Using human stem cell embryo models, we find that AFF3 is
enriched in neuroectoderm and is required for neural lineage specification in
female, but not male embryos. Mechanistically, AFF3 prevents premature
upregulation of the X-inactivation mediator XIST during neural specification in
both embryo models and cerebral organoids. Loss of AFF3 leads to excessive
XIST expression and ectopic repression of autosomal neural programs.
Pharmacogenetic attenuation of XIST restores neural lineage in KINSSHIP-
associated AFF3 variants. These findings identify AFF3 as a safeguard of early
neural lineage specification and reveal XIST dysregulation as a targetable axis in
sex-biased neurodevelopmental disorders.


### Download link for processed files

```text
goes here
```

### Analysis script structure

## Figure-2

```text
Figure-2/
├── Fig2B_scatter_female.R
├── Fig2B_scatter_male.R
├── Fig2C_venn_down.R
├── Fig2C_venn_up.R
├── Fig2D_GO_female.R
├── Fig2D_GO_male.R
├── Fig2E_heatmap_cardiac.R
├── Fig2E_heatmap_germlayer.R
├── Fig2E_heatmap_hematopoiesis.R
└── Fig2E_heatmap_neural.R
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
└── AFF3_M_WT_2_mm_counts_revised.txt
```

</details>

---

## Figure-3

```text
Figure-3/
├── Fig3C_3E_lineage_trace.ipynb
├── Fig3D_cluster_percentage.R
└── Fig3F_fate_percentage.R
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-3/source_files/
├── celltype_v2_top50_markers.csv
├── d8_ko_1_1_mm_counts_revised.txt
├── d8_ko_1_2_mm_counts_revised.txt
├── d8_ko_2_1_mm_counts_revised.txt
├── d8_ko_2_2_mm_counts_revised.txt
├── d8_wt_1_mm_counts_revised.txt
└── d8_wt_2_mm_counts_revised.txt
```

</details>

---

## Figure-4

```text
Figure-4/
├── Fig4A_public_dataset_violin.R
├── Fig4B_AFF3_expression_timepoint_comparison.R 
├── Fig4C_4D_scRNA.ipynb (need source file with scRNA raw QC and processing)
├── Fig4F_RNA_FISH.R (need file or code to show how to get intensity data)
├── Fig4K_heatmap.R (need to clarify where the gene list in the code is coming from)
└── Fig4L_GO.R (need to clarify where the gene_list and gene_list_GO object come from)
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
└── AFF3_F_WT_2_mm_counts_revised.txt
```

</details>

---

## Figure-5

```text
Figure-5/
├── Fig5C_5F_brain_organoid_UMAP.R
├── Fig5D_neural_lineage_percent.R (need code to show how to get the percentage, can be either R or python, may be incorporated in previous Fig5C_5F code)
└── Fig5E_5G_Venn_GO.R (need to show how to get the gene list)
```

<details>
<summary><strong>source_files</strong></summary>

```text
Figure-5/source_files/
```

</details>

---
