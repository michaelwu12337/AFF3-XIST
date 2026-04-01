# Code for AFF3 paper

### Download link for processed files
```text
goes here
```

### Analysis script structure

<details>
<summary><strong>Figure-2</strong></summary>

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

</details>

<details>
<summary><strong>Figure-3</strong></summary>

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

</details>

<details>
<summary><strong>Figure-4</strong></summary>

```text
Figure-4/
├── Fig4A_public_dataset_violin.R
├── Fig4B_AFF3_expression_timepoint_comparison.R
├── Fig4C_4D_scRNA.ipynb
├── Fig4F_RNA_FISH.R
├── Fig4K_heatmap.R
└── Fig4L_GO.R
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

</details>

<details>
<summary><strong>Figure-5</strong></summary>

```text
Figure-5/
├── Fig3C_3E_lineage_trace.ipynb
├── 02_integration.ipynb
├── 03_annotation.ipynb
├── 04_blood_trajectory.ipynb
└── 05_endo_subtype.ipynb
```

</details>
│   │   ├── AFF3_F_WT_1_mm_counts_revised.txt
│   │   └── AFF3_F_WT_1_mm_counts_revised.txt
│   ├── Fig4A_public_dataset_violin.R 
│   ├── Fig4B_AFF3_expression_timepoint_comparison.R **(change source file format to txt, currently is csv)**
│   ├── Fig4C_4D_scRNA.ipynb **(need source file with scRNA raw QC and processing)**
│   ├── Fig4F_RNA_FISH.R
│   ├── Fig4K_Heatmap **(need to clarify where the gene list in the code is coming from)**
│   └── Fig4L_GO **(need to clarify where the gene_list and gene_list_GO object)**
│
├── Figure-5/
│   ├── Fig3C_3E_lineage_trace.ipynb
│   ├── 02_integration.ipynb
│   ├── 03_annotation.ipynb
│   ├── 04_blood_trajectory.ipynb
│   └── 05_endo_subtype.ipynb
```
