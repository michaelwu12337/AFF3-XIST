# Supplementary Figure 1: dot plots and heatmaps

### Mouse embryo dataset dotplot

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import pandas as pd
import os
import matplotlib
import numpy as np


###
from bbknn import bbknn

import palantir
import scanpy as sc
import pandas as pd
import os
import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import seaborn as sns

adata_vv1 = sc.read_h5ad('GSM9046243_Embryo_E7.5_stereo_rep1.h5ad')
adata_vv2 = sc.read_h5ad('GSM9046244_Embryo_E7.5_stereo_rep2.h5ad')
adata_vv3 = sc.read_h5ad('GSM9046245_Embryo_E7.75_stereo_rep1.h5ad')
adata_vv4 = sc.read_h5ad('GSM9046246_Embryo_E7.75_stereo_rep2.h5ad')
adata_vv5 = sc.read_h5ad('GSM9046247_Embryo_E8.0_stereo_rep1.h5ad')
adata_vv6 = sc.read_h5ad('GSM9046248_Embryo_E8.0_stereo_rep2.h5ad')


# run the same preprocessing over all six embryos
adata_list = [adata_vv1, adata_vv2, adata_vv3, adata_vv4, adata_vv5, adata_vv6]

for i, adata in enumerate(adata_list, start=1):
    print(f"Processing adata_vv{i}...")
    
    # normalize only when SCT is not already there
    if 'SCT' not in adata.layers:
        print(f"Performing basic normalization for adata_vv{i}...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    
    # pick variable genes if this object does not have them yet
    if 'highly_variable' not in adata.var.columns:
        print(f"Identifying highly variable genes for adata_vv{i}...")
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]
    
    # scaling before PCA
    print(f"Scaling data for adata_vv{i}...")
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    print(f"Running PCA for adata_vv{i}...")
    sc.tl.pca(adata, svd_solver='arpack')
    
    # neighbors for the UMAP
    print(f"Computing neighborhood graph for adata_vv{i}...")
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
    
    # UMAP embedding
    print(f"Generating UMAP embedding for adata_vv{i}...")
    sc.tl.umap(adata)
    
    print(f"Finished processing adata_vv{i}\n")






# Aff3/Ell3 dotplots for each embryo
adata_list = [adata_vv1, adata_vv2, adata_vv3, adata_vv4, adata_vv5, adata_vv6]
genes_of_interest = ["Aff3", "Ell3"]

for i, adata in enumerate(adata_list, start=1):
    print(f"Creating dotplot for adata_vv{i}...")
    
    # dotplot for this embryo
    sc.pl.dotplot(adata, var_names=genes_of_interest, groupby='germ_layer', show=False, figsize=(4, 4))
    
    # keep a seperate file for each embryo
    plt.savefig(f"/Users/yangxiang/Desktop/Dotplot_germ_layer_vv{i}.png", dpi=300, bbox_inches='tight')
    
    # close it before the next one, otherwise memory piles up
    plt.close()
    
    print(f"Dotplot saved for adata_vv{i}")

print("All dotplots have been generated successfully!")





#For mouse datasets
#Mouse digital embryo dataset #vv1 male #vv4 male, others are female ##v3
#Compare vv3-female with vv4-male

adata_list = [adata_vv1, adata_vv2, adata_vv3, adata_vv4, adata_vv5, adata_vv6]
genes_of_interest = ["Aff3","Mllt3","Ell3","Paf1","Tcea2"]

for i, adata in enumerate(adata_list, start=1):
    print(f"Creating dotplot for adata_vv{i}...")
    
    # dotplot for this embryo
    sc.pl.dotplot(adata, var_names=genes_of_interest, groupby='germ_layer', show=False, figsize=(4, 4))
    
    # one output per embryo
    plt.savefig(f"/Users/yangxiang/Desktop/Dotplot_EF_vv{i}.png", dpi=300, bbox_inches='tight')
    
    # close before moving on
    plt.close()
    
    print(f"Dotplot saved for adata_vv{i}")

print("All dotplots have been generated successfully!")






#### Gastruloid Dotplot


genes_of_interest = ["AFF3","MLLT3","ELL3","PAF1","TCEA2","SUPT6H"]
sc.pl.dotplot(adata_combined, var_names=genes_of_interest, groupby='celltype_v2', show=False)
ax.set_xlabel('UMAP1', fontsize=24)
ax.set_ylabel('UMAP2', fontsize=24)
plt.savefig('EF_dotplot_v2.tiff', 
            format='tiff', 
            dpi=600, 
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none')
plt.tight_layout()
plt.show()





#### Human Embryo Dotplot

genes_of_interest = ["AFF3","MLLT3","ELL3","PAF1","TCEA2","SUPT6H"]

#Dotplot
sc.pl.dotplot(adata_w345w, genes_of_interest, groupby="celltype_w345w", standard_scale="var",show=False)

# grab the current figure and tighten it before saving
plt.gcf().tight_layout()

# save this one as TIFF
plt.savefig('adata_w345w_EF_Dot_dotplot_v1.tiff', 
            format='tiff', 
            dpi=600, 
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none')
plt.show()





#### Human embryo & Gastruloids Heatmaps (Pearson correlation)
# combined is the integrated object from Figure 1

#Revised correlation analysis heatmap-version2

combined_Amnion = combined[combined.obs['celltype_w345w'] == 'Amnion'].copy()
combined_Neural_Progenitors = combined[combined.obs['celltype_w345w'] == 'Neural Progenitors'].copy()
combined_VEC = combined[combined.obs['celltype_w345w'] == 'VEC'].copy()
#combined_Cardiac_Mesoderm = combined[combined.obs['celltype_w345w'] == 'Cardiac Mesoderm'].copy()
combined_Blood_Progenitors = combined[combined.obs['celltype_w345w'] == 'Blood Progenitors'].copy()

combined_Early_Neural_Progenitors = combined[combined.obs['celltype_w345w'] == 'Early Neural Progenitors'].copy()
combined_Vascular_Endothelial_cells = combined[combined.obs['celltype_w345w'] == 'Vascular Endothelial cells'].copy()
combined_Amniotic_Epithelial_Cells = combined[combined.obs['celltype_w345w'] == 'Amniotic Epithelial Cells'].copy()
#combined_Lateral_plate_mesoderm_Cardiac_mesoderm= combined[combined.obs['celltype_w345w'] == 'Lateral plate mesoderm/Cardiac mesoderm'].copy()
combined_Immune_cells = combined[combined.obs['celltype_w345w'] == 'Immune cells'].copy()

# put the selected cell types back together
combined_merged = ad.concat([combined_Amnion, combined_Neural_Progenitors, combined_VEC, combined_Blood_Progenitors, combined_Early_Neural_Progenitors, combined_Vascular_Endothelial_cells, combined_Amniotic_Epithelial_Cells,combined_Immune_cells], join="outer")



# average expression for each cluster in the merged subset
# this second pass uses the merged object above
cluster_profiles = pd.DataFrame()
for cluster in combined_merged.obs['celltype_w345w'].unique():
    # mean expression within this cluster
    cluster_cells = combined_merged[combined_merged.obs['celltype_w345w'] == cluster]
    mean_expression = cluster_cells.to_df().mean(axis=0) # .to_df() gets the expression matrix
    cluster_profiles[cluster] = mean_expression

# rows are genes and columns are clusters now

# correlation between the cluster profiles
correlation_matrix = cluster_profiles.corr(method='pearson') # spearman was another option here

# heatmap
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# colors used for this heatmap
colors = ['#3170AB', '#527FB7', '#E3E9F2', '#F7DED9', '#BA3330']
positions = [0, 0.30, 0.50, 0.60, 1.0]
custom_cmap = LinearSegmentedColormap.from_list('custom_blue_white_red', 
                                               list(zip(positions, colors)))
#BA3330
#F7DED9
#EDF1F8
#8FA8CD
#3170AB


# revised heatmap version
plt.figure(figsize=(10, 8))
heatmap = sns.heatmap(correlation_matrix,
                      vmin=0, vmax=1,  # 0-1 is useful here; an older try used -1 to 1
                      cmap=custom_cmap,
                      annot=False,
                      fmt=".2f",
                      square=True,
                      annot_kws={"size": 15},
                      cbar_kws={
                          'ticks': [0,1],
                          'shrink': 0.5,
                          'aspect': 8,
                          'anchor': (0.0, 0.0),
                          'panchor': (1.0, 0.0)
                      })

# larger labels to match the slide
heatmap.set_xticklabels(heatmap.get_xticklabels(), fontsize=28)
heatmap.set_yticklabels(heatmap.get_yticklabels(), fontsize=28)

# colorbar sizing was done seperately
# this bit is only for the colorbar
cbar = heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)

# move the colorbar by hand, it was easier this way
cbar.ax.set_position([0.80, 0.15, 0.15, 0.3])  # [left, bottom, width, height]

plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/Pearson_correlation_heatmap_v20260307_notitle.png", 
            dpi=300, bbox_inches='tight')
plt.show()

# heatmap part ends here
# thats all for this file
