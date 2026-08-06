# Figure 1C-1D: scRNA-seq integration and population highlighting

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse

adata_w345w = sc.read_h5ad('/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/adata_w345w_v3_aff3.h5ad')

sc.pl.umap(adata_w345w, color=['leiden'], show=False, title='Week stage', size=10)
plt.show()


import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from bbknn import bbknn

import pandas as pd
import os
import matplotlib

import numpy as np
import scipy.stats

import skmisc
from skmisc.loess import loess
from matplotlib.colors import LinearSegmentedColormap
adata_v5= sc.read_h5ad('/home/davidxiang/Larry_scRNA_seq_Barcode_analysis_h5ad/adata_combined.h5ad')
sc.pl.umap(adata_v5, color=["celltype_v2"], legend_loc="on data")
adata_combined = adata_v5
adata_combined_larry = adata_combined

# older integration pass kept here
# the two datasets used below
adata_v2 = adata_combined_larry
adata_public_v3 = adata_w345w

adata11 = adata_v2
adata22 = adata_public_v3

# tag where each cell came from
adata11.obs['sample_ident'] = "Gastruloids"
adata22.obs['sample_ident'] = "GW345"

# merge them before BBKNN
combined = adata11.concatenate(adata22, batch_key='sample_ident')

# zero-count cells were a problem in an older run, so check them here
print(f"Initial dataset shape: {combined.shape}")

# permissive filtering here so we arent dropping cells by accident
sc.pp.filter_cells(combined, min_counts=0)
sc.pp.filter_genes(combined, min_counts=0)
print(f"After filtering: {combined.shape}")

# normalize, with a fallback for the occasional zero-count issue
try:
    sc.pp.normalize_total(combined, target_sum=1e4)
except Exception as e:
    print(f"Normalization error: {e}")
    # fallback if normalization complains about an empty cell
    combined = combined[combined.obs['n_counts'] > 0].copy()
    if combined.n_obs == 0:
        raise ValueError("No cells remaining after filtering zero counts")
    sc.pp.normalize_total(combined, target_sum=1e4)

# log1p after dealing with a possible sparse matrix
if hasattr(combined.X, 'toarray'):
    # sparse matrix first
    combined.X = combined.X.toarray()
combined.X = np.nan_to_num(np.log1p(combined.X))

# PCA before BBKNN
print("Performing PCA...")
sc.tl.pca(combined, n_comps=100, svd_solver='arpack')

# BBKNN settings used for this integration
print("Running BBKNN integration...")
sc.external.pp.bbknn(
    combined,
    batch_key='sample_ident',
    n_pcs=50,  # 50 PCs worked better than 100 here
    neighbors_within_batch=3,
    metric='angular',  # angular distance gave the cleaner integration
    trim=0 # no trimming, these datasets are fairly small
)

# BBKNN already gives the neighbors; now make the UMAP
print("Computing UMAP...")
sc.tl.umap(combined)


sc.pl.umap(combined, color='AFF3') # quick AFF3 check; use_raw was tried before
print("In adata11:", 'AFF3' in adata11.var_names)
print("In adata22:", 'AFF3' in adata22.var_names)


# readable dataset names for plotting
combined.obs['sample_ident_new'] = combined.obs['sample_ident'].astype(str).replace({
    '0': 'hGastruloids',
    '1': 'GW345'
}).astype('category')

adata_combined.var



# two quick views of the same embedding
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# dataset source
sc.pl.embedding(combined, basis='X_umap', color='sample_ident', 
                ax=ax1, show=False, title='By Dataset')
ax1.set_xlabel('UMAP1', fontsize=12)
ax1.set_ylabel('UMAP2', fontsize=12)

# renamed dataset groups
sc.pl.embedding(combined, basis='X_umap', color='sample_ident_new', 
                ax=ax2, show=False, title='By Cell Type')
ax2.set_xlabel('UMAP1', fontsize=12)
ax2.set_ylabel('UMAP2', fontsize=12)

plt.tight_layout()
plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/Integration_plot_v202603012_v1.png", dpi=300, bbox_inches='tight')
plt.show()


# fill missing public-dataset labels from the gastruloid labels
print("Performing label transfer...")

# this has to be categorical before adding any new labels
combined.obs['celltype_w345w'] = combined.obs['celltype_w345w'].astype('category')

# labels that only occur in the gastruloid data
current_cats = combined.obs['celltype_w345w'].cat.categories
new_cats = combined.obs['celltype_v2'].dropna().unique()
to_add = [cat for cat in new_cats if cat not in current_cats]

# add any missing labels first
if to_add:
    combined.obs['celltype_w345w'] = combined.obs['celltype_w345w'].cat.add_categories(to_add)

# only fill rows without a public-dataset annotation
mask = combined.obs['celltype_w345w'].isna() & combined.obs['celltype_v2'].notna()

# categorical assignment gets unhappy otherwise
if mask.any():
    combined.obs.loc[mask, 'celltype_w345w'] = (
        combined.obs.loc[mask, 'celltype_v2']
        .astype(combined.obs['celltype_w345w'].dtype)
    )

# a few sanity checks before saving
print("\n=== Integration Summary ===")
print(f"Final dataset shape: {combined.shape}")
print(f"Number of cells from Gastruloids: {(combined.obs['sample_ident_new'] == 'hGastruloids').sum()}")
print(f"Number of cells from GW345: {(combined.obs['sample_ident_new'] == 'GW345').sum()}")
print(f"Unique cell types in sample_ident_new: {combined.obs['sample_ident_new'].nunique()}")
print(f"Unique labels: {combined.obs['celltype_w345w'].nunique()}")

# keep the integrated object for the later panels
print("Saving integrated dataset...")
combined.write_h5ad("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/integrated_dataset_v20260312.h5ad")

print("Integration pipeline completed successfully!")



## first version of the dataset comparison UMAP

# comparison UMAP: gastruloids red, CS7/CS8 Stereo-seq blue
plt.figure(figsize=(10, 8))

# colors used in the slide version
#custom_palette = {'0': '#D62728', '1': '#1F77B4'}
custom_palette = {'0': '#E41A1C', '1': '#377EB8'}


# draw the comparison UMAP
ax = sc.pl.embedding(combined, basis='X_umap', color='sample_ident', 
                     palette=custom_palette, show=False, title='Dataset Integration: Gastruloids (Red) vs GW345 (Blue)')
ax.set_xlabel('UMAP1', fontsize=12)
ax.set_ylabel('UMAP2', fontsize=12)

# legend uses the same two colors
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#E41A1C', markersize=18, label='Gastruloids'),
           plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#377EB8', markersize=18, label='GW345')]
ax.legend(handles=handles, loc='best')

legend = ax.legend(handles=handles, loc='best', 
                   fontsize=16,
                   frameon=False,      # no frame
                   handletextpad=0.5,  # keep marker close to text
                   borderpad=0,        # no extra border space
                   labelspacing=0.5)   # labels a little closer together
plt.tight_layout()
plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/Gastruloids_vs_GW345_UMAP_v20260312.png", dpi=300, bbox_inches='tight')
plt.show()

# some integration counts, mostly a sanity check
print("\n=== Integration Statistics ===")
print(f"Total cells: {combined.n_obs}")
print(f"Gastruloids cells: {(combined.obs['sample_ident'] == '0').sum()} ({(combined.obs['sample_ident'] == '0').sum()/combined.n_obs*100:.1f}%)")
print(f"GW345 cells: {(combined.obs['sample_ident'] == '1').sum()} ({(combined.obs['sample_ident'] == '1').sum()/combined.n_obs*100:.1f}%)")




### AFF3 expression on the combined UMAP

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# expression range used for color scaling
aff3_expression = combined.obs['AFF3'] if 'AFF3' in combined.obs else combined[:, 'AFF3'].X.toarray().flatten()

# keep weak expression visible
custom_palette = sns.color_palette("Reds", as_cmap=True)
ax = sc.pl.embedding(combined, basis='X_umap', color='AFF3', 
                    cmap=custom_palette, 
                    vmin=0,  # start from zero
                    vmax=np.percentile(aff3_expression, 90),  # trim the very high tail
                    show=False)

ax.set_xlabel('UMAP1', fontsize=12)
ax.set_ylabel('UMAP2', fontsize=12)
plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/AFF3_combined_umap_plot_v20260312.png", dpi=300, bbox_inches='tight')





### neural populations highlighted for the figure

# same UMAP, but only the two neural groups get color
plt.figure(figsize=(12, 8))

# the two groups used in this panel
highlight_populations = ["Neural Progenitors", "Early Neural Progenitors"]

# everything starts light gray
# all labels use that color at first
default_color = '#D3D3D3'  # Light gray

# color dictionary for all labels
all_populations = combined.obs['celltype_w345w'].cat.categories
custom_palette = {pop: default_color for pop in all_populations}

# overwrite just the two groups we need
highlight_colors = {
    "Neural Progenitors": '#E41A1C',      # Red
    "Early Neural Progenitors": '#377EB8',        # Blue  
}

# add the highlight colors
custom_palette.update(highlight_colors)

# final highlighted UMAP
ax = sc.pl.embedding(combined, basis='X_umap', color='celltype_w345w', 
                     palette=custom_palette, show=False, 
                     )

ax.set_xlabel('UMAP1', fontsize=14)
ax.set_ylabel('UMAP2', fontsize=14)



plt.tight_layout()
plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/Highlighted_population_UMAP_20260312.png", 
            dpi=300, bbox_inches='tight')
plt.show()

# counts for one last quick check
print("\n=== Highlighted Population Statistics ===")
for pop in highlight_populations:
    count = (combined.obs['celltype_w345w'] == pop).sum()
    percentage = (count / combined.n_obs) * 100
    print(f"{pop}: {count} cells ({percentage:.1f}%)")

# quick checks end here
