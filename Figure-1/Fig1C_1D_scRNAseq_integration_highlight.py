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

# Old version
# 1. Load both datasets
adata_v2 = adata_combined_larry
adata_public_v3 = adata_w345w

adata11 = adata_v2
adata22 = adata_public_v3

# 2. Add dataset identifier
adata11.obs['sample_ident'] = "Gastruloids"
adata22.obs['sample_ident'] = "GW345"

# 3. Merge datasets
combined = adata11.concatenate(adata22, batch_key='sample_ident')

# 4. Enhanced preprocessing with better zero-count handling
print(f"Initial dataset shape: {combined.shape}")

# Remove cells with zero counts and genes with zero counts #avoid cell loss
sc.pp.filter_cells(combined, min_counts=0)
sc.pp.filter_genes(combined, min_counts=0)
print(f"After filtering: {combined.shape}")

# Normalization with better error handling
try:
    sc.pp.normalize_total(combined, target_sum=1e4)
except Exception as e:
    print(f"Normalization error: {e}")
    # More robust handling of potential issues
    combined = combined[combined.obs['n_counts'] > 0].copy()
    if combined.n_obs == 0:
        raise ValueError("No cells remaining after filtering zero counts")
    sc.pp.normalize_total(combined, target_sum=1e4)

# Log transformation with proper sparse matrix handling
if hasattr(combined.X, 'toarray'):
    # For sparse matrices
    combined.X = combined.X.toarray()
combined.X = np.nan_to_num(np.log1p(combined.X))

# 5. PCA for dimensionality reduction before BBKNN
print("Performing PCA...")
sc.tl.pca(combined, n_comps=100, svd_solver='arpack')

# 6. BBKNN integration with optimized parameters, #3
print("Running BBKNN integration...")
sc.external.pp.bbknn(
    combined,
    batch_key='sample_ident',
    n_pcs=50,  # Reduced from 100 to 50 for better performance
    neighbors_within_batch=3,
    metric='angular',  # Added metric for better integration
    trim=0 # No trimming for small datasets
)

# 7. Compute neighborhood graph and UMAP
print("Computing UMAP...")
sc.tl.umap(combined)


sc.pl.umap(combined, color='AFF3') #use_raw= True
print("In adata11:", 'AFF3' in adata11.var_names)
print("In adata22:", 'AFF3' in adata22.var_names)


# Method 1: Create a new column with descriptive names
combined.obs['sample_ident_new'] = combined.obs['sample_ident'].astype(str).replace({
    '0': 'hGastruloids',
    '1': 'GW345'
}).astype('category')

adata_combined.var



# Create a figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot by dataset origin
sc.pl.embedding(combined, basis='X_umap', color='sample_ident', 
                ax=ax1, show=False, title='By Dataset')
ax1.set_xlabel('UMAP1', fontsize=12)
ax1.set_ylabel('UMAP2', fontsize=12)

# Plot by cell identity
sc.pl.embedding(combined, basis='X_umap', color='sample_ident_new', 
                ax=ax2, show=False, title='By Cell Type')
ax2.set_xlabel('UMAP1', fontsize=12)
ax2.set_ylabel('UMAP2', fontsize=12)

plt.tight_layout()
plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/Integration_plot_v202603012_v1.png", dpi=300, bbox_inches='tight')
plt.show()


# 10. Enhanced label transfer for cell annotation
print("Performing label transfer...")

# Step 1: Convert 'label' to categorical
combined.obs['celltype_w345w'] = combined.obs['celltype_w345w'].astype('category')

# Step 2: Get unique non-NaN values from Tyser_label not already in categories
current_cats = combined.obs['celltype_w345w'].cat.categories
new_cats = combined.obs['celltype_v2'].dropna().unique()
to_add = [cat for cat in new_cats if cat not in current_cats]

# Step 3: Add new categories (if any exist)
if to_add:
    combined.obs['celltype_w345w'] = combined.obs['celltype_w345w'].cat.add_categories(to_add)

# Step 4: Create mask for rows to update
mask = combined.obs['celltype_w345w'].isna() & combined.obs['celltype_v2'].notna()

# Step 5: Convert values to categorical BEFORE assignment
if mask.any():
    combined.obs.loc[mask, 'celltype_w345w'] = (
        combined.obs.loc[mask, 'celltype_v2']
        .astype(combined.obs['celltype_w345w'].dtype)
    )

# 11. Additional quality checks
print("\n=== Integration Summary ===")
print(f"Final dataset shape: {combined.shape}")
print(f"Number of cells from Gastruloids: {(combined.obs['sample_ident_new'] == 'hGastruloids').sum()}")
print(f"Number of cells from GW345: {(combined.obs['sample_ident_new'] == 'GW345').sum()}")
print(f"Unique cell types in sample_ident_new: {combined.obs['sample_ident_new'].nunique()}")
print(f"Unique labels: {combined.obs['celltype_w345w'].nunique()}")

# 12. Save the integrated dataset
print("Saving integrated dataset...")
combined.write_h5ad("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/integrated_dataset_v20260312.h5ad")

print("Integration pipeline completed successfully!")



##Version#1

# Create a new UMAP plot with custom colors for Gastruloids (Red) and CS7/CS8_Stereo-seq (Blue)
plt.figure(figsize=(10, 8))

# Create custom color mapping
#custom_palette = {'0': '#D62728', '1': '#1F77B4'}
custom_palette = {'0': '#E41A1C', '1': '#377EB8'}


# Plot the UMAP with custom colors
ax = sc.pl.embedding(combined, basis='X_umap', color='sample_ident', 
                     palette=custom_palette, show=False, title='Dataset Integration: Gastruloids (Red) vs GW345 (Blue)')
ax.set_xlabel('UMAP1', fontsize=12)
ax.set_ylabel('UMAP2', fontsize=12)

# Add legend with custom colors
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#E41A1C', markersize=18, label='Gastruloids'),
           plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#377EB8', markersize=18, label='GW345')]
ax.legend(handles=handles, loc='best')

legend = ax.legend(handles=handles, loc='best', 
                   fontsize=16,
                   frameon=False,      # Remove frame
                   handletextpad=0.5,  # Reduce space between marker and text
                   borderpad=0,        # Remove border padding
                   labelspacing=0.5)   # Reduce space between labels
plt.tight_layout()
plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/Gastruloids_vs_GW345_UMAP_v20260312.png", dpi=300, bbox_inches='tight')
plt.show()

# Additional: Print some statistics about the integration
print("\n=== Integration Statistics ===")
print(f"Total cells: {combined.n_obs}")
print(f"Gastruloids cells: {(combined.obs['sample_ident'] == '0').sum()} ({(combined.obs['sample_ident'] == '0').sum()/combined.n_obs*100:.1f}%)")
print(f"GW345 cells: {(combined.obs['sample_ident'] == '1').sum()} ({(combined.obs['sample_ident'] == '1').sum()/combined.n_obs*100:.1f}%)")




###Check AFF3 expression in the combined dataset

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Get expression values to set appropriate limits
aff3_expression = combined.obs['AFF3'] if 'AFF3' in combined.obs else combined[:, 'AFF3'].X.toarray().flatten()

# Set vmin to emphasize lower expression values
custom_palette = sns.color_palette("Reds", as_cmap=True)
ax = sc.pl.embedding(combined, basis='X_umap', color='AFF3', 
                    cmap=custom_palette, 
                    vmin=0,  # Emphasize from zero
                    vmax=np.percentile(aff3_expression, 90),  # Cap at 90th percentile
                    show=False)

ax.set_xlabel('UMAP1', fontsize=12)
ax.set_ylabel('UMAP2', fontsize=12)
plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/AFF3_combined_umap_plot_v20260312.png", dpi=300, bbox_inches='tight')





###Highlight specific populations

# Create a new UMAP plot highlighting specific Tyser_label populations
plt.figure(figsize=(12, 8))

# Define the populations to highlight
highlight_populations = ["Neural Progenitors", "Early Neural Progenitors"]

# Create a custom color mapping
# First, set all cells to light gray
default_color = '#D3D3D3'  # Light gray

# Create a color map for all Tyser_label categories
all_populations = combined.obs['celltype_w345w'].cat.categories
custom_palette = {pop: default_color for pop in all_populations}

# Now assign specific colors to our highlighted populations
highlight_colors = {
    "Neural Progenitors": '#E41A1C',      # Red
    "Early Neural Progenitors": '#377EB8',        # Blue  
}

# Update the palette with highlight colors
custom_palette.update(highlight_colors)

# Plot the UMAP with custom colors
ax = sc.pl.embedding(combined, basis='X_umap', color='celltype_w345w', 
                     palette=custom_palette, show=False, 
                     )

ax.set_xlabel('UMAP1', fontsize=14)
ax.set_ylabel('UMAP2', fontsize=14)



plt.tight_layout()
plt.savefig("/mnt/davidxiang/nfs_share2/Jupyter_lab_remote/Highlighted_population_UMAP_20260312.png", 
            dpi=300, bbox_inches='tight')
plt.show()

# Print statistics about the highlighted populations
print("\n=== Highlighted Population Statistics ===")
for pop in highlight_populations:
    count = (combined.obs['celltype_w345w'] == pop).sum()
    percentage = (count / combined.n_obs) * 100
    print(f"{pop}: {count} cells ({percentage:.1f}%)")


