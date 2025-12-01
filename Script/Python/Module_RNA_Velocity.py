# Single-cell RNA velocity analysis workflow
# Core libraries
import scanpy as sc
import anndata
from scipy import io
import numpy as np
import pandas as pd
import os


# === Color ===

color_palette = {
    "Astrocytes Like Cells": "#56017b",
    "NPC Like Cells": "#039cb1",
    "Oligodendrocytes Like Cells": "#FFA500",
    "Unknown": "#5B8C1C",
    "OPC Like Cells": "#fae831"
}


# === Data Loading ===

# Load sparse count matrix (.mtx format)
X = io.mmread("new_counts.mtx")

# Create AnnData object, transpose matrix if needed
adata = anndata.AnnData(X=X.transpose().tocsr())

# Load cell metadata (CSV)
cell_meta = pd.read_csv("new_metadata.csv")

# Load gene names
with open("new_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# Assign AnnData observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# Load and assign dimensional reduction coordinates
pca = pd.read_csv("new_pca.csv")
pca.index = adata.obs.index
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((
    adata.obs['UMAP_1'].to_numpy(),
    adata.obs['UMAP_2'].to_numpy()
)).T

# Initial UMAP plot (colored by classification)
sc.pl.umap(adata, color=['scType_Annotations'], frameon=False, save=True)

# Save AnnData dataset
adata.write('my_data.h5ad')

# === RNA Velocity Analysis ===

import scvelo as scv
import cellrank as cr

# Settings for figures and output
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

# Reload AnnData for velocity workflow
adata = sc.read_h5ad('my_data.h5ad')

# Load loom file for spliced/unspliced counts
ldata1 = scv.read('LGG275_diff.loom', cache=True)
adata = scv.utils.merge(adata, ldata1)

# Checking UMAP after merge
sc.pl.umap(adata, color='scType_Annotations', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')

# Plot: Spliced/Unspliced proportions by cell type
scv.pl.proportions(adata, groupby='scType_Annotations', save='Spliced_Unspliced_Proportion.pdf')

# Preprocess for velocity
scv.pp.filter_and_normalize(adata, enforce=True)
sc.pp.neighbors(adata, n_neighbors=30)
scv.pp.moments(adata, n_pcs=30)

# Compute RNA velocity (dynamical model)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

# Save processed dataset
adata.write('Dynamic_my_adata.h5ad')
adata = scv.read('Dynamic_my_adata.h5ad')

# Plot velocity embedding
scv.pl.velocity_embedding_stream(adata, basis='umap')
scv.pl.velocity_embedding_stream(
    adata, basis='umap',
    color=['scType_Annotations'],
    save='Dynamycal_embedding_stream.pdf',
    title='Embedding Stream'
)

# Gene-level velocity plots
genes_of_interest = ['OLIG1','OLIG2','DLL1','DLL3','CD44','APOE','CRYAB','GABBR2','KCN3','MKI67','PCNA']
scv.pl.velocity(adata, var_names=genes_of_interest, color='scType_Annotations', save='Dynamical_Genes_RNA_Velocity.png')
scv.pl.velocity(adata, var_names=genes_of_interest, color='scType_Annotations', save='Dynamical_Genes_RNA_Velocity.svg')
scv.pl.velocity(adata, var_names=genes_of_interest, color='scType_Annotations', save='Dynamical_Genes_RNA_Velocity.pdf')

# Terminal state analysis and PAGA trajectory
scv.tl.terminal_states(adata)
scv.pl.scatter(adata, color=["root_cells", "end_points"], save='Dynamical_Terminal_Points.pdf')
scv.pl.scatter(adata, color=["root_cells", "end_points"], save='Dynamical_Terminal_Points.png')

scv.tl.paga(adata, groups='scType_Annotations')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save='Dynamical_PAGA.pdf')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save='Dynamical_PAGA.png')

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save='Dynamical_Trajectory.png')
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save='Dynamical_Trajectory.pdf')

# Final save
adata.write('Dynamic_my_adata.h5ad')
