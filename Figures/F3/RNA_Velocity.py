# Single-cell RNA velocity analysis workflow
# Core libraries
import scvelo as scv
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

# Reload AnnData for velocity workflow
adata = scv.read('Dynamic_my_adata.h5ad')

# Save AnnData dataset
adata.write('my_data.h5ad')

# === RNA Velocity Analysis ===

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


# Terminal state analysis and PAGA trajectory

scv.tl.paga(adata, groups='scType_Annotations')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save='Dynamical_PAGA.png')




