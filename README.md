# LGG275 Glioma Single-Cell Analysis Suite

This repository contains reproducible scripts and workflows for single-cell RNA-seq analysis of low-grade glioma (LGG275). The suite covers quality control, cell type annotation, cell cycle analysis, cell–cell communication mapping, and RNA velocity modeling in both **R** and **Python**.

---

## Repository Contents

- `Module_QC.R`  
  Quality control and basic preprocessing with Seurat.
- `Module_Cell_Annotation.R`  
  Cell type annotation using SCType and differential expression analysis.
- `Module_Cell_cell_communication.R`  
  Cell–cell communication analysis using CellChat.
- `Module_CellCycle_ccAFvé.R`  
  Cell cycle annotation with ccAFv2 and visualization.
- `Module_RNA_Velocity.py`  
  RNA velocity workflow using Scanpy, scVelo, and CellRank.

---

## Prerequisites

- **R** (≥ 4.0) and **Python** (≥ 3.8)
- Data files: count matrix, metadata, gene names, reductions, loom files (see scripts for expected filenames)
- Recommended package installation commands:
  - R:  
    ```
    install.packages(c("Seurat", "Matrix", "dplyr", "ggplot2", "openxlsx", "pheatmap", "MAST", "cluster", "CellChat", "tidyr", "patchwork", "scales", "NMF", "ggalluvial", "ccAFv2", "vcd", "gplots", "corrplot"))
    ```
  - Python:  
    ```
    pip install scanpy anndata scipy numpy pandas scvelo cellrank
    ```

---

## Usage

1. **Prepare Data**  
   Place all raw data files in the `/data/` directory or update script paths accordingly.

2. **Run Scripts Sequentially**  
   - Start with `Module_QC.R` for initial QC and Seurat object creation.
   - Use outputs in subsequent scripts for annotation, cell cycle, communication, and velocity analysis.
   - Intermediate results and figures are saved automatically in output folders.

3. **Output**  
   - Figures, processed data objects, and results tables (CSV, RDS, H5AD, PDF/PNG/SVG).

---

## Example Directory Structure


---

## Contact

Dassou SIKA  
Bioinformatics Engineer, IGF CNRS Montpellier  
dassou.sika@igf.cnrs.fr

---

## Citation

If you use these workflows, please cite this repository and the relevant primary literature for Seurat, SCType, CellChat, ccAFv2, Scanpy, scVelo, and CellRank.

