# Dassou SIKA
# dassou.sika@igf.cnrs.fr
# Bioinformatics Engenier
# Plasticité cérébrale, cellules souches et tumeurs gliales de bas grade / IGF
# CNRS | Montpelier-France
#
# ------------------------------------------------------------------
# LGG275_diff: DATA INTEGRATION AND PROCESSING
# ------------------------------------------------------------------


# Libraries
library(Seurat)
library(Matrix)


# Data: Diff condition only
diff.data <- Read10X("raw_feature_bc_matrix")

diff.data <- CreateSeuratObject(counts = diff.data, project = "LGG275.diff", min.cells = 3, min.features = 200)
diff.data$stim <- "diff"
diff.data$cell_line <- "LGG275"
diff.data[["percent.mt"]] <- PercentageFeatureSet(diff.data, pattern = "^MT-")
diff.data[["percent.ribo"]] <- PercentageFeatureSet(diff.data, pattern = "^RP")
diff.data@meta.data[["cell_name"]] <- paste(rownames(diff.data@meta.data), sep = "")
diff.data <- subset(diff.data, subset = nFeature_RNA > 200 & percent.mt < 20 & nCount_RNA < 60000 & nFeature_RNA < 9000)
diff.data <- NormalizeData(diff.data, verbose = FALSE)
diff.data <- FindVariableFeatures(diff.data, selection.method = "vst", nfeatures = 2000)

# PCA
diff.data <- ScaleData(diff.data, verbose = TRUE)
diff.data <- RunPCA(diff.data, npcs = 30, verbose = TRUE)

# PCA Plots
VizDimLoadings(diff.data, dims = 1:2, reduction = "pca")
DimPlot(diff.data, reduction = "pca") + NoLegend()
DimHeatmap(diff.data, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(diff.data)

# UMAP reduction and clustering
diff.data <- RunUMAP(diff.data, reduction = "pca", dims = 1:15)
diff.data <- FindNeighbors(diff.data, reduction = "pca", dims = 1:15)
diff.data <- FindClusters(diff.data, resolution = 0.3)


# Save Data
saveRDS(diff.data, "LGG275_diff.rds")






