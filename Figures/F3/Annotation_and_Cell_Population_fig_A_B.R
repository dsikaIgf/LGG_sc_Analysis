# Dassou SIKA
# dassou.sika@igf.cnrs.fr
# Bioinformatics Engenier
# Plasticité cérébrale, cellules souches et tumeurs gliales de bas grade / IGF
# CNRS | Montpelier-France
#
# ------------------------------------------------------------
# Rennotation: Cell Type Annotation and DEG Analysis
# ------------------------------------------------------------

# Libraries ==================================================================
library(knitr)        # opts_chunk$set, chunk options
library(Seurat)       # Core analysis: objects, plots, marker analysis
library(dplyr)        # Data manipulation (group_by, slice_max, etc.)
library(Matrix)       # Sparse matrix handling for Seurat
library(ggplot2)      # Plots and visualizations
library(openxlsx)     # Reading in Excel-based annotation databases
library(pheatmap)     # Heatmap visualization
library(MAST)         # Differential expression analysis in FindAllMarkers
library(cluster)      # Clustering support (if needed)

# Data Loading ==============================================================

Object.data <- readRDS("LGG275_diff.rds") # Example of Data

# SCType dependencies & gene sets ============================================

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ <- paste(wrk_path, Sc_customDB, sep = "")
tissue <- "Cells_DB"
gs_list <- gene_sets_prepare(db_, tissue)

# Seurat object version check ===============================================

seurat_package_v5 <- isFALSE('counts' %in% names(attributes(Object.data[["RNA"]])))
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(Object.data[["RNA"]]$scale.data) else as.matrix(Object.data[["RNA"]]@scale.data)

es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = NULL)

cL_resutls <- do.call("rbind", lapply(unique(Object.data@meta.data$RNA_snn_res.1.5), function(cl){
  es.max.cl = sort(rowSums(es.max[ , rownames(Object.data@meta.data[Object.data@meta.data$RNA_snn_res.1.5 == cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Object.data@meta.data$RNA_snn_res.1.5 == cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < 0] <- "Unknown"
print(sctype_scores[,1:3])

# Save scores ==============================================================

dir.create(paste( "Cells_DB", sep = ""), showWarnings = FALSE)
write.csv(cL_resutls, paste( "Cells_DB/", "cL_resutls.tsv", sep = ""))
write.csv(sctype_scores, paste( "Cells_DB/", "sctype_scores.tsv", sep = ""))


# Update metadata with Cluster Classification ===============================

Object.data@meta.data$scType_Annotations = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster == j, ]
  Object.data@meta.data$scType_Annotations[Object.data@meta.data$RNA_snn_res.1.5 == j] = as.character(cl_type$type[1])
}


# ========================================================================
# Cell Population
# ========================================================================

# Colors (reusing your annotation_colors)
annotation_colors <- list(
  cell_type = c(
    "Astrocytes Like Cells"       = "#56017b",
    "Astro Like Cells"            = "#56017b",
    "NPC Like Cells"              = "#039cb1",
    "Oligodendrocytes Like Cells" = "#FFA500",
    "Oligo Like Cells"            = "#FFA500",
    "Unknown"                     = "#5B8C1C",
    "OPC Like Cells"              = "#fae831"
  )
)

# Example 1: Proportions by stim and scType_Annotations (pie charts)
# Make sure you have a 'stim' column in meta.data as in your other script.
ct_table <- table(Object.data$stim, Object.data$scType_Annotations)

# Optionally rename rows, e.g. "diff" -> "-GF" if applicable
rownames(ct_table)[rownames(ct_table) == "diff"] <- "-GF"

# Row-wise proportions (per stim)
ct_prop_row <- prop.table(ct_table, margin = 1)

# Base R pie charts per stim
par(mfrow = c(1, nrow(ct_prop_row)))
for (i in 1:nrow(ct_prop_row)) {
  proportions <- ct_prop_row[i, ]
  cell_types  <- names(proportions)
  colors      <- annotation_colors$cell_type[cell_types]
  pie(
    proportions,
    main = rownames(ct_prop_row)[i],
    col  = colors
  )
}
par(mfrow = c(1, 1))

# Example 2: Save ggplot2 pie charts by stim
ct_prop_df <- as.data.frame(as.table(ct_prop_row))
colnames(ct_prop_df) <- c("Stim", "CellType", "Proportion")

dir.create(paste("scType_plots", sep = ""), showWarnings = FALSE)

for (stim in unique(ct_prop_df$Stim)) {
  subset_df <- ct_prop_df %>% dplyr::filter(Stim == stim)
  
  p <- ggplot(subset_df, aes(x = "", y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = annotation_colors$cell_type) +
    theme_void() +
    ggtitle(stim)
  ggsave(filename = paste("scType_plots/",
                          paste0("Piechart_scType_Proportion_", stim, ".pdf")), plot = p)
  ggsave(filename = paste("scType_plots/",
                          paste0("Piechart_scType_Proportion_", stim, ".png")), plot = p)
  
}


# Save Data
saveRDS(diff.data, "LGG275_diff.rds")


