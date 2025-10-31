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

Object.data <- readRDS("LGG275_diff.rds")

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

# DEG Analysis ==============================================================

Idents(Object.data) <- "scType_Annotations"
all_markers <- FindAllMarkers(Object.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
For_all_markers <- FindAllMarkers(Object.data, only.pos = FALSE, test.use = "MAST")

# Overlap Visualization and Quantification ===================================

For_all_markers_filtered <- For_all_markers[For_all_markers$p_val_adj < 0.05, ]

calculate_proportion <- function(gs_list, markers_filtered, cell_type) {
  genes_gs <- gs_list[["gs_positive"]][[cell_type]]
  markers <- markers_filtered[markers_filtered$cluster == cell_type, ]
  common_genes <- intersect(genes_gs, markers$gene)
  proportion <- length(common_genes) / length(genes_gs)
  return(list(common_genes = length(common_genes),
              common_genes_list = common_genes,
              total_genes = length(genes_gs), 
              proportion = proportion))
}

cell_types <- names(gs_list[["gs_positive"]])
results <- lapply(cell_types, function(ct) calculate_proportion(gs_list, For_all_markers_filtered, ct))
names(results) <- cell_types

# Print overlap results to console
for (ct in cell_types) {
  cat(ct, ":\n")
  cat("  Common genes:", results[[ct]]$common_genes, "\n")
  cat("  Total genes:", results[[ct]]$total_genes, "\n")
  cat("  Proportion:", results[[ct]]$proportion, "\n\n")
}

# Barplot annotation overlap
plot_data <- do.call(rbind, lapply(names(results), function(ct) {
  data.frame(
    CellType = ct,
    CommonGenes = results[[ct]]$common_genes,
    TotalGenes = results[[ct]]$total_genes,
    Proportion = results[[ct]]$proportion
  )
}))

my_plot <- ggplot(plot_data, aes(x = CellType, y = TotalGenes)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_bar(aes(y = CommonGenes), stat = "identity", fill = "darkblue") +
  geom_text(aes(label = paste0(round(Proportion * 100, 1), "%"), y = CommonGenes), vjust = -0.5, color = "black", size = 3.5) +
  labs(title = "Overlap between Annotation genes and DEG Markers genes", x = "Cell Type", y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

for(ext in c("bmp","eps","pdf","png","svg","tif")) {
  ggsave(paste0( "Cells_DB/", "gene_overlap_barplot.", ext), plot = my_plot, width = 10, height = 6, units = "in", dpi = 300)
}

# Save summary overlap table
common_genes_lists <- lapply(results, function(x) x$common_genes_list)
max_length <- max(sapply(common_genes_lists, length))
padded_lists <- lapply(common_genes_lists, function(x) c(x, rep(NA, max_length - length(x))))
summary_df <- as.data.frame(padded_lists)
colnames(summary_df) <- names(results)
summary_df$Gene <- rownames(summary_df)
summary_df <- summary_df[, c("Gene", setdiff(names(summary_df), "Gene"))]
summary_df <- summary_df[rowSums(is.na(summary_df[, -1])) != ncol(summary_df) - 1, ]
write.csv(summary_df, paste( "Cells_DB/", "common_genes_summary.csv"), row.names = FALSE)

# Per-cell-type DE marker CSVs
create_cell_type_df <- function(cell_type, common_genes, markers_filtered) {
  cell_type_markers <- markers_filtered[markers_filtered$cluster == cell_type & markers_filtered$gene %in% common_genes, ]
  return(cell_type_markers)
}
for(cell_type in names(results)) {
  common_genes <- results[[cell_type]]$common_genes_list
  cell_type_df <- create_cell_type_df(cell_type, common_genes, For_all_markers_filtered)
  filename <- paste( "Cells_DB/", paste0(gsub(" ", "_", cell_type), "_differential_expression.csv"))
  write.csv(cell_type_df, filename, row.names = FALSE)
}

# End of script ==============================================================



