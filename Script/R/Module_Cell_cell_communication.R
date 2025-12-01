# Dassou SIKA
# dassou.sika@igf.cnrs.fr
# Bioinformatics Engenier
# Plasticité cérébrale, cellules souches et tumeurs gliales de bas grade / IGF
# CNRS | Montpelier-France
#
# ------------------------------------------------------------
# CellChat: CellMarkersV2 Classification Workflow
# ------------------------------------------------------------

# Required Libraries =========================================================
library(Seurat)       # Seurat object operations
library(dplyr)        # Data manipulation
library(Matrix)       # Sparse matrix support for Seurat
library(CellChat)     # Cell-cell communication analysis
library(ggplot2)      # Visualization
library(tidyr)        # pivot_longer, data reshaping
library(patchwork)    # Combining ggplots
library(scales)       # Visualization scaling
library(NMF)          # Non-negative matrix factorization for CellChat
library(ggalluvial)   # River plot for CellChat analysis

# Paths ======================================================================

cellchat_path <- paste0("CellChat/")

# Data Loading ==============================================================

Object.data <- readRDS("LGG275_diff.rds")

# CellChat Database =========================================================

CellChatDB <- CellChatDB.human # If mouse, set to CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB)

# Helper Function: Prepare Seurat Metadata for CellChat ======================

prepareSeuratData <- function(seurat_obj, primary_ident, secondary_ident) {
  Idents(seurat_obj) <- primary_ident
  data_input <- seurat_obj[["RNA"]]@layers[["data"]]
  rownames(data_input) <- rownames(seurat_obj[["RNA"]])
  labels <- Idents(seurat_obj)
  Idents(seurat_obj) <- secondary_ident
  samples <- Idents(seurat_obj)
  meta <- data.frame(labels = labels, samples = samples, row.names = names(labels))
  return(list(data_input = data_input, meta = meta))
}

# Prepare CellChat Input =====================================================

result <- prepareSeuratData(
  seurat_obj = Object.data,
  primary_ident = "scType_Annotations",
  secondary_ident = "cell_line"
)
data.input <- result$data_input
meta <- result$meta

# Initialize and Analyze CellChat ============================================

cellchat_CM_V2 <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat_CM_V2@DB <- CellChatDB.use

cellchat_CM_V2 <- subsetData(cellchat_CM_V2)
future::plan("multisession", workers = 4)
cellchat_CM_V2 <- identifyOverExpressedGenes(cellchat_CM_V2)
cellchat_CM_V2 <- identifyOverExpressedInteractions(cellchat_CM_V2)
cellchat_CM_V2 <- projectData(cellchat_CM_V2, PPI.human)
cellchat_CM_V2 <- computeCommunProb(cellchat_CM_V2, type = "triMean")
cellchat_CM_V2 <- filterCommunication(cellchat_CM_V2, min.cells = 10)
cellchat_CM_V2 <- computeCommunProbPathway(cellchat_CM_V2)
cellchat_CM_V2 <- aggregateNet(cellchat_CM_V2)
cellchat_CM_V2 <- netAnalysis_computeCentrality(cellchat_CM_V2, slot.name = "netP")

saveRDS(cellchat_CM_V2, paste0(cellchat_path, "Cellchat_CM_V2.rds"))


# Visualization: Annotation Colors ===========================================

annotation_colors <- list(
  cell_type = c(
    "Astrocytes Like Cells" = "#56017b",
    "NPC Like Cells" = "#039cb1",
    "Oligodendrocytes Like Cells" = "#FFA500",
    "Unknown" = "#5B8C1C",
    "OPC Like Cells" = "#fae831"
  )
)

groupSize <- as.numeric(table(cellchat_CM_V2@idents))

# Number/Weight of Interactions Heatmap ======================================

count_df <- as.data.frame(cellchat_CM_V2@net$count)
count_df$source <- rownames(count_df)
count_df <- count_df %>% pivot_longer(cols = -source, names_to = "target", values_to = "count")

weight_df <- as.data.frame(cellchat_CM_V2@net$weight)
weight_df$source <- rownames(weight_df)
weight_df <- weight_df %>% pivot_longer(cols = -source, names_to = "target", values_to = "weight")

ggplot_count <- ggplot(count_df, aes(x = source, y = target, fill = count)) +
  geom_tile() +
  geom_text(aes(label = count), color = "#8e44ad", size = 5) +
  scale_fill_gradient(low = "#eaecee", high = "#283747") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid = element_blank()) +
  labs(fill = "Number of Interactions") +
  coord_fixed() +
  scale_x_discrete(limits = names(sort(groupSize, decreasing = TRUE))) +
  scale_y_discrete(limits = names(sort(groupSize, decreasing = TRUE)))

ggplot_weight <- ggplot(weight_df, aes(x = source, y = target, fill = weight)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", weight)), color = "#8e44ad", size = 5) +
  scale_fill_gradient(low = "#eaecee", high = "#283747") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid = element_blank()) +
  labs(fill = "Weight of Interactions") +
  coord_fixed() +
  scale_x_discrete(limits = names(sort(groupSize, decreasing = TRUE))) +
  scale_y_discrete(limits = names(sort(groupSize, decreasing = TRUE)))

ggsave(paste0(cellchat_path, "Heatmap_of_interactions.png"),
       plot = ggplot_count + ggplot_weight, width = 20, height = 8)

# Circle Network Visualization ===============================================

Cp1 <- netVisual_circle(cellchat_CM_V2@net$count, color.use = annotation_colors$cell_type, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
Cp2 <- netVisual_circle(cellchat_CM_V2@net$weight, color.use = annotation_colors$cell_type, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")

png(paste0(cellchat_path, "Number_of_interactions.png"), width = 2000, height = 800)
print(Cp1)
print(Cp2)
dev.off()

# Pathways Analysis ==========================================================

pathways.show.all <- cellchat_CM_V2@netP$pathways
write.csv(pathways.show.all, paste0(cellchat_path, "Pathways.csv"))
Cp3 <- netAnalysis_contribution(cellchat_CM_V2, signaling = pathways.show.all)
saveRDS(Cp3, paste0(cellchat_path, "Ligand_Receptor_Contribution.rds"))

# Signaling Role Heatmaps ====================================================

png(paste0(cellchat_path, "Signaling_Roles_Out.png"), width = 800, height = 800)
netAnalysis_signalingRole_heatmap(cellchat_CM_V2, pattern = "outgoing", color.use = annotation_colors$cell_type)
dev.off()

png(paste0(cellchat_path, "Signaling_Roles_In.png"), width = 800, height = 800)
netAnalysis_signalingRole_heatmap(cellchat_CM_V2, pattern = "incoming", color.use = annotation_colors$cell_type)
dev.off()

# Communication Patterns =====================================================

nPatterns <- 3
cellchat_CM_V2 <- identifyCommunicationPatterns(cellchat_CM_V2, pattern = "outgoing", k = nPatterns)
cellchat_CM_V2 <- identifyCommunicationPatterns(cellchat_CM_V2, pattern = "incoming", k = nPatterns)

Cp6 <- netAnalysis_river(cellchat_CM_V2, pattern = "outgoing", color.use = annotation_colors$cell_type, font.size = 7, font.size.title = 40)
ggsave(paste0(cellchat_path, "Net_Analysis_River_Out.png"), Cp6, width = 20, height = 20)

Cp7 <- netAnalysis_river(cellchat_CM_V2, pattern = "incoming", color.use = annotation_colors$cell_type, font.size = 7, font.size.title = 40)
ggsave(paste0(cellchat_path, "Net_Analysis_River_In.png"), Cp7, width = 20, height = 20)

# Ligand-Receptor Interactions ===============================================

pdf(paste0(cellchat_path, "ligand_receptor_interactions_net.pdf"), width = 17, height = 17)
netVisual_chord_gene(cellchat_CM_V2, slot.name = "net", color.use = annotation_colors$cell_type, title.name = "Ligand-Receptor Interactions")
dev.off()

for (cell in c("Astrocytes Like Cells", "Oligodendrocytes Like Cells", "NPC Like Cells")) {
  pdf(paste0(cellchat_path, "ligand_receptor_interactions_net_", gsub(" ", "", cell), ".pdf"), width = 17, height = 17)
  netVisual_chord_gene(cellchat_CM_V2, slot.name = "net",
                       sources.use = c(cell),
                       targets.use = c("Astrocytes Like Cells", "NPC Like Cells", "Oligodendrocytes Like Cells", "Unknown"),
                       color.use = annotation_colors$cell_type,
                       title.name = "Ligand-Receptor Interactions")
  dev.off()
}

# End of script ==============================================================





