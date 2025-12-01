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

# ========================================================================
# DEG Analysis
# ========================================================================

deg_path <- paste(out_path, "DEGs/", sep = "")

# Part 1: Run DEG with MAST per scType_Annotations ------------------------

Idents(Object.data) <- "scType_Annotations"

# Positive markers per cell type
all_markers <- FindAllMarkers(
  Object.data,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "MAST"
)  

# All markers (both up and down) per cell type
For_all_markers <- FindAllMarkers(
  Object.data,
  only.pos = FALSE,
  test.use = "MAST"
)  

# Optional: save raw marker tables
dir.create(paste0(deg_path), showWarnings = FALSE)
saveRDS(all_markers, file = file.path(deg_path, "All_Positive_Markers_MAST.rds"))
saveRDS(For_all_markers, file = file.path(deg_path, "All_Markers_MAST.rds"))

# Optional: filter significant & remove unwanted gene patterns [web:27]
significant_markers <- For_all_markers[For_all_markers$p_val_adj < 0.05, ]
filtered_markers <- significant_markers[
  !grepl("^(LINC|ENSG|MIR)", significant_markers$gene) &
    !grepl("-AS1$", significant_markers$gene),
]


saveRDS(filtered_markers, file = file.path(deg_path, "All_filtered_Markers.rds"))


# Part 2: Curated markers + DotPlot ---------------------------------------

# If you already have filtered_markers in memory from above, you can skip this line
# filtered_markers <- readRDS(file.path(deg_path, "All_filtered_Markers.rds"))

filtered_markers$color <- "#1707AB"

AC_df      <- filtered_markers[filtered_markers$cluster == "Astrocytes Like Cells", ]
OC_df      <- filtered_markers[filtered_markers$cluster == "Oligodendrocytes Like Cells", ]
NPC_df     <- filtered_markers[filtered_markers$cluster == "NPC Like Cells", ]
Unknown_df <- filtered_markers[filtered_markers$cluster == "Unknown", ]

# Curated gene lists and colors (adapt or extend as you wish)
Curated_Astro <- c("ALDOC", "APOE", "AQP4", "CD44", "CRYAB", "GABBR2", "KCNN3", "LPL", "SLC1A3")
Curated_Astro_df <- data.frame(
  gene  = Curated_Astro,
  color = "#751505"
)

Curated_Oligo <- c("ASCL1", "COL20A1", "DLL1", "DLL3", "EGFR", "OLIG1", "OLIG2", "OPCML", "PDGFRA")
Curated_Oligo_df <- data.frame(
  gene  = Curated_Oligo,
  color = "#751505"
)

Curated_NPC <- c("AURKB", "BUB1", "E2F7", "E2F8", "EZH2", "FOXM1", "MKI67", "PCNA")
Curated_NPC_df <- data.frame(
  gene  = Curated_NPC,
  color = "#751505"
)

# Helper to force curated genes into top20 / top5 lists
add_curated_to_markers <- function(curated_genes, curated_df, top20_markers, top5_markers) {
  for (gene in curated_genes) {
    if (!gene %in% top20_markers$gene) {
      gene_row <- curated_df[curated_df$gene == gene, ]
      top20_markers <- rbind(top20_markers, gene_row)
    }
    gene_row_20 <- top20_markers[top20_markers$gene == gene, ]
    if (!gene %in% top5_markers$gene) {
      top5_markers <- rbind(top5_markers, gene_row_20)
    }
  }
  return(list(top20 = top20_markers, top5 = top5_markers))
}

# Top markers per cell type (from filtered_markers)
top20.markers_Astro <- AC_df %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) %>%
  as.data.frame()
top20.markers_Astro <- top20.markers_Astro[c("gene", "color")]

top5.markers_Astro <- AC_df %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  as.data.frame()
top5.markers_Astro <- top5.markers_Astro[c("gene", "color")]

top20.markers_Oligo <- OC_df %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) %>%
  as.data.frame()
top20.markers_Oligo <- top20.markers_Oligo[c("gene", "color")]

top5.markers_Oligo <- OC_df %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  as.data.frame()
top5.markers_Oligo <- top5.markers_Oligo[c("gene", "color")]

top20.markers_NPC <- NPC_df %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) %>%
  as.data.frame()
top20.markers_NPC <- top20.markers_NPC[c("gene", "color")]

top5.markers_NPC <- NPC_df %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  as.data.frame()
top5.markers_NPC <- top5.markers_NPC[c("gene", "color")]

top20.markers_Unknown <- Unknown_df %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) %>%
  as.data.frame()
top20.markers_Unknown <- top20.markers_Unknown[c("gene", "color")]

top5.markers_Unknown <- Unknown_df %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  as.data.frame()
top5.markers_Unknown <- top5.markers_Unknown[c("gene", "color")]

# Enforce curated genes into top lists
astro_result <- add_curated_to_markers(Curated_Astro, Curated_Astro_df, top20.markers_Astro, top5.markers_Astro)
top20.markers_Astro <- astro_result$top20
top5.markers_Astro  <- astro_result$top5

oligo_result <- add_curated_to_markers(Curated_Oligo, Curated_Oligo_df, top20.markers_Oligo, top5.markers_Oligo)
top20.markers_Oligo <- oligo_result$top20
top5.markers_Oligo  <- oligo_result$top5

npc_result <- add_curated_to_markers(Curated_NPC, Curated_NPC_df, top20.markers_NPC, top5.markers_NPC)
top20.markers_NPC <- npc_result$top20
top5.markers_NPC  <- npc_result$top5

# Combine the 5 top markers from each type (you can also use top20 if needed)
All_interest_df <- rbind(
  top5.markers_Astro,
  top5.markers_NPC,
  top5.markers_Oligo,
  top5.markers_Unknown
)

# Build gene-specific colors for x-axis labels
gene_colors  <- setNames(All_interest_df$color, All_interest_df$gene)
unique_genes <- unique(All_interest_df$gene)
unique_colors <- gene_colors[unique_genes]

# Set factor order for scType_Annotations if you want a specific order in DotPlot
Object.data$scType_Annotations <- factor(
  Object.data$scType_Annotations,
  levels = c("Astrocytes Like Cells",
             "NPC Like Cells",
             "Oligodendrocytes Like Cells",
             "Unknown")
)

desired_order <- rev(levels(Object.data$scType_Annotations))

p <- DotPlot(
  Object.data,
  features = unique_genes,
  group.by = "scType_Annotations"
) +  # [web:33][web:11][web:38]
  scale_color_gradient(low = "white", high = "#4B0082") +
  scale_y_discrete(limits = desired_order) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1,
      color = unique_colors,
      size = 13
    ),
    axis.text.y = element_text(color = "black", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = "black"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 0.90) +
  scale_size_continuous(range = c(0, 9)) +
  guides(
    size  = guide_legend(title = "Fraction of cells\nin group (%)"),
    color = guide_colorbar(title = "Mean expression\nin group")
  )

ggsave(file.path(deg_path, "dotplot_scType.pdf"), p, width = 20, height = 5, dpi = 300)
ggsave(file.path(deg_path, "dotplot_scType.png"), p, width = 20, height = 5, dpi = 300)



