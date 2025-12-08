# Dassou SIKA
# dassou.sika@igf.cnrs.fr
# Bioinformatics Engenier
# Plasticité cérébrale, cellules souches et tumeurs gliales de bas grade / IGF
# CNRS | Montpelier-France
#
# ------------------------------------------------------------
# Cell Cycle Annotation via ccAFv2
# ------------------------------------------------------------

# Required Libraries =========================================================
library(Seurat)         # Core Seurat object manipulation
library(dplyr)          # Data manipulation
library(Matrix)         # Sparse matrices
library(ggplot2)        # Plots
library(ccAFv2)         # Cell cycle annotation with ccAFv2
library(patchwork)      # ggplot composition
library(tidyr)          # Data pivoting
library(vcd)            # Assoc plot and mosaicplot
library(gplots)         # Balloonplot
library(corrplot)       # Correlation heatmaps

# Paths ======================================================================

ccAFv2_path <- paste0("ccAFV2/")

# Data Loading ==============================================================

Object.data <- readRDS("LGG275_diff.rds")

# Cell Cycle Prediction with ccAFv2 ==========================================
Object.data <- PredictCellCycle(
  Object.data,
  cutoff = 0.5,
  do_sctransform = TRUE,
  assay = 'SCT',
  species = 'human',
  gene_id = 'symbol',
  spatial = FALSE
)

# Color Definitions ==========================================================
annotation_colors <- list(
  Cell_Cycle_Group = c(
    "G1" = "brown",
    "G2M" = "#039cb1",
    "G2/M" = "#229954",
    "S" = "orange",
    "Neural G0" = "black",
    "Late G1" = "#56017b",
    "S/G2" = "#DAF7A6",
    "M/Early G1" = "#ba4a00",
    "Unknown" = "#5d6d7e"
  )
)

Dim2 <- DimPlot.ccAFv2(seur_obj, reduction = "umap", split.by = "stim", label=FALSE,label.size = 8,pt.size =6)+ xlim(c(-10, 15)) + ylim(c(-10, 10))


p <- Dim2 +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 40),
    strip.text.x = element_text(size = 40),
    plot.title = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1, bycol = TRUE, override.aes = list(size = 8), title.size = 20)) +
  facet_wrap(~ stim, labeller = labeller(stim = stim_labels))



# Cell Cycle correlation by Cell Type =======================================
phase_table <- table(Object.data$scType_Annotations, Object.data$ccAFv2)
phase_df <- as.data.frame.matrix(phase_table)
phase_df$Cell_Type <- rownames(phase_df)

dt <- as.table(as.matrix(phase_df[, setdiff(colnames(phase_df), "Cell_Type")]))
chisq <- chisq.test(dt)
cat("\nChi-squared Test\n")
print(chisq)
cat("\nResiduals:\n")
print(round(chisq$residuals, 3))

contrib <- 100 * chisq$residuals^2 / chisq$statistic
print(round(contrib, 3))

pdf(paste0(ccAFv2_path, "Correlation_Cell_Type.pdf"), width = 10, height = 8)
corrplot(chisq$residuals, is.cor = FALSE)
dev.off()

