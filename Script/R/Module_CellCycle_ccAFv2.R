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

# Save Cell Cycle Plots ======================================================
for (ext in c("bmp", "eps", "pdf", "png", "svg", "tiff")) {
  ggsave(
    p,
    width = 25, height = 15,
    filename = paste0(ccAFv2_path, "CcAFv2_Annotation.", ext)
  )
}

# Cell Cycle Distribution by Cell Type =======================================
phase_table <- table(Object.data$CellMarkersV2_classification_V3, Object.data$ccAFv2)
phase_df <- as.data.frame.matrix(phase_table)
phase_df$Cell_Type <- rownames(phase_df)
long_data <- pivot_longer(phase_df, cols = -Cell_Type, names_to = "Cell_Cycle_Phase", values_to = "Count")

# Stacked bar plot: absolute counts
ggplot(long_data, aes(x = Cell_Type, y = Count, fill = Cell_Cycle_Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = annotation_colors$Cell_Cycle_Group) +
  theme_minimal() +
  labs(title = "Cell Type Distribution by Cell Cycle Phase", x = "Cell Type", y = "Count", fill = "Cell Cycle Phase")

# Frequency by Cell Type
total_counts <- long_data %>%
  group_by(Cell_Type) %>%
  summarise(Total = sum(Count))

frequency_data <- long_data %>%
  left_join(total_counts, by = "Cell_Type") %>%
  mutate(Frequency = Count / Total) %>%
  select(-Count, -Total)

frequency_table <- frequency_data %>%
  pivot_wider(names_from = Cell_Cycle_Phase, values_from = Frequency)


# Chi-squared Test & Association Plots =======================================

dt <- as.table(as.matrix(phase_df[, setdiff(colnames(phase_df), "Cell_Type")]))
chisq <- chisq.test(dt)
cat("\nChi-squared Test\n")
print(chisq)
cat("\nResiduals:\n")
print(round(chisq$residuals, 3))

contrib <- 100 * chisq$residuals^2 / chisq$statistic
print(round(contrib, 3))

pdf(paste0(ccAFv2_path, "Correlation_Cell_Type.pdf"), width = 10, height = 8)
balloonplot(t(dt), main = "Contingency table", xlab = "", ylab = "", label = FALSE, show.margins = FALSE)
mosaicplot(dt, shade = TRUE, las = 2, main = "Cell type tasks")
assoc(head(dt, 5), shade = TRUE, las = 3)
corrplot(chisq$residuals, is.cor = FALSE)
corrplot(contrib, is.cor = FALSE)
dev.off()


# End of script ==============================================================


