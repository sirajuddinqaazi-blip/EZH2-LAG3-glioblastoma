# =====================================================================
# 06_single_cell_TISCH2.R
# Cell-type-resolved expression of EZH2, LAG3 and PDCD1 in adult GBM
# single-cell data (TISCH2 heatmap exports; Neftel et al. 2019,
# GSE131928 10X and Smart-seq2). Values: TISCH2 mean log(TPM/10+1).
# Input: share/TISCH_EZH2_heatmap.csv, TISCH_LAG3_heatmap.csv, TISCH_PDCD1_heatmap.csv
# Output: figures/Figure4_single_cell.*, results/Table_S4_single_cell.csv
# NOTE: restricted to cell types annotated in GSE131928 (verify on the
# TISCH2 dataset page before submission).
# =====================================================================
source("00_config_and_functions.R")
read_tisch <- function(g) {
  x <- read.csv(file.path(DATA_DIR, paste0("TISCH_", g, "_heatmap.csv")), check.names = FALSE, fileEncoding = "UTF-8-BOM")
  colnames(x) <- c("cell_type", "dataset", "value"); x$gene <- g; x
}
sc <- do.call(rbind, lapply(c("EZH2", "LAG3", "PDCD1"), read_tisch))
DATASETS <- c("Glioma_GSE131928_10X" = "Neftel 10X", "Glioma_GSE131928_Smartseq2" = "Neftel Smart-seq2")
CELLS <- c("CD8Tex" = "Exhausted CD8 T", "Mono/Macro" = "Monocyte/macrophage",
           "AC-like Malignant" = "Malignant AC-like", "MES-like Malignant" = "Malignant MES-like",
           "NPC-like Malignant" = "Malignant NPC-like", "OPC-like Malignant" = "Malignant OPC-like")
sc <- subset(sc, dataset %in% names(DATASETS) & cell_type %in% names(CELLS))
sc$dataset   <- factor(DATASETS[sc$dataset], levels = DATASETS)
sc$cell_type <- factor(CELLS[sc$cell_type], levels = CELLS)
sc$gene      <- factor(sc$gene, levels = c("PDCD1", "LAG3", "EZH2"))
write.csv(sc, file.path(RES_DIR, "Table_S4_single_cell_TISCH2.csv"), row.names = FALSE)
print(reshape(sc[, c("dataset", "gene", "cell_type", "value")], idvar = c("dataset", "gene"), timevar = "cell_type", direction = "wide"))

p4 <- ggplot(sc, aes(cell_type, gene, fill = value)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3.2) + facet_wrap(~dataset, ncol = 1) +
  scale_fill_gradient(low = "grey97", high = "#B03A2E", name = "Mean\nlog(TPM/10+1)") +
  labs(x = NULL, y = NULL, title = "Cell-type expression in adult GBM single-cell data (GSE131928)") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(face = "italic"))
save_fig(p4, "Figure4_single_cell", 7.5, 5.5)
cat("✓ Figure 4 saved\n")
