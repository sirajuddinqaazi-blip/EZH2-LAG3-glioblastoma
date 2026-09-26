# =====================================================================
# 05_immune_landscape.R
# EZH2 vs validated immune/stromal estimates (ESTIMATE, MCP-counter),
# unadjusted and adjusted for proliferation, 3 cohorts + meta-analysis.
# Replaces the custom ssGSEA signatures of the earlier version.
# Outputs: results/Table_S3_immune_landscape.csv, figures/FigureS1_*
# =====================================================================
source("00_config_and_functions.R")
L <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts

CELL_VARS <- c(ESTIMATE_immune = "ESTIMATE immune score", ESTIMATE_stromal = "ESTIMATE stromal score",
  MCP_T_cell = "T cells", MCP_T_cell_CD8_ = "CD8 T cells", MCP_cytotoxicity_score = "Cytotoxic lymphocytes",
  MCP_NK_cell = "NK cells", MCP_B_cell = "B lineage", MCP_Monocyte = "Monocytic lineage",
  MCP_Myeloid_dendritic_cell = "Myeloid dendritic cells", MCP_Neutrophil = "Neutrophils",
  MCP_Endothelial_cell = "Endothelial cells", MCP_Cancer_associated_fibroblast = "Fibroblasts")
ADJ <- list(Unadjusted = character(0), `Adjusted for proliferation` = "prolif")

rows <- list()
for (a in names(ADJ)) for (k in names(L)) for (v in names(CELL_VARS)) {
  s <- pspear("EZH2", v, ADJ[[a]], L[[k]])
  rows[[length(rows) + 1]] <- data.frame(adjustment = a, cohort = k, variable = v, label = CELL_VARS[[v]],
                                         n = s[["n"]], r = s[["r"]], p = s[["p"]])
}
tab <- do.call(rbind, rows)
tab$q <- ave(tab$p, tab$adjustment, tab$cohort, FUN = function(p) p.adjust(p, "BH"))
meta <- do.call(rbind, lapply(split(tab, list(tab$adjustment, tab$variable)), function(x) {
  m <- meta_r(x$r, x$n, length(ADJ[[x$adjustment[1]]]))
  data.frame(adjustment = x$adjustment[1], variable = x$variable[1], label = x$label[1],
             TCGA = x$r[x$cohort == "TCGA"], CGGA_325 = x$r[x$cohort == "CGGA_325"], CGGA_693 = x$r[x$cohort == "CGGA_693"],
             r_pooled = m[["r"]], p_pooled = m[["p"]], I2 = m[["I2"]], same_direction = length(unique(sign(x$r))) == 1)
}))
meta$q_pooled <- ave(meta$p_pooled, meta$adjustment, FUN = function(p) p.adjust(p, "BH"))
meta$robust <- meta$q_pooled < 0.05 & meta$I2 < 0.5 & meta$same_direction
write.csv(tab,  file.path(RES_DIR, "Table_S3a_immune_by_cohort.csv"), row.names = FALSE)
write.csv(meta, file.path(RES_DIR, "Table_S3b_immune_meta.csv"), row.names = FALSE)
print(transform(meta[order(meta$adjustment, meta$r_pooled), c("adjustment", "label", "TCGA", "CGGA_325", "CGGA_693", "r_pooled", "q_pooled", "I2", "robust")],
                TCGA = round(TCGA, 2), CGGA_325 = round(CGGA_325, 2), CGGA_693 = round(CGGA_693, 2),
                r_pooled = round(r_pooled, 2), q_pooled = fmt_p(q_pooled), I2 = round(I2, 2)), row.names = FALSE)

tab$label <- factor(tab$label, levels = rev(unname(CELL_VARS)))
tab$star  <- ifelse(tab$q < 0.05, "*", "")
pS <- ggplot(tab, aes(cohort, label, fill = r)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f%s", r, star)), size = 2.8) +
  facet_wrap(~adjustment) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C0392B", limits = c(-0.7, 0.7)) +
  labs(x = NULL, y = NULL, fill = "Spearman r\nwith EZH2",
       title = "EZH2 and the tumor microenvironment (ESTIMATE / MCP-counter)",
       caption = "* BH q < 0.05 within cohort") + theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(pS, "FigureS1_immune_landscape", 9, 6)
cat("✓ Figure S1 and Table S3 saved\n")
