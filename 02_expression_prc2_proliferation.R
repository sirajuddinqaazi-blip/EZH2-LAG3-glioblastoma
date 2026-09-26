# =====================================================================
# 02_expression_prc2_proliferation.R
# EZH2 tumor vs normal, PRC2 co-expression, clinical correlates,
# and EZH2-proliferation coupling in all three cohorts.
# Outputs: results/Table_S1_expression.csv, figures/Figure1_*
# =====================================================================
source("00_config_and_functions.R")
obj <- readRDS(file.path(RES_DIR, "cohorts.rds")); L <- obj$cohorts; N <- obj$tcga_normals
tcga <- L$TCGA

# ---- A. Tumor vs normal (TCGA) ----
tn <- rbind(data.frame(group = "GBM (IDH-wt)", EZH2 = tcga$EZH2),
            data.frame(group = "Normal brain", EZH2 = N$EZH2))
tn$group <- factor(tn$group, levels = c("Normal brain", "GBM (IDH-wt)"))
w_tn <- wilcox.test(EZH2 ~ group, data = tn)
fc   <- median(tcga$EZH2) / median(N$EZH2)
cat(sprintf("Tumor vs normal: median %.2f vs %.2f TPM, fold = %.1f, Wilcoxon p = %s\n",
            median(tcga$EZH2), median(N$EZH2), fc, fmt_p(w_tn$p.value)))
pA <- ggplot(tn, aes(group, log2(EZH2 + 1), fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  scale_fill_manual(values = c("grey70", "#C0392B"), guide = "none") +
  annotate("text", x = 1.5, y = max(log2(tn$EZH2 + 1)) * 1.08,
           label = paste0("Wilcoxon p = ", fmt_p(w_tn$p.value)), size = 3.2) +
  labs(x = NULL, y = "EZH2 log2(TPM+1)", title = "A  TCGA: tumor vs normal") + theme_pub

# ---- B. PRC2 + proliferation correlations, all cohorts ----
rows <- list()
for (k in names(L)) for (g in c(PRC2_GENES, "prolif")) {
  s <- pspear("EZH2", g, character(0), L[[k]])
  rows[[length(rows) + 1]] <- data.frame(cohort = k, variable = g, n = s[["n"]], rho = s[["r"]], p = s[["p"]])
}
cor_tab <- do.call(rbind, rows); cor_tab$q <- p.adjust(cor_tab$p, "BH")
print(transform(cor_tab, rho = round(rho, 3), p = fmt_p(p), q = fmt_p(q)))

bdat <- do.call(rbind, lapply(L, function(d) d[, c("cohort", "log2EZH2", "prolif")])); bdat$cohort <- factor(bdat$cohort, levels = names(L))
cor_tab$cohort <- factor(cor_tab$cohort, levels = names(L))
pB <- ggplot(bdat,
             aes(prolif, log2EZH2)) +
  geom_point(size = 0.8, alpha = 0.6) + geom_smooth(method = "lm", se = TRUE, colour = "#C0392B", linewidth = 0.7) +
  facet_wrap(~cohort, scales = "free") +
  geom_text(data = subset(cor_tab, variable == "prolif"),
            aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f\np %s", rho, ifelse(p < 1e-4, "< 0.0001", paste("=", fmt_p(p))))),
            hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
  labs(x = "Proliferation score (10 genes)", y = "EZH2 log2(expression+1)",
       title = "C  EZH2 tracks tumor proliferation") + theme_pub

pC <- ggplot(subset(cor_tab, variable != "prolif"), aes(variable, cohort, fill = rho)) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", rho)), size = 3.2) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C0392B", limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman\nrho", title = "B  EZH2 vs PRC2 genes") + theme_pub

# ---- D. Clinical correlates (TCGA) ----
kw  <- kruskal.test(EZH2 ~ subtype, data = tcga)
pw  <- pairwise.wilcox.test(tcga$EZH2, tcga$subtype, p.adjust.method = "BH")
w_m <- wilcox.test(EZH2 ~ MGMT, data = tcga); w_s <- wilcox.test(EZH2 ~ sex, data = tcga)
cat(sprintf("Subtype KW p = %s | MGMT p = %s | sex p = %s\n", fmt_p(kw$p.value), fmt_p(w_m$p.value), fmt_p(w_s$p.value)))
print(pw$p.value)
pD <- ggplot(subset(tcga, !is.na(subtype)), aes(subtype, log2EZH2, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = "EZH2 log2(TPM+1)", title = paste0("D  TCGA subtype (Kruskal-Wallis p = ", fmt_p(kw$p.value), ")")) + theme_pub

fig1 <- (pA | pC) / pB / pD + plot_layout(heights = c(1, 1, 1))
save_fig(fig1, "Figure1_expression_PRC2_proliferation", 10, 11)

clin_tab <- data.frame(test = c("Tumor vs normal (TCGA)", "MGMT unmeth vs meth (TCGA)", "Sex (TCGA)", "Subtype KW (TCGA)"),
                       p = c(w_tn$p.value, w_m$p.value, w_s$p.value, kw$p.value))
write.csv(cor_tab, file.path(RES_DIR, "Table_S1a_EZH2_PRC2_proliferation_correlations.csv"), row.names = FALSE)
write.csv(clin_tab, file.path(RES_DIR, "Table_S1b_EZH2_clinical_tests.csv"), row.names = FALSE)
cat("✓ Figure 1 and Table S1 saved\n")
