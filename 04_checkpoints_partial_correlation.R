# =====================================================================
# 04_checkpoints_partial_correlation.R
# Core analysis: EZH2 vs 8 immune-checkpoint genes, raw and adjusted
# for confounders, in 3 cohorts + fixed-effect meta-analysis.
#   M0 raw | M1 + proliferation | M2 + ESTIMATE immune & stromal
#   M3 + proliferation + ESTIMATE | M4 = M3 + MCP-counter T cells
#   TCGA only: M3 + ABSOLUTE (DNA-based) purity
# BH correction across the 8 genes within each cohort x model.
# Outputs: results/Table3_*.csv, figures/Figure3_*
# =====================================================================
source("00_config_and_functions.R")
L <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts

MODELS <- list(
  M0_raw                    = character(0),
  M1_proliferation          = "prolif",
  M2_ESTIMATE               = c("ESTIMATE_immune", "ESTIMATE_stromal"),
  M3_prolif_ESTIMATE        = c("prolif", "ESTIMATE_immune", "ESTIMATE_stromal"),
  M4_prolif_ESTIMATE_Tcells = c("prolif", "ESTIMATE_immune", "ESTIMATE_stromal", "MCP_T_cell"))

rows <- list()
for (m in names(MODELS)) for (k in names(L)) for (g in CHECKPOINTS) {
  s <- pspear("EZH2", g, MODELS[[m]], L[[k]])
  rows[[length(rows) + 1]] <- data.frame(model = m, cohort = k, gene = g, n = s[["n"]], r = s[["r"]], p = s[["p"]])
}
tab <- do.call(rbind, rows)
tab$q <- ave(tab$p, tab$model, tab$cohort, FUN = function(p) p.adjust(p, "BH"))

# TCGA-only DNA purity model
tp <- do.call(rbind, lapply(CHECKPOINTS, function(g) {
  s <- pspear("EZH2", g, c(MODELS$M3_prolif_ESTIMATE, "ABSOLUTE_purity"), L$TCGA)
  data.frame(model = "M3 + ABSOLUTE purity (TCGA)", cohort = "TCGA", gene = g, n = s[["n"]], r = s[["r"]], p = s[["p"]]) }))
tp$q <- p.adjust(tp$p, "BH")

# meta-analysis per gene x model
meta_tab <- do.call(rbind, lapply(split(tab, list(tab$model, tab$gene)), function(x) {
  m <- meta_r(x$r, x$n, k = length(MODELS[[x$model[1]]]))
  data.frame(model = x$model[1], gene = x$gene[1], r_pooled = m[["r"]], lo = m[["lo"]], hi = m[["hi"]],
             p_pooled = m[["p"]], I2 = m[["I2"]], Q_p = m[["Q_p"]],
             n_cohorts_sig = sum(x$q < 0.05 & sign(x$r) == sign(m[["r"]])),
             same_direction_all = length(unique(sign(x$r))) == 1)
}))
meta_tab$q_pooled <- ave(meta_tab$p_pooled, meta_tab$model, FUN = function(p) p.adjust(p, "BH"))
meta_tab <- meta_tab[order(meta_tab$model, -meta_tab$r_pooled), ]

write.csv(rbind(tab, tp), file.path(RES_DIR, "Table3a_checkpoint_partial_correlations_by_cohort.csv"), row.names = FALSE)
write.csv(meta_tab, file.path(RES_DIR, "Table3b_checkpoint_meta_analysis.csv"), row.names = FALSE)

cat("\n=== Pooled partial correlations (EZH2 vs checkpoint) ===\n")
print(transform(meta_tab[, c("model", "gene", "r_pooled", "lo", "hi", "p_pooled", "q_pooled", "I2", "n_cohorts_sig")],
                r_pooled = round(r_pooled, 3), lo = round(lo, 3), hi = round(hi, 3), I2 = round(I2, 2),
                p_pooled = fmt_p(p_pooled), q_pooled = fmt_p(q_pooled)), row.names = FALSE)
cat("\n=== LAG3 by cohort ===\n")
print(transform(subset(rbind(tab, tp), gene == "LAG3"), r = round(r, 3), p = fmt_p(p), q = fmt_p(q)), row.names = FALSE)

# ---- supporting analyses reported in the text ----
supp <- do.call(rbind, lapply(names(L), function(k) {
  a1 <- pspear("LAG3", "MCP_T_cell", character(0), L[[k]])        # LAG3 tracks T-cell abundance
  a2 <- pspear("prolif", "LAG3", "EZH2", L[[k]])                  # proliferation beyond EZH2
  a3 <- pspear("EZH2", "LAG3", "prolif", L[[k]])                  # EZH2 beyond proliferation
  data.frame(cohort = k, LAG3_vs_Tcells_rho = a1[["r"]], LAG3_vs_Tcells_p = a1[["p"]],
             prolif_LAG3_given_EZH2_r = a2[["r"]], prolif_LAG3_given_EZH2_p = a2[["p"]],
             EZH2_LAG3_given_prolif_r = a3[["r"]], EZH2_LAG3_given_prolif_p = a3[["p"]]) }))
write.csv(supp, file.path(RES_DIR, "Table3c_LAG3_supporting_analyses.csv"), row.names = FALSE)
cat("\n=== LAG3 supporting analyses ===\n"); print(transform(supp, LAG3_vs_Tcells_rho = round(LAG3_vs_Tcells_rho, 3),
  prolif_LAG3_given_EZH2_r = round(prolif_LAG3_given_EZH2_r, 3), EZH2_LAG3_given_prolif_r = round(EZH2_LAG3_given_prolif_r, 3)))

# ------------------------------ Figure 3 ------------------------------
lab_models <- c(M0_raw = "Unadjusted", M1_proliferation = "+ Proliferation", M2_ESTIMATE = "+ Immune/stromal (ESTIMATE)",
                M3_prolif_ESTIMATE = "+ Proliferation + ESTIMATE", M4_prolif_ESTIMATE_Tcells = "+ Prolif. + ESTIMATE + T cells")
# A: pooled r per gene, unadjusted vs fully adjusted (M3)
fa <- subset(meta_tab, model %in% c("M0_raw", "M3_prolif_ESTIMATE"))
fa$model <- factor(lab_models[fa$model], levels = lab_models[c("M0_raw", "M3_prolif_ESTIMATE")])
fa$gene  <- factor(fa$gene, levels = rev(subset(meta_tab, model == "M3_prolif_ESTIMATE")$gene))
fa$sig   <- ifelse(fa$q_pooled < 0.05 & fa$I2 < 0.5, "q<0.05 & I2<50%", "not robust")
pA <- ggplot(fa, aes(r_pooled, gene, colour = model, shape = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.35) +
  scale_colour_manual(values = c("grey55", "#C0392B")) + scale_shape_manual(values = c(`not robust` = 1, `q<0.05 & I2<50%` = 16)) +
  labs(x = "Pooled partial Spearman r with EZH2 (95% CI)", y = NULL, colour = NULL, shape = NULL,
       title = "A  EZH2-checkpoint associations, 3 cohorts pooled") + theme_pub + theme(legend.position = "bottom", legend.box = "vertical")

# B: LAG3 and CD276 per cohort across models
fb <- subset(tab, gene %in% c("LAG3", "CD276")); fb$cohort <- factor(fb$cohort, levels = names(L))
fb$model <- factor(lab_models[fb$model], levels = rev(lab_models))
pB <- ggplot(fb, aes(r, model, colour = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 2.2, position = position_dodge(width = 0.5)) + facet_wrap(~gene) +
  scale_colour_manual(values = c(TCGA = "#2C3E50", CGGA_325 = "#E67E22", CGGA_693 = "#16A085")) +
  labs(x = "Partial Spearman r with EZH2", y = NULL, colour = NULL, title = "B  Effect of confounder adjustment") +
  theme_pub + theme(legend.position = "bottom")

# C: residual plots (LAG3 vs EZH2 after M3 adjustment)
resid_df <- do.call(rbind, lapply(names(L), function(k) {
  dd <- na.omit(L[[k]][, c("EZH2", "LAG3", MODELS$M3_prolif_ESTIMATE)]); r <- as.data.frame(lapply(dd, rank))
  Z <- r[, MODELS$M3_prolif_ESTIMATE]
  data.frame(cohort = k, EZH2_res = resid(lm(r$EZH2 ~ ., data = Z)), LAG3_res = resid(lm(r$LAG3 ~ ., data = Z)))
}))
lab_c <- subset(tab, gene == "LAG3" & model == "M3_prolif_ESTIMATE")
resid_df$cohort <- factor(resid_df$cohort, levels = names(L)); lab_c$cohort <- factor(lab_c$cohort, levels = names(L))
pC <- ggplot(resid_df, aes(EZH2_res, LAG3_res)) + geom_point(size = 0.8, alpha = 0.6) +
  geom_smooth(method = "lm", colour = "#C0392B", linewidth = 0.7) + facet_wrap(~cohort, scales = "free") +
  geom_text(data = lab_c, aes(x = -Inf, y = Inf, label = sprintf("partial r = %.2f\np %s", r, ifelse(p < 1e-4, "< 0.0001", paste("=", fmt_p(p))))),
            hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
  labs(x = "EZH2 (rank residual)", y = "LAG3 (rank residual)",
       title = "C  EZH2-LAG3 after adjusting for proliferation, immune and stromal content") + theme_pub

fig3 <- (pA | pB) / pC + plot_layout(heights = c(1.3, 1))
save_fig(fig3, "Figure3_checkpoints_partial_correlation", 12, 10)
cat("✓ Figure 3 and Table 3 saved\n")
