# =====================================================================
# 03_survival.R
# Prognostic value of EZH2 in three cohorts.
#  - Kaplan-Meier (median split), log-rank
#  - Cox: EZH2 as continuous (HR per 1 SD of log2 expression),
#    univariable and adjusted (TCGA: age, KPS, MGMT; CGGA: age, MGMT,
#    radiotherapy, chemotherapy). Age is modelled continuously.
#  - Reverse Kaplan-Meier median follow-up
#  - Fixed-effect pooled HR across cohorts
# Outputs: results/Table2_Cox.csv, figures/Figure2_survival.*
# =====================================================================
source("00_config_and_functions.R")
L <- readRDS(file.path(RES_DIR, "cohorts.rds"))$cohorts

km_plots <- list(); cox_rows <- list(); ph_rows <- list(); sens_rows <- list()
for (k in names(L)) {
  d <- subset(L[[k]], !is.na(OS_months) & !is.na(event))
  d$EZH2_z <- as.numeric(scale(d$log2EZH2))
  d$age10  <- d$age / 10

  fit <- survfit(Surv(OS_months, event) ~ EZH2_group, data = d)
  lr  <- survdiff(Surv(OS_months, event) ~ EZH2_group, data = d)
  p_lr <- pchisq(lr$chisq, 1, lower.tail = FALSE)
  fu  <- summary(survfit(Surv(OS_months, 1 - event) ~ 1, data = d))$table["median"]
  med <- summary(fit)$table[, "median"]
  cat(sprintf("\n%s: n=%d, events=%d, median follow-up (reverse KM) = %.1f mo; median OS Low/High = %.1f/%.1f; log-rank p = %s\n",
              k, nrow(d), sum(d$event), fu, med[1], med[2], fmt_p(p_lr)))

  g <- ggsurvplot(fit, data = d, palette = unname(COL_LOWHIGH), legend.labs = c("EZH2-low", "EZH2-high"),
                  legend.title = "", xlab = "Months", ylab = "Overall survival", risk.table = FALSE,
                  censor.size = 2, ggtheme = theme_pub)
  km_plots[[k]] <- g$plot + labs(title = sprintf("%s  %s (n = %d)", LETTERS[match(k, names(L))], k, nrow(d))) +
    annotate("text", x = max(d$OS_months) * 0.6, y = 0.9, size = 3.2, label = paste0("log-rank p = ", fmt_p(p_lr)))

  covs <- if (k == "TCGA") c("age10", "KPS", "MGMT") else c("age10", "MGMT", "radio", "chemo")
  for (model in c("Univariable", "Adjusted")) {
    f  <- as.formula(paste("Surv(OS_months, event) ~ EZH2_z", if (model == "Adjusted") paste("+", paste(covs, collapse = " + ")) else ""))
    cx <- coxph(f, data = d); s <- summary(cx)
    zt <- cox.zph(cx)$table; ph <- zt["GLOBAL", "p"]
    ph_rows[[length(ph_rows) + 1]] <- data.frame(cohort = k, model = model, term = rownames(zt), chisq = zt[, 1], p = zt[, 3])
    for (v in rownames(s$coefficients)) cox_rows[[length(cox_rows) + 1]] <- data.frame(
      cohort = k, model = model, n = s$n, events = s$nevent, term = v,
      HR = s$conf.int[v, 1], lo = s$conf.int[v, 3], hi = s$conf.int[v, 4], p = s$coefficients[v, 5],
      C_index = s$concordance[1], PH_global_p = ph)
  }
  # ---- sensitivity analyses for proportional-hazards violations ----
  if (k == "TCGA") {   # age violates PH -> stratify baseline hazard by age tertile
    d$age_strat <- cut(d$age, quantile(d$age, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE)
    s2 <- summary(coxph(Surv(OS_months, event) ~ EZH2_z + KPS + MGMT + strata(age_strat), data = d))
    sens_rows[[length(sens_rows) + 1]] <- data.frame(cohort = k, analysis = "Adjusted, stratified by age tertile",
      HR = s2$conf.int["EZH2_z", 1], lo = s2$conf.int["EZH2_z", 3], hi = s2$conf.int["EZH2_z", 4], p = s2$coefficients["EZH2_z", 5])
  }
  # time-split EZH2 effect (0-24 months vs >24 months), all cohorts
  ds <- survSplit(Surv(OS_months, event) ~ ., data = d, cut = 24, episode = "period")
  s3 <- summary(coxph(as.formula(paste("Surv(tstart, OS_months, event) ~ EZH2_z:strata(period) +", paste(covs, collapse = " + "))), data = ds))
  for (i in grep("EZH2_z", rownames(s3$coefficients))) sens_rows[[length(sens_rows) + 1]] <- data.frame(
    cohort = k, analysis = paste0("Adjusted, EZH2 effect ", ifelse(grepl("=1", rownames(s3$coefficients)[i]), "0-24 mo", ">24 mo")),
    HR = s3$conf.int[i, 1], lo = s3$conf.int[i, 3], hi = s3$conf.int[i, 4], p = s3$coefficients[i, 5])
}
cox_tab <- do.call(rbind, cox_rows)
print(transform(cox_tab, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2), p = fmt_p(p), C_index = round(C_index, 3), PH_global_p = round(PH_global_p, 3)))
write.csv(cox_tab, file.path(RES_DIR, "Table2_Cox_models.csv"), row.names = FALSE)
ph_tab <- do.call(rbind, ph_rows); write.csv(ph_tab, file.path(RES_DIR, "Table_S2a_PH_assumption.csv"), row.names = FALSE)
sens_tab <- do.call(rbind, sens_rows); write.csv(sens_tab, file.path(RES_DIR, "Table_S2b_Cox_sensitivity.csv"), row.names = FALSE)
cat("\nPH violations (p<0.05):\n"); print(subset(ph_tab, p < 0.05 & term != "GLOBAL"))
cat("\nSensitivity analyses:\n"); print(transform(sens_tab, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2), p = fmt_p(p)))

# ---- pooled HR for EZH2 (inverse-variance, log scale) ----
ez <- subset(cox_tab, term == "EZH2_z")
pool <- do.call(rbind, lapply(split(ez, ez$model), function(x) {
  b <- log(x$HR); se <- (log(x$hi) - log(x$lo)) / (2 * 1.96); w <- 1 / se^2
  bp <- sum(w * b) / sum(w); sp <- 1 / sqrt(sum(w))
  data.frame(cohort = "Pooled", model = x$model[1], HR = exp(bp), lo = exp(bp - 1.96 * sp), hi = exp(bp + 1.96 * sp),
             p = 2 * pnorm(-abs(bp / sp)))
}))
print(pool)
forest <- rbind(ez[, c("cohort", "model", "HR", "lo", "hi", "p")], pool)
forest$cohort <- factor(forest$cohort, levels = rev(c(names(L), "Pooled")))
pF <- ggplot(forest, aes(HR, cohort, colour = model)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), size = 0.35) +
  scale_colour_manual(values = c(Adjusted = "#C0392B", Univariable = "grey40")) +
  scale_x_log10() + labs(x = "HR per 1 SD of log2 EZH2 (95% CI)", y = NULL, colour = NULL,
                         title = "D  EZH2 and overall survival") + theme_pub + theme(legend.position = "bottom")

fig2 <- (km_plots[[1]] | km_plots[[2]]) / (km_plots[[3]] | pF)
save_fig(fig2, "Figure2_survival", 10, 8)
cat("✓ Figure 2 and Table 2 saved\n")
