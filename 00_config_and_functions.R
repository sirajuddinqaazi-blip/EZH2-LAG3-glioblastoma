# =====================================================================
# 00_config_and_functions.R
# EZH2 in IDH-wildtype GBM - final analysis pipeline
# Shared settings, packages and helper functions (sourced by all scripts)
# Run everything from the project folder containing "share/".
# =====================================================================

suppressPackageStartupMessages({
  library(survival); library(survminer)
  library(ggplot2); library(patchwork); library(ggrepel)
})

set.seed(2026)
DATA_DIR <- "share"
RES_DIR  <- "results"
FIG_DIR  <- "figures"
dir.create(RES_DIR, showWarnings = FALSE); dir.create(FIG_DIR, showWarnings = FALSE)

# ---- gene sets used throughout ----
PROLIF_GENES <- c("MKI67","TOP2A","CCNB1","CDK1","BUB1","CENPF","AURKA","PCNA","MCM2","TPX2")
PRC2_GENES   <- c("EZH1","EED","SUZ12")
CHECKPOINTS  <- c("LAG3","CD276","PDCD1","CD274","CTLA4","HAVCR2","TIGIT","BTLA")
GENES_NEEDED <- unique(c("EZH2", PRC2_GENES, CHECKPOINTS, PROLIF_GENES))

# ---- partial Spearman correlation -----------------------------------
# Ranks all variables, regresses covariates Z out of x and y, correlates
# the residuals; t-test with n - 2 - k degrees of freedom.
pspear <- function(x, y, Z = character(0), data) {
  dd <- na.omit(data[, c(x, y, Z), drop = FALSE])
  r  <- as.data.frame(lapply(dd, rank))
  n  <- nrow(dd); k <- length(Z)
  if (k == 0) {
    rr <- cor(r[[x]], r[[y]])
  } else {
    rx <- resid(lm(r[[x]] ~ ., data = r[, Z, drop = FALSE]))
    ry <- resid(lm(r[[y]] ~ ., data = r[, Z, drop = FALSE]))
    rr <- cor(rx, ry)
  }
  tt <- rr * sqrt((n - 2 - k) / (1 - rr^2))
  c(n = n, r = rr, p = 2 * pt(-abs(tt), n - 2 - k))
}

# ---- fixed-effect meta-analysis of correlations (Fisher z) -----------
meta_r <- function(r, n, k = 0) {
  z <- atanh(r); w <- n - 3 - k
  zp <- sum(w * z) / sum(w); se <- 1 / sqrt(sum(w))
  Q  <- sum(w * (z - zp)^2); df <- length(r) - 1
  c(r = tanh(zp), lo = tanh(zp - 1.96 * se), hi = tanh(zp + 1.96 * se),
    p = 2 * pnorm(-abs(zp / se)), Q_p = pchisq(Q, df, lower.tail = FALSE),
    I2 = ifelse(Q > 0, max(0, (Q - df) / Q), 0))
}

# ---- formatting ----
fmt_p <- function(p) ifelse(p < 0.0001, "<0.0001", formatC(p, format = "g", digits = 2))
theme_pub <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
COL_LOWHIGH <- c(Low = "#3B6FB6", High = "#C0392B")

save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
