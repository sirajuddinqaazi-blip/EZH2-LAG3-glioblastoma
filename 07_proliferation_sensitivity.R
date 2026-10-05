# =====================================================================
# 07_proliferation_sensitivity.R
# Sensitivity analysis for the EZH2-checkpoint associations.
# EZH2 correlates with the 10-gene proliferation score at rho ~ 0.9, so
# the "independent of proliferation" result could reflect residual
# confounding if the 10-gene score measures proliferation imprecisely.
# Here the adjustment is repeated with a much larger, published cell-cycle
# signature: MSigDB Hallmark E2F targets + G2M checkpoint (union; EZH2
# removed), scored by ssGSEA and, alternatively, by mean z-score.
# Run from the project folder AFTER run_all.R (needs results/cohorts.rds).
# Output: results/Table_S5_proliferation_sensitivity.csv
# =====================================================================
suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("GSVA", quietly = TRUE)) stop("Package GSVA needed: install.packages('BiocManager'); BiocManager::install('GSVA')")
say <- function(...) cat(sprintf(...), "\n", sep = "")
fmt_p <- function(p) ifelse(p < 1e-4, "<0.0001", formatC(p, format = "g", digits = 2))

# ---- helpers (same partial Spearman and Fisher-z pooling as the main analysis; + random effects) ----
pspear <- function(x, y, Z = character(0), data) {
  dd <- na.omit(data[, c(x, y, Z), drop = FALSE]); r <- as.data.frame(lapply(dd, rank))
  n <- nrow(dd); k <- length(Z)
  rr <- if (k == 0) cor(r[[x]], r[[y]]) else
    cor(resid(lm(r[[x]] ~ ., data = r[, Z, drop = FALSE])), resid(lm(r[[y]] ~ ., data = r[, Z, drop = FALSE])))
  tt <- rr * sqrt((n - 2 - k) / (1 - rr^2))
  c(n = n, r = rr, p = 2 * pt(-abs(tt), n - 2 - k))
}
meta_r <- function(r, n, k = 0) {
  z <- atanh(r); w <- n - 3 - k; zp <- sum(w * z) / sum(w); se <- 1 / sqrt(sum(w))
  Q <- sum(w * (z - zp)^2); df <- length(r) - 1
  tau2 <- if (df > 0) max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w))) else 0
  wr <- 1 / (1 / w + tau2); zr <- sum(wr * z) / sum(wr); ser <- 1 / sqrt(sum(wr))
  c(r = tanh(zp), lo = tanh(zp - 1.96 * se), hi = tanh(zp + 1.96 * se), p = 2 * pnorm(-abs(zp / se)),
    I2 = ifelse(Q > 0, max(0, (Q - df) / Q), 0),
    r_RE = tanh(zr), lo_RE = tanh(zr - 1.96 * ser), hi_RE = tanh(zr + 1.96 * ser))
}
collapse_symbols <- function(mat, sym) {
  ok <- !is.na(sym) & sym != ""; mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
  o <- order(rowMeans(mat), decreasing = TRUE); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  keep <- !duplicated(sym); mat <- mat[keep, , drop = FALSE]; rownames(mat) <- sym[keep]; mat
}

# ---- data ----
co <- readRDS(file.path("results", "cohorts.rds"))$cohorts
tpm <- readRDS(file.path("share", "tcga_tpm.rds")); gi <- readRDS(file.path("share", "tcga_gene_info.rds"))
mats <- list(TCGA = collapse_symbols(as.matrix(tpm), gi$gene_name[match(rownames(tpm), rownames(gi))]))
rm(tpm); invisible(gc())
for (id in c("325", "693")) {
  ex <- readRDS(file.path("share", paste0("cgga", id, "_expr_gbm.rds")))
  m <- as.matrix(ex[, -1]); storage.mode(m) <- "numeric"
  mats[[paste0("CGGA_", id)]] <- collapse_symbols(m, as.character(ex[[1]]))
}

# ---- published cell-cycle signature (EZH2 removed) ----
if (!exists("get_hallmark")) get_hallmark <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Package msigdbr needed: install.packages('msigdbr')")
  h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  unique(h$gene_symbol[h$gs_name %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")])
}
cc_genes <- setdiff(get_hallmark(), c("EZH2", "LAG3", "PDCD1", "BTLA", "CD276"))
say("Cell-cycle signature: Hallmark E2F targets + G2M checkpoint, %d genes after removing EZH2 (and target genes)", length(cc_genes))

for (k in names(co)) {
  d <- co[[k]]; m <- mats[[k]][, d$sample]; lm2 <- log2(m + 1); lm2 <- lm2[apply(lm2, 1, var) > 0, ]
  g <- intersect(cc_genes, rownames(lm2))
  es <- suppressMessages(GSVA::gsva(GSVA::ssgseaParam(lm2, list(cc = g)), verbose = FALSE))
  d$cc_ssgsea <- as.numeric(es["cc", ])
  d$cc_meanz  <- colMeans(t(scale(t(lm2[g, , drop = FALSE]))))
  co[[k]] <- d
  say("  [%s] %d/%d signature genes found | rho(EZH2, cell-cycle ssGSEA) = %.2f | rho(EZH2, 10-gene score) = %.2f | rho(10-gene, cell-cycle) = %.2f",
      k, length(g), length(cc_genes), cor(d$EZH2, d$cc_ssgsea, method = "spearman"),
      cor(d$EZH2, d$prolif, method = "spearman"), cor(d$prolif, d$cc_ssgsea, method = "spearman"))
}

# ---- models ----
EST <- c("ESTIMATE_immune", "ESTIMATE_stromal")
MODELS <- list(
  "Original: 10-gene + ESTIMATE"                 = c("prolif", EST),
  "Hallmark cell-cycle (ssGSEA) + ESTIMATE"      = c("cc_ssgsea", EST),
  "Hallmark cell-cycle (mean z) + ESTIMATE"      = c("cc_meanz", EST),
  "10-gene + Hallmark (ssGSEA) + ESTIMATE"       = c("prolif", "cc_ssgsea", EST),
  "MKI67 alone + ESTIMATE"                       = c("MKI67", EST),
  "10-gene + Hallmark + ESTIMATE + T cells"      = c("prolif", "cc_ssgsea", EST, "MCP_T_cell"))
GENES <- c("LAG3", "PDCD1", "BTLA", "CD276")

rows <- list()
for (g in GENES) for (mn in names(MODELS)) {
  Z <- MODELS[[mn]]
  per <- rbindlist(lapply(names(co), function(k) { v <- pspear("EZH2", g, Z, co[[k]])
    data.table(cohort = k, n = v[["n"]], r = v[["r"]], p = v[["p"]]) }))
  m <- meta_r(per$r, per$n, length(Z))
  rows[[length(rows) + 1]] <- data.table(gene = g, model = mn,
    r_TCGA = per$r[1], r_CGGA_325 = per$r[2], r_CGGA_693 = per$r[3],
    p_TCGA = per$p[1], p_CGGA_325 = per$p[2], p_CGGA_693 = per$p[3],
    r_pooled = m[["r"]], lo = m[["lo"]], hi = m[["hi"]], p_pooled = m[["p"]], I2 = m[["I2"]],
    r_RE = m[["r_RE"]], lo_RE = m[["lo_RE"]], hi_RE = m[["hi_RE"]],
    same_direction = length(unique(sign(per$r))) == 1)
}
res <- rbindlist(rows)
res[, q_pooled := p.adjust(p_pooled, "BH"), by = model]

cat("\n---- EZH2 vs checkpoints: partial Spearman r (per cohort: TCGA CGGA_325 CGGA_693) ----\n")
for (i in seq_len(nrow(res))) { t <- res[i]
  say("%-6s %-42s | %5.2f %5.2f %5.2f | pooled r = %5.2f [%5.2f, %5.2f] q = %-7s I2 = %3.0f%% | RE r = %5.2f [%5.2f, %5.2f]",
      t$gene, t$model, t$r_TCGA, t$r_CGGA_325, t$r_CGGA_693, t$r_pooled, t$lo, t$hi, fmt_p(t$q_pooled),
      100 * t$I2, t$r_RE, t$lo_RE, t$hi_RE) }

# reverse check: does the Hallmark score itself track LAG3 once EZH2 is accounted for?
cat("\n---- Reverse check: Hallmark cell-cycle score vs LAG3, adjusted for EZH2 + ESTIMATE ----\n")
for (k in names(co)) { v <- pspear("cc_ssgsea", "LAG3", c("EZH2", EST), co[[k]])
  say("  %-9s r = %5.2f  P = %s", k, v[["r"]], fmt_p(v[["p"]])) }

fwrite(res, file.path("results", "Table_S5_proliferation_sensitivity.csv"))
cat("\nOK results/Table_S5_proliferation_sensitivity.csv saved\n")
