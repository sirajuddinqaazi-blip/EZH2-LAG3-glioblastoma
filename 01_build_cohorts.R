# =====================================================================
# 01_build_cohorts.R
# Builds three analysis-ready cohorts (IDH-wildtype, primary GBM):
#   TCGA-GBM (n=160), CGGA_325 (n=74), CGGA_693 (n=109)
# Output: results/cohorts.rds  (list of 3 data frames + TCGA normals)
# =====================================================================
source("00_config_and_functions.R")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ------------------------- TCGA -------------------------
tpm <- readRDS(file.path(DATA_DIR, "tcga_tpm.rds"))
gi  <- readRDS(file.path(DATA_DIR, "tcga_gene_info.rds"))
cl  <- readRDS(file.path(DATA_DIR, "tcga_cohort_clinical.rds"))
nc  <- readRDS(file.path(DATA_DIR, "tcga_normals_clinical.rds"))
sym <- gi$gene_name[match(rownames(tpm), rownames(gi))]

get_tcga <- function(genes, samples) {
  sapply(genes, function(g) {
    i <- which(sym == g)
    if (length(i) == 0) stop("Gene not found in TCGA: ", g)
    i <- i[which.max(rowMeans(tpm[i, , drop = FALSE]))]   # highest-expressed Ensembl ID
    as.numeric(tpm[i, samples])
  })
}

tcga <- data.frame(sample = rownames(cl), get_tcga(GENES_NEEDED, rownames(cl)), check.names = FALSE)
os_days <- ifelse(cl$vital_status == "Dead", cl$days_to_death, cl$days_to_last_follow_up)
tcga <- within(tcga, {
  OS_months <- os_days / 30.44
  event     <- as.integer(cl$vital_status == "Dead")
  age       <- num(cl$age_at_index)
  sex       <- factor(tolower(cl$paper_Gender), levels = c("female", "male"))
  MGMT      <- factor(cl$paper_MGMT.promoter.status, levels = c("Methylated", "Unmethylated"))
  KPS       <- num(cl$paper_Karnofsky.Performance.Score)
  subtype   <- factor(cl$paper_Transcriptome.Subtype, levels = c("CL", "ME", "NE", "PN"),
                      labels = c("Classical", "Mesenchymal", "Neural", "Proneural"))
  ABSOLUTE_purity <- num(cl$paper_ABSOLUTE.purity)
})
tcga$cohort <- "TCGA"

normals <- data.frame(sample = rownames(nc), get_tcga("EZH2", rownames(nc)))
colnames(normals)[2] <- "EZH2"

# ------------------------- CGGA -------------------------
build_cgga <- function(id) {
  ex <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_expr_gbm.rds")))
  c2 <- readRDS(file.path(DATA_DIR, paste0("cgga", id, "_clinical.rds")))
  for (cc in names(c2)) if (is.character(c2[[cc]])) c2[[cc]] <- trimws(c2[[cc]])
  keep <- c2$PRS_type == "Primary" & c2$Histology == "GBM" & c2$Grade == "WHO IV" &
          c2$IDH_mutation_status %in% "Wildtype"
  ids <- intersect(c2$CGGA_ID[keep], colnames(ex))
  miss <- setdiff(GENES_NEEDED, ex[[1]]); if (length(miss)) stop("Missing in CGGA_", id, ": ", paste(miss, collapse = ","))
  m <- ex[match(GENES_NEEDED, ex[[1]]), ids]
  d <- data.frame(sample = ids, t(apply(m, 2, as.numeric)), check.names = FALSE)
  colnames(d)[-1] <- GENES_NEEDED
  c2 <- c2[match(ids, c2$CGGA_ID), ]
  d$OS_months <- c2$OS / 30.44
  d$event     <- c2$Censor               # CGGA: 1 = dead, 0 = alive
  d$age       <- c2$Age
  d$sex       <- factor(tolower(c2$Gender), levels = c("female", "male"))
  d$MGMT      <- factor(c2$MGMTp_methylation_status, levels = c("methylated", "un-methylated"),
                        labels = c("Methylated", "Unmethylated"))
  d$radio     <- c2$Radio_status; d$chemo <- c2$Chemo_status
  d$cohort    <- paste0("CGGA_", id)
  d
}
cohorts <- list(TCGA = tcga, CGGA_325 = build_cgga("325"), CGGA_693 = build_cgga("693"))

# ------------------- scores, derived variables -------------------
score_files <- c(TCGA = "scores_TCGA.csv", CGGA_325 = "scores_CGGA325.csv", CGGA_693 = "scores_CGGA693.csv")
for (k in names(cohorts)) {
  d  <- cohorts[[k]]
  sc <- read.csv(file.path(DATA_DIR, score_files[[k]]), check.names = FALSE)
  d  <- merge(d, sc, by = "sample", all.x = TRUE, sort = FALSE)
  if (any(is.na(d$ESTIMATE_immune))) stop("Missing ESTIMATE scores in ", k)
  d$prolif     <- rowMeans(scale(log2(d[, PROLIF_GENES] + 1)))   # 10-gene proliferation score
  d$log2EZH2   <- log2(d$EZH2 + 1)
  d$EZH2_group <- factor(ifelse(d$EZH2 >= median(d$EZH2), "High", "Low"), levels = c("Low", "High"))
  cohorts[[k]] <- d
}

saveRDS(list(cohorts = cohorts, tcga_normals = normals), file.path(RES_DIR, "cohorts.rds"))
for (k in names(cohorts)) cat(sprintf("%-9s n = %3d | events = %3d | EZH2 median = %.2f\n",
  k, nrow(cohorts[[k]]), sum(cohorts[[k]]$event, na.rm = TRUE), median(cohorts[[k]]$EZH2)))
cat("TCGA normal brain n =", nrow(normals), "\n✓ results/cohorts.rds saved\n")
