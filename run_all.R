# =====================================================================
# run_all.R  -  runs the full EZH2 GBM pipeline in order.
# Requirements: the "share/" folder (exported data + scores_*.csv) and
# these scripts in the SAME project folder. Then: source("run_all.R")
# =====================================================================
scripts <- c("01_build_cohorts.R", "02_expression_prc2_proliferation.R", "03_survival.R",
             "04_checkpoints_partial_correlation.R", "05_immune_landscape.R", "06_single_cell_TISCH2.R")
log_file <- "results/run_log.txt"; dir.create("results", showWarnings = FALSE)
sink(log_file, split = TRUE)
for (s in scripts) { cat("\n\n##########", s, "##########\n"); source(s, echo = FALSE, local = new.env()) }
cat("\n\n########## sessionInfo ##########\n"); print(sessionInfo())
sink()
cat("\n✓ Pipeline complete. Results in results/, figures in figures/, log in", log_file, "\n")
