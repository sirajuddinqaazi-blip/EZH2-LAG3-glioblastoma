# EZH2, LAG3 and PDCD1 in IDH-wildtype glioblastoma — analysis code

Code for: *"EZH2 expression is associated with LAG3 and PDCD1 independently of tumor proliferation in IDH-wildtype glioblastoma: a three-cohort analysis"* (submitted).

## Requirements
R ≥ 4.3 (analysis run in R 4.6.0), internet access for the first step, ~3 GB disk space.
Packages are installed automatically by the scripts (TCGAbiolinks, SummarizedExperiment, estimate, immunedeconv, survival, survminer, ggplot2, patchwork, ggrepel, dplyr, stringr).

## How to reproduce (run from the project folder, in this order)

| Step | Script | What it does |
|---|---|---|
| 1 | `00a_download_and_prepare_data.R` | Downloads TCGA-GBM RNA-seq (GDC, STAR – Counts), builds the IDH-wildtype primary GBM cohort (n = 160), and exports compact input files to `share/`. CGGA files must be downloaded manually first (see below). |
| 2 | `11_estimate_mcp_scores.R` | Computes ESTIMATE and MCP-counter scores for all cohorts → `share/scores_*.csv` |
| 3 | TISCH2 files | Download the LAG3, EZH2 and PDCD1 cell-type heatmap tables (Celltype major-lineage, glioma datasets) from http://tisch.compbio.cn and save as `share/TISCH_<GENE>_heatmap.csv` |
| 4 | `run_all.R` | Runs scripts 01–06 and writes all figures (`figures/`), tables (`results/`) and a log with sessionInfo (`results/run_log.txt`) |

### CGGA data (manual download)
From http://www.cgga.org.cn download and unzip into `Data/geo/`:
- `CGGA.mRNAseq_325.RSEM-genes.20200506.txt`, `CGGA.mRNAseq_325_clinical.20200506.txt`
- `CGGA.mRNAseq_693.RSEM-genes.20200506.txt`, `CGGA.mRNAseq_693_clinical.20200506.txt`

## Analysis scripts
| Script | Analysis | Outputs |
|---|---|---|
| 00_config_and_functions.R | Shared settings and functions (partial Spearman correlation, fixed-effect meta-analysis) | – |
| 01_build_cohorts.R | TCGA (n = 160), CGGA_325 (n = 74), CGGA_693 (n = 109) IDH-wildtype primary GBM; proliferation score | results/cohorts.rds |
| 02_expression_prc2_proliferation.R | Tumor vs normal, PRC2 co-expression, proliferation, clinical correlates | Figure 1, Table S1 |
| 03_survival.R | Kaplan–Meier, Cox models, proportional-hazards checks, sensitivity analyses, pooled HR | Figure 2, Table 2, Table S2 |
| 04_checkpoints_partial_correlation.R | EZH2 vs 8 immune checkpoints: raw and confounder-adjusted, meta-analysis | Figure 3, Table S3 |
| 05_immune_landscape.R | EZH2 vs ESTIMATE / MCP-counter populations | Figure S1, Table S3 |
| 06_single_cell_TISCH2.R | Cell-type expression in GSE131928 (TISCH2) | Figure 4, Table S4 |

## Data sources
TCGA: Genomic Data Commons (https://portal.gdc.cancer.gov). CGGA: http://www.cgga.org.cn. Single-cell: TISCH2 (GEO GSE131928).
No patient-level data are redistributed in this repository.

## Note on data versions
The cohort sizes printed by step 1 (289 → 185 → 162 → 160 patients) correspond to the GDC release used for the manuscript. Later GDC releases may differ slightly; the script warns if so and records the release number in `Data/processed/GDC_release.txt`.

## License
MIT
