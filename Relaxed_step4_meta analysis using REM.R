################################################################################
# Script: 04_PFC_Meta_Analysis_REM.R
# Description: Random-effects model (REML) meta-analysis of study-wise
#              astrocyte pseudo-bulk DEA results from Step 3.
#
# Input:  DEA_Results_limma_voom/DEA_GSE######.csv
# Output: Meta-analysis results with combined beta, SE, FDR, I2, tau2.
#
# Notes:
# - Step 3 DEA files are not pre-filtered by FDR. They contain all genes passing
#   study-level expression filtering.
# - This script applies the final FDR < 0.05 filter after REM meta-analysis.
# - GSE243292 is intentionally excluded because it was excluded from Step 2/3.
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(metafor)
})

step2_dir <- Sys.getenv(
  "PFC_STEP2_DIR",
  unset = file.path("~", "data", "meta_analysis", "results", "PFC", "2026.04")
)
step2_dir <- normalizePath(path.expand(step2_dir), winslash = "/", mustWork = FALSE)

input_dir <- file.path(step2_dir, "DEA_Results_limma_voom")
output_dir <- file.path(step2_dir, "Meta_Analysis_REM")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Keep TRUE to match the previous workflow: meta-analysis only on genes present
# in every loaded study-level DEA result. Set FALSE for a sensitivity analysis
# allowing genes present in at least MIN_STUDIES_FOR_META studies.
REQUIRE_ALL_STUDIES_FOR_META <- TRUE
MIN_STUDIES_FOR_META <- 2

studies <- tibble::tribble(
  ~study,      ~dataset_id,  ~label,       ~file,
  "GSE157827", "gse157827",  "Lau",        "DEA_GSE157827.csv",
  "GSE167490", "gse167490",  "SadickMain", "DEA_GSE167490.csv",
  "GSE167492", "gse167492",  "SadickPilot","DEA_GSE167492.csv",
  "GSE214979", "gse214979",  "Anderson",   "DEA_GSE214979.csv",
  "GSE263468", "gse263468",  "Cobos",      "DEA_GSE263468.csv",
  "GSE268599", "gse268599",  "Serrano",    "DEA_GSE268599.csv",
  "GSE303823", "gse303823",  "Chia",       "DEA_GSE303823.csv"
) %>%
  mutate(path = file.path(input_dir, file))

# ==============================================================================
# 1. Load Step 3 DEA Results
# ==============================================================================

message("Loading study-level DEA results...")

missing_files <- studies %>% filter(!file.exists(path))
if (nrow(missing_files) > 0L) {
  warning(
    "Missing DEA file(s); these studies will be skipped:\n  ",
    paste(missing_files$path, collapse = "\n  ")
  )
}

loaded_studies <- studies %>% filter(file.exists(path))
if (nrow(loaded_studies) < MIN_STUDIES_FOR_META) {
  stop("Fewer than ", MIN_STUDIES_FOR_META,
       " study-level DEA files are available for meta-analysis.")
}

read_dea_file <- function(path, study, dataset_id, label) {
  df <- read_csv(path, show_col_types = FALSE)
  required_cols <- c("gene", "logFC", "SE", "PValue", "FDR")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop("[", study, "] Missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  
  df %>%
    transmute(
      gene = as.character(gene),
      study = study,
      dataset_id = dataset_id,
      label = label,
      logFC = as.numeric(logFC),
      SE = as.numeric(SE),
      PValue = as.numeric(PValue),
      FDR = as.numeric(FDR)
    ) %>%
    filter(
      !is.na(gene),
      gene != "",
      is.finite(logFC),
      is.finite(SE),
      SE > 0,
      is.finite(PValue)
    ) %>%
    distinct(gene, study, .keep_all = TRUE)
}

dea_long <- pmap_dfr(
  loaded_studies %>% select(path, study, dataset_id, label),
  read_dea_file
)

write_csv(dea_long, file.path(output_dir, "PFC_meta_input_long.csv"))

study_summary <- dea_long %>%
  count(study, dataset_id, label, name = "n_genes") %>%
  arrange(study)
write_csv(study_summary, file.path(output_dir, "PFC_meta_input_study_summary.csv"))

message(sprintf("  Loaded studies: %d", n_distinct(dea_long$study)))
print(study_summary)

gene_coverage <- dea_long %>%
  count(gene, name = "n_studies") %>%
  arrange(desc(n_studies), gene)
write_csv(gene_coverage, file.path(output_dir, "PFC_meta_gene_coverage.csv"))

if (REQUIRE_ALL_STUDIES_FOR_META) {
  required_n <- n_distinct(dea_long$study)
} else {
  required_n <- MIN_STUDIES_FOR_META
}

eligible_genes <- gene_coverage %>%
  filter(n_studies >= required_n) %>%
  pull(gene)

message(sprintf(
  "  Genes eligible for REM meta-analysis: %d (required studies per gene: %d)",
  length(eligible_genes), required_n
))

if (length(eligible_genes) < 100L) {
  stop("Too few eligible genes. Check gene identifiers and Step 3 DEA outputs.")
}

meta_input <- dea_long %>%
  filter(gene %in% eligible_genes) %>%
  arrange(gene, study)

meta_input_wide <- meta_input %>%
  select(gene, study, logFC, SE, PValue, FDR) %>%
  tidyr::pivot_wider(
    names_from = study,
    values_from = c(logFC, SE, PValue, FDR),
    names_sep = "_"
  )
write_csv(meta_input_wide, file.path(output_dir, "PFC_meta_input_wide.csv"))

# ==============================================================================
# 2. Random-Effects Meta-Analysis (REML)
# ==============================================================================

message("Running REML random-effects meta-analysis...")

run_gene_meta <- function(df_gene) {
  gene_id <- df_gene$gene[1]
  df_gene <- df_gene %>%
    filter(is.finite(logFC), is.finite(SE), SE > 0)
  
  if (nrow(df_gene) < MIN_STUDIES_FOR_META) {
    return(NULL)
  }
  
  fit <- tryCatch(
    metafor::rma(yi = logFC, sei = SE, data = df_gene, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NULL)
  }
  
  tibble(
    gene = gene_id,
    beta_meta = as.numeric(fit$beta),
    se_meta = fit$se,
    z_meta = fit$zval,
    p_meta = fit$pval,
    ci_lb = fit$ci.lb,
    ci_ub = fit$ci.ub,
    I2 = fit$I2,
    tau2 = fit$tau2,
    QE = fit$QE,
    QEp = fit$QEp,
    n_studies = fit$k,
    studies = paste(df_gene$study, collapse = ";")
  )
}

meta_results <- meta_input %>%
  group_split(gene) %>%
  map_dfr(run_gene_meta)

if (nrow(meta_results) == 0L) {
  stop("No genes successfully completed REM meta-analysis.")
}

meta_final <- meta_results %>%
  mutate(
    padj_meta = p.adjust(p_meta, method = "fdr"),
    direction = case_when(
      beta_meta > 0 ~ "Up",
      beta_meta < 0 ~ "Down",
      TRUE ~ "No_change"
    )
  ) %>%
  arrange(p_meta)

write_csv(meta_final, file.path(output_dir, "PFC_metaDEGs_All_Results.csv"))

message(sprintf(
  "  REM meta-analysis complete: %d genes", nrow(meta_final)
))

# ==============================================================================
# 3. Significant Meta-DEG Tables
# ==============================================================================

message("Filtering significant meta-DEGs (FDR < 0.05)...")

meta_degs <- meta_final %>% filter(padj_meta < 0.05)
meta_up <- meta_degs %>% filter(beta_meta > 0)
meta_down <- meta_degs %>% filter(beta_meta < 0)

message(sprintf(
  "  Significant meta-DEGs: %d (Up=%d / Down=%d)",
  nrow(meta_degs), nrow(meta_up), nrow(meta_down)
))

write_csv(meta_degs, file.path(output_dir, "PFC_metaDEGs_Significant.csv"))
write_csv(meta_up, file.path(output_dir, "PFC_metaDEGs_Upregulated.csv"))
write_csv(meta_down, file.path(output_dir, "PFC_metaDEGs_Downregulated.csv"))

# ==============================================================================
# 4. Optional SAG Overlap Validation
# ==============================================================================

message("Checking optional SAG overlap list...")

sags_path <- Sys.getenv(
  "PFC_SAG_LIST",
  unset = file.path("~", "data", "meta_analysis", "NewSAGs_Unique.csv")
)
sags_path <- normalizePath(path.expand(sags_path), winslash = "/", mustWork = FALSE)

if (file.exists(sags_path)) {
  sags_df <- read_csv(sags_path, show_col_types = FALSE)
  sags_col <- colnames(sags_df)[grep("gene", colnames(sags_df), ignore.case = TRUE)][1]
  
  if (!is.na(sags_col)) {
    sag_genes <- toupper(trimws(unique(as.character(sags_df[[sags_col]]))))
    meta_genes <- toupper(trimws(meta_final$gene))
    sig_genes <- toupper(trimws(meta_degs$gene))
    
    sag_overlap_all <- intersect(meta_genes, sag_genes)
    sag_overlap_sig <- intersect(sig_genes, sag_genes)
    
    sag_summary <- tibble(
      sag_file = sags_path,
      sag_gene_column = sags_col,
      n_sag_genes = length(sag_genes),
      n_meta_genes = length(meta_genes),
      n_significant_meta_degs = length(sig_genes),
      n_sag_overlap_all_meta_genes = length(sag_overlap_all),
      n_sag_overlap_significant_meta_degs = length(sag_overlap_sig)
    )
    
    write_csv(sag_summary, file.path(output_dir, "SAG_overlap_summary.csv"))
    write_csv(
      tibble(gene = sag_overlap_sig),
      file.path(output_dir, "SAG_overlap_significant_metaDEGs.csv")
    )
    
    message(sprintf(
      "  SAG overlap among significant meta-DEGs: %d",
      length(sag_overlap_sig)
    ))
  } else {
    warning("SAG overlap skipped: no gene-like column found in SAG list.")
  }
} else {
  warning("SAG overlap skipped: SAG file not found at ", sags_path)
}

message("\nStep 4 REM meta-analysis complete.")
