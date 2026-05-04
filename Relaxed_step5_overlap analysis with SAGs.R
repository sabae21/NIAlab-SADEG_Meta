################################################################################
# Script: 05_PFC_Overlap_Analysis_with_SAGs.R
# Description: Overlap analysis between PFC astrocyte meta-DEGs and
#              senescence-associated genes (SAGs) to define SADEGs.
#
# Input:
# - Step 4 REM meta-analysis outputs in Meta_Analysis_REM/
# - SAG reference list, default: ~/data/meta_analysis/NewSAGs_Unique.csv
#
# Output:
# - SADEG gene tables and hypergeometric enrichment summaries.
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

step2_dir <- Sys.getenv(
  "PFC_STEP2_DIR",
  unset = file.path("~", "data", "meta_analysis", "results", "PFC", "2026.04")
)
step2_dir <- normalizePath(path.expand(step2_dir), winslash = "/", mustWork = FALSE)

meta_dir <- file.path(step2_dir, "Meta_Analysis_REM")
output_dir <- file.path(step2_dir, "SADEG_Overlap_Results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sags_path <- Sys.getenv(
  "PFC_SAG_LIST",
  unset = file.path("~", "data", "meta_analysis", "NewSAGs_Unique.csv")
)
sags_path <- normalizePath(path.expand(sags_path), winslash = "/", mustWork = FALSE)

# ==============================================================================
# 1. Helper Functions
# ==============================================================================

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

find_gene_column <- function(df, file_label) {
  candidates <- grep("^(gene|genes|symbol|gene_symbol|hgnc_symbol)$",
                     colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(candidates) == 0L) {
    candidates <- grep("gene|symbol", colnames(df), ignore.case = TRUE, value = TRUE)
  }
  if (length(candidates) == 0L) {
    stop("[", file_label, "] Could not identify a gene/symbol column.")
  }
  candidates[1L]
}

load_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop("[", label, "] File not found:\n  ", path)
  }
  read_csv(path, show_col_types = FALSE)
}

run_hypergeo_test <- function(target_genes, universe_genes, sag_genes, label) {
  universe_genes <- unique(normalize_gene(universe_genes))
  target_genes <- intersect(unique(normalize_gene(target_genes)), universe_genes)
  sag_in_universe <- intersect(unique(normalize_gene(sag_genes)), universe_genes)
  
  N <- length(universe_genes)
  M <- length(sag_in_universe)
  n <- length(target_genes)
  overlap <- intersect(target_genes, sag_in_universe)
  k <- length(overlap)
  
  p_value <- if (N > 0L && M > 0L && n > 0L) {
    phyper(k - 1L, M, N - M, n, lower.tail = FALSE)
  } else {
    NA_real_
  }
  
  expected_k <- if (N > 0L) n * M / N else NA_real_
  fold_enrichment <- if (!is.na(expected_k) && expected_k > 0) {
    k / expected_k
  } else {
    NA_real_
  }
  
  list(
    summary = tibble(
      category = label,
      universe_N = N,
      SAGs_in_universe_M = M,
      DEG_set_n = n,
      overlap_k = k,
      expected_overlap = expected_k,
      fold_enrichment = fold_enrichment,
      p_value = p_value
    ),
    overlap_genes = overlap
  )
}

# ==============================================================================
# 2. Load Meta-DEGs and SAG Reference
# ==============================================================================

message("Loading Step 4 meta-analysis results and SAG reference list...")

meta_all <- load_required_csv(
  file.path(meta_dir, "PFC_metaDEGs_All_Results.csv"),
  "All REM meta-analysis results"
)
meta_sig <- load_required_csv(
  file.path(meta_dir, "PFC_metaDEGs_Significant.csv"),
  "Significant meta-DEGs"
)
meta_up <- load_required_csv(
  file.path(meta_dir, "PFC_metaDEGs_Upregulated.csv"),
  "Upregulated significant meta-DEGs"
)
meta_down <- load_required_csv(
  file.path(meta_dir, "PFC_metaDEGs_Downregulated.csv"),
  "Downregulated significant meta-DEGs"
)

required_meta_cols <- c("gene", "beta_meta", "padj_meta")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_all))
if (length(missing_meta_cols) > 0L) {
  stop("Meta-analysis all-results file is missing required column(s): ",
       paste(missing_meta_cols, collapse = ", "))
}

sags_raw <- load_required_csv(sags_path, "SAG reference list")
sags_col <- find_gene_column(sags_raw, "SAG reference list")

sags_tbl <- sags_raw %>%
  mutate(
    sag_gene_original = as.character(.data[[sags_col]]),
    gene_key = normalize_gene(.data[[sags_col]])
  ) %>%
  filter(!is.na(gene_key), gene_key != "") %>%
  distinct(gene_key, .keep_all = TRUE)

sag_genes <- sags_tbl$gene_key

meta_all <- meta_all %>%
  mutate(gene_key = normalize_gene(gene)) %>%
  filter(!is.na(gene_key), gene_key != "") %>%
  distinct(gene_key, .keep_all = TRUE)

meta_sig <- meta_sig %>%
  mutate(gene_key = normalize_gene(gene)) %>%
  filter(!is.na(gene_key), gene_key != "") %>%
  distinct(gene_key, .keep_all = TRUE)
meta_up <- meta_up %>%
  mutate(gene_key = normalize_gene(gene)) %>%
  filter(!is.na(gene_key), gene_key != "") %>%
  distinct(gene_key, .keep_all = TRUE)
meta_down <- meta_down %>%
  mutate(gene_key = normalize_gene(gene)) %>%
  filter(!is.na(gene_key), gene_key != "") %>%
  distinct(gene_key, .keep_all = TRUE)

universe_genes <- meta_all$gene_key

message(sprintf("  Meta-analysis universe genes: %d", length(universe_genes)))
message(sprintf("  Significant meta-DEGs: %d", nrow(meta_sig)))
message(sprintf("  SAG reference genes: %d", length(sag_genes)))
message(sprintf(
  "  SAGs in meta-analysis universe: %d",
  length(intersect(sag_genes, universe_genes))
))

# ==============================================================================
# 3. Define SADEGs
# ==============================================================================

message("Defining SADEGs as significant meta-DEGs overlapping SAGs...")

sadeg_all <- meta_sig %>%
  filter(gene_key %in% sag_genes) %>%
  mutate(
    SADEG = TRUE,
    direction = case_when(
      beta_meta > 0 ~ "Up",
      beta_meta < 0 ~ "Down",
      TRUE ~ "No_change"
    )
  ) %>%
  left_join(
    sags_tbl %>% select(gene_key, sag_gene_original),
    by = "gene_key"
  ) %>%
  arrange(padj_meta, desc(abs(beta_meta)))

sadeg_up <- sadeg_all %>% filter(direction == "Up")
sadeg_down <- sadeg_all %>% filter(direction == "Down")

message(sprintf(
  "  SADEGs: %d (Up=%d / Down=%d)",
  nrow(sadeg_all), nrow(sadeg_up), nrow(sadeg_down)
))

write_csv(sadeg_all, file.path(output_dir, "PFC_SADEGs_All.csv"))
write_csv(sadeg_up, file.path(output_dir, "PFC_SADEGs_Upregulated.csv"))
write_csv(sadeg_down, file.path(output_dir, "PFC_SADEGs_Downregulated.csv"))

# Also export the same table under the older descriptive naming convention.
write_csv(sadeg_all, file.path(output_dir, "PFC_SAGs_Overlap_Detailed.csv"))

# ==============================================================================
# 4. Hypergeometric Enrichment Tests
# ==============================================================================

message("Running hypergeometric enrichment tests...")

res_all <- run_hypergeo_test(
  target_genes = meta_sig$gene_key,
  universe_genes = universe_genes,
  sag_genes = sag_genes,
  label = "All significant meta-DEGs"
)
res_up <- run_hypergeo_test(
  target_genes = meta_up$gene_key,
  universe_genes = universe_genes,
  sag_genes = sag_genes,
  label = "Upregulated significant meta-DEGs"
)
res_down <- run_hypergeo_test(
  target_genes = meta_down$gene_key,
  universe_genes = universe_genes,
  sag_genes = sag_genes,
  label = "Downregulated significant meta-DEGs"
)

hypergeo_summary <- bind_rows(
  res_all$summary,
  res_up$summary,
  res_down$summary
) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

write_csv(
  hypergeo_summary,
  file.path(output_dir, "PFC_SAGs_Hypergeometric_Summary.csv")
)

# ==============================================================================
# 5. Background and Reporting Tables
# ==============================================================================

sag_universe_tbl <- meta_all %>%
  filter(gene_key %in% sag_genes) %>%
  left_join(
    sags_tbl %>% select(gene_key, sag_gene_original),
    by = "gene_key"
  ) %>%
  arrange(gene)

write_csv(
  sag_universe_tbl,
  file.path(output_dir, "PFC_SAGs_in_MetaAnalysis_Universe.csv")
)

overlap_sets <- tibble(
  category = c("All significant meta-DEGs",
               "Upregulated significant meta-DEGs",
               "Downregulated significant meta-DEGs"),
  overlap_genes = c(
    paste(sort(res_all$overlap_genes), collapse = ";"),
    paste(sort(res_up$overlap_genes), collapse = ";"),
    paste(sort(res_down$overlap_genes), collapse = ";")
  )
)

write_csv(
  overlap_sets,
  file.path(output_dir, "PFC_SAGs_Overlap_GeneSets.csv")
)

message("\nStep 5 SAG overlap analysis complete.")
