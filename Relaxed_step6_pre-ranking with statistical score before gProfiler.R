################################################################################
# Script: 06_PFC_PreRanking_for_gProfiler.R
# Description: Generate the minimal pre-ranked files needed for g:Profiler,
#              GSEA, and EnrichmentMap from Step 4/5 outputs.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ==============================================================================
# 0. Setup
# ==============================================================================

step2_dir <- Sys.getenv(
  "PFC_STEP2_DIR",
  unset = file.path("~", "data", "meta_analysis", "results", "PFC", "2026.04")
)
step2_dir <- normalizePath(path.expand(step2_dir), winslash = "/", mustWork = FALSE)

input_meta_dir <- file.path(step2_dir, "Meta_Analysis_REM")
input_overlap_dir <- file.path(step2_dir, "SADEG_Overlap_Results")
output_dir <- file.path(step2_dir, "PreRanked_gProfiler_EnrichmentMap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

CLEAN_OUTPUT_DIR <- TRUE
if (CLEAN_OUTPUT_DIR) {
  old_outputs <- list.files(output_dir, full.names = TRUE)
  if (length(old_outputs) > 0L) {
    unlink(old_outputs, recursive = TRUE, force = TRUE)
  }
}

message("Generating minimal pre-ranked lists...")

load_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop("[", label, "] File not found:\n  ", path)
  }
  read_csv(path, show_col_types = FALSE)
}

check_required_cols <- function(df, label) {
  missing_cols <- setdiff(c("gene", "beta_meta"), colnames(df))
  if (length(missing_cols) > 0L) {
    stop("[", label, "] Missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
}

write_gene_list <- function(genes, path) {
  genes <- genes[!is.na(genes) & genes != ""]
  writeLines(genes, con = path)
}

# ==============================================================================
# 1. All Meta-Analysis Genes Ranked by Absolute log2FC
# ==============================================================================

message("  [1] Ranking all REM meta-analysis genes by absolute log2FC...")

df_all <- load_required_csv(
  file.path(input_meta_dir, "PFC_metaDEGs_All_Results.csv"),
  "All REM meta-analysis results"
)
check_required_cols(df_all, "All REM meta-analysis results")

df_all_ranked <- df_all %>%
  mutate(
    gene = as.character(gene),
    abs_log2fc = abs(as.numeric(beta_meta))
  ) %>%
  filter(!is.na(gene), gene != "", is.finite(abs_log2fc)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(abs_log2fc)) %>%
  select(gene, abs_log2fc, beta_meta, p_meta, padj_meta, everything())

write_csv(
  df_all_ranked,
  file.path(output_dir, "PFC_metaDEGs_All_Results_log2FC_Ranked.csv")
)
write_gene_list(
  df_all_ranked$gene,
  file.path(output_dir, "PFC_All_Genes_log2FC_Ranked.txt")
)

message(sprintf("      Saved all-gene ranked list: %d genes", nrow(df_all_ranked)))

# ==============================================================================
# 2. SADEG Up/Down Gene Lists for g:Profiler
# ==============================================================================

message("  [2] Ranking SADEGs by signed log2FC...")

df_overlap <- load_required_csv(
  file.path(input_overlap_dir, "PFC_SAGs_Overlap_Detailed.csv"),
  "SADEG overlap table"
)
check_required_cols(df_overlap, "SADEG overlap table")

if (!"direction" %in% colnames(df_overlap)) {
  df_overlap$direction <- ifelse(df_overlap$beta_meta > 0, "Up", "Down")
}

df_overlap <- df_overlap %>%
  mutate(
    gene = as.character(gene),
    beta_meta = as.numeric(beta_meta),
    direction = as.character(direction)
  ) %>%
  filter(!is.na(gene), gene != "", is.finite(beta_meta)) %>%
  distinct(gene, .keep_all = TRUE)

df_overlap_up <- df_overlap %>%
  filter(direction == "Up" | beta_meta > 0) %>%
  arrange(desc(beta_meta)) %>%
  select(gene, direction, beta_meta, padj_meta, everything())

df_overlap_down <- df_overlap %>%
  filter(direction == "Down" | beta_meta < 0) %>%
  arrange(beta_meta) %>%
  select(gene, direction, beta_meta, padj_meta, everything())

write_csv(
  df_overlap_up,
  file.path(output_dir, "PFC_MetaDEGs_SAGs_Overlap_UP_log2FC_Ranked.csv")
)
write_gene_list(
  df_overlap_up$gene,
  file.path(output_dir, "PFC_MetaDEGs_SAGs_Overlap_UP_Genes.txt")
)

write_csv(
  df_overlap_down,
  file.path(output_dir, "PFC_MetaDEGs_SAGs_Overlap_DOWN_log2FC_Ranked.csv")
)
write_gene_list(
  df_overlap_down$gene,
  file.path(output_dir, "PFC_MetaDEGs_SAGs_Overlap_DOWN_Genes.txt")
)

message(sprintf(
  "      Saved SADEG gene lists: Up=%d | Down=%d",
  nrow(df_overlap_up), nrow(df_overlap_down)
))

# ==============================================================================
# 3. Signed log2FC .rnk for GSEA / EnrichmentMap
# ==============================================================================

message("  [3] Writing signed log2FC .rnk for all SADEGs...")

df_overlap_signed <- df_overlap %>%
  arrange(desc(beta_meta)) %>%
  select(gene, beta_meta)

write.table(
  df_overlap_signed,
  file.path(output_dir, "PFC_MetaDEGs_SAGs_Overlap_Signed_log2FC.rnk"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message(sprintf("      Saved signed .rnk: %d genes", nrow(df_overlap_signed)))
message("Step 6 pre-ranking complete.")
