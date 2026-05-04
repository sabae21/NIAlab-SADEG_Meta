################################################################################
# Script: 01_PFC_Data_Loading_and_Preprocessing.R
# Description: Optionally downloads GEO supplementary files, loads raw
#              snRNA-seq data for 8 AD PFC datasets, performs
#              metadata curation using relaxed criteria
#              (AD/control diagnosis labels + age >= 65, no Braak requirement),
#              QC filtering, and saves Seurat objects for downstream analysis.
# Datasets:    GSE157827, GSE167490, GSE167492, GSE214979,
#              GSE263468, GSE268599, GSE303823
################################################################################

# ==============================================================================
# 0. Global configuration and optional data acquisition
# ==============================================================================

# Keep installation outside analysis scripts for reproducibility. Install these
# packages once before running:
# install.packages(c("dplyr", "readr", "stringr", "openxlsx", "Seurat", "Matrix"))
# install.packages("R.utils")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GEOquery", "zellkonverter", "SingleCellExperiment"))

# Set RUN_DOWNLOADS <- TRUE when running a fresh clone or adding a new dataset.
# DOWNLOAD_ALL_AT_START should usually stay FALSE so each section downloads only
# the dataset it needs.
RUN_DOWNLOADS <- FALSE
DOWNLOAD_ALL_AT_START <- FALSE
RUN_PREP_LARGE_FILES <- FALSE

control_terms <- c("control", "unaffected", "ns", "non-symptomatic",
                   "non symptomatic", "nc", "nci", "neurologically normal",
                   "normal")
ad_terms <- c("ad", "alzheimer's disease", "alzheimer disease",
              "alzheimer’s disease", "add", "probable ad", "possible ad")

data_root <- Sys.getenv(
  "PFC_RAW_DATA_DIR",
  unset = file.path("~", "data", "meta_analysis", "raw_data", "PFC")
)
output_dir <- Sys.getenv(
  "PFC_RESULTS_DIR",
  unset = file.path("~", "data", "meta_analysis", "results", "PFC", "2026.04")
)
data_root <- normalizePath(path.expand(data_root), winslash = "/", mustWork = FALSE)
output_dir <- normalizePath(path.expand(output_dir), winslash = "/", mustWork = FALSE)
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

join_layers_if_available <- function(obj) {
  if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
    Seurat::JoinLayers(obj)
  } else {
    obj
  }
}

merge_seurat_list <- function(seurat_list, project) {
  if (length(seurat_list) == 0L) {
    stop(sprintf("[%s] No samples loaded.", project))
  }
  if (length(seurat_list) == 1L) {
    return(seurat_list[[1L]])
  }
  merge(
    x = seurat_list[[1L]],
    y = seurat_list[-1L],
    add.cell.ids = names(seurat_list),
    project = project
  )
}

not_reported_if_missing <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Not reported"
  x
}

normalize_diagnosis <- function(x) {
  x <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
  x <- stringr::str_replace_all(x, "\u2019", "'")
  x
}

assign_relaxed_condition <- function(diagnosis) {
  dx <- normalize_diagnosis(diagnosis)
  dplyr::case_when(
    dx %in% ad_terms |
      stringr::str_detect(dx, "alzheimer") ~ "AD",
    dx %in% control_terms ~ "Control",
    TRUE ~ NA_character_
  )
}

first_existing_column <- function(df, candidates, default = NA_character_) {
  hit <- intersect(candidates, colnames(df))[1]
  if (is.na(hit)) {
    rep(default, nrow(df))
  } else {
    as.character(df[[hit]])
  }
}

count_unique_ids <- function(df, id_col) {
  if (nrow(df) == 0L || is.null(id_col) || !id_col %in% colnames(df)) {
    return(0L)
  }
  ids <- as.character(df[[id_col]])
  length(unique(ids[!is.na(ids) & trimws(ids) != ""]))
}

count_condition_rows <- function(df, condition_value) {
  if (nrow(df) == 0L || !"condition" %in% colnames(df)) {
    return(0L)
  }
  sum(df$condition == condition_value, na.rm = TRUE)
}

count_condition_donors <- function(df, donor_col, condition_value) {
  if (nrow(df) == 0L || !donor_col %in% colnames(df) ||
      !"condition" %in% colnames(df)) {
    return(0L)
  }
  count_unique_ids(df[df$condition == condition_value, , drop = FALSE], donor_col)
}

count_excluded_donors <- function(df, donor_col) {
  if (nrow(df) == 0L || !donor_col %in% colnames(df) ||
      !"included_in_analysis" %in% colnames(df)) {
    return(0L)
  }
  count_unique_ids(df[df$included_in_analysis == "NO", , drop = FALSE], donor_col)
}

download_geo_supp <- function(gse_id, filter_regex = NULL, untar_archives = TRUE) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Install GEOquery before running downloads: BiocManager::install('GEOquery')")
  }
  
  args <- list(GEO = gse_id, baseDir = data_root, makeDirectory = TRUE)
  if (!is.null(filter_regex)) args$filter_regex <- filter_regex
  do.call(GEOquery::getGEOSuppFiles, args)
  
  gse_dir <- normalizePath(file.path(data_root, gse_id), winslash = "/", mustWork = FALSE)
  dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)
  if (untar_archives && dir.exists(gse_dir)) {
    tar_files <- list.files(gse_dir, pattern = "\\.tar$", full.names = TRUE)
    for (tar_file in tar_files) {
      utils::untar(tar_file, exdir = gse_dir, tar = "internal")
    }
  }
  invisible(gse_dir)
}

has_files_matching <- function(data_dir, patterns) {
  all(vapply(
    patterns,
    function(pattern) {
      length(list.files(data_dir, pattern = pattern, recursive = TRUE)) > 0L
    },
    logical(1)
  ))
}

ensure_geo_supp_files <- function(gse_id, required_patterns,
                                  filter_regex = NULL,
                                  untar_archives = TRUE) {
  gse_dir <- file.path(data_root, gse_id)
  dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (has_files_matching(gse_dir, required_patterns)) {
    return(invisible(gse_dir))
  }
  
  if (!RUN_DOWNLOADS) {
    stop(
      "[", gse_id, "] Required raw files are missing under:\n  ", gse_dir,
      "\nSet RUN_DOWNLOADS <- TRUE and rerun the data acquisition section, ",
      "or download GEO supplementary files manually."
    )
  }
  
  message("[", gse_id, "] Raw files missing; downloading GEO supplementary files...")
  download_geo_supp(gse_id, filter_regex = filter_regex, untar_archives = untar_archives)
  
  if (!has_files_matching(gse_dir, required_patterns)) {
    stop(
      "[", gse_id, "] Download finished, but required files are still missing. ",
      "Check GEO supplementary file names and extraction under: ", gse_dir
    )
  }
  
  invisible(gse_dir)
}

find_series_matrix <- function(gse_id, data_dir) {
  candidates <- list.files(
    data_dir,
    pattern = paste0(gse_id, ".*series_matrix.*\\.txt(\\.gz)?$"),
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(candidates) == 0L) NA_character_ else candidates[1]
}

load_geo_series <- function(gse_id, data_dir) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Install GEOquery before reading GEO series matrices: BiocManager::install('GEOquery')")
  }
  
  series_file <- find_series_matrix(gse_id, data_dir)
  if (!is.na(series_file) && file.exists(series_file)) {
    gse <- GEOquery::getGEO(filename = series_file, GSEMatrix = TRUE, getGPL = FALSE)
  } else {
    if (!RUN_DOWNLOADS) {
      stop(
        "[", gse_id, "] Series matrix file is missing under:\n  ", data_dir,
        "\nSet RUN_DOWNLOADS <- TRUE and rerun, or download the series matrix manually."
      )
    }
    message("[", gse_id, "] Series matrix missing; downloading via GEOquery::getGEO()...")
    gse <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
  }
  
  if (is.list(gse)) gse[[1]] else gse
}

if (RUN_DOWNLOADS && DOWNLOAD_ALL_AT_START) {
  download_geo_supp("GSE157827")
  download_geo_supp("GSE167490")
  download_geo_supp("GSE167492")
  download_geo_supp(
    "GSE214979",
    filter_regex = "cell_metadata\\.csv\\.gz|filtered_feature_bc_matrix\\.h5$",
    untar_archives = FALSE
  )
  download_geo_supp(
    "GSE243292",
    filter_regex = "ADsnRNAseq_GEO.h5ad.gz|processed_data_readme",
    untar_archives = FALSE
  )
  download_geo_supp("GSE263468", untar_archives = FALSE)
  download_geo_supp("GSE268599")
  download_geo_supp("GSE303823", untar_archives = FALSE)
}

if (RUN_PREP_LARGE_FILES) {
  suppressPackageStartupMessages(library(Seurat))
  
  gse303823_big_file <- file.path(data_root, "GSE303823", "GSE303823_rna_annotated.rds")
  gse303823_slim_file <- file.path(data_root, "GSE303823", "GSE303823_raw_counts.rds")
  
  if (!file.exists(gse303823_big_file)) {
    stop("Missing file: ", gse303823_big_file)
  }
  
  temp_obj <- readRDS(gse303823_big_file)
  raw_counts <- GetAssayData(temp_obj, layer = "counts")
  metadata <- temp_obj@meta.data
  stopifnot(identical(colnames(raw_counts), rownames(metadata)))
  saveRDS(list(counts = raw_counts, meta = metadata), file = gse303823_slim_file)
  
  rm(temp_obj, raw_counts, metadata)
  gc()
  
  message(
    "For GSE263468, convert the h5ad file before this R step, for example:\n",
    "python scripts/convert_h5ad_to_mtx.py --input ",
    file.path(data_root, "GSE263468", "GSE263468_processed_data.h5ad"),
    " --output-dir ", file.path(data_root, "GSE263468")
  )
}

# ==============================================================================
# 1. GSE157827 (Lau et al.)
# ==============================================================================

message("Processing GSE157827...")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(openxlsx)
  library(Seurat)
})

data_dir <- file.path(data_root, "GSE157827")

series_matrix_path <- file.path(data_dir, "GSE157827_series_matrix.txt")

# ------------------------------------------------------------------------------
# 1.1  Parse GEO series matrix
# ------------------------------------------------------------------------------

parse_gse157827_series_matrix <- function(path) {
  lines        <- readLines(path, warn = FALSE)
  sample_lines <- lines[grepl("^!Sample_", lines)]
  
  get_field <- function(prefix) {
    x <- sample_lines[grepl(paste0("^", prefix), sample_lines)]
    if (length(x) == 0L) return(NULL)
    x[1L]
  }
  split_line <- function(x) {
    if (is.null(x)) return(character(0))
    gsub('^"|"$', "", strsplit(x, "\\t")[[1L]][-1L])
  }
  
  char_lines <- sample_lines[grepl("^!Sample_characteristics_ch1", sample_lines)]
  char_map   <- list()
  for (ln in char_lines) {
    vals          <- gsub('^"|"$', "", strsplit(ln, "\\t")[[1L]][-1L])
    key           <- str_trim(strsplit(vals[1L], ":")[[1L]][1L])
    char_map[[key]] <- str_trim(sub(paste0("^", key, ":\\s*"), "", vals))
  }
  
  tibble(
    sample_name      = split_line(get_field("!Sample_title")),
    gsm_id           = split_line(get_field("!Sample_geo_accession")),
    source_name      = split_line(get_field("!Sample_source_name_ch1")),
    tissue           = char_map[["tissue"]],
    diagnosis_series = char_map[["diagnosis"]],
    instrument_model = split_line(get_field("!Sample_instrument_model")),
    library_strategy = split_line(get_field("!Sample_library_strategy"))
  )
}

lookup_series <- parse_gse157827_series_matrix(series_matrix_path)

# ------------------------------------------------------------------------------
# 1.2  Merge with donor metadata
# ------------------------------------------------------------------------------

braak_info <- read_csv(
  file.path(data_dir, "Patient_braak_info.csv"),
  show_col_types = FALSE
) %>%
  rename(
    sample_name          = ID,
    condition_csv        = CONDITION,
    age                  = AGE,
    sex                  = SEX,
    pmd_hr               = PMD,
    hist_diagnosis       = `HIST DIAGNOSIS`,
    diag_1               = `DIAG 1`,
    diag_2               = `DIAG 2`,
    diag_3               = `DIAG 3`,
    amyloid_plaque_score = `AGE-RELATED PLAQUE SCORE`,
    braak_numeric        = `Braak tangle stage`,
    apoe                 = APOE
  ) %>%
  mutate(age = suppressWarnings(as.numeric(age)))

lookup_table <- lookup_series %>%
  left_join(braak_info, by = "sample_name") %>%
  mutate(
    region    = "PFC",
    study     = "GSE157827",
    age_pass  = !is.na(age) & age >= 65,
    condition = assign_relaxed_condition(condition_csv),
    included_in_analysis = if_else(!is.na(condition) & age_pass, "YES", "NO"),
    exclusion_reason = case_when(
      included_in_analysis == "YES"        ~ NA_character_,
      !age_pass & !is.na(condition)        ~ paste0("Age < 65 (age=", age, ")"),
      TRUE                                 ~ paste0("Diagnosis not AD/control: diag=", condition_csv)
    )
  ) %>%
  arrange(
    factor(condition, levels = c("AD", "Control"), exclude = NULL),
    desc(braak_numeric), sample_name
  )

write_csv(lookup_table,  file.path(output_dir, "gse157827_lookup_table.csv"))
write.xlsx(lookup_table, file.path(output_dir, "gse157827_lookup_table.xlsx"),
           asTable = TRUE)

message(sprintf(
  "  [GSE157827] total rows=%d / donors=%d | AD rows=%d / donors=%d | Control rows=%d / donors=%d | excluded rows=%d / donors=%d",
  nrow(lookup_table),
  count_unique_ids(lookup_table, "sample_name"),
  count_condition_rows(lookup_table, "AD"),
  count_condition_donors(lookup_table, "sample_name", "AD"),
  count_condition_rows(lookup_table, "Control"),
  count_condition_donors(lookup_table, "sample_name", "Control"),
  sum(lookup_table$included_in_analysis == "NO"),
  count_excluded_donors(lookup_table, "sample_name")
))

# ------------------------------------------------------------------------------
# 1.3  Sample selection
# ------------------------------------------------------------------------------

meta_ref <- lookup_table %>%
  filter(included_in_analysis == "YES") %>%
  select(gsm_id, sample_name, condition, braak_numeric, region,
         age, sex, pmd_hr, apoe, hist_diagnosis, amyloid_plaque_score)

message(sprintf(
  "  [GSE157827] Selected donors: AD=%d / Control=%d (rows=%d)",
  count_condition_donors(meta_ref, "sample_name", "AD"),
  count_condition_donors(meta_ref, "sample_name", "Control"),
  nrow(meta_ref)
))

# ------------------------------------------------------------------------------
# 1.4  Data loading
# ------------------------------------------------------------------------------

seurat_list <- list()

for (i in seq_len(nrow(meta_ref))) {
  gsm_id <- meta_ref$gsm_id[i]
  sid    <- meta_ref$sample_name[i]
  prefix <- paste0(gsm_id, "_", sid)
  
  mtx_path  <- list.files(data_dir, pattern = paste0("^", prefix, ".*matrix.mtx.gz$"),
                          full.names = TRUE)[1]
  bar_path  <- list.files(data_dir, pattern = paste0("^", prefix, ".*barcodes.tsv.gz$"),
                          full.names = TRUE)[1]
  feat_path <- list.files(data_dir, pattern = paste0("^", prefix, ".*features.tsv.gz$"),
                          full.names = TRUE)[1]
  
  if (any(is.na(c(mtx_path, bar_path, feat_path)))) {
    mtx_path  <- list.files(data_dir, pattern = paste0("^", gsm_id, ".*matrix.mtx.gz$"),
                            full.names = TRUE)[1]
    bar_path  <- list.files(data_dir, pattern = paste0("^", gsm_id, ".*barcodes.tsv.gz$"),
                            full.names = TRUE)[1]
    feat_path <- list.files(data_dir, pattern = paste0("^", gsm_id, ".*features.tsv.gz$"),
                            full.names = TRUE)[1]
  }
  
  stopifnot(file.exists(mtx_path), file.exists(bar_path), file.exists(feat_path))
  
  counts     <- ReadMtx(mtx = mtx_path, cells = bar_path, features = feat_path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = gsm_id,
                                   min.cells = 10, min.features = 200)
  
  seurat_obj$study                <- "GSE157827"
  seurat_obj$sample_id            <- sid
  seurat_obj$donor_id             <- sid
  seurat_obj$gsm_id               <- gsm_id
  seurat_obj$condition            <- meta_ref$condition[i]
  seurat_obj$braak_numeric        <- meta_ref$braak_numeric[i]
  seurat_obj$region               <- meta_ref$region[i]
  seurat_obj$age                  <- meta_ref$age[i]
  seurat_obj$sex                  <- meta_ref$sex[i]
  seurat_obj$pmd_hr               <- meta_ref$pmd_hr[i]
  seurat_obj$apoe                 <- meta_ref$apoe[i]
  seurat_obj$hist_diagnosis       <- meta_ref$hist_diagnosis[i]
  seurat_obj$amyloid_plaque_score <- meta_ref$amyloid_plaque_score[i]
  
  seurat_list[[gsm_id]] <- seurat_obj
}

# ------------------------------------------------------------------------------
# 1.5  Merge & QC
# ------------------------------------------------------------------------------

merged_obj <- merge_seurat_list(seurat_list, project = "GSE157827")
merged_obj <- join_layers_if_available(merged_obj)
merged_obj[["percent.mt"]] <- PercentageFeatureSet(merged_obj, pattern = "^MT-")
merged_obj <- subset(merged_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 5000 &
                       percent.mt  < 10)

message(sprintf("  [GSE157827] Post-QC cells: %d", ncol(merged_obj)))
saveRDS(merged_obj, file.path(output_dir, "gse157827_azimuth_input.rds"))

rm(seurat_list, merged_obj); gc()



# ==============================================================================
# 2. GSE167490 (Sadick et al. main cohort)
# ==============================================================================
#
# Reference : Sadick et al., Neuron 2022 (PMID: 35381194)
# Design    : snRNA-seq, astrocytes & oligodendrocytes, PFC
# Criteria  : Relaxed criteria: AD/control diagnosis labels + age >= 65
# ==============================================================================

message("Processing GSE167490 (Sadick et al. main cohort)...")

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
  library(stringr)
  library(openxlsx)
  library(Seurat)
})

data_dir_167 <- file.path(data_root, "GSE167490")
ensure_geo_supp_files(
  "GSE167490",
  required_patterns = c("matrix\\.mtx\\.gz$", "barcodes\\.tsv\\.gz$", "features\\.tsv\\.gz$")
)

# ------------------------------------------------------------------------------
# 2.1  Parse GEO series matrix
# ------------------------------------------------------------------------------

gse <- load_geo_series("GSE167490", data_dir_167)

lookup_geo <- pData(gse) %>%
  select(-matches("extract_protocol|data_processing|supplementary|relation|description")) %>%
  rename(gsm_id    = geo_accession,
         geo_title = title,
         tissue    = source_name_ch1) %>%
  filter(grepl("^Donor", geo_title)) %>%
  mutate(
    donor_key = gsub("^Donor", "D", geo_title),
    age_series = suppressWarnings(as.numeric(first_existing_column(
      cur_data_all(),
      c("age:ch1", "Age:ch1", "age", "Age")
    ))),
    sorting = first_existing_column(
      cur_data_all(),
      c("sorting_condition:ch1", "sorting condition:ch1",
        "sorting_condition", "sorting condition"),
      default = "Not reported"
    ),
    disease_state = first_existing_column(
      cur_data_all(),
      c("disease state:ch1", "disease_state:ch1",
        "disease state", "disease_state")
    ),
    sex = first_existing_column(cur_data_all(), c("Sex:ch1", "sex:ch1", "Sex", "sex")),
    sorting = if_else(is.na(sorting) | sorting == "", "Not reported", sorting)
  )

# ------------------------------------------------------------------------------
# 2.2  Load donor metadata
# ------------------------------------------------------------------------------

donor_info_path <- c(
  file.path(data_dir_167, "gse167490_donor_info.csv"),
  file.path(data_dir_167, "donor_braak_stage_info.csv")
)
donor_info_path <- donor_info_path[file.exists(donor_info_path)][1]

if (!is.na(donor_info_path)) {
  donor_info <- read.csv(donor_info_path, stringsAsFactors = FALSE)
} else {
  table_s2_path <- c(
    file.path(data_dir_167, "NIHMS1795138-supplement-4.xlsx"),
    file.path(data_root, "NIHMS1795138-supplement-4.xlsx")
  )
  table_s2_path <- table_s2_path[file.exists(table_s2_path)][1]
  if (is.na(table_s2_path)) {
    stop(
      "[GSE167490] Missing donor metadata. Provide gse167490_donor_info.csv, ",
      "donor_braak_stage_info.csv, or NIHMS1795138-supplement-4.xlsx."
    )
  }
  donor_info <- openxlsx::read.xlsx(table_s2_path, sheet = "ADRC_Dx_donor_information")
}

braak_csv <- donor_info %>%
  filter(!grepl("\\(IF\\)", DONOR_NUMBER)) %>%
  select(donor_key      = DONOR_NUMBER,
         braak_numeric  = BRAAK,
         pmi            = PMI,
         apoe           = APOE,
         rin            = RIN,
         sample_bank_id = SAMPLE_ID) %>%
  mutate(braak_numeric = suppressWarnings(as.numeric(braak_numeric)))

# ------------------------------------------------------------------------------
# 2.3  Build lookup table
# ------------------------------------------------------------------------------

lookup_table <- lookup_geo %>%
  left_join(braak_csv, by = "donor_key") %>%
  mutate(
    age      = age_series,
    age_pass = !is.na(age) & age >= 65,
    removed_by_publication_qc = donor_key %in% c("D5", "D9"),
    publication_qc_note = if_else(
      removed_by_publication_qc,
      "Original publication removed this donor from final analyses after sequencing analysis; retained for independent reprocessing and flagged here.",
      NA_character_
    ),
    diag_raw = str_to_lower(str_trim(disease_state)),
    condition = assign_relaxed_condition(diag_raw),
    included_in_analysis = if_else(!is.na(condition) & age_pass, "YES", "NO"),
    exclusion_reason = case_when(
      included_in_analysis == "YES"  ~ NA_character_,
      !age_pass & !is.na(condition)  ~ paste0("Age < 65 (age=", age, ")"),
      TRUE                           ~ paste0("Diagnosis not AD/control: diag=", disease_state)
    ),
    region = "PFC",
    study  = "GSE167490"
  ) %>%
  select(gsm_id, geo_title, donor_key, sample_bank_id, tissue,
         disease_state,
         sex,
         age,
         sorting,
         braak_numeric, condition, pmi, apoe, rin,
         removed_by_publication_qc, publication_qc_note,
         region, study, included_in_analysis, exclusion_reason) %>%
  arrange(
    factor(condition, levels = c("AD", "Control"), exclude = NULL),
    desc(braak_numeric), geo_title
  )

write_csv(lookup_table,  file.path(output_dir, "gse167490_lookup_table.csv"))
write.xlsx(lookup_table, file.path(output_dir, "gse167490_lookup_table.xlsx"),
           asTable = TRUE)

message(sprintf(
  "  [GSE167490] total rows=%d / donors=%d | AD rows=%d / donors=%d | Control rows=%d / donors=%d | excluded rows=%d / donors=%d",
  nrow(lookup_table),
  count_unique_ids(lookup_table, "donor_key"),
  count_condition_rows(lookup_table, "AD"),
  count_condition_donors(lookup_table, "donor_key", "AD"),
  count_condition_rows(lookup_table, "Control"),
  count_condition_donors(lookup_table, "donor_key", "Control"),
  sum(lookup_table$included_in_analysis == "NO"),
  count_excluded_donors(lookup_table, "donor_key")
))

# ------------------------------------------------------------------------------
# 2.4  Sample selection
# ------------------------------------------------------------------------------

meta_ref <- lookup_table %>%
  filter(included_in_analysis == "YES") %>%
  select(gsm_id, donor_key, condition, braak_numeric,
         age, sex, pmi, apoe, rin, region)

if (nrow(meta_ref) == 0L) {
  stop(
    "[GSE167490] No samples passed relaxed criteria. Check lookup table at: ",
    file.path(output_dir, "gse167490_lookup_table.csv")
  )
}

message(sprintf(
  "  [GSE167490] Selected donors: AD=%d / Control=%d (rows=%d)",
  count_condition_donors(meta_ref, "donor_key", "AD"),
  count_condition_donors(meta_ref, "donor_key", "Control"),
  nrow(meta_ref)
))

# ------------------------------------------------------------------------------
# 2.5  Data loading
# ------------------------------------------------------------------------------

seurat_list <- list()

for (i in seq_len(nrow(meta_ref))) {
  gsm_id <- meta_ref$gsm_id[i]
  
  mtx_path  <- list.files(data_dir_167,
                          pattern    = paste0(gsm_id, ".*matrix.mtx.gz"),
                          full.names = TRUE, recursive = TRUE)[1]
  bar_path  <- list.files(data_dir_167,
                          pattern    = paste0(gsm_id, ".*barcodes.tsv.gz"),
                          full.names = TRUE, recursive = TRUE)[1]
  feat_path <- list.files(data_dir_167,
                          pattern    = paste0(gsm_id, ".*features.tsv.gz"),
                          full.names = TRUE, recursive = TRUE)[1]
  
  if (any(is.na(c(mtx_path, bar_path, feat_path)))) {
    warning(sprintf("  [GSE167490] Files missing for %s — skipping.", gsm_id))
    next
  }
  
  counts     <- ReadMtx(mtx = mtx_path, cells = bar_path, features = feat_path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = gsm_id,
                                   min.cells = 10, min.features = 200)
  
  seurat_obj$study         <- "GSE167490"
  seurat_obj$sample_id     <- gsm_id
  seurat_obj$gsm_id        <- gsm_id
  seurat_obj$donor_id      <- meta_ref$donor_key[i]
  seurat_obj$donor_key     <- meta_ref$donor_key[i]
  seurat_obj$condition     <- meta_ref$condition[i]
  seurat_obj$braak_numeric <- meta_ref$braak_numeric[i]
  seurat_obj$region        <- meta_ref$region[i]
  seurat_obj$age           <- meta_ref$age[i]
  seurat_obj$sex           <- meta_ref$sex[i]
  seurat_obj$pmi           <- meta_ref$pmi[i]
  seurat_obj$apoe          <- meta_ref$apoe[i]
  seurat_obj$rin           <- meta_ref$rin[i]
  
  seurat_list[[gsm_id]] <- seurat_obj
}

# ------------------------------------------------------------------------------
# 2.6  Merge & QC
# ------------------------------------------------------------------------------

merged_obj <- merge_seurat_list(seurat_list, project = "GSE167490")
merged_obj <- join_layers_if_available(merged_obj)
merged_obj[["percent.mt"]] <- PercentageFeatureSet(merged_obj, pattern = "^MT-")
merged_obj <- subset(merged_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 5000 &
                       percent.mt  < 10)

message(sprintf("  [GSE167490] Post-QC cells: %d", ncol(merged_obj)))
saveRDS(merged_obj, file.path(output_dir, "gse167490_azimuth_input.rds"))

rm(gse, lookup_geo, donor_info, braak_csv, lookup_table,
   meta_ref, seurat_list, merged_obj)
gc()



# ==============================================================================
# 3. GSE167492 (Sadick et al. pilot cohort)
# ==============================================================================
#
# Reference : Sadick et al., Neuron 2022 (PMID: 35381189)
# Design    : Pilot snRNA-seq cohort, SOX9-positive/negative sorted PFC nuclei
# Criteria  : Relaxed criteria: AD/control diagnosis labels + age >= 65
# ==============================================================================

message("Processing GSE167492 (Sadick et al. pilot cohort)...")

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
  library(stringr)
  library(openxlsx)
  library(Seurat)
})

data_dir_492 <- file.path(data_root, "GSE167492")
ensure_geo_supp_files(
  "GSE167492",
  required_patterns = c("matrix\\.mtx\\.gz$", "barcodes\\.tsv\\.gz$", "features\\.tsv\\.gz$")
)

# ------------------------------------------------------------------------------
# 3.1 Parse GEO series matrix
# ------------------------------------------------------------------------------

gse <- load_geo_series("GSE167492", data_dir_492)

lookup_geo <- pData(gse) %>%
  select(-matches("extract_protocol|data_processing|supplementary|relation|description")) %>%
  rename(
    gsm_id = geo_accession,
    geo_title = title,
    tissue = source_name_ch1
  ) %>%
  mutate(
    donor_key = case_when(
      str_detect(geo_title, "^NS1") ~ "Pilot_NS1",
      str_detect(geo_title, "^NS2") ~ "Pilot_NS2",
      str_detect(geo_title, "^NS3") ~ "Pilot_NS3",
      str_detect(geo_title, "^AD1") ~ "Pilot_AD1",
      str_detect(geo_title, "^AD2") ~ "Pilot_AD2",
      TRUE ~ NA_character_
    ),
    sorting = first_existing_column(
      cur_data_all(),
      c("sorting_condition:ch1", "sorting condition:ch1",
        "sorting_condition", "sorting condition"),
      default = "Not reported"
    ),
    diag_raw = first_existing_column(
      cur_data_all(),
      c("disease state:ch1", "disease_state:ch1",
        "disease state", "disease_state")
    ),
    age_series = suppressWarnings(as.numeric(first_existing_column(
      cur_data_all(),
      c("age:ch1", "Age:ch1", "age", "Age")
    ))),
    sex_series = first_existing_column(cur_data_all(), c("Sex:ch1", "sex:ch1", "Sex", "sex"))
  )

# ------------------------------------------------------------------------------
# 3.2 Load pilot donor metadata
# ------------------------------------------------------------------------------

pilot_info_path <- file.path(data_dir_492, "gse167492_donor_info.csv")
if (file.exists(pilot_info_path)) {
  pilot_info <- read_csv(pilot_info_path, show_col_types = FALSE)
} else {
  warning("[GSE167492] gse167492_donor_info.csv not found; using embedded Table S2 pilot metadata.")
  pilot_info <- tibble::tribble(
    ~DONOR_NUMBER, ~SAMPLE_ID, ~RIN, ~PMI, ~AGE, ~SEX, ~APOE, ~BRAAK,
    "Pilot_NS1", "BU_NS3", 7, "<24", 74, "F", 33, NA,
    "Pilot_NS2", "BU_NS4", NA, "<24", 60, "M", 33, NA,
    "Pilot_NS3", "NYUTN11-48", 8.3, "12", 90, "M", 23, NA,
    "Pilot_AD1", "BU_AD1", 9.3, "<24", 80, "M", 34, 6,
    "Pilot_AD2", "BU_AD2", 9.6, "<24", 90, "M", 34, 6
  )
}

pilot_info <- pilot_info %>%
  transmute(
    donor_key = DONOR_NUMBER,
    sample_bank_id = SAMPLE_ID,
    rin = suppressWarnings(as.numeric(RIN)),
    pmi = as.character(PMI),
    age_donor = suppressWarnings(as.numeric(AGE)),
    sex_donor = SEX,
    apoe = as.character(APOE),
    braak_numeric = suppressWarnings(as.numeric(BRAAK))
  )

# ------------------------------------------------------------------------------
# 3.3 Build lookup table
# ------------------------------------------------------------------------------

lookup_table <- lookup_geo %>%
  left_join(pilot_info, by = "donor_key") %>%
  mutate(
    age = coalesce(age_donor, age_series),
    sex = coalesce(sex_donor, sex_series),
    age_pass = !is.na(age) & age >= 65,
    condition = assign_relaxed_condition(diag_raw),
    included_in_analysis = if_else(!is.na(condition) & age_pass, "YES", "NO"),
    exclusion_reason = case_when(
      included_in_analysis == "YES" ~ NA_character_,
      !age_pass & !is.na(condition) ~ paste0("Age < 65 (age=", age, ")"),
      TRUE ~ paste0("Diagnosis not AD/control: diag=", diag_raw)
    ),
    region = "PFC",
    study = "GSE167492"
  ) %>%
  select(
    gsm_id, geo_title, donor_key, sample_bank_id, tissue, diag_raw,
    sorting, age, sex, condition, braak_numeric, pmi, apoe, rin,
    region, study, included_in_analysis, exclusion_reason
  ) %>%
  arrange(
    factor(condition, levels = c("AD", "Control"), exclude = NULL),
    donor_key, sorting, gsm_id
  )

write_csv(lookup_table, file.path(output_dir, "gse167492_lookup_table.csv"))
write.xlsx(lookup_table, file.path(output_dir, "gse167492_lookup_table.xlsx"),
           asTable = TRUE)

message(sprintf(
  "  [GSE167492] total libraries=%d / donors=%d | AD libraries=%d / donors=%d | Control libraries=%d / donors=%d | excluded libraries=%d / donors=%d",
  nrow(lookup_table),
  count_unique_ids(lookup_table, "donor_key"),
  sum(lookup_table$condition == "AD" & lookup_table$included_in_analysis == "YES", na.rm = TRUE),
  count_condition_donors(lookup_table[lookup_table$included_in_analysis == "YES", , drop = FALSE], "donor_key", "AD"),
  sum(lookup_table$condition == "Control" & lookup_table$included_in_analysis == "YES", na.rm = TRUE),
  count_condition_donors(lookup_table[lookup_table$included_in_analysis == "YES", , drop = FALSE], "donor_key", "Control"),
  sum(lookup_table$included_in_analysis == "NO"),
  count_excluded_donors(lookup_table, "donor_key")
))

# ------------------------------------------------------------------------------
# 3.4 Sample selection
# ------------------------------------------------------------------------------

meta_ref <- lookup_table %>%
  filter(included_in_analysis == "YES") %>%
  select(gsm_id, geo_title, donor_key, sample_bank_id, condition,
         braak_numeric, age, sex, pmi, apoe, rin, sorting, region)

message(sprintf(
  "  [GSE167492] Selected donors: AD=%d / Control=%d (libraries=%d)",
  count_condition_donors(meta_ref, "donor_key", "AD"),
  count_condition_donors(meta_ref, "donor_key", "Control"),
  nrow(meta_ref)
))

# ------------------------------------------------------------------------------
# 3.5 Data loading
# ------------------------------------------------------------------------------

seurat_list <- list()

for (i in seq_len(nrow(meta_ref))) {
  gsm_id <- meta_ref$gsm_id[i]
  
  mtx_path <- list.files(data_dir_492,
                         pattern = paste0(gsm_id, ".*matrix.mtx.gz"),
                         full.names = TRUE, recursive = TRUE)[1]
  bar_path <- list.files(data_dir_492,
                         pattern = paste0(gsm_id, ".*barcodes.tsv.gz"),
                         full.names = TRUE, recursive = TRUE)[1]
  feat_path <- list.files(data_dir_492,
                          pattern = paste0(gsm_id, ".*features.tsv.gz"),
                          full.names = TRUE, recursive = TRUE)[1]
  
  if (any(is.na(c(mtx_path, bar_path, feat_path)))) {
    warning(sprintf("  [GSE167492] Files missing for %s — skipping.", gsm_id))
    next
  }
  
  counts <- ReadMtx(mtx = mtx_path, cells = bar_path, features = feat_path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = gsm_id,
                                   min.cells = 10, min.features = 200)
  
  seurat_obj$study <- "GSE167492"
  seurat_obj$sample_id <- gsm_id
  seurat_obj$gsm_id <- gsm_id
  seurat_obj$library_id <- meta_ref$geo_title[i]
  seurat_obj$donor_id <- meta_ref$donor_key[i]
  seurat_obj$donor_key <- meta_ref$donor_key[i]
  seurat_obj$sample_bank_id <- meta_ref$sample_bank_id[i]
  seurat_obj$condition <- meta_ref$condition[i]
  seurat_obj$braak_numeric <- meta_ref$braak_numeric[i]
  seurat_obj$region <- meta_ref$region[i]
  seurat_obj$age <- meta_ref$age[i]
  seurat_obj$sex <- meta_ref$sex[i]
  seurat_obj$pmi <- meta_ref$pmi[i]
  seurat_obj$apoe <- meta_ref$apoe[i]
  seurat_obj$rin <- meta_ref$rin[i]
  seurat_obj$sorting <- meta_ref$sorting[i]
  
  seurat_list[[gsm_id]] <- seurat_obj
}

# ------------------------------------------------------------------------------
# 3.6 Merge & QC
# ------------------------------------------------------------------------------

merged_obj <- merge_seurat_list(seurat_list, project = "GSE167492")
merged_obj <- join_layers_if_available(merged_obj)
merged_obj[["percent.mt"]] <- PercentageFeatureSet(merged_obj, pattern = "^MT-")
merged_obj <- subset(merged_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 5000 &
                       percent.mt < 10)

message(sprintf("  [GSE167492] Post-QC cells: %d", ncol(merged_obj)))
saveRDS(merged_obj, file.path(output_dir, "gse167492_azimuth_input.rds"))

rm(gse, lookup_geo, pilot_info, lookup_table, meta_ref, seurat_list, merged_obj)
gc()



# ==============================================================================
# 4. GSE214979 (Anderson et al.)
# ==============================================================================
#
# Reference : Anderson et al., Nat Neurosci 2023 (PMID: 36950385)
# Design    : single-nucleus multiome, DLPFC, AD and unaffected controls
# File      : filtered_feature_bc_matrix.h5 + cell_metadata.csv.gz
# Criteria  : Relaxed criteria: AD/unaffected diagnosis labels + age >= 65.
#             Braak stage is retained as metadata but is not used for inclusion.
# ==============================================================================

message("Processing GSE214979 (Anderson et al.)...")

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
  library(stringr)
  library(openxlsx)
  library(Seurat)
})

data_dir_214 <- file.path(data_root, "GSE214979")
ensure_geo_supp_files(
  "GSE214979",
  required_patterns = c(
    "GSE214979_filtered_feature_bc_matrix\\.h5$",
    "GSE214979_cell_metadata\\.csv\\.gz$"
  ),
  filter_regex = "cell_metadata\\.csv\\.gz|filtered_feature_bc_matrix\\.h5$",
  untar_archives = FALSE
)

donor_info_path <- file.path(data_dir_214, "GSE214979_braak_stage_info.csv")
if (!file.exists(donor_info_path)) {
  stop(
    "[GSE214979] Missing donor metadata file: ", donor_info_path, "\n",
    "Place GSE214979_braak_stage_info.csv under the GSE214979 raw-data directory."
  )
}

# ------------------------------------------------------------------------------
# 4.1 Parse GEO series matrix and keep RNA libraries only
# ------------------------------------------------------------------------------

gse <- load_geo_series("GSE214979", data_dir_214)
geo_pd <- pData(gse)

lookup_geo <- tibble::tibble(
  sample_id = as.character(geo_pd$geo_accession),
  geo_title = as.character(geo_pd$title),
  tissue = first_existing_column(
    geo_pd,
    c("source_name_ch1", "tissue:ch1", "source name:ch1"),
    default = "Not reported"
  ),
  molecule = first_existing_column(
    geo_pd,
    c("molecule_ch1", "molecule:ch1", "molecule"),
    default = "Not reported"
  )
)

lookup_rna <- lookup_geo %>%
  filter(
    str_to_lower(molecule) == "polya rna" |
      (str_detect(str_to_lower(geo_title), "scrnaseq") &
         !str_detect(str_to_lower(geo_title), "scatac"))
  ) %>%
  mutate(
    donor_id = geo_title %>%
      str_replace(regex(",?\\s*scRNAseq.*$", ignore_case = TRUE), "") %>%
      str_replace(regex("\\s*,?\\s*rep\\s*[0-9]+$", ignore_case = TRUE), "") %>%
      str_squish()
  )

if (nrow(lookup_rna) == 0L) {
  stop("[GSE214979] No scRNA-seq/polyA RNA samples were found in the series matrix.")
}

# ------------------------------------------------------------------------------
# 4.2 Load donor metadata
# ------------------------------------------------------------------------------

donor_raw <- read.csv(donor_info_path, skip = 1, stringsAsFactors = FALSE,
                      check.names = FALSE)

donor_meta <- tibble::tibble(
  donor_id = first_existing_column(donor_raw, c("ID", "Donor ID", "donor_id")),
  age = first_existing_column(donor_raw, c("Age", "age")),
  ethnicity = first_existing_column(donor_raw, c("Ethnicity", "ethnicity", "Race", "race")),
  sex = first_existing_column(donor_raw, c("Gender", "Sex", "sex")),
  apoe = first_existing_column(donor_raw, c("APOE_Status", "APOE Status", "APOE", "apoe")),
  diagnosis_raw = first_existing_column(donor_raw, c("Diagnosis", "diagnosis")),
  braak_numeric = first_existing_column(donor_raw, c("Braak", "braak", "Braak stage")),
  pmi = first_existing_column(donor_raw, c("PMI.h.", "PMI.h", "PMI (h)", "PMI", "pmi")),
  rin = first_existing_column(donor_raw, c("RIN", "rin")),
  structure = first_existing_column(donor_raw, c("structure", "Structure", "Region", "region"))
) %>%
  mutate(
    donor_id = str_squish(as.character(donor_id)),
    age = suppressWarnings(as.numeric(age)),
    braak_numeric = suppressWarnings(as.numeric(braak_numeric)),
    pmi = suppressWarnings(as.numeric(pmi)),
    rin = suppressWarnings(as.numeric(rin))
  ) %>%
  filter(!is.na(donor_id), donor_id != "")

if (nrow(donor_meta) == 0L) {
  stop(
    "[GSE214979] Donor metadata could not be parsed from: ", donor_info_path,
    "\nCheck that the file has columns such as ID, Age, Diagnosis, Braak, APOE_Status."
  )
}

# ------------------------------------------------------------------------------
# 4.3 Build lookup table
# ------------------------------------------------------------------------------

missing_donor_meta <- setdiff(unique(lookup_rna$donor_id), unique(donor_meta$donor_id))
if (length(missing_donor_meta) > 0L) {
  warning(
    "[GSE214979] GEO RNA donor IDs missing from donor metadata: ",
    paste(missing_donor_meta, collapse = ", ")
  )
}

lookup_table <- lookup_rna %>%
  left_join(donor_meta, by = "donor_id") %>%
  mutate(
    tissue_norm = str_to_lower(str_squish(paste(tissue, structure))),
    region_pass = str_detect(tissue_norm, "dlpfc|ba9|ba46|prefrontal|pfc"),
    age_pass = !is.na(age) & age >= 65,
    condition = assign_relaxed_condition(diagnosis_raw),
    included_in_analysis = if_else(!is.na(condition) & age_pass & region_pass,
                                   "YES", "NO"),
    exclusion_reason = case_when(
      included_in_analysis == "YES" ~ NA_character_,
      is.na(diagnosis_raw) | diagnosis_raw == "" ~ "Missing diagnosis metadata",
      is.na(condition) ~ paste0("Diagnosis not AD/control: diag=", diagnosis_raw),
      !age_pass ~ paste0("Age < 65 or missing (age=", age, ")"),
      !region_pass ~ paste0("Not a PFC/DLPFC sample: tissue=", tissue,
                            ", structure=", structure),
      TRUE ~ "Not selected"
    ),
    region = "PFC",
    study = "GSE214979"
  ) %>%
  select(
    sample_id, geo_title, donor_id, tissue, structure,
    diagnosis_raw, sex, age, ethnicity, apoe, braak_numeric, pmi, rin,
    condition, region, study, included_in_analysis, exclusion_reason
  ) %>%
  arrange(
    factor(condition, levels = c("AD", "Control"), exclude = NULL),
    donor_id, sample_id
  )

write_csv(lookup_table, file.path(output_dir, "gse214979_lookup_table.csv"))
write.xlsx(lookup_table, file.path(output_dir, "gse214979_lookup_table.xlsx"),
           asTable = TRUE)

message(sprintf(
  "  [GSE214979] total RNA libraries=%d / donors=%d | AD libraries=%d / donors=%d | Control libraries=%d / donors=%d | excluded libraries=%d / donors=%d",
  nrow(lookup_table),
  count_unique_ids(lookup_table, "donor_id"),
  count_condition_rows(lookup_table, "AD"),
  count_condition_donors(lookup_table, "donor_id", "AD"),
  count_condition_rows(lookup_table, "Control"),
  count_condition_donors(lookup_table, "donor_id", "Control"),
  sum(lookup_table$included_in_analysis == "NO"),
  count_excluded_donors(lookup_table, "donor_id")
))

# ------------------------------------------------------------------------------
# 4.4 Sample selection
# ------------------------------------------------------------------------------

meta_ref <- lookup_table %>%
  filter(included_in_analysis == "YES") %>%
  select(sample_id, donor_id, condition, braak_numeric, region,
         age, sex, ethnicity, apoe, pmi, rin, structure, diagnosis_raw)

if (nrow(meta_ref) == 0L) {
  stop(
    "[GSE214979] No samples passed relaxed criteria. Check lookup table at: ",
    file.path(output_dir, "gse214979_lookup_table.csv")
  )
}

message(sprintf(
  "  [GSE214979] Selected donors: AD=%d / Control=%d (RNA libraries=%d)",
  count_condition_donors(meta_ref, "donor_id", "AD"),
  count_condition_donors(meta_ref, "donor_id", "Control"),
  nrow(meta_ref)
))

donor_meta_lut <- meta_ref %>%
  select(donor_id, condition, braak_numeric, region, age, sex, ethnicity,
         apoe, pmi, rin, structure, diagnosis_raw) %>%
  distinct()

duplicated_donors <- unique(donor_meta_lut$donor_id[duplicated(donor_meta_lut$donor_id)])
if (length(duplicated_donors) > 0L) {
  stop(
    "[GSE214979] Conflicting donor-level metadata for: ",
    paste(duplicated_donors, collapse = ", ")
  )
}
rownames(donor_meta_lut) <- donor_meta_lut$donor_id

make_donor_lookup_214 <- function(val_col) {
  v <- donor_meta_lut[[val_col]]
  names(v) <- donor_meta_lut$donor_id
  v
}

lut_condition <- make_donor_lookup_214("condition")
lut_braak <- make_donor_lookup_214("braak_numeric")
lut_region <- make_donor_lookup_214("region")
lut_age <- make_donor_lookup_214("age")
lut_sex <- make_donor_lookup_214("sex")
lut_ethnicity <- make_donor_lookup_214("ethnicity")
lut_apoe <- make_donor_lookup_214("apoe")
lut_pmi <- make_donor_lookup_214("pmi")
lut_rin <- make_donor_lookup_214("rin")
lut_dx <- make_donor_lookup_214("diagnosis_raw")

# ------------------------------------------------------------------------------
# 4.5 Load combined 10X h5 counts and cell metadata
# ------------------------------------------------------------------------------

h5_file <- file.path(data_dir_214, "GSE214979_filtered_feature_bc_matrix.h5")
meta_file <- file.path(data_dir_214, "GSE214979_cell_metadata.csv.gz")

counts <- Read10X_h5(h5_file)
if (is.list(counts)) {
  assay_name <- intersect(names(counts), "Gene Expression")[1]
  if (is.na(assay_name)) assay_name <- names(counts)[1]
  counts <- counts[[assay_name]]
}

cell_meta <- read.csv(meta_file, row.names = 1, check.names = FALSE)
common_cells <- intersect(colnames(counts), rownames(cell_meta))
if (length(common_cells) == 0L) {
  stop("[GSE214979] No overlapping cell barcodes between h5 counts and cell metadata.")
}

seurat_obj <- CreateSeuratObject(
  counts = counts[, common_cells, drop = FALSE],
  meta.data = cell_meta[common_cells, , drop = FALSE],
  project = "GSE214979",
  min.cells = 10,
  min.features = 200
)
rm(counts); gc()

donor_col <- intersect(
  c("donor_id", "id", "donor", "sample", "sample_id", "orig.ident"),
  colnames(seurat_obj@meta.data)
)[1]
if (is.na(donor_col)) {
  stop(
    "[GSE214979] Could not find donor ID column in cell metadata. Available columns: ",
    paste(colnames(seurat_obj@meta.data), collapse = ", ")
  )
}
seurat_obj$donor_id <- as.character(seurat_obj@meta.data[[donor_col]])

selected_cells <- colnames(seurat_obj)[seurat_obj$donor_id %in% rownames(donor_meta_lut)]
if (length(selected_cells) == 0L) {
  stop("[GSE214979] No cells matched selected donor IDs.")
}
seurat_obj <- subset(seurat_obj, cells = selected_cells)

sid <- seurat_obj$donor_id
seurat_obj$study <- "GSE214979"
seurat_obj$sample_id <- sid
seurat_obj$condition <- unname(lut_condition[sid])
seurat_obj$braak_numeric <- unname(lut_braak[sid])
seurat_obj$region <- unname(lut_region[sid])
seurat_obj$age <- unname(lut_age[sid])
seurat_obj$sex <- unname(lut_sex[sid])
seurat_obj$ethnicity <- unname(lut_ethnicity[sid])
seurat_obj$apoe <- unname(lut_apoe[sid])
seurat_obj$pmi <- unname(lut_pmi[sid])
seurat_obj$rin <- unname(lut_rin[sid])
seurat_obj$diagnosis_raw <- unname(lut_dx[sid])

if (any(is.na(seurat_obj$condition))) {
  stop("[GSE214979] Condition metadata injection failed for one or more selected cells.")
}

message(sprintf("  [GSE214979] Post-donor-filter cells: %d", ncol(seurat_obj)))

# ------------------------------------------------------------------------------
# 4.6 QC filtering
# ------------------------------------------------------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 5000 &
                       percent.mt < 10)

message(sprintf("  [GSE214979] Post-QC cells: %d", ncol(seurat_obj)))
saveRDS(seurat_obj, file.path(output_dir, "gse214979_azimuth_input.rds"))

rm(gse, geo_pd, lookup_geo, lookup_rna, donor_raw, donor_meta, lookup_table,
   meta_ref, donor_meta_lut, cell_meta, common_cells, selected_cells, sid,
   seurat_obj)
gc()




# ==============================================================================
# 5. GSE263468 (Cobos et al.)
# ==============================================================================
#
# Reference : Cobos et al. (PMID: TBD)
# Design    : snRNA-seq, PFC, two platforms (GPL18573 + GPL24676)
#             Counts provided as single integrated h5ad → counts.mtx + metadata.csv
# Criteria  : Relaxed criteria: NP diagnosis AD/control + age >= 65
#             Excluded= PART, other diagnoses, age < 65
# Selection : BA9 (Prefrontal Cortex), sort == "All Nuclei"
# ==============================================================================

message("Processing GSE263468 (Cobos et al.)...")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(openxlsx)
  library(Seurat)
  library(Matrix)
})

data_dir_263 <- file.path(data_root, "GSE263468")

gse263468_required_files <- file.path(
  data_dir_263,
  c("counts.mtx", "genes.csv", "metadata.csv")
)
if (!all(file.exists(gse263468_required_files))) {
  stop(
    "[GSE263468] Missing h5ad-derived files: ",
    paste(basename(gse263468_required_files[!file.exists(gse263468_required_files)]),
          collapse = ", "),
    "\nGenerate them before running this section, for example:\n",
    "python scripts/convert_h5ad_to_mtx.py --input ",
    file.path(data_dir_263, "GSE263468_processed_data.h5ad"),
    " --output-dir ", data_dir_263
  )
}

# Braak roman numeral → integer
braak_map <- c("0"=0L, "I"=1L, "II"=2L, "III"=3L, "IV"=4L, "V"=5L, "VI"=6L)

# ------------------------------------------------------------------------------
# 5.1  Parse GSM accession maps from both platform series matrices
# ------------------------------------------------------------------------------

parse_gsm_runid_map <- function(path, platform_label) {
  lines        <- readLines(path, warn = FALSE)
  sample_lines <- lines[grepl("^!Sample_", lines)]
  
  get_field <- function(prefix) {
    x <- sample_lines[grepl(paste0("^", prefix), sample_lines)]
    if (length(x) == 0L) return(NULL)
    x[1L]
  }
  split_line <- function(x) {
    if (is.null(x)) return(character(0))
    gsub('^"|"$', "", strsplit(x, "\\t")[[1L]][-1L])
  }
  
  tibble(
    GSM_accession = split_line(get_field("!Sample_geo_accession")),
    geo_title     = split_line(get_field("!Sample_title")),
    platform      = platform_label
  ) %>%
    mutate(run_id = str_extract(geo_title, "(?<=ID:)[CD]\\d+"))
}

gsm_map <- bind_rows(
  parse_gsm_runid_map(
    file.path(data_dir_263, "GSE263468-GPL18573_series_matrix.txt"), "GPL18573"
  ),
  parse_gsm_runid_map(
    file.path(data_dir_263, "GSE263468-GPL24676_series_matrix.txt"), "GPL24676"
  )
)

# ------------------------------------------------------------------------------
# 5.2  Load supplementary metadata
# ------------------------------------------------------------------------------

supp_meta <- read_csv(
  file.path(data_dir_263, "GSE263468_braak_stage_info.csv"),
  show_col_types = FALSE
) %>%
  rename(
    run_id        = `Unique identifier`,
    assay         = Assay,
    donor_id      = `Donor ID`,
    apoe          = APOE,
    brain_region  = `Brain region`,
    sort          = `FANS SORT`,
    np_diagnosis  = `NP Diagnosis`,
    amyloid       = `Amyloid in sample`,
    nfts          = `NFTs in sample`,
    braak_stage   = `Braak stage`,
    cerad         = CERAD,
    disease_group = `Disease group`,
    brain_weight  = `Brain weight`,
    pmi_hr        = `PMI (hrs)`,
    race          = Race,
    age           = Age,
    sex           = Sex,
    rin           = RIN
  ) %>%
  mutate(
    np_diagnosis = str_trim(np_diagnosis),
    age          = suppressWarnings(as.numeric(age)),
    pmi_hr       = suppressWarnings(as.numeric(pmi_hr))
  )

# ------------------------------------------------------------------------------
# 5.3  Build lookup table
# ------------------------------------------------------------------------------

lookup_table <- supp_meta %>%
  left_join(gsm_map %>% select(run_id, GSM_accession, platform), by = "run_id") %>%
  mutate(
    brain_region_clean = str_squish(as.character(brain_region)),
    sort_clean = str_squish(as.character(sort)),
    braak_numeric = unname(braak_map[braak_stage]),
    age_pass  = !is.na(age) & age >= 65,
    condition = assign_relaxed_condition(np_diagnosis),
    included_in_analysis = if_else(!is.na(condition) & age_pass, "YES", "NO"),
    is_pfc_all_nuclei = brain_region_clean == "Prefrontal Cortex (BA9)" &
      str_to_lower(sort_clean) == "all nuclei",
    selected_for_pfc_meta_analysis = if_else(
      included_in_analysis == "YES" & is_pfc_all_nuclei,
      "YES", "NO"
    ),
    exclusion_reason = case_when(
      included_in_analysis == "YES"   ~ NA_character_,
      !age_pass & !is.na(condition)   ~ paste0("Age < 65 (age=", age, ")"),
      np_diagnosis == "PART"          ~ paste0("PART; Braak ", braak_stage),
      is.na(braak_numeric)            ~ paste0("Invalid Braak stage ('", braak_stage, "')"),
      TRUE                            ~ paste0("Diagnosis not AD/control: diag=", np_diagnosis)
    ),
    region = "PFC",
    study  = "GSE263468"
  ) %>%
  arrange(
    factor(condition, levels = c("AD", "Control"), exclude = NULL),
    desc(braak_numeric), brain_region, run_id
  )

write_csv(lookup_table,  file.path(output_dir, "gse263468_lookup_table.csv"))
write.xlsx(lookup_table, file.path(output_dir, "gse263468_lookup_table.xlsx"),
           asTable = TRUE)

pfc_all_nuclei_lookup <- lookup_table %>%
  filter(selected_for_pfc_meta_analysis == "YES")

message(sprintf(
  "  [GSE263468] lookup total across all regions/sorts: rows=%d / donors=%d | AD rows=%d / donors=%d | Control rows=%d / donors=%d | excluded rows=%d / donors=%d (PART rows=%d / donors=%d)",
  nrow(lookup_table),
  count_unique_ids(lookup_table, "donor_id"),
  count_condition_rows(lookup_table, "AD"),
  count_condition_donors(lookup_table, "donor_id", "AD"),
  count_condition_rows(lookup_table, "Control"),
  count_condition_donors(lookup_table, "donor_id", "Control"),
  sum(lookup_table$included_in_analysis == "NO"),
  count_excluded_donors(lookup_table, "donor_id"),
  sum(lookup_table$np_diagnosis %in% "PART"),
  count_unique_ids(lookup_table[lookup_table$np_diagnosis %in% "PART", , drop = FALSE], "donor_id")
))

message(sprintf(
  "  [GSE263468] PFC All Nuclei candidates only: rows=%d / donors=%d | AD donors=%d | Control donors=%d",
  nrow(pfc_all_nuclei_lookup),
  count_unique_ids(pfc_all_nuclei_lookup, "donor_id"),
  count_condition_donors(pfc_all_nuclei_lookup, "donor_id", "AD"),
  count_condition_donors(pfc_all_nuclei_lookup, "donor_id", "Control")
))

# ------------------------------------------------------------------------------
# 5.4  Sample selection: BA9 / All Nuclei / included_in_analysis == "YES"
# ------------------------------------------------------------------------------

meta_ref <- lookup_table %>%
  filter(selected_for_pfc_meta_analysis == "YES") %>%
  select(run_id, GSM_accession, donor_id, platform, assay,
         condition, braak_numeric, braak_stage,
         np_diagnosis, disease_group, brain_region, sort,
         age, sex, apoe, pmi_hr, amyloid, nfts, cerad, rin,
         region, study)

if (nrow(meta_ref) == 0L) {
  stop("[GSE263468] No BA9 / All Nuclei samples passed the relaxed criteria.")
}
if (any(str_to_lower(str_squish(meta_ref$sort)) != "all nuclei", na.rm = TRUE)) {
  stop("[GSE263468] Non-All Nuclei rows entered meta_ref; check sort filtering.")
}

message(sprintf(
  "  [GSE263468] Selected donors (BA9 / All Nuclei): AD=%d / Control=%d (rows=%d)",
  count_condition_donors(meta_ref, "donor_id", "AD"),
  count_condition_donors(meta_ref, "donor_id", "Control"),
  nrow(meta_ref)
))

# donor-level lookup for metadata injection (indexed by run_id)
duplicated_run_ids <- unique(meta_ref$run_id[duplicated(meta_ref$run_id)])
if (length(duplicated_run_ids) > 0L) {
  stop(
    "[GSE263468] Duplicated run_id values in selected metadata: ",
    paste(duplicated_run_ids, collapse = ", ")
  )
}

donor_lut <- meta_ref %>%
  select(run_id, donor_id, condition, braak_numeric,
         age, sex, apoe, pmi_hr)
rownames(donor_lut) <- donor_lut$run_id

# ------------------------------------------------------------------------------
# 5.5  Load integrated counts matrix
# ------------------------------------------------------------------------------

cell_meta  <- read.csv(file.path(data_dir_263, "metadata.csv"),
                       row.names = 1, check.names = FALSE)
gene_names <- read.csv(file.path(data_dir_263, "genes.csv"),
                       header = FALSE, stringsAsFactors = FALSE)[[1]]
counts_mat <- readMM(file.path(data_dir_263, "counts.mtx"))

n_genes <- length(gene_names)
n_cells <- nrow(cell_meta)

if (nrow(counts_mat) == n_genes && ncol(counts_mat) == n_cells) {
  rownames(counts_mat) <- gene_names
  colnames(counts_mat) <- rownames(cell_meta)
} else if (nrow(counts_mat) == n_cells && ncol(counts_mat) == n_genes) {
  counts_mat           <- t(counts_mat)
  rownames(counts_mat) <- gene_names
  colnames(counts_mat) <- rownames(cell_meta)
} else {
  stop(sprintf(
    "Matrix dimension mismatch: matrix %d×%d, genes %d, cells %d",
    nrow(counts_mat), ncol(counts_mat), n_genes, n_cells
  ))
}

message(sprintf("  [GSE263468] Matrix loaded: %d genes × %d cells",
                nrow(counts_mat), ncol(counts_mat)))

# Verify sample column (verified from adata.obs structure)
sample_id_col <- "sample"
stopifnot(
  "sample column not found in metadata" =
    sample_id_col %in% colnames(cell_meta)
)

# ------------------------------------------------------------------------------
# 5.6  Filter cells to target run_ids before creating Seurat object
# ------------------------------------------------------------------------------

target_ids  <- rownames(donor_lut)
matched_ids <- intersect(unique(cell_meta[[sample_id_col]]), target_ids)
unmatched   <- setdiff(target_ids, matched_ids)

message(sprintf("  [GSE263468] Sample match: %d / %d",
                length(matched_ids), length(target_ids)))
if (length(unmatched) > 0L)
  message("  [GSE263468] Unmatched run_ids: ",
          paste(head(unmatched, 10L), collapse = ", "))
if (length(matched_ids) == 0L) {
  stop("[GSE263468] No selected run_ids matched the metadata sample column.")
}

selected_cells <- rownames(cell_meta)[cell_meta[[sample_id_col]] %in% matched_ids]
missing_count_cells <- setdiff(selected_cells, colnames(counts_mat))
if (length(missing_count_cells) > 0L) {
  stop("[GSE263468] Selected metadata cells are missing from the count matrix.")
}

# Create the object after BA9 / All Nuclei selection so min.cells/min.features
# are applied consistently to the analysis subset.
seurat_obj <- CreateSeuratObject(
  counts = counts_mat[, selected_cells, drop = FALSE],
  meta.data = cell_meta[selected_cells, , drop = FALSE],
  project = "GSE263468",
  min.cells = 10,
  min.features = 200
)
rm(counts_mat); gc()

# ------------------------------------------------------------------------------
# 5.7  Inject sample-level metadata
# ------------------------------------------------------------------------------

rid <- seurat_obj@meta.data[[sample_id_col]]

make_lookup <- function(df, key_col, val_col) {
  v <- df[[val_col]]
  names(v) <- df[[key_col]]
  v
}

lut_condition     <- make_lookup(donor_lut, "run_id", "condition")
lut_donor_id      <- make_lookup(donor_lut, "run_id", "donor_id")
lut_braak         <- make_lookup(donor_lut, "run_id", "braak_numeric")
lut_age           <- make_lookup(donor_lut, "run_id", "age")
lut_sex           <- make_lookup(donor_lut, "run_id", "sex")
lut_apoe          <- make_lookup(donor_lut, "run_id", "apoe")
lut_pmi           <- make_lookup(donor_lut, "run_id", "pmi_hr")

seurat_obj$study         <- "GSE263468"
seurat_obj$sample_id     <- rid
seurat_obj$donor_id      <- as.character(lut_donor_id[rid])
seurat_obj$condition     <- unname(lut_condition[rid])
seurat_obj$braak_numeric <- unname(lut_braak[rid])
seurat_obj$age           <- unname(lut_age[rid])
seurat_obj$sex           <- unname(lut_sex[rid])
seurat_obj$apoe          <- unname(lut_apoe[rid])
seurat_obj$pmi_hr        <- unname(lut_pmi[rid])
seurat_obj$region        <- "PFC"

# ------------------------------------------------------------------------------
# 5.8  QC filtering
# ------------------------------------------------------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 5000 &
                       percent.mt  < 10)

message(sprintf("  [GSE263468] Post-QC cells: %d", ncol(seurat_obj)))
saveRDS(seurat_obj, file.path(output_dir, "gse263468_azimuth_input.rds"))

rm(cell_meta, gene_names, gsm_map, supp_meta, lookup_table, pfc_all_nuclei_lookup,
   meta_ref, donor_lut, matched_ids, target_ids, selected_cells, rid, seurat_obj)
gc()



# ==============================================================================
# 6. GSE268599 (Serrano-Pozo et al.)
# ==============================================================================
#
# Reference : Serrano-Pozo et al., Nat Neuroscience 2024 (PMID: 39528672)
# Design    : snRNA-seq, 5 brain regions, 32 donors, PathStage 1–4
# Region    : BA46 (dorsolateral PFC) only
# Criteria  : No clinical diagnosis field is available in the local metadata;
#             use study-provided PathStage grouping with age >= 65.
#             Control = PathStage 1, AD = PathStage 3-4, exclude PathStage 2.
# ==============================================================================

message("Processing GSE268599 (Serrano-Pozo et al.)...")

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
  library(stringr)
  library(openxlsx)
  library(Seurat)
})

data_dir_268 <- file.path(data_root, "GSE268599")

# ------------------------------------------------------------------------------
# 6.1  Parse GEO series matrix
# ------------------------------------------------------------------------------

gse <- getGEO(
  filename  = file.path(data_dir_268, "GSE268599_series_matrix.txt"),
  GSEMatrix = TRUE,
  getGPL    = FALSE
)

lookup_geo <- pData(gse) %>%
  select(
    gsm_id             = geo_accession,
    sample_name        = title,
    supplementary_file = supplementary_file_1,
    tissue             = `tissue:ch1`,
    age                = `age:ch1`,
    sex                = `Sex:ch1`
  ) %>%
  mutate(
    SampleName = str_extract(supplementary_file, "6289-MW-[0-9]+"),
    sex = case_when(
      str_to_lower(sex) == "male"   ~ "M",
      str_to_lower(sex) == "female" ~ "F",
      TRUE                          ~ sex
    )
  )

stopifnot(
  "MW number extraction failed for one or more rows" =
    !any(is.na(lookup_geo$SampleName))
)

# ------------------------------------------------------------------------------
# 6.2  Load donor metadata
# ------------------------------------------------------------------------------

# Braak roman numeral → integer
braak_map <- c("0"=0L, "I"=1L, "II"=2L, "III"=3L, "IV"=4L, "V"=5L, "VI"=6L)

sample_csv <- read.csv(
  file.path(data_dir_268, "DonorID_SampleName.csv"),
  stringsAsFactors = FALSE
) %>%
  transmute(
    SampleName           = SampleName,
    tissue_csv           = region,
    number_nuclei_sorted = number_nuclei_sorted,
    rin                  = suppressWarnings(as.numeric(RIN)),
    donor_id             = as.integer(Donor.ID)
  )

braak_csv <- read.csv(
  file.path(data_dir_268, "BraakStage_info.csv"),
  stringsAsFactors = FALSE
) %>%
  # NOTE: "Pathology Stage" (no dot) → R auto-converts to Pathology.Stage
  rename(
    donor_id    = Donor.ID,
    braak_stage = Braak,
    cerad       = CERAD,
    path_stage  = Pathology.Stage,
    apoe        = APOE
  ) %>%
  select(donor_id, braak_stage, cerad, path_stage, apoe) %>%
  mutate(
    donor_id      = as.integer(donor_id),
    braak_numeric = unname(braak_map[braak_stage]),
    # PathStage classification per Serrano-Pozo et al. 2024:
    #   PathStage 1 → no/low ADNC, Braak 0/I/II  → Control
    #   PathStage 2 → intermediate ADNC           → Excluded
    #   PathStage 3 → high ADNC, Braak V          → AD
    #   PathStage 4 → high ADNC, Braak VI         → AD
    condition_braak = case_when(
      path_stage == "PathStage 1"                      ~ "Control",
      path_stage %in% c("PathStage 3", "PathStage 4") ~ "AD",
      TRUE                                             ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------
# 6.3  Build lookup table
# ------------------------------------------------------------------------------

lookup_table <- lookup_geo %>%
  left_join(sample_csv, by = "SampleName") %>%
  left_join(braak_csv,  by = "donor_id") %>%
  mutate(
    # "90+" → 90 before numeric coercion
    age      = suppressWarnings(as.numeric(gsub("\\+$", "", age))),
    age_pass = !is.na(age) & age >= 65,
    condition = condition_braak,
    included_in_analysis = if_else(!is.na(condition) & age_pass, "YES", "NO"),
    exclusion_reason = case_when(
      included_in_analysis == "YES" ~ NA_character_,
      path_stage == "PathStage 2"   ~ "Intermediate pathology (PathStage 2)",
      !age_pass                     ~ paste0("Age < 65 (age=", age, ")"),
      TRUE                          ~ paste0(
        "Unknown: path_stage=", path_stage, ", braak=", braak_stage
      )
    ),
    region = "PFC",
    study  = "GSE268599"
  ) %>%
  select(
    gsm_id, sample_name, SampleName, tissue, donor_id,
    path_stage, braak_stage, braak_numeric, cerad, condition, apoe,
    age, sex, rin, number_nuclei_sorted,
    region, study, included_in_analysis, exclusion_reason
  ) %>%
  arrange(
    factor(condition, levels = c("AD", "Control"), exclude = NULL),
    desc(braak_numeric), sample_name
  )

n_join_fail <- sum(is.na(lookup_table$donor_id))
if (n_join_fail > 0L)
  warning(sprintf(
    "  [GSE268599] DonorID join failed for %d row(s) — check DonorID_SampleName.csv.",
    n_join_fail
  ))

write_csv(lookup_table,  file.path(output_dir, "gse268599_lookup_table.csv"))
write.xlsx(lookup_table, file.path(output_dir, "gse268599_lookup_table.xlsx"),
           asTable = TRUE)

message(sprintf(
  "  [GSE268599] total rows=%d / donors=%d | AD rows=%d / donors=%d | Control rows=%d / donors=%d | excluded rows=%d / donors=%d",
  nrow(lookup_table),
  count_unique_ids(lookup_table, "donor_id"),
  count_condition_rows(lookup_table, "AD"),
  count_condition_donors(lookup_table, "donor_id", "AD"),
  count_condition_rows(lookup_table, "Control"),
  count_condition_donors(lookup_table, "donor_id", "Control"),
  sum(lookup_table$included_in_analysis == "NO"),
  count_excluded_donors(lookup_table, "donor_id")
))

# ------------------------------------------------------------------------------
# 6.4  Sample selection: BA46, included_in_analysis == "YES"
# ------------------------------------------------------------------------------

meta_ref <- lookup_table %>%
  filter(tissue == "BA46", included_in_analysis == "YES") %>%
  select(gsm_id, sample_name, SampleName, donor_id,
         condition, braak_numeric, age, sex, apoe, rin, cerad, region)

message(sprintf(
  "  [GSE268599] Selected donors (BA46): AD=%d / Control=%d (rows=%d)",
  count_condition_donors(meta_ref, "donor_id", "AD"),
  count_condition_donors(meta_ref, "donor_id", "Control"),
  nrow(meta_ref)
))

# ------------------------------------------------------------------------------
# 6.5  Data loading
# ------------------------------------------------------------------------------

seurat_list <- list()

for (i in seq_len(nrow(meta_ref))) {
  gsm_id <- meta_ref$gsm_id[i]
  sid    <- meta_ref$sample_name[i]
  mw_val <- meta_ref$SampleName[i]
  prefix <- file.path(data_dir_268, paste0(gsm_id, "_", mw_val))
  
  mtx_path  <- paste0(prefix, "_matrix.mtx.gz")
  bar_path  <- paste0(prefix, "_barcodes.tsv.gz")
  feat_path <- paste0(prefix, "_features.tsv.gz")
  
  if (!file.exists(mtx_path)) {
    warning(sprintf("  [GSE268599] File not found: %s — skipping.", basename(mtx_path)))
    next
  }
  
  counts     <- ReadMtx(mtx = mtx_path, cells = bar_path, features = feat_path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = gsm_id,
                                   min.cells = 10, min.features = 200)
  
  seurat_obj$study         <- "GSE268599"
  seurat_obj$sample_id     <- sid
  seurat_obj$gsm_id        <- gsm_id
  seurat_obj$donor_id      <- meta_ref$donor_id[i]
  seurat_obj$condition     <- meta_ref$condition[i]
  seurat_obj$braak_numeric <- meta_ref$braak_numeric[i]
  seurat_obj$region        <- meta_ref$region[i]
  seurat_obj$age           <- meta_ref$age[i]
  seurat_obj$sex           <- meta_ref$sex[i]
  seurat_obj$apoe          <- meta_ref$apoe[i]
  seurat_obj$rin           <- meta_ref$rin[i]
  seurat_obj$cerad         <- meta_ref$cerad[i]
  
  seurat_list[[gsm_id]] <- seurat_obj
}

# ------------------------------------------------------------------------------
# 6.6  Merge & QC
# ------------------------------------------------------------------------------

merged_obj <- merge_seurat_list(seurat_list, project = "GSE268599")
merged_obj <- join_layers_if_available(merged_obj)
merged_obj[["percent.mt"]] <- PercentageFeatureSet(merged_obj, pattern = "^MT-")
merged_obj <- subset(merged_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 5000 &
                       percent.mt  < 10)

message(sprintf("  [GSE268599] Post-QC cells: %d", ncol(merged_obj)))
saveRDS(merged_obj, file.path(output_dir, "gse268599_azimuth_input.rds"))

rm(gse, lookup_geo, sample_csv, braak_csv, lookup_table,
   meta_ref, seurat_list, merged_obj)
gc()



# ==============================================================================
# 7. GSE303823 (Chia et al.)
# ==============================================================================
#
# Reference : Chia et al., Brain Commun. 2025 (PMC12359983)
# Design    : snRNA-seq, PFC (BA9), counts provided as RDS (list: counts + meta)
# Criteria  : Relaxed criteria: diagnosis AD/NCI + age >= 65.
#             Age     : all donors 73-92 (paper-confirmed; no filter needed)
#             Excluded= DLB, PDD, other non-AD/non-control groups
# ==============================================================================

message("Processing GSE303823 (Chia et al.)...")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(openxlsx)
  library(Seurat)
})

data_dir_303 <- file.path(data_root, "GSE303823")

# ------------------------------------------------------------------------------
# 7.1  Load donor metadata
# ------------------------------------------------------------------------------

braak_raw           <- read.csv(file.path(data_dir_303, "GSE303823_braak_stage_info.csv"),
                                header = FALSE, stringsAsFactors = FALSE)
colnames(braak_raw) <- braak_raw[3, ]
braak_raw           <- braak_raw[4:nrow(braak_raw), 2:ncol(braak_raw)]

supp_meta <- braak_raw %>%
  rename(
    paper_id      = Samples,
    diagnosis_raw = Diagnosis,
    braak_stage   = `AD grade`
  ) %>%
  filter(paper_id != "") %>%
  mutate(
    braak_numeric = case_when(
      grepl("stage I$",   braak_stage) ~ 1L, grepl("stage II$",  braak_stage) ~ 2L,
      grepl("stage III$", braak_stage) ~ 3L, grepl("stage IV$",  braak_stage) ~ 4L,
      grepl("stage V$",   braak_stage) ~ 5L, grepl("stage VI$",  braak_stage) ~ 6L,
      TRUE                             ~ NA_integer_
    ),
    # Map paper_id → sample_id used in RDS meta$sample
    sample_id = case_when(
      paper_id == "NCI_1" ~ "CTRL_4", paper_id == "NCI_2" ~ "CTRL_5",
      paper_id == "NCI_3" ~ "CTRL_6", paper_id == "NCI_4" ~ "CTRL_7",
      paper_id == "AD_1"  ~ "ADD_4",  paper_id == "AD_2"  ~ "ADD_5",
      paper_id == "AD_3"  ~ "ADD_6",  paper_id == "AD_4"  ~ "ADD_7",
      TRUE                ~ paper_id
    ),
    group = case_when(
      grepl("^NCI", paper_id) ~ "NCI",
      grepl("^AD_", paper_id) ~ "AD",
      grepl("^DLB", paper_id) ~ "DLB",
      grepl("^PDD", paper_id) ~ "PDD",
      TRUE                    ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------
# 7.2  Build lookup table
# ------------------------------------------------------------------------------

lookup_table <- supp_meta %>%
  mutate(
    condition = assign_relaxed_condition(group),
    included_in_analysis = if_else(!is.na(condition), "YES", "NO"),
    exclusion_reason = case_when(
      included_in_analysis == "YES"       ~ NA_character_,
      group %in% c("DLB", "PDD")         ~ paste0("Non-AD group (", group, ")"),
      is.na(braak_numeric)               ~ "Missing/unknown Braak stage",
      TRUE                               ~ paste0("Diagnosis not AD/control: diag=", group)
    ),
    region = "PFC",
    study  = "GSE303823"
  ) %>%
  select(sample_id, paper_id, group, diagnosis_raw, braak_stage, braak_numeric,
         condition, region, study, included_in_analysis, exclusion_reason) %>%
  arrange(
    factor(condition, levels = c("AD", "Control"), exclude = NULL),
    desc(braak_numeric), sample_id
  )

write_csv(lookup_table,  file.path(output_dir, "gse303823_lookup_table.csv"))
write.xlsx(lookup_table, file.path(output_dir, "gse303823_lookup_table.xlsx"),
           asTable = TRUE)

message(sprintf(
  "  [GSE303823] total rows=%d / donors=%d | AD rows=%d / donors=%d | Control rows=%d / donors=%d | excluded rows=%d / donors=%d",
  nrow(lookup_table),
  count_unique_ids(lookup_table, "paper_id"),
  count_condition_rows(lookup_table, "AD"),
  count_condition_donors(lookup_table, "paper_id", "AD"),
  count_condition_rows(lookup_table, "Control"),
  count_condition_donors(lookup_table, "paper_id", "Control"),
  sum(lookup_table$included_in_analysis == "NO"),
  count_excluded_donors(lookup_table, "paper_id")
))

# ------------------------------------------------------------------------------
# 7.3  Sample selection
# ------------------------------------------------------------------------------

meta_ref <- lookup_table %>%
  filter(included_in_analysis == "YES") %>%
  select(sample_id, paper_id, group, condition, braak_numeric, region, study)

message(sprintf(
  "  [GSE303823] Selected donors: AD=%d / Control=%d (rows=%d)",
  count_condition_donors(meta_ref, "paper_id", "AD"),
  count_condition_donors(meta_ref, "paper_id", "Control"),
  nrow(meta_ref)
))

# Named vector lookup (avoids C stack overflow on large cell vectors)
make_lookup <- function(df, key_col, val_col) {
  v <- df[[val_col]]
  names(v) <- df[[key_col]]
  v
}

lut_condition <- make_lookup(meta_ref, "sample_id", "condition")
lut_braak     <- make_lookup(meta_ref, "sample_id", "braak_numeric")
lut_paper_id  <- make_lookup(meta_ref, "sample_id", "paper_id")

# ------------------------------------------------------------------------------
# 7.4  Load counts & create Seurat object
# ------------------------------------------------------------------------------

raw_data     <- readRDS(file.path(data_dir_303, "GSE303823_raw_counts.rds"))
meta_full    <- raw_data$meta
target_ids   <- meta_ref$sample_id
common_cells <- rownames(meta_full)[meta_full$sample %in% target_ids]
missing_sample_ids <- setdiff(target_ids, unique(meta_full$sample))
if (length(missing_sample_ids) > 0L) {
  warning(
    "[GSE303823] Selected sample_id values not found in raw metadata: ",
    paste(missing_sample_ids, collapse = ", ")
  )
}
if (length(common_cells) == 0L) {
  stop("[GSE303823] No cells matched selected sample IDs.")
}
missing_count_cells <- setdiff(common_cells, colnames(raw_data$counts))
if (length(missing_count_cells) > 0L) {
  stop("[GSE303823] Metadata cells missing from count matrix.")
}

seurat_obj <- CreateSeuratObject(
  counts    = raw_data$counts[, common_cells, drop = FALSE],
  meta.data = meta_full[common_cells, ],
  project   = "GSE303823",
  min.cells = 10, min.features = 200
)
rm(raw_data); gc()

sid <- seurat_obj$sample

seurat_obj$study         <- "GSE303823"
seurat_obj$sample_id     <- sid
seurat_obj$donor_id      <- unname(lut_paper_id[sid])
seurat_obj$condition     <- unname(lut_condition[sid])
seurat_obj$braak_numeric <- unname(lut_braak[sid])
seurat_obj$region        <- "PFC"

# ------------------------------------------------------------------------------
# 7.5  QC filtering
# ------------------------------------------------------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 5000 &
                       percent.mt  < 10)

message(sprintf("  [GSE303823] Post-QC cells: %d", ncol(seurat_obj)))
saveRDS(seurat_obj, file.path(output_dir, "gse303823_azimuth_input.rds"))

rm(braak_raw, supp_meta, lookup_table, meta_ref, meta_full,
   common_cells, target_ids, sid, seurat_obj)
gc()

