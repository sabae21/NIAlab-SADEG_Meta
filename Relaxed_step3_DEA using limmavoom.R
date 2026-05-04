################################################################################
# Script: 03_PFC_DEA_using_limma_voom.R
# Description: Study-wise astrocyte pseudo-bulk DEA using edgeR/limma-voom.
#
# Input:  Step 2 dual-evidence astrocyte objects:
#         step2_astrocyte_objects/gse######_astrocytes_step2.rds
#
# Output: One DEA table per study plus donor-level pseudo-bulk metadata.
#
# Notes:
# - DEA uses raw RNA counts aggregated at donor level.
# - The relaxed Step 1 criteria are carried through via obj$condition.
# - GSE243292 is intentionally excluded from Step 2/3 because only one control
#   donor passed the Step 1 criteria.
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(readr)
  library(tibble)
})

step2_dir <- Sys.getenv(
  "PFC_STEP2_DIR",
  unset = file.path("~", "data", "meta_analysis", "results", "PFC", "2026.04")
)
step2_dir <- normalizePath(path.expand(step2_dir), winslash = "/", mustWork = FALSE)

astro_dir <- file.path(step2_dir, "step2_astrocyte_objects")
output_dir <- file.path(step2_dir, "DEA_Results_limma_voom")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Primary analysis does not exclude pseudo-bulk samples by astrocyte cell count.
# Set this to a positive value only for a prespecified sensitivity analysis.
MIN_CELLS_PER_PSEUDOBULK <- 0
MIN_DONORS_PER_GROUP <- 2
CPM_CUTOFF <- 0.5
MIN_SAMPLES_EXPRESSED <- 2

datasets <- tibble::tribble(
  ~dataset_id, ~study_id,   ~input_file,
  "gse157827", "GSE157827", "gse157827_astrocytes_step2.rds",
  "gse167490", "GSE167490", "gse167490_astrocytes_step2.rds",
  "gse167492", "GSE167492", "gse167492_astrocytes_step2.rds",
  "gse214979", "GSE214979", "gse214979_astrocytes_step2.rds",
  "gse263468", "GSE263468", "gse263468_astrocytes_step2.rds",
  "gse268599", "GSE268599", "gse268599_astrocytes_step2.rds",
  "gse303823", "GSE303823", "gse303823_astrocytes_step2.rds"
) %>%
  mutate(input_path = file.path(astro_dir, input_file))

missing_inputs <- datasets$input_path[!file.exists(datasets$input_path)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Step 2 astrocyte object(s):\n  ",
    paste(missing_inputs, collapse = "\n  "),
    "\nFinish Step 2 before running Step 3."
  )
}

# ==============================================================================
# 1. Helpers
# ==============================================================================

safe_id <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", ".", x)
  x <- gsub("^\\.|\\.$", "", x)
  paste0("S.", x)
}

first_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (is.character(x)) {
    x <- x[x != ""]
  }
  if (length(x) == 0L) {
    return(NA)
  }
  x[1L]
}

choose_pseudobulk_variable <- function(meta, study_id) {
  candidates <- c("donor_id", "sample_id", "orig.ident")
  for (candidate in candidates) {
    if (candidate %in% colnames(meta)) {
      x <- as.character(meta[[candidate]])
      if (!all(is.na(x)) && length(unique(na.omit(x))) > 1L) {
        return(candidate)
      }
    }
  }
  stop("[", study_id, "] No usable pseudo-bulk grouping variable found.")
}

get_counts_matrix <- function(obj, assay = "RNA", study_id = NULL) {
  if (!assay %in% names(obj@assays)) {
    stop("[", study_id, "] Assay not found: ", assay)
  }
  
  assay_obj <- obj[[assay]]
  layers <- if ("Layers" %in% getNamespaceExports("SeuratObject")) {
    SeuratObject::Layers(assay_obj)
  } else {
    character(0)
  }
  
  if ("counts" %in% layers) {
    counts <- GetAssayData(obj, assay = assay, layer = "counts")
    message(sprintf(
      "    [%s] Using RNA counts layer: %d genes x %d cells",
      study_id, nrow(counts), ncol(counts)
    ))
    return(counts)
  }
  
  counts_layers <- grep("^counts", layers, value = TRUE)
  if (length(counts_layers) == 0L) {
    stop(
      "[", study_id, "] No RNA counts layer found. Available RNA layers: ",
      paste(layers, collapse = ", ")
    )
  }
  
  message(sprintf(
    "    [%s] Found %d split RNA counts layers; combining them.",
    study_id, length(counts_layers)
  ))
  
  layer_mats <- lapply(counts_layers, function(layer_name) {
    mat <- SeuratObject::LayerData(obj, assay = assay, layer = layer_name)
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop("[", study_id, "] Layer lacks rownames or colnames: ", layer_name)
    }
    mat
  })
  
  all_features <- Reduce(union, lapply(layer_mats, rownames))
  all_cells <- unlist(lapply(layer_mats, colnames), use.names = FALSE)
  if (anyDuplicated(all_cells) > 0L) {
    stop("[", study_id, "] Duplicated cell barcodes across split counts layers.")
  }
  
  layer_mats <- lapply(layer_mats, function(mat) {
    missing_features <- setdiff(all_features, rownames(mat))
    if (length(missing_features) > 0L) {
      zero_mat <- Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        dims = c(length(missing_features), ncol(mat)),
        dimnames = list(missing_features, colnames(mat))
      )
      mat <- rbind(mat, zero_mat)
    }
    mat[all_features, , drop = FALSE]
  })
  
  counts <- do.call(cbind, layer_mats)
  counts <- counts[, all_cells, drop = FALSE]
  message(sprintf(
    "    [%s] Combined RNA counts layers: %d genes x %d cells",
    study_id, nrow(counts), ncol(counts)
  ))
  counts
}

make_pseudobulk <- function(obj, study_id) {
  DefaultAssay(obj) <- "RNA"
  counts <- get_counts_matrix(obj, assay = "RNA", study_id = study_id)
  meta <- obj@meta.data
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  if (length(common_cells) == 0L) {
    stop("[", study_id, "] No shared cell names between counts and metadata.")
  }
  
  counts <- counts[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  
  if (!"condition" %in% colnames(meta)) {
    stop("[", study_id, "] condition column is missing from metadata.")
  }
  
  pb_var <- choose_pseudobulk_variable(meta, study_id)
  meta$pb_source_id <- as.character(meta[[pb_var]])
  meta$pb_id <- safe_id(meta$pb_source_id)
  meta$condition <- as.character(meta$condition)
  
  keep_cells <- !is.na(meta$pb_id) &
    meta$pb_id != "S." &
    meta$condition %in% c("AD", "Control")
  
  counts <- counts[, keep_cells, drop = FALSE]
  meta <- meta[keep_cells, , drop = FALSE]
  
  if (ncol(counts) == 0L) {
    stop("[", study_id, "] No AD/Control astrocytes available for DEA.")
  }
  
  donor_condition_check <- meta %>%
    select(pb_id, condition) %>%
    distinct() %>%
    count(pb_id, name = "n_conditions") %>%
    filter(n_conditions > 1L)
  if (nrow(donor_condition_check) > 0L) {
    stop("[", study_id, "] A pseudo-bulk ID maps to multiple conditions: ",
         paste(donor_condition_check$pb_id, collapse = ", "))
  }
  
  pb_meta <- meta %>%
    mutate(
      donor_id = if ("donor_id" %in% colnames(meta)) as.character(donor_id) else NA_character_,
      sample_id = if ("sample_id" %in% colnames(meta)) as.character(sample_id) else NA_character_,
      age = if ("age" %in% colnames(meta)) suppressWarnings(as.numeric(age)) else NA_real_,
      sex = if ("sex" %in% colnames(meta)) as.character(sex) else NA_character_,
      braak_numeric = if ("braak_numeric" %in% colnames(meta)) {
        suppressWarnings(as.numeric(braak_numeric))
      } else {
        NA_real_
      }
    ) %>%
    group_by(pb_id, pb_source_id, condition) %>%
    summarize(
      donor_id = first_or_na(donor_id),
      sample_id = first_or_na(sample_id),
      age = first_or_na(age),
      sex = first_or_na(sex),
      braak_numeric = first_or_na(braak_numeric),
      n_astrocytes = n(),
      .groups = "drop"
    ) %>%
    mutate(
      study = study_id,
      donor_id = ifelse(is.na(donor_id) | donor_id == "", pb_source_id, donor_id),
      sample_id = ifelse(is.na(sample_id) | sample_id == "", pb_source_id, sample_id)
    )
  
  if (MIN_CELLS_PER_PSEUDOBULK > 0L) {
    low_cell_pb <- pb_meta %>%
      filter(n_astrocytes < MIN_CELLS_PER_PSEUDOBULK) %>%
      pull(pb_id)
    if (length(low_cell_pb) > 0L) {
      message(sprintf(
        "    [%s] Excluding %d pseudo-bulk sample(s) with <%d astrocytes.",
        study_id, length(low_cell_pb), MIN_CELLS_PER_PSEUDOBULK
      ))
    }
    
    keep_pb <- setdiff(pb_meta$pb_id, low_cell_pb)
    keep_cells <- meta$pb_id %in% keep_pb
    counts <- counts[, keep_cells, drop = FALSE]
    meta <- meta[keep_cells, , drop = FALSE]
    pb_meta <- pb_meta %>%
      filter(pb_id %in% keep_pb)
  }
  
  pb_meta <- pb_meta %>% arrange(pb_id)
  
  group <- factor(meta$pb_id, levels = pb_meta$pb_id)
  design_pb <- Matrix::sparse.model.matrix(~ 0 + group)
  colnames(design_pb) <- levels(group)
  
  pb_counts <- counts %*% design_pb
  pb_counts <- as.matrix(pb_counts)
  storage.mode(pb_counts) <- "integer"
  pb_counts <- pb_counts[, pb_meta$pb_id, drop = FALSE]
  
  list(counts = pb_counts, metadata = pb_meta, pb_var = pb_var)
}

run_limma_voom <- function(pb_counts, pb_meta, study_id) {
  rownames(pb_meta) <- pb_meta$pb_id
  pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
  
  n_ad <- sum(pb_meta$condition == "AD")
  n_control <- sum(pb_meta$condition == "Control")
  message(sprintf(
    "    [%s] Pseudo-bulk donors: %d total (AD=%d / Control=%d)",
    study_id, ncol(pb_counts), n_ad, n_control
  ))
  
  if (n_ad < MIN_DONORS_PER_GROUP || n_control < MIN_DONORS_PER_GROUP) {
    warning(sprintf(
      "[%s] Skipping DEA: need >=%d donors per group (AD=%d, Control=%d).",
      study_id, MIN_DONORS_PER_GROUP, n_ad, n_control
    ))
    return(NULL)
  }
  
  y <- DGEList(counts = pb_counts, samples = pb_meta)
  group <- factor(y$samples$condition, levels = c("Control", "AD"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  keep <- rowSums(cpm(y) > CPM_CUTOFF) >= MIN_SAMPLES_EXPRESSED
  message(sprintf(
    "    [%s] Genes retained after CPM filtering: %d / %d",
    study_id, sum(keep), nrow(y)
  ))
  if (sum(keep) < 100L) {
    warning(sprintf("[%s] Skipping DEA: too few genes passed CPM filtering.", study_id))
    return(NULL)
  }
  
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  contrast <- makeContrasts(AD_vs_Control = AD - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  res <- topTable(fit2, coef = "AD_vs_Control", number = Inf, sort.by = "none")
  res$SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, "AD_vs_Control"]
  res$gene <- rownames(res)
  res$study <- study_id
  res$n_AD <- n_ad
  res$n_Control <- n_control
  
  res %>%
    select(gene, study, logFC, SE, AveExpr, t, P.Value, adj.P.Val,
           B, n_AD, n_Control) %>%
    rename(PValue = P.Value, FDR = adj.P.Val)
}

# ==============================================================================
# 2. Run study-wise DEA
# ==============================================================================

message("Starting study-wise astrocyte pseudo-bulk DEA...")

dea_results <- list()
dataset_summary <- list()

for (i in seq_len(nrow(datasets))) {
  study_id <- datasets$study_id[i]
  dataset_id <- datasets$dataset_id[i]
  input_path <- datasets$input_path[i]
  
  message(sprintf("\n>>> Processing %s (%s)", study_id, basename(input_path)))
  obj <- readRDS(input_path)
  
  pb <- make_pseudobulk(obj, study_id = study_id)
  
  write_csv(
    pb$metadata,
    file.path(output_dir, paste0("pseudobulk_metadata_", study_id, ".csv"))
  )
  
  dea_res <- run_limma_voom(pb$counts, pb$metadata, study_id = study_id)
  if (!is.null(dea_res)) {
    out_csv <- file.path(output_dir, paste0("DEA_", study_id, ".csv"))
    write_csv(dea_res, out_csv)
    dea_results[[study_id]] <- dea_res
    message(sprintf("    [%s] Saved: %s", study_id, basename(out_csv)))
  }
  
  dataset_summary[[study_id]] <- tibble(
    study = study_id,
    dataset_id = dataset_id,
    input_cells = ncol(obj),
    pseudobulk_variable = pb$pb_var,
    pseudobulk_samples = nrow(pb$metadata),
    AD_pseudobulk_samples = sum(pb$metadata$condition == "AD"),
    Control_pseudobulk_samples = sum(pb$metadata$condition == "Control"),
    total_astrocytes_in_pseudobulk = sum(pb$metadata$n_astrocytes),
    dea_completed = !is.null(dea_res)
  )
  
  rm(obj, pb, dea_res)
  gc()
}

dataset_summary <- bind_rows(dataset_summary)
write_csv(dataset_summary, file.path(output_dir, "DEA_dataset_summary.csv"))

if (length(dea_results) > 0L) {
  dea_all <- bind_rows(dea_results)
  write_csv(dea_all, file.path(output_dir, "DEA_all_studies_long.csv"))
}

message("\nDEA pipeline complete.")
