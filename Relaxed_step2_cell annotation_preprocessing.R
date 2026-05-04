################################################################################
# Script: 02_PFC_Cell_Annotation_and_Preprocessing.R
# Description: Azimuth cell type annotation, BRETIGEA hybrid detection,
#              astrocyte subsetting, and Harmony batch correction.
# Datasets:    GSE157827, GSE167490, GSE167492, GSE214979, GSE263468,
#              GSE268599, GSE303823
#
# Note: GSE243292 is intentionally excluded because only one control donor
#       passed the Step 1 criteria.
################################################################################

# ==============================================================================
# 0. Setup
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
  library(SeuratObject)
  library(BRETIGEA)
  library(dplyr)
  library(readr)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

input_dir <- Sys.getenv(
  "PFC_STEP1_DIR",
  unset = file.path("~", "data", "meta_analysis", "results", "PFC", "2026.04", "Azimuth_input")
)
output_dir <- Sys.getenv(
  "PFC_STEP2_DIR",
  unset = file.path("~", "data", "meta_analysis", "results", "PFC", "2026.04")
)
ref_dir <- Sys.getenv(
  "AZIMUTH_REF_DIR",
  unset = file.path("~", "data", "meta_analysis", "azimuth_references")
)

input_dir <- normalizePath(path.expand(input_dir), winslash = "/", mustWork = FALSE)
output_dir <- normalizePath(path.expand(output_dir), winslash = "/", mustWork = FALSE)
ref_dir <- normalizePath(path.expand(ref_dir), winslash = "/", mustWork = FALSE)

annotated_dir <- file.path(output_dir, "step2_azimuth_annotated")
astro_dir <- file.path(output_dir, "step2_astrocyte_objects")
plot_dir <- file.path(output_dir, "step2_plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotated_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(astro_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Reuse existing per-dataset outputs only after confirming they were generated
# from the same Step 1 input RDS files. Keep FALSE for clean paper runs.
REUSE_ANNOTATED_OBJECTS <- FALSE

# Set TRUE only when the reference folder may contain files from different
# downloads. This forces a fresh matching ref.Rds / idx.annoy pair.
FORCE_REFERENCE_DOWNLOAD <- FALSE

# Keep FALSE for the main paper workflow. If TRUE, a rescue path attempts
# Azimuth label transfer without the final reference UMAP coordinate projection
# when full RunAzimuth fails. This is useful for troubleshooting, but the
# primary analysis should use full RunAzimuth whenever possible.
USE_AZIMUTH_LABEL_ONLY_FALLBACK <- FALSE

# Rebuild each Step 1 Seurat object as a simple RNA-counts-only query before
# RunAzimuth. This removes saved reductions, graphs, and complex Seurat v5 layer
# structure that can trigger matrix-dimension errors inside Azimuth projection.
RECREATE_QUERY_FOR_AZIMUTH <- TRUE

datasets <- tibble::tribble(
  ~dataset_id, ~study_id,    ~input_file,
  "gse157827", "GSE157827",  "gse157827_azimuth_input.rds",
  "gse167490", "GSE167490",  "gse167490_azimuth_input.rds",
  "gse167492", "GSE167492",  "gse167492_azimuth_input.rds",
  "gse214979", "GSE214979",  "gse214979_azimuth_input.rds",
  "gse263468", "GSE263468",  "gse263468_azimuth_input.rds",
  "gse268599", "GSE268599",  "gse268599_azimuth_input.rds",
  "gse303823", "GSE303823",  "gse303823_azimuth_input.rds"
) %>%
  mutate(
    input_path = file.path(input_dir, input_file),
    annotated_path = file.path(annotated_dir,
                               paste0(dataset_id, "_azimuth_bretigea.rds")),
    astro_path = file.path(astro_dir,
                           paste0(dataset_id, "_astrocytes_step2.rds"))
  )

missing_inputs <- datasets$input_path[!file.exists(datasets$input_path)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Step 1 input RDS file(s):\n  ",
    paste(missing_inputs, collapse = "\n  ")
  )
}

join_layers_if_available <- function(obj) {
  if ("JoinLayers" %in% getNamespaceExports("Seurat") &&
      "RNA" %in% names(obj@assays)) {
    obj[["RNA"]] <- Seurat::JoinLayers(obj[["RNA"]])
  }
  DefaultAssay(obj) <- "RNA"
  obj
}

get_counts_matrix <- function(obj, assay = "RNA", project = NULL) {
  assay_obj <- obj[[assay]]
  layers <- if ("Layers" %in% getNamespaceExports("SeuratObject")) {
    SeuratObject::Layers(assay_obj)
  } else {
    character(0)
  }
  
  if ("counts" %in% layers) {
    counts <- GetAssayData(obj, assay = assay, layer = "counts")
    if (!is.null(project)) {
      message(sprintf(
        "  [%s] Using joined RNA counts layer: %d genes x %d cells",
        project, nrow(counts), ncol(counts)
      ))
    }
    return(counts)
  }
  
  counts_layers <- grep("^counts", layers, value = TRUE)
  if (length(counts_layers) == 0L) {
    stop(
      "No RNA counts layer found. Available RNA layers: ",
      paste(layers, collapse = ", ")
    )
  }
  
  if (!is.null(project)) {
    message(sprintf(
      "  [%s] Found %d split RNA counts layers; combining them before Azimuth.",
      project, length(counts_layers)
    ))
  }
  
  layer_mats <- lapply(counts_layers, function(layer_name) {
    mat <- tryCatch(
      SeuratObject::LayerData(obj, assay = assay, layer = layer_name),
      error = function(e) GetAssayData(obj, assay = assay, layer = layer_name)
    )
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop("Layer lacks rownames or colnames: ", layer_name)
    }
    mat
  })
  names(layer_mats) <- counts_layers
  
  all_features <- Reduce(union, lapply(layer_mats, rownames))
  all_cells <- unlist(lapply(layer_mats, colnames), use.names = FALSE)
  if (anyDuplicated(all_cells) > 0L) {
    stop(
      "Duplicated cell barcodes across split RNA counts layers. ",
      "Layer-specific prefixes should be preserved before combining."
    )
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
  if (!is.null(project)) {
    message(sprintf(
      "  [%s] Combined split RNA counts layers: %d genes x %d cells",
      project, nrow(counts), ncol(counts)
    ))
  }
  counts
}

prepare_azimuth_query <- function(obj, project) {
  obj <- join_layers_if_available(obj)
  counts <- get_counts_matrix(obj, assay = "RNA", project = project)
  meta <- obj@meta.data
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  if (length(common_cells) == 0L) {
    stop("[", project, "] No shared cell names between RNA counts and metadata.")
  }
  counts <- counts[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  
  if (anyDuplicated(rownames(counts)) > 0L) {
    stop("[", project, "] Duplicated gene names found in RNA counts.")
  }
  if (anyDuplicated(colnames(counts)) > 0L) {
    stop("[", project, "] Duplicated cell barcodes found in RNA counts.")
  }
  
  query <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    project = project,
    min.cells = 0,
    min.features = 0
  )
  DefaultAssay(query) <- "RNA"
  query
}

drop_heavy_slots <- function(obj, keep_data = FALSE) {
  DietSeurat(
    obj,
    counts = TRUE,
    data = keep_data,
    scale.data = FALSE,
    dimreducs = NULL,
    graphs = NULL
  )
}

make_lookup <- function(df, key_col, val_col) {
  v <- df[[val_col]]
  names(v) <- df[[key_col]]
  v
}

extract_azimuth_umap <- function(obj, study_id, dataset_id) {
  azimuth_umap_name <- intersect(
    c("ref.umap", "umap", "wnn.umap", "predicted.umap"),
    Reductions(obj)
  )[1]
  if (is.na(azimuth_umap_name)) {
    return(NULL)
  }
  
  label_col <- intersect(
    c("azimuth_subclass", "predicted.subclass", "final_label"),
    colnames(obj@meta.data)
  )[1]
  if (is.na(label_col)) {
    return(NULL)
  }
  
  umap_df <- as.data.frame(Embeddings(obj, azimuth_umap_name)[, 1:2, drop = FALSE])
  colnames(umap_df) <- c("umap_1", "umap_2")
  umap_df$cell_barcode <- rownames(umap_df)
  umap_df$study <- study_id
  umap_df$dataset_id <- dataset_id
  umap_df$azimuth_subclass <- obj@meta.data[rownames(umap_df), label_col]
  umap_df
}

run_azimuth_labels_only <- function(query, reference_dir, assay = "RNA",
                                    annotation.levels = NULL,
                                    k.weight = 50,
                                    n.trees = 20,
                                    mapping.score.k = 100,
                                    verbose = TRUE) {
  if (!dir.exists(reference_dir)) {
    stop("Azimuth reference directory not found: ", reference_dir)
  }
  
  reference <- Azimuth:::LoadReference(reference_dir)$map
  if (is.null(reference)) {
    stop("Could not load Azimuth map reference from: ", reference_dir)
  }
  
  if (!"num_precomputed_nns" %in% names(Misc(reference[["refUMAP"]])$model)) {
    Misc(reference[["refUMAP"]], slot = "model")$num_precomputed_nns <- 1
  }
  
  key.pattern <- "^[^_]*_"
  new.colnames <- gsub(
    pattern = key.pattern,
    replacement = Key(reference[["refDR"]]),
    x = colnames(Loadings(object = reference[["refDR"]], projected = FALSE))
  )
  colnames(Loadings(object = reference[["refDR"]], projected = FALSE)) <- new.colnames
  
  dims <- as.double(slot(reference, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)
  if (is.null(annotation.levels)) {
    annotation.levels <- names(slot(object = reference, name = "meta.data"))
    annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
  }
  
  query <- Azimuth:::ConvertGeneNames(
    object = query,
    reference.names = rownames(reference),
    homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds"
  )
  
  if (!all(c("nCount_RNA", "nFeature_RNA") %in% colnames(query[[]]))) {
    calcn <- as.data.frame(Seurat:::CalcN(object = query[[assay]]))
    colnames(calcn) <- paste(colnames(calcn), assay, sep = "_")
    query <- AddMetaData(object = query, metadata = calcn)
    rm(calcn)
  }
  
  if (any(grepl(pattern = "^MT-", x = rownames(query)))) {
    query <- PercentageFeatureSet(
      object = query,
      pattern = "^MT-",
      col.name = "percent.mt",
      assay = assay
    )
  }
  
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = assay,
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = rownames(Loadings(reference[["refDR"]])),
    dims = 1:dims,
    n.trees = n.trees,
    mapping.score.k = mapping.score.k,
    verbose = verbose
  )
  
  refdata <- lapply(annotation.levels, function(x) {
    reference[[x, drop = TRUE]]
  })
  names(refdata) <- annotation.levels
  
  query <- TransferData(
    reference = reference,
    query = query,
    query.assay = assay,
    dims = 1:dims,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE,
    k.weight = k.weight,
    verbose = verbose
  )
  
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference,
    query = query,
    query.assay = assay,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    verbose = verbose
  )
  
  query <- AddMetaData(
    object = query,
    metadata = Azimuth:::MappingScore(anchors = anchors, ndim = dims),
    col.name = "mapping.score"
  )
  
  query
}

# Azimuth Human Motor Cortex Reference (Bakken et al. 2021; Allen Institute)
# https://azimuth.hubmapconsortium.org/references/human_motorcortex/
reference_files <- file.path(ref_dir, c("idx.annoy", "ref.Rds"))
if (FORCE_REFERENCE_DOWNLOAD && any(file.exists(reference_files))) {
  unlink(reference_files[file.exists(reference_files)])
}

if (!all(file.exists(reference_files))) {
  message("Downloading Azimuth Human Motor Cortex reference...")
  options(timeout = 3600)
  download.file(
    "https://zenodo.org/records/4546932/files/idx.annoy?download=1",
    destfile = file.path(ref_dir, "idx.annoy"),
    mode = "wb"
  )
  download.file(
    "https://zenodo.org/records/4546932/files/ref.Rds?download=1",
    destfile = file.path(ref_dir, "ref.Rds"),
    mode = "wb"
  )
}

missing_reference <- reference_files[!file.exists(reference_files)]
if (length(missing_reference) > 0L) {
  stop("Azimuth reference download failed; missing:\n  ",
       paste(missing_reference, collapse = "\n  "))
}
message(sprintf("Using Azimuth reference directory: %s", ref_dir))
message(sprintf(
  "  ref.Rds=%s bytes | idx.annoy=%s bytes",
  file.info(file.path(ref_dir, "ref.Rds"))$size,
  file.info(file.path(ref_dir, "idx.annoy"))$size
))

# ==============================================================================
# 1. BRETIGEA marker setup
# ==============================================================================
message("Preparing BRETIGEA marker sets...")

data("markers_df_human_brain", package = "BRETIGEA")
cell_types <- unique(markers_df_human_brain$cell)
brain_markers <- lapply(cell_types, function(ct) {
  head(markers_df_human_brain$markers[markers_df_human_brain$cell == ct], 50)
})
names(brain_markers) <- cell_types

if (!"ast" %in% cell_types) {
  stop("BRETIGEA marker object does not contain the expected 'ast' cell type.")
}

# ==============================================================================
# 2. Per-dataset Azimuth annotation and astrocyte extraction
# ==============================================================================
message("Running per-dataset Azimuth annotation and BRETIGEA filtering...")

annotation_summary_list <- list()
astro_funnel_list <- list()
all_cell_umap_list <- list()

for (i in seq_len(nrow(datasets))) {
  ds <- datasets$dataset_id[i]
  study_id <- datasets$study_id[i]
  in_path <- datasets$input_path[i]
  annotated_path <- datasets$annotated_path[i]
  astro_path <- datasets$astro_path[i]
  
  message(sprintf("  [%s] Processing %s", study_id, basename(in_path)))
  
  if (REUSE_ANNOTATED_OBJECTS && file.exists(astro_path)) {
    message(sprintf("  [%s] Existing astrocyte object found; reusing.", study_id))
    astro_obj <- readRDS(astro_path)
    annotation_summary_list[[study_id]] <- read_csv(
      file.path(annotated_dir, paste0(ds, "_annotation_summary.csv")),
      show_col_types = FALSE
    )
    astro_funnel_list[[study_id]] <- read_csv(
      file.path(annotated_dir, paste0(ds, "_astrocyte_filtering_funnel.csv")),
      show_col_types = FALSE
    )
    if (file.exists(annotated_path)) {
      obj_slim <- readRDS(annotated_path)
      umap_df <- extract_azimuth_umap(obj_slim, study_id, ds)
      if (!is.null(umap_df)) all_cell_umap_list[[study_id]] <- umap_df
      rm(obj_slim, umap_df)
    }
    rm(astro_obj); gc()
    next
  }
  
  obj <- readRDS(in_path)
  if (RECREATE_QUERY_FOR_AZIMUTH) {
    obj <- prepare_azimuth_query(obj, project = study_id)
  } else {
    obj <- join_layers_if_available(obj)
  }
  
  obj$study <- study_id
  obj$dataset_id <- ds
  message(sprintf(
    "  [%s] Query prepared for Azimuth: %d genes x %d cells",
    study_id, nrow(obj), ncol(obj)
  ))
  
  obj <- tryCatch(
    RunAzimuth(obj, reference = ref_dir),
    error = function(e) {
      err_msg <- conditionMessage(e)
      if (USE_AZIMUTH_LABEL_ONLY_FALLBACK &&
          grepl("non-conformable arguments", err_msg, fixed = TRUE)) {
        warning(
          "[", study_id, "] Full RunAzimuth failed. Falling back to ",
          "Azimuth label transfer without the final reference UMAP coordinate ",
          "projection. Use this only as a troubleshooting/rescue mode.\n",
          "Original error: ", err_msg,
          call. = FALSE
        )
        return(run_azimuth_labels_only(
          query = obj,
          reference_dir = ref_dir,
          assay = "RNA",
          verbose = TRUE
        ))
      }
      stop(
        "[", study_id, "] Full RunAzimuth failed.\n",
        "Original error: ", err_msg, "\n",
        "The script is intentionally not using label-only fallback for the ",
        "main workflow. First check package/reference compatibility and the ",
        "RNA layer join step.",
        call. = FALSE
      )
    }
  )
  if (!"predicted.subclass" %in% colnames(obj@meta.data)) {
    stop(sprintf("[%s] Azimuth did not return predicted.subclass.", study_id))
  }
  obj$azimuth_subclass <- obj$predicted.subclass
  obj$final_label <- obj$azimuth_subclass
  
  umap_df <- extract_azimuth_umap(obj, study_id, ds)
  if (!is.null(umap_df)) {
    all_cell_umap_list[[study_id]] <- umap_df
  } else {
    warning(sprintf(
      "[%s] No Azimuth UMAP reduction found; all-cell overview UMAP will skip this study.",
      study_id
    ))
  }
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- AddModuleScore(obj, features = brain_markers, name = "BRET_")
  
  score_cols <- grep("^BRET_", colnames(obj@meta.data), value = TRUE)
  if (length(score_cols) != length(cell_types)) {
    stop(sprintf("[%s] BRETIGEA score column mismatch: expected %d, found %d.",
                 study_id, length(cell_types), length(score_cols)))
  }
  colnames(obj@meta.data)[colnames(obj@meta.data) %in% score_cols] <-
    paste0("Score_", cell_types)
  
  score_mat <- as.matrix(obj@meta.data[, paste0("Score_", cell_types)])
  obj$max_score <- apply(score_mat, 1, max)
  obj$sec_score <- apply(score_mat, 1, function(x) sort(x, decreasing = TRUE)[2])
  obj$hybrid_ratio <- ifelse(
    obj$max_score > 0,
    (obj$max_score - obj$sec_score) / obj$max_score,
    1
  )
  obj$is_hybrid <- obj$hybrid_ratio < 0.2
  obj$bret_top_cell <- cell_types[apply(score_mat, 1, which.max)]
  obj$final_label <- ifelse(obj$is_hybrid, "Hybrid", obj$azimuth_subclass)
  
  total_cells <- ncol(obj)
  n_azimuth_astro <- sum(obj$azimuth_subclass == "Astro", na.rm = TRUE)
  n_final_astro <- sum(
    obj$azimuth_subclass == "Astro" &
      obj$bret_top_cell == "ast" &
      !obj$is_hybrid,
    na.rm = TRUE
  )
  
  annotation_summary <- as.data.frame(table(final_label = obj$final_label)) %>%
    rename(cell_type = final_label, n_cells = Freq) %>%
    mutate(
      study = study_id,
      dataset_id = ds,
      percentage = round(n_cells / total_cells * 100, 2)
    ) %>%
    arrange(desc(n_cells))
  
  astro_funnel <- tibble::tibble(
    study = study_id,
    dataset_id = ds,
    step = c(
      "1. All cells",
      "2. Azimuth: Astro",
      "3. Dual-evidence astrocytes",
      "4. Azimuth Astro not retained"
    ),
    n_cells = c(
      total_cells,
      n_azimuth_astro,
      n_final_astro,
      n_azimuth_astro - n_final_astro
    ),
    pct_of_total = round(c(
      100,
      n_azimuth_astro / total_cells * 100,
      n_final_astro / total_cells * 100,
      (n_azimuth_astro - n_final_astro) / total_cells * 100
    ), 2)
  )
  
  write_csv(annotation_summary,
            file.path(annotated_dir, paste0(ds, "_annotation_summary.csv")))
  write_csv(astro_funnel,
            file.path(annotated_dir, paste0(ds, "_astrocyte_filtering_funnel.csv")))
  
  annotation_summary_list[[study_id]] <- annotation_summary
  astro_funnel_list[[study_id]] <- astro_funnel
  
  message(sprintf(
    "  [%s] cells=%d | Azimuth Astro=%d | dual-evidence astrocytes=%d | hybrids=%d",
    study_id, total_cells, n_azimuth_astro, n_final_astro,
    sum(obj$is_hybrid, na.rm = TRUE)
  ))
  
  astro_obj <- subset(
    obj,
    subset = azimuth_subclass == "Astro" &
      bret_top_cell == "ast" &
      is_hybrid == FALSE
  )
  if (ncol(astro_obj) < 100) {
    stop(sprintf("[%s] Fewer than 100 dual-evidence astrocytes retained.", study_id))
  }
  
  astro_obj <- RenameCells(astro_obj, add.cell.id = study_id)
  astro_obj <- drop_heavy_slots(astro_obj, keep_data = FALSE)
  saveRDS(astro_obj, astro_path)
  
  obj_slim <- drop_heavy_slots(obj, keep_data = FALSE)
  saveRDS(obj_slim, annotated_path)
  
  rm(obj, obj_slim, astro_obj, score_mat, annotation_summary, astro_funnel)
  gc()
}

annotation_summary_all <- bind_rows(annotation_summary_list)
astro_funnel_all <- bind_rows(astro_funnel_list)
all_cell_umap <- bind_rows(all_cell_umap_list)

write_csv(annotation_summary_all,
          file.path(output_dir, "cell_annotation_summary_by_study.csv"))
write_csv(astro_funnel_all,
          file.path(output_dir, "astrocyte_filtering_funnel_by_study.csv"))

annotation_summary_overall <- annotation_summary_all %>%
  group_by(cell_type) %>%
  summarize(n_cells = sum(n_cells), .groups = "drop") %>%
  mutate(percentage = round(n_cells / sum(n_cells) * 100, 2)) %>%
  arrange(desc(n_cells))

write_csv(annotation_summary_overall,
          file.path(output_dir, "cell_annotation_summary_overall.csv"))

total_cells <- sum(annotation_summary_overall$n_cells)
n_azimuth_astro <- annotation_summary_overall$n_cells[
  annotation_summary_overall$cell_type == "Astro"
]
if (length(n_azimuth_astro) == 0L) n_azimuth_astro <- 0L

annotation_summary <- bind_rows(
  annotation_summary_overall %>% filter(cell_type == "Astro"),
  annotation_summary_overall %>% filter(cell_type != "Astro"),
  tibble::tibble(cell_type = "Total", n_cells = total_cells, percentage = 100)
)
write_csv(annotation_summary,
          file.path(output_dir, "cell_annotation_summary.csv"))

astro_funnel <- astro_funnel_all %>%
  group_by(step) %>%
  summarize(n_cells = sum(n_cells), .groups = "drop") %>%
  mutate(
    pct_of_total = round(n_cells / total_cells * 100, 2),
    step = factor(
      step,
      levels = c(
        "1. All cells",
        "2. Azimuth: Astro",
        "3. Dual-evidence astrocytes",
        "4. Azimuth Astro not retained"
      )
    )
  ) %>%
  arrange(step) %>%
  mutate(step = as.character(step))

write_csv(astro_funnel,
          file.path(output_dir, "astrocyte_filtering_funnel.csv"))

astro_fraction_diagnostics <- astro_funnel_all %>%
  select(study, dataset_id, step, n_cells) %>%
  tidyr::pivot_wider(names_from = step, values_from = n_cells) %>%
  mutate(
    pct_azimuth_astro = round(`2. Azimuth: Astro` / `1. All cells` * 100, 2),
    pct_dual_evidence_astro = round(
      `3. Dual-evidence astrocytes` / `1. All cells` * 100, 2
    ),
    pct_of_total_cells = round(`1. All cells` / sum(`1. All cells`) * 100, 2),
    pct_of_total_dual_evidence_astro = round(
      `3. Dual-evidence astrocytes` /
        sum(`3. Dual-evidence astrocytes`) * 100,
      2
    )
  ) %>%
  arrange(pct_azimuth_astro)

write_csv(
  astro_fraction_diagnostics,
  file.path(output_dir, "astrocyte_fraction_diagnostics_by_study.csv")
)

n_final_astro <- astro_funnel$n_cells[
  astro_funnel$step == "3. Dual-evidence astrocytes"
]
if (length(n_final_astro) == 0L) n_final_astro <- 0L

message(sprintf(
  "  Total cells: %d | Azimuth Astro: %d (%.1f%%) | Dual-evidence astrocytes: %d (%.1f%%)",
  total_cells,
  n_azimuth_astro,
  n_azimuth_astro / total_cells * 100,
  n_final_astro,
  n_final_astro / total_cells * 100
))

if (nrow(all_cell_umap) > 0L) {
  write_csv(all_cell_umap,
            file.path(output_dir, "all_cells_azimuth_umap_coordinates.csv"))
  
  cell_type_levels <- sort(unique(all_cell_umap$azimuth_subclass))
  cell_type_cols <- setNames(scales::hue_pal()(length(cell_type_levels)),
                             cell_type_levels)
  cell_type_label_df <- all_cell_umap %>%
    group_by(azimuth_subclass) %>%
    summarize(
      umap_1 = median(umap_1, na.rm = TRUE),
      umap_2 = median(umap_2, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
  
  p_all <- ggplot(
    all_cell_umap,
    aes(x = umap_1, y = umap_2, color = azimuth_subclass)
  ) +
    geom_point(size = 0.08, alpha = 0.45) +
    scale_color_manual(values = cell_type_cols, breaks = cell_type_levels) +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    theme_classic(base_size = 11) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_all <- p_all +
      ggrepel::geom_label_repel(
        data = cell_type_label_df,
        aes(x = umap_1, y = umap_2, label = azimuth_subclass),
        inherit.aes = FALSE,
        size = 3,
        label.size = 0.15,
        label.padding = ggplot2::unit(0.15, "lines"),
        fill = "white",
        alpha = 0.85,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 1
      )
  } else {
    p_all <- p_all +
      geom_label(
        data = cell_type_label_df,
        aes(x = umap_1, y = umap_2, label = azimuth_subclass),
        inherit.aes = FALSE,
        size = 3,
        label.size = 0.15,
        fill = "white",
        alpha = 0.85
      )
  }
  
  ggsave(file.path(plot_dir, "Figure1_AllCellTypes_Azimuth_UMAP.png"),
         p_all, width = 12, height = 9, dpi = 300)
}

# ==============================================================================
# 3. Merge astrocytes only
# ==============================================================================
message("Merging dual-evidence astrocytes only...")

astro_files <- datasets$astro_path
names(astro_files) <- datasets$study_id

missing_astro <- astro_files[!file.exists(astro_files)]
if (length(missing_astro) > 0L) {
  stop("Missing astrocyte object(s):\n  ", paste(missing_astro, collapse = "\n  "))
}

seurat_astro <- NULL
for (study_id in names(astro_files)) {
  message(sprintf("  Loading astrocytes: %s", study_id))
  obj <- readRDS(astro_files[[study_id]])
  if (is.null(seurat_astro)) {
    seurat_astro <- obj
  } else {
    seurat_astro <- merge(seurat_astro, y = obj,
                          project = "Meta_Analysis_PFC_Astrocytes")
  }
  rm(obj)
  gc()
}

seurat_astro <- join_layers_if_available(seurat_astro)
message(sprintf("  Astrocytes retained after merge: %d", ncol(seurat_astro)))
print(table(seurat_astro$study))
print(table(seurat_astro$condition, useNA = "ifany"))

# ==============================================================================
# 4. Astrocyte preprocessing and Harmony integration
# ==============================================================================
message("Preprocessing merged astrocytes...")

seurat_astro <- NormalizeData(seurat_astro, verbose = FALSE)
seurat_astro <- FindVariableFeatures(seurat_astro, selection.method = "vst",
                                     nfeatures = 2000)
seurat_astro <- ScaleData(seurat_astro, features = VariableFeatures(seurat_astro),
                          verbose = FALSE)
seurat_astro <- RunPCA(seurat_astro, npcs = 30, verbose = FALSE)

message("Running Harmony batch correction...")
seurat_astro <- RunHarmony(seurat_astro, group.by.vars = "study",
                           plot_convergence = FALSE)
seurat_astro <- RunUMAP(seurat_astro, reduction = "harmony",
                        dims = 1:30, verbose = FALSE)
seurat_astro <- FindNeighbors(seurat_astro, reduction = "harmony", dims = 1:30)
seurat_astro <- FindClusters(seurat_astro, resolution = 0.4)

# ==============================================================================
# 5. Plots and outputs
# ==============================================================================
message("Saving plots and processed astrocyte objects...")

study_colors <- c(
  GSE157827 = "#E07070",
  GSE167490 = "#B08B00",
  GSE167492 = "#9A6F2D",
  GSE214979 = "#7B61B3",
  GSE263468 = "#2CA02C",
  GSE268599 = "#1F9FD5",
  GSE303823 = "#D962D9"
)

p_study <- DimPlot(
  seurat_astro,
  reduction = "umap",
  group.by = "study",
  raster = TRUE,
  alpha = 0.6,
  cols = study_colors
) + labs(title = NULL)

p_cluster <- DimPlot(
  seurat_astro,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  raster = TRUE
) + labs(title = NULL)

ggsave(file.path(plot_dir, "Figure1_Astrocyte_Integration_by_Study.png"),
       p_study, width = 10, height = 7, dpi = 300)
ggsave(file.path(plot_dir, "Figure1_Astrocyte_Subclusters.png"),
       p_cluster, width = 10, height = 7, dpi = 300)

umap_coords <- as.data.frame(Embeddings(seurat_astro, "umap"))
colnames(umap_coords) <- c("umap_1", "umap_2")
umap_coords$study <- seurat_astro$study

plots_split <- lapply(names(study_colors), function(s) {
  ggplot() +
    geom_point(
      data = umap_coords,
      aes(x = umap_1, y = umap_2),
      color = "grey88",
      size = 0.05,
      alpha = 0.3
    ) +
    geom_point(
      data = umap_coords[umap_coords$study == s, ],
      aes(x = umap_1, y = umap_2),
      color = study_colors[s],
      size = 0.1,
      alpha = 0.6
    ) +
    coord_fixed() +
    labs(title = s, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(linewidth = 0.4)
    )
})

p_split <- wrap_plots(plots_split, ncol = 3)
ggsave(file.path(plot_dir, "Figure1_Astrocyte_Integration_Split.png"),
       p_split, width = 15, height = 10, dpi = 300)

saveRDS(seurat_astro,
        file.path(output_dir, "seurat_pfc_astrocytes_preprocessed.rds"))

obj_list_final <- SplitObject(seurat_astro, split.by = "study")
for (study_id in names(obj_list_final)) {
  out_path <- file.path(output_dir,
                        sprintf("seurat_%s_preprocessed.rds",
                                tolower(study_id)))
  message(sprintf("  Saving: %s", basename(out_path)))
  saveRDS(obj_list_final[[study_id]], out_path)
}

rm(seurat_astro, obj_list_final, umap_coords, plots_split)
gc()

message("Script 02 complete.")
