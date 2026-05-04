################################################################################
# Script: 08_PFC_Visualization.R
# Description: Publication-oriented figures for the relaxed PFC astrocyte
#              meta-analysis workflow.
#
# Figures:
# 1. Donut charts: SADEG proportion among significant meta-DEGs
# 2. Volcano plot: REM meta-analysis overview with top SADEG labels
# 3. Forest plots: top SADEGs by absolute meta log2FC
# 4. Study-consistency heatmap: study-level log2FC for top SADEGs
# 5. Pseudobulk hierarchical heatmap: donor-level astrocyte expression
# 6. BRETIGEA astrocyte marker expression UMAP feature plots
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(BRETIGEA)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(metafor)
  library(pheatmap)
})

# ==============================================================================
# 0. Setup
# ==============================================================================

step2_dir <- Sys.getenv(
  "PFC_STEP2_DIR",
  unset = file.path("~", "data", "meta_analysis", "results", "PFC", "2026.04")
)
step2_dir <- normalizePath(path.expand(step2_dir), winslash = "/", mustWork = FALSE)

meta_dir <- file.path(step2_dir, "Meta_Analysis_REM")
sadeg_dir <- file.path(step2_dir, "SADEG_Overlap_Results")
prerank_dir <- file.path(step2_dir, "PreRanked_gProfiler_EnrichmentMap")
astro_dir <- file.path(step2_dir, "step2_astrocyte_objects")
annotated_dir <- file.path(step2_dir, "step2_azimuth_annotated")
output_dir <- file.path(step2_dir, "Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Heatmap genes are all SADEGs ordered by signed REM meta-analysis log2FC.
# Rows are not hierarchically clustered so that upregulated SADEGs appear at the
# top and downregulated SADEGs appear at the bottom. Columns are still clustered.
HEATMAP_RANK_METRIC <- "signed_log2FC"

# BRETIGEA marker feature plots. Step 2 used the top 50 markers per cell type
# from BRETIGEA::markers_df_human_brain for AddModuleScore().
N_BRETIGEA_MARKERS_PER_CELLTYPE <- 50
N_ASTRO_MARKERS_FOR_FEATURE_PLOT <- 6

datasets <- tibble::tribble(
  ~study,      ~dataset_id,  ~label,            ~astro_file,
  "GSE157827", "gse157827",  "Lau",             "gse157827_astrocytes_step2.rds",
  "GSE167490", "gse167490",  "Sadick main",     "gse167490_astrocytes_step2.rds",
  "GSE167492", "gse167492",  "Sadick pilot",    "gse167492_astrocytes_step2.rds",
  "GSE214979", "gse214979",  "Anderson",        "gse214979_astrocytes_step2.rds",
  "GSE263468", "gse263468",  "Cobos",           "gse263468_astrocytes_step2.rds",
  "GSE268599", "gse268599",  "Serrano-Pozo",    "gse268599_astrocytes_step2.rds",
  "GSE303823", "gse303823",  "Chia",            "gse303823_astrocytes_step2.rds"
) %>%
  mutate(
    astro_path = file.path(astro_dir, astro_file),
    annotated_path = file.path(annotated_dir, paste0(dataset_id, "_azimuth_bretigea.rds"))
  )

# ==============================================================================
# 1. Helpers
# ==============================================================================

load_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop("[", label, "] File not found:\n  ", path)
  }
  read_csv(path, show_col_types = FALSE)
}

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

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

choose_pseudobulk_variable <- function(meta) {
  for (candidate in c("donor_id", "sample_id", "orig.ident")) {
    if (candidate %in% colnames(meta)) {
      x <- as.character(meta[[candidate]])
      if (!all(is.na(x)) && length(unique(na.omit(x))) > 1L) {
        return(candidate)
      }
    }
  }
  NA_character_
}

get_counts_for_genes <- function(obj, genes, assay = "RNA", study = NULL) {
  genes <- unique(as.character(genes))
  assay_obj <- obj[[assay]]
  layers <- if ("Layers" %in% getNamespaceExports("SeuratObject")) {
    SeuratObject::Layers(assay_obj)
  } else {
    character(0)
  }
  
  read_layer <- function(layer_name) {
    mat <- SeuratObject::LayerData(obj, assay = assay, layer = layer_name)
    available <- intersect(genes, rownames(mat))
    mat <- mat[available, , drop = FALSE]
    missing_genes <- setdiff(genes, rownames(mat))
    if (length(missing_genes) > 0L) {
      zero_mat <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0),
        dims = c(length(missing_genes), ncol(mat)),
        dimnames = list(missing_genes, colnames(mat))
      )
      mat <- rbind(mat, zero_mat)
    }
    mat[genes, , drop = FALSE]
  }
  
  if ("counts" %in% layers) {
    mat <- read_layer("counts")
    message(sprintf(
      "  [%s] Heatmap counts layer: %d genes x %d cells",
      study, nrow(mat), ncol(mat)
    ))
    return(mat)
  }
  
  counts_layers <- grep("^counts", layers, value = TRUE)
  if (length(counts_layers) == 0L) {
    stop("[", study, "] No counts layer found for heatmap.")
  }
  
  mats <- lapply(counts_layers, read_layer)
  all_cells <- unlist(lapply(mats, colnames), use.names = FALSE)
  counts <- do.call(cbind, mats)
  counts <- counts[, all_cells, drop = FALSE]
  message(sprintf(
    "  [%s] Combined heatmap split counts layers: %d genes x %d cells",
    study, nrow(counts), ncol(counts)
  ))
  counts
}

get_expression_for_feature_plot <- function(obj, genes, assay = "RNA") {
  DefaultAssay(obj) <- assay
  assay_obj <- obj[[assay]]
  layers <- if ("Layers" %in% getNamespaceExports("SeuratObject")) {
    SeuratObject::Layers(assay_obj)
  } else {
    character(0)
  }
  
  data_layers <- grep("^data", layers, value = TRUE)
  if (length(data_layers) == 0L) {
    message("  No RNA data layer found for marker feature plots; running NormalizeData().")
    obj <- NormalizeData(obj, verbose = FALSE)
    assay_obj <- obj[[assay]]
    layers <- SeuratObject::Layers(assay_obj)
    data_layers <- grep("^data", layers, value = TRUE)
  }
  
  layer_type <- "data"
  target_layers <- data_layers
  if (length(target_layers) == 0L) {
    target_layers <- grep("^counts", layers, value = TRUE)
    layer_type <- "counts"
  }
  if (length(target_layers) == 0L) {
    stop("No data or counts layer found for marker feature plots.")
  }
  
  read_feature_layer <- function(layer_name) {
    mat <- SeuratObject::LayerData(obj, assay = assay, layer = layer_name)
    present <- intersect(genes, rownames(mat))
    mat <- mat[present, , drop = FALSE]
    missing_genes <- setdiff(genes, rownames(mat))
    if (length(missing_genes) > 0L) {
      zero_mat <- Matrix::sparseMatrix(
        i = integer(0), j = integer(0),
        dims = c(length(missing_genes), ncol(mat)),
        dimnames = list(missing_genes, colnames(mat))
      )
      mat <- rbind(mat, zero_mat)
    }
    mat[genes, , drop = FALSE]
  }
  
  expr <- do.call(cbind, lapply(target_layers, read_feature_layer))
  expr <- expr[, colnames(obj), drop = FALSE]
  attr(expr, "layer_type") <- layer_type
  expr
}

make_marker_feature_plots <- function(obj, marker_genes, output_dir,
                                      file_prefix,
                                      ncol = 3,
                                      point_size = 0.05) {
  if (!"umap" %in% Reductions(obj)) {
    warning("Skipping marker feature plots; the object has no UMAP reduction.")
    return(invisible(NULL))
  }
  
  present_genes <- marker_genes[marker_genes %in% rownames(obj)]
  if (length(present_genes) == 0L) {
    warning("Skipping marker feature plots; none of the requested markers are present.")
    return(invisible(NULL))
  }
  
  expr <- get_expression_for_feature_plot(obj, present_genes, assay = "RNA")
  layer_type <- attr(expr, "layer_type")
  
  umap_df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap_df)[seq_len(2)] <- c("umap_1", "umap_2")
  umap_df$cell <- rownames(umap_df)
  umap_df$cluster <- if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    as.character(obj$seurat_clusters)
  } else {
    NA_character_
  }
  
  cluster_label_df <- umap_df %>%
    filter(!is.na(cluster)) %>%
    group_by(cluster) %>%
    summarize(
      umap_1 = median(umap_1, na.rm = TRUE),
      umap_2 = median(umap_2, na.rm = TRUE),
      .groups = "drop"
    )
  
  plot_one_gene <- function(gene) {
    plot_df <- umap_df
    plot_df$expression <- as.numeric(expr[gene, plot_df$cell])
    upper <- stats::quantile(plot_df$expression, probs = 0.99, na.rm = TRUE)
    if (!is.finite(upper) || upper <= 0) {
      upper <- max(plot_df$expression, na.rm = TRUE)
    }
    
    ggplot(plot_df, aes(x = umap_1, y = umap_2, color = expression)) +
      geom_point(size = point_size, alpha = 0.75) +
      scale_color_gradientn(
        colors = c("grey85", "#BDA7FF", "#2D00FF"),
        limits = c(0, upper),
        oob = scales::squish,
        name = layer_type
      ) +
      coord_fixed() +
      labs(title = gene, x = "UMAP_1", y = "UMAP_2") +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      ) +
      geom_text(
        data = cluster_label_df,
        aes(x = umap_1, y = umap_2, label = cluster),
        inherit.aes = FALSE,
        size = 2.4,
        color = "black"
      )
  }
  
  plots <- lapply(present_genes, plot_one_gene)
  p <- patchwork::wrap_plots(plots, ncol = ncol)
  height <- max(4, ceiling(length(plots) / ncol) * 3.1)
  width <- min(16, ncol * 3.8)
  
  ggsave(
    file.path(output_dir, paste0(file_prefix, ".png")),
    p,
    width = width,
    height = height,
    dpi = 300
  )
  
  invisible(p)
}

make_all_cell_marker_feature_plots <- function(plot_df, marker_genes, output_dir,
                                               file_prefix,
                                               ncol = 3,
                                               point_size = 0.04) {
  marker_genes <- marker_genes[marker_genes %in% colnames(plot_df)]
  if (length(marker_genes) == 0L) {
    warning("Skipping all-cell marker feature plots; no marker expression columns found.")
    return(invisible(NULL))
  }
  
  cell_type_label_df <- plot_df %>%
    group_by(azimuth_subclass) %>%
    summarize(
      umap_1 = median(umap_1, na.rm = TRUE),
      umap_2 = median(umap_2, na.rm = TRUE),
      .groups = "drop"
    )
  
  plot_one_gene <- function(gene) {
    gene_df <- plot_df %>%
      select(umap_1, umap_2, azimuth_subclass, expression = all_of(gene))
    upper <- stats::quantile(gene_df$expression, probs = 0.99, na.rm = TRUE)
    if (!is.finite(upper) || upper <= 0) {
      upper <- max(gene_df$expression, na.rm = TRUE)
    }
    
    ggplot(gene_df, aes(x = umap_1, y = umap_2, color = expression)) +
      geom_point(size = point_size, alpha = 0.65) +
      scale_color_gradientn(
        colors = c("grey88", "#BDA7FF", "#2D00FF"),
        limits = c(0, upper),
        oob = scales::squish,
        name = "data"
      ) +
      coord_fixed() +
      labs(title = gene, x = "UMAP_1", y = "UMAP_2") +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      ) +
      ggrepel::geom_label_repel(
        data = cell_type_label_df,
        aes(x = umap_1, y = umap_2, label = azimuth_subclass),
        inherit.aes = FALSE,
        size = 2.4,
        label.size = 0.12,
        label.padding = ggplot2::unit(0.12, "lines"),
        fill = "white",
        alpha = 0.75,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 1
      )
  }
  
  plots <- lapply(marker_genes, plot_one_gene)
  p <- patchwork::wrap_plots(plots, ncol = ncol)
  height <- max(4, ceiling(length(plots) / ncol) * 3.5)
  width <- ncol * 4.4
  
  ggsave(
    file.path(output_dir, paste0(file_prefix, ".png")),
    p,
    width = width,
    height = height,
    dpi = 300
  )
  
  invisible(p)
}

make_pseudobulk_counts_for_heatmap <- function(obj, genes, study) {
  DefaultAssay(obj) <- "RNA"
  counts <- get_counts_for_genes(obj, genes = genes, assay = "RNA", study = study)
  meta <- obj@meta.data
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  counts <- counts[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  
  pb_var <- choose_pseudobulk_variable(meta)
  if (is.na(pb_var)) {
    stop("[", study, "] No usable donor/sample ID column for heatmap.")
  }
  
  meta$pb_source_id <- as.character(meta[[pb_var]])
  meta$pb_id <- safe_id(meta$pb_source_id)
  meta$condition <- as.character(meta$condition)
  
  keep_cells <- !is.na(meta$pb_id) & meta$pb_id != "S." &
    meta$condition %in% c("AD", "Control")
  counts <- counts[, keep_cells, drop = FALSE]
  meta <- meta[keep_cells, , drop = FALSE]
  
  pb_meta <- meta %>%
    mutate(
      age = if ("age" %in% colnames(meta)) suppressWarnings(as.numeric(age)) else NA_real_,
      sex = if ("sex" %in% colnames(meta)) as.character(sex) else NA_character_
    ) %>%
    group_by(pb_id, pb_source_id, condition) %>%
    summarize(
      age = first_or_na(age),
      sex = first_or_na(sex),
      n_astrocytes = n(),
      .groups = "drop"
    ) %>%
    mutate(
      study = study,
      sample_label = paste(study, pb_source_id, sep = "__")
    ) %>%
    arrange(sample_label)
  
  group <- factor(meta$pb_id, levels = pb_meta$pb_id)
  design <- Matrix::sparse.model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  pb_counts <- counts %*% design
  pb_counts <- as.matrix(pb_counts)
  colnames(pb_counts) <- paste(study, pb_meta$pb_source_id, sep = "__")
  pb_counts <- pb_counts[, pb_meta$sample_label, drop = FALSE]
  
  list(counts = pb_counts, metadata = pb_meta)
}

generate_donut_chart <- function(overlap_count, total_count, title,
                                 color_sadeg, color_other) {
  other_count <- total_count - overlap_count
  df <- tibble(
    Category = factor(c("SADEGs", "Non-SADEGs"),
                      levels = c("Non-SADEGs", "SADEGs")),
    Count = c(overlap_count, other_count)
  ) %>%
    mutate(
      Fraction = Count / sum(Count),
      ymax = cumsum(Fraction),
      ymin = c(0, head(ymax, n = -1)),
      LabelPosition = (ymax + ymin) / 2,
      Label = sprintf("%d\n(%.1f%%)", Count, Fraction * 100)
    )
  
  ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Category)) +
    geom_rect(color = "white", linewidth = 0.5) +
    geom_text(aes(x = 3.5, y = LabelPosition, label = Label),
              size = 4.5, color = "white", fontface = "bold") +
    scale_fill_manual(values = c("SADEGs" = color_sadeg,
                                 "Non-SADEGs" = color_other)) +
    coord_polar(theta = "y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
    ) +
    ggtitle(title)
}

# ==============================================================================
# 2. Load Data
# ==============================================================================

message("Loading visualization inputs...")

meta_res <- load_required_csv(
  file.path(meta_dir, "PFC_metaDEGs_All_Results.csv"),
  "REM meta-analysis all results"
)
meta_input_long <- load_required_csv(
  file.path(meta_dir, "PFC_meta_input_long.csv"),
  "REM meta-analysis long input"
)
meta_sig <- load_required_csv(
  file.path(meta_dir, "PFC_metaDEGs_Significant.csv"),
  "Significant meta-DEGs"
)
sadeg_all <- load_required_csv(
  file.path(sadeg_dir, "PFC_SADEGs_All.csv"),
  "SADEGs"
)
sadeg_up <- load_required_csv(
  file.path(sadeg_dir, "PFC_SADEGs_Upregulated.csv"),
  "Upregulated SADEGs"
)
sadeg_down <- load_required_csv(
  file.path(sadeg_dir, "PFC_SADEGs_Downregulated.csv"),
  "Downregulated SADEGs"
)

rnk_path <- file.path(prerank_dir, "PFC_MetaDEGs_SAGs_Overlap_Signed_log2FC.rnk")
if (file.exists(rnk_path)) {
  rnk_df <- read.table(rnk_path, header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE)
  colnames(rnk_df) <- c("gene", "signed_log2fc")
} else {
  rnk_df <- sadeg_all %>%
    transmute(gene = gene, signed_log2fc = beta_meta)
}

rnk_df <- rnk_df %>%
  mutate(gene = as.character(gene),
         signed_log2fc = as.numeric(signed_log2fc),
         abs_log2fc = abs(signed_log2fc)) %>%
  filter(!is.na(gene), gene != "", is.finite(signed_log2fc)) %>%
  distinct(gene, .keep_all = TRUE)

top_label_genes <- bind_rows(
  rnk_df %>% filter(signed_log2fc > 0) %>%
    arrange(desc(signed_log2fc)) %>% head(5),
  rnk_df %>% filter(signed_log2fc < 0) %>%
    arrange(signed_log2fc) %>% head(5)
) %>%
  pull(gene) %>%
  unique()

heatmap_gene_table <- rnk_df %>%
  arrange(desc(signed_log2fc)) %>%
  mutate(
    heatmap_rank = row_number(),
    selection_rule = paste0(
      "All SADEGs ordered by ", HEATMAP_RANK_METRIC,
      " from most positive to most negative"
    )
  ) %>%
  select(heatmap_rank, gene, signed_log2fc, abs_log2fc, selection_rule)

heatmap_genes <- heatmap_gene_table$gene

write_csv(
  heatmap_gene_table,
  file.path(output_dir, "Figure5_Selected_SADEGs_for_Heatmap_by_signed_log2FC.csv")
)

message(sprintf(
  "Selected %d SADEGs for heatmaps using %s.",
  length(heatmap_genes), HEATMAP_RANK_METRIC
))

# ==============================================================================
# 3. Figure 1: Donut Charts
# ==============================================================================

message("Drawing donut charts...")

n_sig <- nrow(meta_sig)
n_sig_up <- sum(meta_sig$beta_meta > 0, na.rm = TRUE)
n_sig_down <- sum(meta_sig$beta_meta < 0, na.rm = TRUE)

png(file.path(output_dir, "Figure1A_Donut_All_MetaDEGs.png"),
    width = 1800, height = 1500, res = 300)
print(generate_donut_chart(nrow(sadeg_all), n_sig, "All significant meta-DEGs",
                           "#EFC000FF", "#868686FF"))
dev.off()

png(file.path(output_dir, "Figure1B_Donut_Up_MetaDEGs.png"),
    width = 1800, height = 1500, res = 300)
print(generate_donut_chart(nrow(sadeg_up), n_sig_up, "Upregulated meta-DEGs",
                           "#CD534CFF", "#868686FF"))
dev.off()

png(file.path(output_dir, "Figure1C_Donut_Down_MetaDEGs.png"),
    width = 1800, height = 1500, res = 300)
print(generate_donut_chart(nrow(sadeg_down), n_sig_down, "Downregulated meta-DEGs",
                           "#0073C2FF", "#868686FF"))
dev.off()

# ==============================================================================
# 4. Figure 2: Volcano Plot
# ==============================================================================

message("Drawing volcano plot...")

volcano_data <- meta_res %>%
  mutate(
    p_safe = ifelse(is.na(p_meta) | p_meta <= 0, .Machine$double.xmin, p_meta),
    regulation = case_when(
      padj_meta < 0.05 & beta_meta > 0 ~ "Up",
      padj_meta < 0.05 & beta_meta < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    is_sadeg = gene %in% sadeg_all$gene
  )

if (any(volcano_data$padj_meta < 0.05, na.rm = TRUE)) {
  fdr_threshold <- max(volcano_data$p_safe[volcano_data$padj_meta < 0.05],
                       na.rm = TRUE)
  y_cutoff <- -log10(fdr_threshold)
} else {
  y_cutoff <- NA_real_
}

label_data <- volcano_data %>%
  filter(gene %in% top_label_genes)

p_volcano <- ggplot(volcano_data,
                    aes(x = beta_meta, y = -log10(p_safe), color = regulation)) +
  geom_point(alpha = 0.45, size = 1.1) +
  geom_point(data = subset(volcano_data, is_sadeg),
             aes(x = beta_meta, y = -log10(p_safe)),
             inherit.aes = FALSE, shape = 21, fill = NA,
             color = "black", stroke = 0.25, size = 1.5) +
  scale_color_manual(values = c("Down" = "#377EB8",
                                "NS" = "grey80",
                                "Up" = "#E41A1C")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Meta-analysis log2 fold change",
       y = "-log10(P-value)",
       color = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

if (is.finite(y_cutoff)) {
  p_volcano <- p_volcano +
    geom_hline(yintercept = y_cutoff, linetype = "dashed", color = "black")
}
if (nrow(label_data) > 0L) {
  p_volcano <- p_volcano +
    geom_text_repel(
      data = label_data,
      aes(label = gene),
      color = "black",
      size = 4,
      fontface = "bold.italic",
      box.padding = 0.8,
      point.padding = 0.4,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 42
    )
}

ggsave(file.path(output_dir, "Figure2_Volcano_MetaAnalysis_SADEGs.png"),
       p_volcano, width = 8, height = 7, dpi = 300)

# ==============================================================================
# 5. Figure 3: Forest Plots
# ==============================================================================

message("Drawing forest plots...")

forest_genes <- sadeg_all %>%
  mutate(abs_beta = abs(beta_meta)) %>%
  arrange(desc(abs_beta)) %>%
  head(3) %>%
  pull(gene)

draw_forest_plot <- function(gene_name) {
  g_data <- meta_input_long %>%
    filter(gene == gene_name) %>%
    arrange(study)
  if (nrow(g_data) < 2L) {
    warning("Skipping forest plot; fewer than 2 studies for gene: ", gene_name)
    return(invisible(NULL))
  }
  
  fit <- tryCatch(
    metafor::rma(yi = logFC, sei = SE, data = g_data, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    warning("Forest model failed for gene: ", gene_name)
    return(invisible(NULL))
  }
  
  labels <- datasets$label[match(g_data$study, datasets$study)]
  labels <- ifelse(is.na(labels), g_data$study, labels)
  
  png(file.path(output_dir,
                paste0("Figure3_Forest_", safe_filename(gene_name), ".png")),
      width = 2400, height = 1500, res = 300)
  metafor::forest(
    fit,
    slab = labels,
    header = "Study",
    xlab = "Log2 fold change",
    col = "royalblue",
    border = "royalblue",
    cex = 0.9
  )
  text(min(fit$ci.lb, na.rm = TRUE), -1.5,
       bquote(paste("Heterogeneity: ", I^2, " = ", .(round(fit$I2, 1)), "%")),
       pos = 4, cex = 0.8)
  dev.off()
}

for (g in forest_genes) {
  draw_forest_plot(g)
}

# ==============================================================================
# 6. Figure 4: Study-Level Consistency Heatmap
# ==============================================================================

message("Drawing study-level consistency heatmap...")

study_labels <- datasets %>%
  transmute(study, display = paste0(study, "\n(", label, ")"))

consistency_mat <- meta_input_long %>%
  filter(gene %in% heatmap_genes) %>%
  select(gene, study, logFC) %>%
  tidyr::pivot_wider(names_from = study, values_from = logFC) %>%
  column_to_rownames("gene") %>%
  as.matrix()

available_studies <- intersect(study_labels$study, colnames(consistency_mat))
consistency_mat <- consistency_mat[heatmap_genes[heatmap_genes %in% rownames(consistency_mat)],
                                   available_studies, drop = FALSE]
colnames(consistency_mat) <- study_labels$display[
  match(colnames(consistency_mat), study_labels$study)
]

if (nrow(consistency_mat) >= 2L && ncol(consistency_mat) >= 2L) {
  png(file.path(output_dir, "Figure4_StudyLevel_LogFC_Consistency_Heatmap.png"),
      width = 2200, height = 3000, res = 300)
  pheatmap(
    consistency_mat,
    color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(50),
    breaks = seq(-1.5, 1.5, length.out = 51),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    display_numbers = FALSE,
    fontsize_row = ifelse(nrow(consistency_mat) > 45, 6, 8),
    fontsize_col = 9,
    angle_col = 45,
    main = ""
  )
  dev.off()
} else {
  warning("Skipping study-level consistency heatmap; matrix is too small.")
}

# ==============================================================================
# 7. Figure 5: Donor-Level Pseudobulk Hierarchical Heatmap
# ==============================================================================

message("Drawing donor-level pseudobulk hierarchical heatmap...")

missing_astro <- datasets %>% filter(!file.exists(astro_path))
if (nrow(missing_astro) > 0L) {
  warning("Skipping pseudobulk heatmap; missing astrocyte object(s):\n  ",
          paste(missing_astro$astro_path, collapse = "\n  "))
} else if (length(heatmap_genes) < 2L) {
  warning("Skipping pseudobulk heatmap; fewer than 2 SADEG heatmap genes.")
} else {
  pb_counts_list <- list()
  pb_meta_list <- list()
  
  for (i in seq_len(nrow(datasets))) {
    study_id <- datasets$study[i]
    message(sprintf("  Loading astrocyte object for heatmap: %s", study_id))
    obj <- readRDS(datasets$astro_path[i])
    pb <- make_pseudobulk_counts_for_heatmap(obj, heatmap_genes, study_id)
    pb_counts_list[[study_id]] <- pb$counts
    pb_meta_list[[study_id]] <- pb$metadata
    rm(obj, pb)
    gc()
  }
  
  all_genes <- Reduce(union, lapply(pb_counts_list, rownames))
  pb_counts_list <- lapply(pb_counts_list, function(mat) {
    missing_genes <- setdiff(all_genes, rownames(mat))
    if (length(missing_genes) > 0L) {
      zero_mat <- matrix(
        0,
        nrow = length(missing_genes),
        ncol = ncol(mat),
        dimnames = list(missing_genes, colnames(mat))
      )
      mat <- rbind(mat, zero_mat)
    }
    mat[all_genes, , drop = FALSE]
  })
  
  pb_counts <- do.call(cbind, pb_counts_list)
  pb_meta <- bind_rows(pb_meta_list) %>%
    distinct(sample_label, .keep_all = TRUE)
  rownames(pb_meta) <- pb_meta$sample_label
  pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
  
  y <- edgeR::DGEList(counts = pb_counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  
  keep_sd <- apply(logcpm, 1, sd, na.rm = TRUE) > 0
  logcpm <- logcpm[keep_sd, , drop = FALSE]
  if (nrow(logcpm) < 2L || ncol(logcpm) < 2L) {
    warning("Skipping pseudobulk heatmap; expression matrix is too small after filtering.")
  } else {
    
    z_mat <- t(scale(t(logcpm)))
    z_mat[is.na(z_mat)] <- 0
    z_mat[z_mat > 2.5] <- 2.5
    z_mat[z_mat < -2.5] <- -2.5
    
    annotation_col <- pb_meta %>%
      transmute(
        Condition = factor(condition, levels = c("Control", "AD")),
        Study = factor(study),
        Sex = ifelse(is.na(sex) | sex == "", "Not reported", sex),
        Age = age,
        Astrocytes = n_astrocytes
      ) %>%
      as.data.frame()
    rownames(annotation_col) <- pb_meta$sample_label
    
    png(file.path(output_dir, "Figure5_SADEG_Pseudobulk_Hierarchical_Heatmap.png"),
        width = 3200, height = 2800, res = 300)
    pheatmap(
      z_mat,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
      breaks = seq(-2.5, 2.5, length.out = 102),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      annotation_col = annotation_col,
      show_colnames = FALSE,
      show_rownames = TRUE,
      fontsize_row = ifelse(nrow(z_mat) > 45, 6, 8),
      main = ""
    )
    dev.off()
    
    write_csv(
      pb_meta,
      file.path(output_dir, "Figure5_Pseudobulk_Heatmap_Sample_Metadata.csv")
    )
    write_csv(
      as.data.frame(z_mat) %>% rownames_to_column("gene"),
      file.path(output_dir, "Figure5_Pseudobulk_Heatmap_RowZScore_Matrix.csv")
    )
  }
}

# ==============================================================================
# 8. BRETIGEA Astrocyte Marker Feature Plots
# ==============================================================================

message("Exporting BRETIGEA markers and plotting astrocyte marker expression...")

data("markers_df_human_brain", package = "BRETIGEA")
required_marker_cols <- c("cell", "markers")
if (!all(required_marker_cols %in% colnames(markers_df_human_brain))) {
  stop(
    "BRETIGEA::markers_df_human_brain must contain columns: ",
    paste(required_marker_cols, collapse = ", ")
  )
}

bretigea_marker_table <- markers_df_human_brain %>%
  group_by(cell) %>%
  mutate(
    marker_rank = row_number(),
    used_for_step2_scoring = marker_rank <= N_BRETIGEA_MARKERS_PER_CELLTYPE
  ) %>%
  ungroup() %>%
  select(cell, marker_rank, marker = markers, used_for_step2_scoring)

bretigea_top50_markers <- bretigea_marker_table %>%
  filter(used_for_step2_scoring)

bretigea_astro_top50 <- bretigea_top50_markers %>%
  filter(cell == "ast")

if (nrow(bretigea_astro_top50) == 0L) {
  stop("BRETIGEA::markers_df_human_brain does not contain astrocyte cell type 'ast'.")
}

write_csv(
  bretigea_top50_markers,
  file.path(output_dir, "BRETIGEA_human_brain_top50_markers_used_for_step2_scoring.csv")
)
write_csv(
  bretigea_astro_top50,
  file.path(output_dir, "BRETIGEA_astrocyte_top50_markers_used_for_step2_scoring.csv")
)

all_cell_umap_path <- file.path(step2_dir, "all_cells_azimuth_umap_coordinates.csv")
missing_annotated <- datasets %>% filter(!file.exists(annotated_path))

if (!file.exists(all_cell_umap_path)) {
  warning(
    "Skipping all-cell BRETIGEA astrocyte marker feature plots; file not found:\n  ",
    all_cell_umap_path
  )
} else if (nrow(missing_annotated) > 0L) {
  warning(
    "Skipping all-cell BRETIGEA astrocyte marker feature plots; missing annotated object(s):\n  ",
    paste(missing_annotated$annotated_path, collapse = "\n  ")
  )
} else {
  marker_candidates <- bretigea_astro_top50$marker
  
  marker_presence <- lapply(seq_len(nrow(datasets)), function(i) {
    study_id <- datasets$study[i]
    message(sprintf("  Checking BRETIGEA astrocyte marker availability: %s", study_id))
    obj <- readRDS(datasets$annotated_path[i])
    tibble(
      study = study_id,
      marker = marker_candidates,
      present = marker_candidates %in% rownames(obj)
    )
  })
  marker_presence <- bind_rows(marker_presence)
  
  astro_marker_coverage <- bretigea_astro_top50 %>%
    left_join(
      marker_presence %>%
        group_by(marker) %>%
        summarize(
          n_studies_present = sum(present),
          present_in_all_cell_objects = any(present),
          .groups = "drop"
        ),
      by = "marker"
    ) %>%
    mutate(
      n_studies_present = ifelse(is.na(n_studies_present), 0L, n_studies_present),
      present_in_all_cell_objects = ifelse(
        is.na(present_in_all_cell_objects),
        FALSE,
        present_in_all_cell_objects
      ),
      selected_for_feature_plot = FALSE
    )
  
  feature_genes <- astro_marker_coverage %>%
    filter(present_in_all_cell_objects) %>%
    arrange(marker_rank) %>%
    slice_head(n = N_ASTRO_MARKERS_FOR_FEATURE_PLOT) %>%
    pull(marker)
  
  astro_marker_coverage <- astro_marker_coverage %>%
    mutate(selected_for_feature_plot = marker %in% feature_genes)
  
  write_csv(
    astro_marker_coverage,
    file.path(output_dir, "BRETIGEA_astrocyte_top50_marker_feature_plot_coverage.csv")
  )
  
  if (length(feature_genes) < N_ASTRO_MARKERS_FOR_FEATURE_PLOT) {
    warning(sprintf(
      "Only %d BRETIGEA astrocyte markers are present in the all-cell objects.",
      length(feature_genes)
    ))
  }
  
  if (length(feature_genes) == 0L) {
    warning("No BRETIGEA astrocyte markers selected for all-cell feature plotting.")
  } else {
    message(sprintf(
      "  Plotting %d BRETIGEA astrocyte markers on all-cell Azimuth UMAP: %s",
      length(feature_genes),
      paste(feature_genes, collapse = ", ")
    ))
    
    all_cell_umap <- read_csv(all_cell_umap_path, show_col_types = FALSE) %>%
      mutate(
        study = as.character(study),
        cell_barcode = as.character(cell_barcode)
      )
    
    expression_list <- lapply(seq_len(nrow(datasets)), function(i) {
      study_id <- datasets$study[i]
      message(sprintf("  Extracting all-cell marker expression: %s", study_id))
      obj <- readRDS(datasets$annotated_path[i])
      expr <- get_expression_for_feature_plot(obj, feature_genes, assay = "RNA")
      expr_df <- as.data.frame(t(as.matrix(expr))) %>%
        rownames_to_column("cell_barcode") %>%
        mutate(study = study_id)
      rm(obj, expr)
      gc()
      expr_df
    })
    
    marker_expression <- bind_rows(expression_list)
    all_cell_marker_df <- all_cell_umap %>%
      left_join(marker_expression, by = c("study", "cell_barcode")) %>%
      mutate(across(all_of(feature_genes), ~ tidyr::replace_na(.x, 0)))
    
    write_csv(
      all_cell_marker_df %>%
        select(study, dataset_id, cell_barcode, azimuth_subclass, all_of(feature_genes)),
      file.path(output_dir, "Figure6_AllCell_BRETIGEA_Astrocyte_Top6Marker_Expression.csv")
    )
    
    make_all_cell_marker_feature_plots(
      plot_df = all_cell_marker_df,
      marker_genes = feature_genes,
      output_dir = output_dir,
      file_prefix = sprintf(
        "Figure6_AllCell_AzimuthUMAP_BRETIGEA_Astrocyte_Top%dMarker_FeaturePlots",
        length(feature_genes)
      ),
      ncol = 3,
      point_size = 0.035
    )
    
    rm(all_cell_umap, marker_expression, all_cell_marker_df, expression_list)
    gc()
  }
}

message("Step 8 visualization complete.")
