# ============================================================
# FragPipeAnalystR - QC plots export (MaxLFQ protein level)
# ============================================================


suppressPackageStartupMessages({
  library(FragPipeAnalystR)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(SummarizedExperiment)
  library(ggplot2)
})


#### Adding the datasets I need to run
datasets <- list(
  Ace_Alex = list(
    combined_file = "../alexandre/resultsalex/combined_protein_Ace_Alex.tsv",
    design_file   = "../alexandre/resultsalex/experiment_annotation_Ace_Alex.tsv",
    outdir        = "FragPipeAnalystR_QC_Ace_Alex"
  ),
  
  Muriformes_Mirelle = list(
    combined_file = "../mirelle/resultsmirelle/combined_protein_Pedrosoi.tsv",
    design_file   = "../mirelle/resultsmirelle/experiment_annotation_Pedrosoi.tsv",
    outdir        = "FragPipeAnalystR_QC_Pedrosoi"
  ),
  
  HapX_Olivia = list(
    combined_file = "../celia/olivia/combined_protein_HapX_Olivia.tsv",
    design_file   = "../celia/olivia/experiment_annotation_HapX_Olivia.tsv",
    outdir        = "FragPipeAnalystR_QC_HapX"
  ),
  
  Zrt_Edu = list(
    combined_file = "../celia/edu/combined_protein_Edu_Zrt.tsv",
    design_file   = "../celia/edu/experiment_annotation_Edu_Zrt.tsv",
    outdir        = "FragPipeAnalystR_QC_Zrt_Edu"
  ),
  
  A549_Melissa = list(
    combined_file = "../celia/melissa/combined_protein_Melissa_A549.tsv",
    design_file   = "../celia/melissa/experiment_annotation_Melissa_A549.tsv",
    outdir        = "FragPipeAnalystR_QC_A549_Melissa"
  )
)


### loop para rodar todos os plots 
for (ds_name in names(datasets)) {
  
  message("Running dataset: ", ds_name)
  
  combined_file <- datasets[[ds_name]]$combined_file
  design_file   <- datasets[[ds_name]]$design_file
  outdir        <- datasets[[ds_name]]$outdir
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- A PARTIR DAQUI SEU CÓDIGO FICA IGUAL ----
  proteins <- readr::read_tsv(combined_file, show_col_types = FALSE)
  exp_design_raw <- readr::read_tsv(design_file, show_col_types = FALSE)
  
  exp_design <- exp_design_raw %>%
    dplyr::mutate(
      label     = dplyr::coalesce(.data$sample_name, .data$sample),
      condition = .data$condition,
      replicate = .data$replicate
    ) %>%
    dplyr::select(label, condition, replicate, dplyr::everything())
  
  lfq_cols <- grep(lfq_suffix_pattern, colnames(proteins), value = TRUE, perl = TRUE)
  if (length(lfq_cols) == 0) stop("No MaxLFQ columns in ", ds_name)
  
  lfq_labels <- stringr::str_replace(lfq_cols, lfq_suffix_pattern, "")
  proteins2 <- proteins
  colnames(proteins2)[match(lfq_cols, colnames(proteins2))] <- lfq_labels
  
  exp_design2 <- exp_design %>%
    dplyr::mutate(
      label = stringr::str_trim(label),
      label = stringr::str_replace(label, lfq_suffix_pattern, ""),
      label = stringr::str_replace_all(label, "\\s+", " ")
    ) %>%
    dplyr::filter(label %in% lfq_labels)
  
  if (nrow(exp_design2) == 0) {
    stop("Label mismatch in dataset: ", ds_name)
  }
  
  tmp_combined <- file.path(outdir, "combined_protein_MaxLFQ_renamed.tsv")
  tmp_design   <- file.path(outdir, "experiment_annotation_cleaned.tsv")
  readr::write_tsv(proteins2, tmp_combined)
  readr::write_tsv(exp_design2, tmp_design)
  
  se <- FragPipeAnalystR::make_se_from_files(
    tmp_combined,
    tmp_design,
    type  = "LFQ",
    level = "protein"
  )
  
  if (isTRUE(log2transform)) {
    mat <- SummarizedExperiment::assay(se)
    SummarizedExperiment::assay(se) <- log2(mat)
  }
  
  # QC
  p_feat <- FragPipeAnalystR::plot_feature_numbers(se)
  save_gg_safe(p_feat, file.path(outdir, "QC_FeatureNumbers.png"), 8, 5)
  
  ht_cor <- FragPipeAnalystR::plot_correlation_heatmap(se)
  save_complexheatmap(ht_cor, "QC_CorrelationHeatmap")
  
  ht_miss <- FragPipeAnalystR::plot_missval_heatmap(se)
  save_complexheatmap(ht_miss, "QC_MissingValueHeatmap")
  
  p_pca <- FragPipeAnalystR::plot_pca(se) +
    ggplot2::labs(title = paste("PCA —", ds_name)) +
    ggplot2::theme_classic(base_size = 12)
  
  save_gg_safe(p_pca, file.path(outdir, "QC_PCA.png"), 9, 8)
  
  message("Finished dataset: ", ds_name)
}
