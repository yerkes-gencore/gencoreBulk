#' Normalize an pre-processed expression matrix for heatmaps
#' 
#' Normalize an RNA expression matrix by samples or groups of samples within the
#'  matrix. The intended use case is to prepare data for plotting with heatmaps.
#'  You define a grouping variable and a baseline level of that variable to normalize
#'  other groups by. If there are more than one sample in the baseline group, the
#'  median expression value is used to normalize other values.
#'  
#'  This is NOT intended to normalized data as part of the upstream processing. 
#'  Data fed to this function should already be log normalized, as the baseline
#'  normalization assumes the expression is on a log scale.
#'
#' @param count_matrix Matrix of processed, log normalized counts data
#' @param metadata Metadata for samples in counts matrix
#' @param group_var Column of metadata to group samples by
#' @param baseline Level of `group_var` to use as baseline for normalization
#' @param remove_baseline Bool, remove baseline samples from output
#'
#' @importFrom matrixStats rowMedians
#' 
#' @returns A numeric matrix
#' @export
#'
#' @examples
#' \dontrun{
#' tmp <- normalizeCounts(assay(assays(obj.pbmc)$rld),
#'         obj.pbmc@colData,
#'         group_var = 'Timepoint', baseline = 'pre', 
#'         remove_baseline = TRUE)
#' heatmapFromGenelist(genes, tmp)
#' }
normalizeCountsForHeatmap <- function(count_matrix, 
                                      metadata, 
                                      group_var, 
                                      baseline, 
                                      remove_baseline=FALSE) {
  if (nrow(metadata) != ncol(count_matrix)) {
    stop("Number of rows in metadata and count_matrix must match.")
  }
  if (!(group_var %in% colnames(metadata))) {
    stop("group_var, individual_var, and baseline must be valid column names in metadata.")
  }
  if (!(baseline %in% metadata[[group_var]])) {
    stop('Baseline should be a level of group_var')
  }
  baseline_indices <- metadata[[group_var]]==baseline
  baseline_counts <- count_matrix[,baseline_indices]
  if (length(baseline_counts) == 0) {
    stop("Baseline not found.")
  }
  if (!is.null(dim(baseline_counts))) {
    message('More than 1 sample in baseline, using median')
    baseline_counts <- matrixStats::rowMedians(baseline_counts)
  }
  groups <- unique(metadata[[group_var]])
  normalized_counts <- count_matrix
  for (group in groups) {
    group_indices <- metadata[[group_var]] == group
    group_metadata <- metadata[group_indices, ]
    group_counts <- count_matrix[,group_indices]
    normalized_counts[,group_indices] <- group_counts - baseline_counts
  }
  if (remove_baseline) {
    normalized_counts <- normalized_counts[,!baseline_indices]
  }
  return(normalized_counts)
}

#' Normalize an pre-processed expression matrix for heatmaps within individuals
#' 
#' Extension of `gencoreBulk::normalizeCounts` to normalize data within individual.
#' Samples are grouped by individual, then a baseline is defined on a per-individual
#' level. Samples for an individual are normalized by that individual's baseline.
#' If there are more than one sample in the individual's baseline, the median expression
#'  value is used as the baseline (and you have a well designed study). 
#' 
#' This is intended for use with time-course studies with repeated measures from
#' the same individual. Rather than normalizing by the group-level statistics,
#' you can normalize within individual to account for different baseline levels 
#' of expression.
#' 
#' This is NOT intended to normalized data as part of the upstream processing. 
#'  Data fed to this function should already be log normalized, as the baseline
#'  normalization assumes the expression is on a log scale.
#'
#' @param count_matrix Matrix of processed, log normalized counts data
#' @param metadata Metadata for samples in counts matrix
#' @param group_var Column of metadata to group samples by
#' @param baseline Level of `group_var` to use as baseline for normalization
#' @param individual_var Column of metadata to group individuals by
#' @param remove_baseline Bool, remove baseline samples from output
#'
#' @returns A numeric matrix
#' @export
#'
#' @examples
#' \dontrun{
#' tmp <- normalizeCountsByIndividual(assay(assays(obj.pbmc)$rld),
#'         obj.pbmc@colData,
#'         group_var = 'Timepoint', baseline = 'pre', 
#'         individual_var = 'Individual',
#'         remove_baseline = TRUE)
#' heatmapFromGenelist(genes, tmp)
#' }
normalizeCountsForHeatmapByIndividual <- function(count_matrix, 
                                                  metadata, 
                                                  group_var, baseline, 
                                                  individual_var, 
                                                  remove_baseline=FALSE) {
  if (nrow(metadata) != ncol(count_matrix)) {
    stop("Number of rows in metadata and count_matrix must match.")
  }
  if (!(group_var %in% colnames(metadata)) | !(individual_var %in% colnames(metadata))) {
    stop("group_var and individual_var must be valid column names in metadata.")
  }
  if (!(baseline %in% metadata[[group_var]])) {
    stop('Baseline should be a level of group_var')
  }
  individuals <- unique(metadata[[individual_var]])
  normalized_counts <- count_matrix
  for (individual in individuals) {
    # Subset metadata and count matrix for current individual
    individual_indices <- metadata[[individual_var]] == individual
    individual_metadata <- metadata[individual_indices, ]
    individual_counts <- count_matrix[,individual_indices]
    
    individual_counts <- normalizeCountsForHeatmap(individual_counts, individual_metadata, 
                                                   group_var, baseline,
                                                   remove_baseline=FALSE)
    
    normalized_counts[, individual_indices] <- individual_counts
  }
  if (remove_baseline) {
    baseline_indices <- metadata[[group_var]]==baseline
    normalized_counts <- normalized_counts[,!baseline_indices]
  }
  return(normalized_counts)
}
