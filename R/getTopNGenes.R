#' Retreive the top differentially expressed genes from a DESeq results obj
#'
#' @param result DESeqResults object
#' @param N Number of genes to return
#' @param min_logFC Minimum absolute value of log2FoldChange
#' @param min_padj Minimum adjusted p value
#' @param exclude_ENS Boolean: exclude genes with ENSEMBL IDs instead of symbols
#' @param direction Which direction log-fold change to return. Accepts values 
#'  c('up', 'down', 'mixed', 'equal'). Mixed returns the top N genes regardless
#'  of direction, while equal returns an equal number of up and down regulated
#'  genes
#' @param ENS_pattern Pattern used for ENSEMBL id recognition, default "^ENS"
#'
#' @returns Character vector of gene names
#' @export
#'
#' @importFrom utils head 
#' @importFrom rlang .data
#' @importFrom tibble rownames_to_column
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'   getTopNGenes(results$Group_FAC_vs_Cont)
#' }
#' 
getTopNGenes <- function(result,
                         N = 30,
                         min_logFC = log2(1.5),
                         min_padj  = 0.05,
                         exclude_ENS = TRUE,
                         ENS_pattern = '^ENS',
                         direction = 'mixed') {
  filtered_results <- data.frame(result) %>%
    dplyr::filter(!is.na(.data$pvalue)) %>%
    dplyr::filter(.data$padj < min_padj) %>%
    dplyr::filter(abs(.data$log2FoldChange) > min_logFC) %>%
    tibble::rownames_to_column('Gene') %>%
    dplyr::arrange(.data$pvalue)
  if (exclude_ENS){
    filtered_results <- filtered_results[grepl(pattern = ENS_pattern,
                                               filtered_results$Gene),]
  }
  if (direction == 'mixed') {
    filtered_results <- filtered_results %>%
      utils::head(N) 
  } else if (direction == 'up'){
    filtered_results <- filtered_results %>%
      dplyr::filter(.data$log2FoldChange > 0) %>%
      utils::head(N) 
  } else if (direction == 'down') {
    filtered_results <- filtered_results %>%
      dplyr::filter(.data$log2FoldChange < 0) %>%
      utils::head(N) 
  } else if (direction == 'equal') {
    filtered_results_up <- filtered_results %>%
      dplyr::filter(.data$log2FoldChange < 0) %>%
      utils::head(round(N/2)) 
    filtered_results_down <- filtered_results %>%
      dplyr::filter(.data$log2FoldChange > 0) %>%
      utils::head(round(N/2)) 
    filtered_results <- rbind(filtered_results_up,filtered_results_down)
  } else {
    stop('Set "direction" to one of "equal", "mixed", "up", or "down"')
  }
  return(filtered_results$Gene)
}