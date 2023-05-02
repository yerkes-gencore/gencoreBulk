#' Title
#'
#' @param result DESeqResults object
#' @param N Number of genes to return
#' @param min_logFC Minimum absolute value of log2FoldChange
#' @param min_padj Minimum adjusted p value
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
#' \dontrun {
#'   getTopNGenes(results$Group_FAC_vs_Cont)
#' }
#' 
getTopNGenes <- function(result,
                         N = 30,
                         min_logFC = log2(1.5),
                         min_padj  = analysis$analysis_config$alpha) {
  filtered_results <- data.frame(result) %>%
    dplyr::filter(!is.na(.data$pvalue)) %>%
    dplyr::filter(.data$padj < min_padj) %>%
    dplyr::filter(abs(.data$log2FoldChange) > min_logFC) %>%
    tibble::rownames_to_column('Gene') %>%
    dplyr::arrange(.data$pvalue) %>%
    utils::head(N) 
  return(filtered_results$Gene)
}