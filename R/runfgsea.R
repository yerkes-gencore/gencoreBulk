#' runfgsea
#'
#' Perform GSEA on a DESeq2 result using the given pathways and fGSEA
#'  The function takes about O(nk^{3/2}) time, where n is number of permutations
#'  and k is a maximal size of the pathways. That means that setting 'maxSize'
#'  parameter with a value of ~500 is strongly recommended.
#'
#' @param result A DESeqResults result object
#' @param breakdown_pathway_names Attempt to split pathway names into source and
#'  pathway for easier reading
#' @inheritParams fgsea::fgseaSimple
#'
#' @inherit fgsea::fgseaSimple return
#' @export
#'
#' @importFrom fgsea fgseaSimple
#' @importFrom stringr str_split
#' @examples
#' \dontrun{
#' runfgsea(results[[1]], gmt.file)
#' }
#'
runfgsea <- function(result,
                     pathways,
                     nperm = 1000,
                     minSize = 10,
                     maxSize = 500,
                     breakdown_pathway_names = FALSE) {
  # result <- result %>% na.omit()
  fgsea_data <- result$stat
  names(fgsea_data) <- rownames(result)
  fgsea_data <- fgsea_data[!is.na(fgsea_data)]
  fgsea_data <- sort(fgsea_data, decreasing = TRUE)
  res <- fgsea::fgseaSimple(
    pathways = pathways,
    stats = fgsea_data,
    nperm = nperm,
    minSize = minSize,
    maxSize = maxSize
  )
  res$source <- unlist(lapply(res$pathway, function(x){stringr::str_split(x, '_')[[1]][1]}))
  if (breakdown_pathway_names) {
    res$pathway_short <- unlist(lapply(gsea_result$pathway, function(x){gsub('^[^_]+_(.+)', replacement = '\\1', x = x)}))
  }
  return(res)
}
