#' runfgsea
#'
#' Perform GSEA on a DESeq2 result using the given pathways and fGSEA
#'  The function takes about O(nk^{3/2}) time, where n is number of permutations
#'  and k is a maximal size of the pathways. That means that setting 'maxSize'
#'  parameter with a value of ~500 is strongly recommended.
#'
#' @param result A DESeqResults result object
#' @inheritParams fgsea::fgseaSimple
#'
#' @inherit fgsea::fgseaSimple return
#' @export
#'
#' @importFrom fgsea fgseaSimple
#' @examples
#' \dontrun{
#' runfgsea(results[[1]], gmt.file)
#' }
#'
runfgsea <- function(result,
                     pathways,
                     nperm = 1000,
                     minSize = 10,
                     maxSize = 500) {
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
  return(res)
}
