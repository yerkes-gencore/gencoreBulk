#' Return a simple ggplot boxplot of gene expression
#'
#' Returns a simple ggplot of gene expression, optionally with jitter or 
#'  boxplots overlayed. These are optional in case you want to do other 
#'  aesthetics. This is meant to be minimal so you can adjust the aesthetics
#'  as necessary, but it can still be used as a quick check. 
#'  
#' @param gene Gene to plot expression for
#' @param grouping Variable to group expression by. Must be a column of `data@colData`
#' @param data A DESeqDataSet object. Use this argument to select specific assays
#'  to plot. The default is `assays(analysis$dds)$rld` (assumed to be the rlog transform).
#' @param boxplot Boolean, plot `geom_boxplot`?
#' @param jitter Boolean, plot `geom_jitter`?
#' @param axis_text_x Option to modify X axis text. Default rotates labels 90
#'  degrees
#'
#' @returns A ggplot object
#' @export
#'
#' @import ggplot2
#' @importFrom SummarizedExperiment assay
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' ## my object is called dds_full, so I have to manually provide data
#' plotGeneExpression('MAMU-A', grouping = 'groupIDs', data = assays(dds_full)$rld)
#' }
#' 
plotGeneExpression <- function(gene,
                               grouping,
                               data = SummarizedExperiment::assays(analysis$dds)$rld,
                               boxplot = TRUE,
                               jitter = TRUE,
                               axis_text_x = element_text(angle = 90, vjust = 0.5, hjust=1)) {
  if (!(gene %in% rownames(data))){
    stop('Gene not found in data')
  }
  plotdata <- as.data.frame(t(SummarizedExperiment::assay(data[gene, ])))
  colnames(plotdata) <- 'Gene'
  plotdata$grouping <- data@colData[[grouping]]
  ggplot2::ggplot(plotdata, aes(x=grouping, y=.data[[Gene]])) + 
    (if (boxplot) { ggplot2::geom_boxplot(outlier.color = if (jitter) {NA} else {'black'}) }) +
    ggplot2::theme_bw() +
    (if (jitter) { ggplot2::geom_jitter() }) +
    ggplot2::labs(x = 'Group', y = 'Expression', main = gene) +
    ggplot2::theme(axis.text.x = axis_text_x)
}