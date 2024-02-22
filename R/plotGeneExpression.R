#' Return a simple ggplot boxplot of gene expression
#'
#' Returns a simple ggplot of gene expression, optionally with jitter or 
#'  boxplots overlayed. These are optional in case you want to do other 
#'  aesthetics. This is meant to be minimal so you can adjust the aesthetics
#'  as necessary, but it can still be used as a quick check. 
#'  
#' @param gene Gene to plot expression for
#' @param counts The expression data to plot
#' @param metadata The metadata to associate with the counts data
#' @param grouping Variable to group expression by. Must be a column of `metadata`
#' @param groups 
#' @param boxplot Boolean, plot `geom_boxplot`?
#' @param jitter Boolean, plot `geom_jitter`?
#' @param axis_text_x Option to modify X axis text. Default rotates labels 90
#'  degrees
#'
#' @returns A ggplot object
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom utils stack
#' @examples
#' \dontrun{
#' 
#'   ## With a limma model fit
#'   plotGeneExpression('MAMU-A', counts = model_fit$EList$E, metadata = model_fit$targets, grouping = 'groupIDs')
#' }
#' 
plotGeneExpression <- function(gene,
                               grouping,
                               groups,
                               counts,
                               metadata,
                               boxplot = TRUE,
                               jitter = TRUE,
                               axis_text_x = element_text(angle = 90, vjust = 0.5, hjust=1)) {
  # if (missing(gene) | missing(grouping) | missing(counts) | missing(metadata)) {
  #   stop('The following arguments are all required: gene, grouping, groups, counts, metadata')
  # }
  if (!(gene %in% rownames(counts))){
    stop('Gene not found in counts data')
  }
  if (!grouping %in% colnames(metadata)) {
    stop('grouping variable not found in metadata')
  }
  plotdata <- utils::stack(counts[gene, ])
  colnames(plotdata) <- 'Gene'
  plotdata$grouping <- metadata[[grouping]]
  if (!missing(groups)) {
    if (!all(groups %in% unique(metadata[[grouping]]))) {
      stop('All groups not found in grouping variable of metadata')
    }
    plotdata$grouping <- factor(plotdata$grouping, levels = groups)
  }
  ggplot2::ggplot(plotdata, aes(x=grouping, y=.data[['Gene']])) + 
    (if (boxplot) { ggplot2::geom_boxplot(outlier.color = if (jitter) {NA} else {'black'}) }) +
    (if (jitter) { ggplot2::geom_jitter() }) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = 'Group', y = 'Expression', main = gene) +
    ggplot2::theme(axis.text.x = axis_text_x)
}