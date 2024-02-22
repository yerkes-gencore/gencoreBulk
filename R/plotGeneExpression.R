#' Return a simple ggplot boxplot of gene expression
#'
#' Returns a simple ggplot of gene expression, optionally with jitter or 
#'  boxplots overlayed. These are optional in case you want to do other 
#'  aesthetics. This is meant to be minimal so you can adjust the aesthetics
#'  as necessary, but it can still be used as a quick check. 
#'  
#' You can subset data using the `subsetting` and `subsets` arguments to focus
#'  the plotting.
#'
#'  
#' @param gene Gene to plot expression for
#' @param counts The expression data to plot
#' @param metadata The metadata to associate with the counts data
#' @param grouping Variable to group expression by. Must be a column of `metadata`
#' @param groups Optional specification of level ordering for grouping var
#' @param subsetting Optional variable to subset expression by. Must be a column of `metadata`
#' @param subsets Specification of levels to subset for
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
#'   plotGeneExpression('MAMU-A', counts = model_fit$EList$E,
#'                       metadata = model_fit$targets, grouping = 'groupIDs')
#'   
#'   ## With a DESeq model
#'   plotGeneExpression('MAMU-A', counts = assay(assays(analysis$dds)$rld), 
#'                       metadata = colData(analysis$dds), grouping = 'groupIDs')
#' }
#' 
plotGeneExpression <- function(gene,
                               grouping,
                               groups,
                               counts,
                               metadata,
                               subsetting,
                               subsets,
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
  if (xor(missing(subsetting), missing(subsets))) {
    stop('Subsetting data requires both subsets and subsetting to be specified')
  } 
  if (any(colnames(counts) != rownames(metadata))) {
    warning('Colnames of counts does not match rownames of metadata. Continuing, but
            this may suggest the data are not properly associated')
  }
  if (!missing(subsetting)) {
    if (!subsetting %in% colnames(metadata)) {
      stop('subsetting variable not found in metadata')
    }
    metadata <- metadata[metadata[[subsetting]] %in% subsets,]
    counts <- counts[,rownames(metadata)]
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




#' Plot the results of a model fit 
#' 
#' Plots the expression, coefficients, and estimated change for a single gene and
#' contrast. You have to pass in a contrasts dataframe (or list) and specify which 
#' contrast to use. See the example for details. The rest of the arguments get 
#' passed to `plotGeneExpression`. 
#' 
#' For example, if the question is:
#'  "the results say JUNB has a significant increase in expression,
#'   in Day 14 compared to day 0, does this seem to be true?"
#' 
#' Your numerator/denominator specification would look something like:
#'  numerator = '(Intercept) + dayDay14', denominator = '(Intercept)'
#'
#' @param gene Gene to plot
#' @param numerator Full numerator of contrast, including all cancelled terms
#' @param denominator Full denominator of contrast, including all cancelled terms
#' @param coefficients Matrix of coefficients fit by model of your choice
#' @inheritParams plotGeneExpression
#'
#' @returns A ggplot object
#' @export
#' 
#' @import ggplot2
#'
#' @examples 
#' \dontrun{ 
#' 
#' contrast_matrix_expanded <- data.frame(
#'   contrast_names = 'P14_grp3vgrp2',
#'   numerators = '(Intercept) + grpgrp3',
#'   denominators = '(Intercept)')
#'   
#' model_plot('JUN', 
#'              counts = model_fit$EList$E, metadata = model_fit$targets,
#'              grouping = 'grp.day',
#'              subsetting = 'grp.day', subsets = c('grp2.P14', 'grp3.P14'),
#'              coefficients = limma_fit$fit$coefficients,
#'              numerator = '(Intercept) + grpgrp3',
#'              denominator = '(Intercept)')
#' }
plotModelCoeffs <- function(gene, 
                            numerator, denominator, 
                            counts, metadata, 
                            grouping, groups, 
                            subsetting, subsets,
                            coefficients) { 
  baseplot <- plotGeneExpression(gene,
                                 counts =  counts,
                                 metadata = metadata,
                                 grouping = grouping,
                                 groups = groups,
                                 subsetting = subsetting,
                                 subsets = subsets)
  contrast_num <- round(.extract_coefficients(gene = gene, 
                                              terms = numerator, 
                                              coefficients = coefficients), 2)
  contrast_den <- round(.extract_coefficients(gene = gene, 
                                              terms = denominator, 
                                              coefficients = coefficients), 2)
  contrast_line_num <- ggplot2::geom_hline(yintercept = contrast_num, linetype = 'dashed') 
  contrast_line_den <- ggplot2::geom_hline(yintercept = contrast_den, linetype = 'dashed') 
  
  arrow <- ggplot2::geom_segment(aes(x = 1.5, xend = 1.5, 
                                     y = contrast_den, 
                                     yend = contrast_num),
                                 arrow = ggplot2::arrow()) 
  baseplot + 
    contrast_line_num +
    contrast_line_den +
    arrow +
    ggplot2::scale_y_continuous(breaks = c(contrast_den, contrast_num))
}

.extract_coefficients <- function(gene,
                                  coefficients,
                                  terms) {
  if (!gene %in% rownames(coefficients)) {
    stop('Gene not found in coefficients')
  }
  coefficients <- coefficients[gene,]
  terms <- unlist(strsplit(terms, '\\s*\\+\\s*', perl = TRUE))
  value <- 0
  for (term in terms) {
    if (!term %in% names(coefficients)) {
      stop(paste0("Error: Term '", term,
                  "' not found in matrix of coefficients"))
    }
    value <- value + unname(coefficients[term])
  }
  return(value)
}
