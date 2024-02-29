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

#' Create arrows for coefficients from model fit
#' 
#' Extracts coefficients from a matrix based on the gene and terms of interest.
#' Creates a dataframe for a geom_segment to be added to a plot like that returned
#' from `gencoreBulk::plotGeneExpression()`
#'
#' @param gene Gene to plot
#' @param expression Character vector of a contrast expression,
#'  consisting of terms present in `coefficients` added or subtracted together
#' @param coefficients Matrix of coefficients from model fit
#' @param data_only Whether to return a dataframe instead of a ggplot object. 
#'  Default FALSE.
#'
#' @returns A ggplot2::geom_segment if data_only = FALSE, else a data.frame
#' @export
#'
#' @examples 
#' \dontrun{
#'  plotGeneExpression(gene = 'JUN', 
#'                    counts = model_fit$EList$E,
#'                    metadata = model_fit$targets,
#'                    grouping = 'grp',
#'                    subsetting = 'day', subsets = 'D28') +
#'  plotModelCoeffs(gene = 'JUN',
#'                  expr = "(Intercept) + dayD28 + grpgrp3 + dayD28:grpgrp3",
#'                  coefficients = model_fit$coefficients)  
#'                  
#'  ## Or you can return the data for manual plotting
#'  arrow_coords <- plotModelCoeffs(gene = 'JUN',
#'                                  "(Intercept) + dayD28 + grpgrp3 + dayD28:grpgrp3",
#'                                  coefficients = model_fit$coefficients, 
#'                                  data_only = TRUE)  
#'  ## Can edit values here if desired, such as X coordinates
#'  ## arrow_coords$x = 0.5
#' 
#'  plotGeneExpression(gene = 'JUN', 
#'                 counts = model_fit$EList$E,
#'                 metadata = model_fit$targets,
#'                 grouping = 'grp',
#'                 subsetting = 'day', subsets = 'D28') +
#'  ## Then manually create the arrows
#'    ggplot2::geom_segment(data = arrow_coords, 
#'                          aes(x = x,
#'                              xend = x,
#'                              y = y, 
#'                              yend = yend,
#'                              color = terms),
#'                          arrow = ggplot2::arrow()) 
#' }
plotModelCoeffs <- function(gene,
                            expression,
                            coefficients,
                            data_only = FALSE) { 
  expression <- .substitute_terms(expression, coefficients[gene,])
  rt <- cumsum(expression$dict)
  arrow_coords <- data.frame(terms = factor(names(rt), levels = names(rt)),
                             y = c(0, rt[-length(rt)]), yend = rt,
                             x = seq(1.4, 1.6, length.out = length(rt)))
  if (data_only) {
    return(arrow_coords)
  }
  ggplot2::geom_segment(data = arrow_coords, 
                        aes(x = .data[['x']],
                            xend = .data[['x']],
                            y = .data[['y']], 
                            yend = .data[['yend']],
                            color = .data[['terms']]),
                        arrow = ggplot2::arrow()) 
}

.substitute_terms <- function(formula, value_list) {
  terms <- unlist(strsplit(as.character(formula), "[+-]"))
  term_dict <- list()
  for (term in terms) {
    term <- trimws(term)
    if (term %in% names(value_list)) {
      value <- value_list[[term]]
      formula <- gsub(term, value, formula, fixed = TRUE)
      term_dict[[term]] <- as.numeric(value)
    }
  }
  return(list(expr = formula, dict = term_dict))
}