#' Running model fitting procedure using voom + LmFit
#'
#' @description
#' Model fitting and contrast extraction using voom + LmFit procedure as implemented in \code{\link[edgeR:voomLmFit]{edgeR::voomLmFit()}}.
#' 
#' @param bulkObj List object with raw counts in `bulkObj$dge$counts` and design matrix in `bulkObj$md$design`.
#' @param contr.matrix Contrast matrix created by `limma::makeContrasts()`
#' @param plotVoom logical, should a plot of the mean-variance trend be displayed?
#' @inheritParams edgeR::voomLmFit
#' @inheritParams limma::contrasts.fit
#' 
#' @details
#' This is a convenience wrapper around \code{\link[edgeR:voomLmFit]{edgeR::voomLmFit()}}, \code{\link[limma:contrasts.fit]{limma::contrasts.fit()}}, \code{\link[limma:eBayes]{limma::eBayes()}}, and \code{\link[limma:decideTests]{limma::decideTests()}}.
#' 
#' # Explanation of workflow:
#' 
#' ## edgeR::voomLmFit() (Adapted from edgeR::voomLmFit() documentation):
#' 
#' This function adapts the limma voom method (Law et al, 2014) to allow for loss of residual degrees of freedom due to exact zero counts (Lun and Smyth, 2017). 
#' 
#' The function is analogous to calling `voom` followed by `duplicateCorrelation` and `lmFit` except for the modified residual df values and residual standard deviation sigma values.
#' 
#' If block is specified, then the intra-block correlation is estimated using `duplicateCorrelation`. In that case, the voom weights and the intra-block correlation are each estimated twice to achieve effective convergence.
#' 
#' Empirical sample quality weights will be estimated if `sample.weights=TRUE` or if `var.design` or `var.group` are non-NULL (Liu et al 2015). In that case, `voomLmFit` is analogous to running `voomWithQualityWeights` followed by `lmFit`.
#' 
#' ## limma::contrasts.fit():
#' 
#' Re-orient the fitted model object from the coefficients of the original design matrix to the set of contrasts defined above in contr.matrix
#' 
#' ## limma::eBayes():
#' 
#' Run empirical Bayes moderation; borrowing information across all the genes to obtain more precise estimates of gene-wise variability
#' 
#' ## plotSA():
#' 
#' Plot the model's residual variances against average expression values; demonstrating that the variance is no longer dependednt on the mean expression level
#' 
#' ## decideTests():
#' 
#' Identify which genes are significantly differentially expressed for each contrast from a fit object containing p-values and test statistics
#' 
#' @export
#' @examples
#' \dontrun{
#' ## Fit a model with repeated measures design with SubjectID treated as a random effect
#' bulk <- runVoomLmFit(bulk, contr.matrix=contr.matrix, 
#'                      sample.weights = TRUE, 
#'                      block = bulk$dge$samples$SubjectID)
#' }
#'
fitVoomLm <- function(bulkObj, contr.matrix, block = NULL, sample.weights = TRUE, var.design = NULL, var.group = NULL, plotVoom = TRUE) {
    bulkObj$fit <- edgeR::voomLmFit(counts = bulkObj$dge$counts, 
                                    design = bulkObj$md$design,
                                    block = block,
                                    sample.weights = sample.weights,
                                    var.design = var.design,
                                    var.group = var.group,
                                    plot = plotVoom)
  bulkObj$contrasts <- limma::contrasts.fit(bulkObj$fit, contrasts = contr.matrix)
  bulkObj$contrasts <- limma::eBayes(bulkObj$contrasts, robust=TRUE)
  limma::plotSA(bulkObj$contrasts, main="Final model: Mean-variance trend")
  bulkObj$de <- limma::decideTests(bulkObj$contrasts, lfc = 0)
  return(bulkObj)
}

#' Running model fitting procedure using edgeR's glmQLFit()
#'
#' @description
#' Model fitting and contrast extraction using quasi-likelihood negative binomial GLM procedure as implemented in \code{\link[edgeR:glmQLFit]{edgeR::glmQLFit()}}.
#' 
#' @param bulkObj List object with raw counts in `bulkObj$dge$counts` and design matrix in `bulkObj$md$design`.
#' @param contr.matrix Contrast matrix created by `limma::makeContrasts()`
#' @export
#' @examples
#' \dontrun{
#' bulk <- fitGlmQL(bulk, contr.matrix=contr.matrix)
#' }
#'
fitGlmQL <- function(bulkObj, contr.matrix) {
  bulkObj$dge <- edgeR::estimateDisp(bulkObj$dge, bulkObj$md$design)
  edgeR::plotBCV(bulkObj$dge)
  bulkObj$fit <- edgeR::glmQLFit(bulkObj$dge, bulkObj$md$design, robust = TRUE, legacy = FALSE)
  bulkObj$res <- list()
  for (i in 1:ncol(contr.matrix)){
    bulkObj$res[[i]]<- edgeR::glmQLFTest(bulkObj$fit, contrast = contr.matrix[,i])
  }
  names(bulkObj$res) <- colnames(contr.matrix)
  return(bulkObj)
}

#' Create results table from `bulk` list object
#'
#' @description
#' Create tibble of top genes from linear model fit for each contrast \code{\link[limma:topTable]{limma::topTable()}}.
#' 
#' @inheritParams limma::topTable
#' @param contr.matrix Contrast matrix created by `limma::makeContrasts()`
#' @export
#' @examples
#' \dontrun{
#' bulk <- runVoomLmFit(bulk, contr.matrix=contr.matrix, 
#'                      sample.weights = TRUE, 
#'                      block = bulk$dge$samples$SubjectID)
#' resultsTable <- createResTable(bulk, contr.matrix)
#' }
#'
createResTable <- function(fit, contr.matrix) {
  resultsTables_list <- list()
  for (contrast in colnames(fit$coefficients)) {
    resultsTables_list[[contrast]] <- limma::topTable(fit, coef = contrast, n = Inf) %>%
      # as_tibble(rownames = "gene") %>%
      dplyr::rename(log2FoldChange = .data$logFC, pvalue = .data$P.Value, padj = .data$adj.P.Val)
  }

  resultsTable <- lapply(resultsTables_list, function(one_tbl) {
    one_tbl %>% tibble::rownames_to_column(var = "gene")
  }) %>% dplyr::bind_rows(.id = "contrast") %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(contrast = forcats::fct(contrast, levels = colnames(contr.matrix)))

  return(resultsTable)
}
