#' Extract results from processed DESeq2 object using Greg's filtering
#'
#' Calls DESeq2::results with the filter specified using Greg's
#' minMaxFilter function. Greg's filter dynamically filters results for
#' minimum supporting reads based on the samples involved in the result calculation.
#' This function is intended to reduce the
#' overall code used to extract all results of interest.
#'
#' @rdname extractFilteredDeseqResults
#' @param comparison Value to pass to [DESeq2::results()] as either a `contrast`
#'  or a `name`
#' @param intgroup Metadata column on which to draw groups for filtering
#' @param filt_groups Character vector of groups to draw from `intgroup`
#' @param test Test to use, default is `Wald`
#' @param dds DESeq2 object to draw results from
#' @param alpha The significance cutoff used for optimizing the independent
#'  filtering (by default 0.05). If the adjusted p-value cutoff (FDR) will be a
#'  value other than 0.05, alpha should be set to that value.
#' @param ... Additional arguments passed to [DESeq2::results()]
#'
#' @returns object of class DESeqResults
#'
#' @import DESeq2
#' @importFrom matrixStats rowMaxs rowMins rowMeans2
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#'
#' @export

extractFilteredDESeqResults <- function(comparison,
                                        intgroup,
                                        filt_groups,
                                        test = "Wald",
                                        dds = analysis$dds,
                                        alpha = 0.05,
                                        listValues = c(1,-1),
                                        ...) {
  if (is(comparison, "character") & length(comparison) == 1) {
    result <- results(dds,
      name = comparison,
      test = test,
      alpha = alpha,
      listValues = listValues,
      ...,
      filter = .maxMinFilter(dds,
        intgroup = intgroup,
        comp = filt_groups
      )
    )
  } else {
    result <- results(dds,
      contrast = comparison,
      test = test,
      alpha = alpha,
      listValues = listValues,
      ...,
      filter = .maxMinFilter(dds,
        intgroup = intgroup,
        comp = filt_groups
      )
    )
  }
  result@metadata$contrast <- paste0(unlist(comparison), collapse = "_")
  return(result)
}


.maxMinFilter <- function(object,
                          intgroup,
                          comp,
                          thresh = 0) {
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])

  # group <- if (length(intgroup) > 1) {
  #   interaction(lapply(intgroup, function(factorname) colData(ddsMoTime20034)[[factorname]]), drop = TRUE)
  #   #    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  # }
  # else {
  group <- SummarizedExperiment::colData(object)[[intgroup]]
  # }
  if (!all(comp %in% levels(group))) {
    stop("the argument 'comp' should specify levels of intgroup")
  }
  if_else(
    matrixStats::rowMaxs(
      sapply(comp, function(lvl) {
        matrixStats::rowMins(counts(object, normalize = TRUE)[, group == lvl])
      })
    ) > thresh, matrixStats::rowMeans2(counts(object, normalized = TRUE)), 0
  )
}
