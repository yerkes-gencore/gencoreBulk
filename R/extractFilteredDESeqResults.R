#' Extract results from processed DESeq2 object using minimum sample detection filtering
#'
#' Calls DESeq2::results with a per-gene filter requiring the numerator and 
#' denominator of the contrast to both have no more than `max_zero_samples` samples
#'  with at least `min_counts` reads.
#' This function is intended to reduce the
#' overall code used to extract all results of interest.
#'
#' @rdname extractFilteredDeseqResults
#'
#' @param contrast Contrast vector to pass to [DESeq2::results()] 
#' @param min_counts Minimum number of counts for a sample to not be ignored
#' @param max_min_samples Maximum number of low-count samples allowed in a group
#'  to still perform the calculation
#' @param dds DESeq2 object to draw results from
#' @param alpha The significance cutoff used for optimizing the independent
#'  filtering (by default 0.05). If the adjusted p-value cutoff (FDR) will be a
#'  value other than 0.05, alpha should be set to that value.
#' @param ... Additional arguments passed to [DESeq2::results()]
#' @param parallel Run DESeq2 commands in parallel
#' @param generate_shrunken_estimates Boolean to run and append `DESeq2::lfcShrink()`
#'  estimates to results
#' @inheritParams DESeq2::lfcShrink
#'
#' @returns object of class DESeqResults
#'
#' @import DESeq2
#' @importFrom matrixStats rowMaxs rowMins rowMeans2
#' @importFrom stats model.matrix
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#' 
#' @examples
#' \dontrun{
#' my_contrast <- limma::makeContrasts(
#' CelltypeBonemarrow_CD19posCD138pos-CelltypeBonemarrow_CD19negCD138pos,
#' (CelltypeBonemarrow_CD19posCD138pos+CelltypeBonemarrow_CD19posCD138neg)/2 -
#'   (CelltypeBonemarrow_CD19negCD138neg-CelltypeBonemarrow_CD19negCD138pos)/2,
#' levels = colnames(analysis$dds$modelMatrix))
#' 
#' tmp <- apply(my_contrast, 2, extractFilteredDESeqResults, dds = analysis$dds)#' 
#' }
#' 
#'
#' @export

extractFilteredDESeqResults <- function(dds, 
                                        contrast,
                                        min_counts = 0,
                                        max_min_samples = 0,
                                        alpha = 0.05,
                                        parallel = FALSE,
                                        generate_shrunken_estimates = TRUE,
                                        type = 'ashr',
                                         ...){
  model_matrix <- stats::model.matrix(design(dds),
                                      SummarizedExperiment::colData(dds))
  samples_used <- .extractSamplesInContrast(model_matrix, contrast)
  deseq_filter <- .maxMinFilter(dds,
                                group_samples = samples_used,
                                min_counts = min_counts,
                                max_min_samples = max_min_samples)
  out <- DESeq2::results(dds, 
                         alpha = alpha,
                         contrast = contrast,
                         filter = deseq_filter,
                         parallel = parallel,
                         ...)
  if (generate_shrunken_estimates){
    shrinks <- DESeq2::lfcShrink(dds, 
                                 res = out,
                                 type = type,
                                 parallel = parallel)
    out$shrunk_LFC <- shrinks$log2FoldChange
    out$shrunk_lfcSE <- shrinks$lfcSD
  }
  out
}


#' Filter genes with low read counts from results calculations
#' 
#' Remove genes with insufficient samples for calculating results. If a 
#'  group of samples has more than `max_min_samples` with fewer than 
#'  `min_counts` reads for a gene, the gene will be filtered out. 
#'
#' @param dds A DESeq object
#' @param group_samples List of character vectors containing sample names used
#'  in comparison groups
#' @param min_counts Minimum number of counts for a sample to not be ignored
#' @param max_min_samples Maximum number of 'zero' samples allowed in a group
#'  to still perform the calculation
#'
#'
#' @return A vector of filtering criteria to pass to `DESeq2::results()`
.maxMinFilter <- function(dds,
                          group_samples,
                          min_counts = 0,
                          max_min_samples = 1) {
  if_else(
    matrixStats::rowMaxs(
      sapply(group_samples, function(samples) {
        ## Rowwise apply
        apply(counts(dds, normalize = TRUE)[,samples], 1, function(x) {
          sum(x <= min_counts)
        })
      })
    ) > max_min_samples , 0, matrixStats::rowMeans2(counts(dds, normalized = TRUE))
  )
}

## A helper function to extact sample names from a design matrix and contrast,
## Used as input to maxMinFilter
.extractSamplesInContrast <- function(model_matrix, contrast){
  group1_terms <- names(contrast[contrast>0])
  group2_terms <- names(contrast[contrast<0])
  
  group1_mat <- model_matrix[,group1_terms]
  group2_mat <- model_matrix[,group2_terms]
  
  extract_samples <- function(mat){
    if (is.matrix(mat)){
      rownames(mat[rowSums(mat) != 0])
    } else {
      names(mat[mat != 0])
    }
  }
  return(list(group1_samples = extract_samples(group1_mat),
              group2_samples = extract_samples(group2_mat)))
}
