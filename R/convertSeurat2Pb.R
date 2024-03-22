#' Create `bulk` list from `SeuratObject`
#'
#' @description
#' Creates a `bulk` list compatible with `gencoreBulk` functions supporting the `sc_template` pseudobulk workflow.
#' 
#' @details
#' If `split_clusters = FALSE` (default), a `bulk` list is create with all clusters included. This is helpful if you intend to draw contrasts between clusters. In practice, this most helpful when few clusters are considered and the study design is relatively simple (few factors with few levels), as large numbers of clusters yield unwieldy sample tables and contrast matrices.
#' 
#' If `split_clusters = TRUE`, a list of `bulk` lists is create with each cluster isolated in its own `bulk` list. This allows you to fit the model on each cluster separately, but make it impossible to draw contrasts between clusters. This can be helpful as an exploratory step when there are many clusters of uncertain importance and/or for complex designs with many factor levels.
#' 
#' @param seurat_obj a Seurat or SeuratObject class data object. Must have counts in its RNA assay
#' @param sample character string specifying the column of `object@meta.data` that contains sample information of the single cell data.
#' @param cluster character string specifying the column of `object@meta.data` that contains single cell cluster labels to split on.
#' @param arrange_by character string specifying the column of `object@meta.data` that contains a factor or character vector to arrange the metadata data table by. Note that counts columns will be in same order as samples table to enable cleaner visualizations.
#' @param sample_md data.frame (or tibble) with all necessary factors for the design and downstream model fitting
#' @param design_str model design string. Default is "~ 0".
#' @param split_clusters whether to isolate clusters in separate objects so that you can fit models and define contrasts on each separately. See `Details` section.
#'
#' @export
#' @examples
#' \dontrun{
#' # Create sample metadata table
#' sample_md <- s.cd8@meta.data %>%
#'   dplyr::distinct(animalID, group, challenge) %>%
#'   as_tibble()
#' 
#' # Create a pseudobulk object with all clusters kept together
#' bulk <- createPBobj(s.cd8, cluster = "s.fullwnn1.0", 
#'                     sample = "animalID", sample_md = sample_md, 
#'                     arrange_by = "challenge", split_clusters = F,
#'                     design_str = "~ 0")
#'                     
#' # Create a pseudobulk object with all clusters kept together
#' bulk <- createPBobj(s.cd8, cluster = "s.fullwnn1.0", 
#'                     sample = "animalID", sample_md = sample_md, 
#'                     arrange_by = "challenge", split_clusters = T,
#'                     design_str = "~ 0")
#' }
convertSeurat2Pb <- function(seurat_obj, cluster = NULL, sample = NULL, arrange_by = NULL,
                             sample_md = NULL, design_str = "~ 0", split_clusters = FALSE) {
  ## Check that Seurat is installed, throw error if not
  rlang::check_installed(c("Seurat"), reason = "to use `convertSeurat2Pb`")
  ## Convert Seurat obj to PB
  pb_obj <- list()
  ## Splitting clusters before pseudobulking
  if (split_clusters) {
    suppressWarnings({ # warns if only one cluster but we want that
      pb_obj <- seurat_obj %>% 
        Seurat::SplitObject(split.by = cluster) %>%
        lapply(., function(x) {
          dge <- x %>%
            edgeR::Seurat2PB(sample = sample, cluster = cluster)
          # # Seurat2PB appends cluster name to colname, but we don't want that
          # colnames(dge) <- dge$samples$sample
          colnames(dge) <- stringr::str_remove(colnames(dge), pattern = "_cluster*$")
          # Add sample metadata
          dge$samples <- dge$samples %>%
            tibble::rownames_to_column(var = "rownames") %>%
            dplyr::select(-dplyr::all_of("group")) %>%
            dplyr::left_join(sample_md, by = c("sample" = sample)) %>%
            dplyr::arrange(.data[[arrange_by]]) %>%
            dplyr::mutate(sample = forcats::fct(sample) %>% forcats::fct_inorder())
          rownames(dge$samples) <- dge$samples$rownames
          # Arrange counts columns in same order as samples table
          dge$counts <- dge$counts[,match(rownames(dge$samples), colnames(dge$counts))]
          # Calculate library normalization factors for use downstream when removing heterscedascity with voom
          dge$samples$norm.factors <- edgeR::calcNormFactors(dge$counts, method = "TMM")
          # Create "bulk" object
          bulk <- list()
          bulk$dge <- dge
          # Add design to bulk$md
          bulk$md$design <- stats::model.matrix(stats::as.formula(design_str), dge$samples)
          return(bulk)
        })
    })
    ## Arrange clusters in alphabetical order
    pb_obj <- pb_obj[order(names(pb_obj))]
    return(pb_obj)
  }
  ## Not splitting clusters before pseudobulking
  else if (!split_clusters) {
    dge <- seurat_obj %>%
      edgeR::Seurat2PB(sample = sample, cluster = cluster)
    # Add sample metadata
    dge$samples <- dge$samples %>%
      tibble::rownames_to_column(var = "rownames") %>%
      dplyr::select(-dplyr::all_of("group")) %>%
      dplyr::left_join(sample_md, by = c("sample" = sample)) %>%
      dplyr::arrange(cluster, .data[[arrange_by]]) %>%
      dplyr::mutate(sample = forcats::fct(sample) %>% forcats::fct_inorder())
    rownames(dge$samples) <- dge$samples$rownames
    # Arrange counts columns in same order as samples table
    dge$counts <- dge$counts[,match(rownames(dge$samples), colnames(dge$counts))]
    # Calculate library normalization factors for use downstream when removing heterscedascity with voom
    dge$samples$norm.factors <- edgeR::calcNormFactors(dge$counts, method = "TMM")
    # Create "bulk" object
    bulk <- list()
    bulk$dge <- dge
    # Add design to bulk$md
    bulk$md$design <- stats::model.matrix(stats::as.formula(design_str), dge$samples)
    return(bulk)
  }
}