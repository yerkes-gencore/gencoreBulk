#' Bar plot of filtered genes
#'
#' Generate a bar plot of genes passing or failing a filtering procedure, eg. \code{\link[edgeR:filterByExpr]{?edgeR::filterByExpr()}}.
#' 
#' @param dge DGEList object
#' @param keep.exprs Named logical vector with gene ids as names. Names must match rownames of `dge$counts`
#'
#' @export
#' @examples
#' \dontrun{
#' keep.exprs <- filterByExpr(bulk$dge, design = bulk$md$design)
#' plotFilterByExpr(bulk$dge, keep.exprs)
#' }
#'
plotFilterByExpr <- function(dge, keep.exprs) {
  dge$counts %>%
    dplyr::as_tibble(rownames = "gene") %>%
    dplyr::left_join(dplyr::tibble(gene = names(keep.exprs), keep = keep.exprs), 
                     by = c("gene")) %>%
    tidyr::pivot_longer(cols = -c("gene", "keep"), names_to = "sampleID", values_to = "counts") %>%
    dplyr::mutate(sampleID = forcats::fct(.data$sampleID, levels = levels(dge$samples$sampleID)),
                  is_zero = .data$counts == 0) %>%
    dplyr::group_by(.data$sampleID) %>%
    dplyr::summarize(keep_zero = sum(.data$is_zero & .data$keep),
                     keep_nonzero = sum(!.data$is_zero & .data$keep),
                     remove_zero = sum(.data$is_zero & !.data$keep),
                     remove_nonzero = sum(!.data$is_zero & !.data$keep)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(cols = c(starts_with("keep_"), starts_with("remove_")), names_to = "cat", values_to = "n_genes") %>%
    dplyr::mutate(cat = forcats::fct(.data$cat, levels = rev(c("keep_nonzero", "keep_zero", "remove_nonzero", "remove_zero")))) %>%
    ggplot(aes(y = forcats::fct_rev(.data$sampleID), x = .data$n_genes, fill = .data$cat)) +
    geom_bar(stat="identity") +
    scale_fill_manual(name = "Filtering",
                      values = c("keep_nonzero"="green2", "keep_zero"="yellow2",
                                 "remove_nonzero"="orange", "remove_zero"="red"),
                      breaks = c("keep_nonzero", "keep_zero",
                                 "remove_nonzero", "remove_zero"),
                      labels = c("Kept non-zero", "Kept zero",
                                 "Removed non-zero", "Removed zero")) +
    xlab("Gene counts") + ylab("Sample ID") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

#' Plot MDS as ggplot
#'
#' Wrapper around \code{\link[limma:plotMDS]{limma::plotMDS()}} to enable ggplot2 formatting.
#' 
#' @param dge DGEList object
#' @param gene.selection `"common"` uses the same genes for all comparisons (i.e. PCA), `"pairwise"` chooses the top genes separately for each pairwise comparison between the samples (i.e MDS); see \code{\link[limma:plotMDS]{?limma::plotMDS()}} for details.
#' @param sampleID Label points by the value of this dge$samples column.
#' @param color Color points by the value of this dge$samples column.
#' @param path Connect points with the same value of this dge$samples column.
#' @param ... Arguments passed to \code{\link[limma:plotMDS]{limma::plotMDS()}}
#' 
#' @export
#' @examples
#' \dontrun{
#' ggplotMDS(bulk$dge, gene.selection = "pairwise", 
#'           sampleID = "sampleID", color = "grp", path = "SubjectID")
#' }
#'
ggplotMDS <- function(dge, sampleID = "sampleID", 
                      gene.selection = "common", 
                      color, path, ...) {
  # Get mds data from edgeR::plotMDS()
  mds_data <- limma::plotMDS(dge, top = 500, plot = FALSE, 
                             gene.selection = gene.selection, ...)
  
  mds_xy <- mds_data[c("x","y")]
  mds_xy[[sampleID]] <- colnames(dge)
  # mds_xy[[sampleID]] <- str_extract(colnames(dge), "(^[A-Za-z0-9]*)")
  mds_xy <- dplyr::as_tibble(mds_xy) %>% dplyr::full_join(dge$samples, by = sampleID)
  
  x_varex <- round(mds_data$var.explained[1]*100, digits = 0)
  y_varex <- round(mds_data$var.explained[2]*100, digits = 0)
  
  mds_xy %>%
    ggplot(aes(x = .data$x, y = .data$y, 
               color = .data[[color]], 
               label = .data[[sampleID]])) +
    geom_text() + 
    geom_path(aes(linetype = .data[[path]])) +
    xlab(paste0(mds_data$axislabel, " 1 (", x_varex, "%)")) +
    ylab(paste0(mds_data$axislabel, " 2 (", y_varex, "%)")) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(clip = "off") +
    ggtitle(ifelse(gene.selection == "common", "PCA", "MDS"))
}