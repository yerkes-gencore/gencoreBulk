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
  mds_xy <- mds_xy %>% dplyr::full_join(dge$samples, by = sampleID)
  
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