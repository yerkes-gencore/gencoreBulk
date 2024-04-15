#' Generate heatmap of gene expression from DESeq2 results
#'
#' Makes a heatmap of the given list of genes, separating samples slightly by
#'  group variable. The heatmap will render samples in the order of colnames
#'  in the provided data, and will be split by the metadata column specified.
#'  Normalization by groups or samples is no longer contained in this function.
#'  See `gencoreBulk::normalizeCounts` for help on normalizing data for heatmaps.
#'
#'  If you want to reorder columns plotted, arrange the data passed into `data`.
#'  See the examples
#'
#' @param geneList Character vector of genes to plot
#' @param data Matrix of processed counts
#' @param slice_labels Optional labels for sliced columns
#' @param colors vector of color values passed to colorRamp2, default is c("blue", "white", "red")
#' @param slice_labels_rot Rotation angle of `slice_labels`
#' @param box_width Width of heatmap cell
#' @param box_height Height of heatmap cell
#' @param width_buffer Additional width buffer to add to final heatmap
#' @param height_buffer Additional height buffer to add to final heatmap
#' @param scale_min Minimum value for scaling colors
#' @param scale_max Maximum value for scaling colors
#' @inheritParams ComplexHeatmap::Heatmap
#' @inheritDotParams ComplexHeatmap::Heatmap
#'
#' @inherit ComplexHeatmap::Heatmap return
#'
#' @import ComplexHeatmap
#' @importFrom SummarizedExperiment assay assays
#' @importFrom plyr mapvalues
#' @importFrom matrixStats rowMedians
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' \dontrun{
#' 
#' ## Normalize data for plotting
#' data_to_plot <- normalizeCountsForHeatmapByIndividual(assay(assays(obj.pbmc)$rld),
#'   obj.pbmc@colData,
#'   group_var = 'Timepoint', baseline = 'pre', 
#'   individual_var = 'Individual',
#'   remove_baseline = FALSE)
#'   
#' ## Simple call
#' heatmapFromGenelist(
#'   geneList = c("Ccl2", "Cxcl1", "Cxcl2", "Postn", "Fn1", "Thbs1"),
#'   data_to_plot,
#'   baseline_grouping = "Group",
#'   baseline = "Cont"
#' )
#'
#' ## run with custom reordering of data and labeled slices
#' heatmapFromGenelist(
#'   geneList = c("Ccl2", "Cxcl1", "Cxcl2", "Postn", "Fn1", "Thbs1"),
#'   data_to_plot,
#'   baseline_grouping = "Group",
#'   baseline = "Cont",
#'   column_split = c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3)),
#'   slice_labels = c("Cont", "Fac", "L", "M"),
#'   data = assays(analysis$dds)$rld[, sort(colnames(assays(analysis$dds)$rld))]
#' )
#' }
#'
#' @export
heatmapFromGenelist <- function(geneList,
                                 data,
                                 column_split = NULL,
                                 slice_labels = NULL,
                                 colors = c("blue", "white", "red"),
                                 column_labels = colnames(data),
                                 slice_labels_rot = 90,
                                 box_width = unit(3.5, "mm"),
                                 box_height = unit(3.5, "mm"),
                                 width_buffer = unit(5, "mm"),
                                 height_buffer = unit(10, "mm"),
                                 column_title = " ",
                                 cluster_rows = FALSE,
                                 cluster_columns = FALSE,
                                 column_gap = unit(2, "mm"),
                                 scale_min = -2,
                                 scale_max = 2,
                                 heatmap_legend_param = list(
                                   at = c(scale_min, 0, scale_max),
                                   labels = c(scale_min,0,scale_max),
                                   title = 'log2 fold\nchange'
                                 ),
                                 ...) {
  if (inherits(data, 'DESeqTransform')) {
    message('Converting DESeqTransform data into matrix')
    data <- SummarizedExperiment::assay(data)
  }
  duds <- geneList[!geneList %in% rownames(data)]
  if (length(duds) > 0){
    geneList <- geneList[geneList %in% rownames(data)]
    if (length(geneList) == 0){
      stop('No data for requested genes')
    } else {
      message(paste0('Genes ', paste0(duds, collapse = ', '), ' not found in data'))
    }
  }
  data <- data[geneList,,drop=FALSE]
  ComplexHeatmap::Heatmap(data,
                          heatmap_legend_param = heatmap_legend_param,
                          #border = "black",
                          width = ncol(data) * box_width + width_buffer,
                          height = nrow(data) * box_height + height_buffer,
                          #rect_gp = grid::gpar(color = "black"),
                          column_title = column_title,
                          cluster_rows = cluster_rows,
                          cluster_columns = cluster_columns,
                          column_split = column_split,
                          top_annotation = (if (!is.null(slice_labels)) {
                            if (is.null(column_split)) {
                              warning("Setting labels requires slices to also be set")
                            }
                            ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(
                              gp = gpar(col = NA),
                              labels = slice_labels,
                              labels_gp = gpar(col = "black", fontsize = 10),
                              labels_rot = slice_labels_rot, height = unit(2, "cm")
                            ))
                          } else {
                            NULL
                          }),
                          column_gap = column_gap,
                          column_labels = column_labels,
                          col = circlize::colorRamp2(c(scale_min, 0, scale_max), colors),
                          ...
  )
}
