#' Title
#'
#' Makes a heatmap of the given list of genes, separating samples slightly by
#'  group variable. The heatmap will render samples in the order of colnames
#'  in the provided data, and will be split by the metadata column specified.
#'  Expression values will be normalized to the median of a specified group.
#'
#'  If you want to reorder columns plotted, arrange the data passed into `data`.
#'
#' @param geneList Character vector of genes to plot
#' @param baseline_grouping Character, level of `colnames(colData(data))`
#' @param baseline level of `baseline_grouping` in `colData(data)`
#' @param data Object of class DESeqTransform
#' @param slice_labels Optional labels for sliced columns
#' @param slice_labels_rot Rotation angle of `slice_labels`
#' @param legend_title Title for legend
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
#' @importFrom SummarizedExperiment assay
#' @importFrom plyr mapvalues
#' @importFrom matrixStats rowMedians
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' \dontrun{
#' ## Simple call
#' heatmapFromGenelist(
#'   geneList = c("Ccl2", "Cxcl1", "Cxcl2", "Postn", "Fn1", "Thbs1"),
#'   baseline_grouping = "Group",
#'   baseline = "Cont"
#' )
#'
#' ## run with custom reordering of data and labeled slices
#' heatmapFromGenelist(
#'   geneList = c("Ccl2", "Cxcl1", "Cxcl2", "Postn", "Fn1", "Thbs1"),
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
                                baseline_grouping,
                                baseline,
                                column_split = NULL,
                                slice_labels = NULL,
                                data = assays(analysis$dds)$rld,
                                column_labels = colnames(data),
                                slice_labels_rot = 90,
                                legend_title = "log2 fold\ndifference\nfrom\nmedian\nbaseline\nexpression",
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
                                ...) {
  if (!baseline_grouping %in% colnames(colData(data))) {
    stop("Argument 'baseline_grouping' should be in colData(data)")
  }
  if (!baseline %in% unique(colData(data)[[baseline_grouping]])) {
    stop("Argument 'baseline' should be a level of 'baseline_grouping' in colData(data)")
  }

  hmap <- data[geneList, ]
  baseline <- matrixStats::rowMedians(assay(hmap[, as.character(hmap@colData[[baseline_grouping]]) %in% baseline]))
  hmap <- assay(hmap) - baseline
  ComplexHeatmap::Heatmap(hmap,
    heatmap_legend_param = list(title = legend_title),
    border = "black",
    width = ncol(hmap) * box_width + width_buffer,
    height = nrow(hmap) * box_height + height_buffer,
    rect_gp = grid::gpar(color = "black"),
    column_title = column_title,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    column_split = column_split,
    top_annotation = (if (!is.null(slice_labels)) {
      if (is.null(column_split)) {
        warning("Setting labels requires slices to also be set")
      }
      HeatmapAnnotation(foo = anno_block(
        gp = gpar(col = NA),
        labels = slice_labels,
        labels_gp = gpar(col = "black", fontsize = 10),
        labels_rot = slice_labels_rot, height = unit(2, "cm")
      ))
    } else {
      NULL
    }),
    column_gap = column_gap,
    col = circlize::colorRamp2(c(scale_min, 0, scale_max), c("blue", "white", "red")),
    ...
  )
}