#' plot PCA from config
#'
#' Returns plot of PCA using aesthetic labels specified in the analysis config
#'
#' @rdname plotPCAFromConfig
#'
#' @param analysis object with config options saved from the template format
#' @param shapes Shape values to pass to `ggplot2::scale_shape_manual()`
#' @param alpha Alpha to use for `ggplot2::geom_point()`
#' @param size Size to use for `ggplot2::geom_point()`
#' @param pc.1 Principal component to plot on X axis
#' @param pc.2 Principal component to plot on Y axis
#' 
#' @returns A ggplot
#'
#' @examples
#' \dontrun{
#' plotPCAFromConfig(analysis)
#' }
#'
#' @import ggplot2
#' @import ggrepel
#' @importFrom SummarizedExperiment assays assay
#' @importFrom rlang .data
#'
#' @export
plotPCAFromConfig <- function(analysis,
                              shapes = 0:25,
                              alpha = 1,
                              size = 5,
                              pc.1 = 1,
                              pc.2 = 2) {
  pcaPlotSimple(counts = assay(assays(analysis$dds)$vst),
    metadata = colData(analysis$dds),
    xpc = pc.1, ypc = pc.2
  ) +
    geom_point(
      aes(
        color = (if (is.null(analysis$qc_config$pcaMapping$color)) {
          NULL
        } else {
          .data[[analysis$qc_config$pcaMapping$color]]
        }),
        shape = (if (is.null(analysis$qc_config$pcaMapping$shape)) {
          NULL
        } else {
          .data[[analysis$qc_config$pcaMapping$shape]]
        }),
        text = (if (is.null(analysis$qc_config$pcaMapping$hover)) {
          NULL
        } else {
          .data[[analysis$qc_config$pcaMapping$hover]]
        })
      ),
      size = size,
      alpha = alpha
    ) +
    labs(
      color = analysis$qc_config$pcaMapping$color,
      shape = analysis$qc_config$pcaMapping$shape
    ) +
    (if (is.null(analysis$qc_config$pcaMapping$label)) {
      NULL
    } else {
      geom_text_repel(aes(label = .data[[analysis$qc_config$pcaMapping$label]]),
        size = 4, hjust = 0.5, vjust = -0.5, alpha = 0.5
      )
    }) +
    scale_x_continuous(expand = c(0.5, 0)) +
    theme_bw() +
    # ggtitle(paste0(analysis$qc_config$analysis," PCA")) +
    (if (!is.null(analysis$qc_config$pcaMapping$path)) {
      geom_path(aes(linetype = .data[[analysis$qc_config$pcaMapping$path]]))
    }) +
    theme(legend.key.width = unit(1.2, "cm")) +
    (if (!is.null(analysis$qc_config$pcaMapping$ellipse)) {
      stat_ellipse(aes(color = .data[[analysis$qc_config$pcaMapping$ellipse]]), type = "norm", level = 0.67)
    }) +
    (if (!is.null(analysis$qc_config$pcaMapping$shape)) {
      scale_shape_manual(values = shapes)
    }) +
    theme(text = element_text(size = 10)) # , arrow=arrow(ends="last", type="closed", length=unit(0.1, "inches")))
}


#' Simple base for PCA plot
#' 
#' Calculates PCA on counts data. Returns a ggplot object with the PCA data
#' and associated metadata. No additional geoms are included, this is meant
#' to be a blank template. You do need to provide metadata here so it can be 
#' accessed with geoms you add to this output. 
#'
#' @param counts Expression data matrix. Colnames of counts should match rownames
#'  of `metadata`
#' @param metadata Dataframe of metadata to associate with counts data. Rownames
#'  of metadata should match colnames of `counts`. 
#' @param xpc Numeric, which PC to plot on the x axis. Default 1
#' @param ypc Numeric, which PC to plot on the y axis. Default 2
#' @param ntop Number of genes to include in PCA. Default 500
#'
#' @returns A ggplot object
#' @export
#' 
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom matrixStats rowVars
#' @importFrom rlang .data
#' 
#'
#' @examples 
#' \dontrun{
#' pcaPlot <- pcaPlotSimple(assay(assays(pbmc_subset)$vst), metadata = colData(pbmc_subset))
#' 
#' ## Manually adding geoms
#' pcaPlot + geom_point(aes(color = Individual, shape=Timepoint)) 
#' 
#' }
#' 
pcaPlotSimple <- function(counts, metadata, xpc = 1, ypc = 2, ntop = 500) {
  if (colnames(counts) != rownames(metadata)) {
    stop('Colnames of counts does not match rownames of metadata')
  }
  rv <- matrixStats::rowVars(counts)
  select <- order(rv, decreasing = TRUE)[seq_len(min(
    ntop,
    length(rv)
  ))]
  pca <- stats::prcomp(t(counts[select,]))
  d <- merge(pca$x, metadata, by.x = 'row.names', by.y = 'row.names')
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  percentVar <- round(100 * percentVar)
  ggplot(d, aes(x = .data[[paste0('PC',xpc)]], y = .data[[paste0('PC',ypc)]])) +
    labs(x = paste0('PC',xpc, ": ", percentVar[xpc], "% variance"), 
         y = paste0('PC',ypc, ": ", percentVar[ypc], "% variance"))
}