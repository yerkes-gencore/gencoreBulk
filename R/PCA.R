#' plot PCA from config
#' 
#' Returns plot of PCA using aesthetic labels specified in the analysis config
#' 
#' @rdname plotPCAFromConfig
#' @param analysis object with config options saved from the template format
#'    
#' @returns A ggplot
#'
#' @examples
#'   plotPCAFromConfig(analysis)
#'   
#' @import ggplot 
#' @import ggrepel
#' 
#' @export
plotPCAFromConfig <- function(analysis){
  .pcaPlotGKT(assays(analysis$dds)$vst,
             intgroup = names(colData(assays(analysis$dds)$vst)),
             xpc = 1, ypc = 2) +
    geom_point(aes(color = (if (is.null(analysis$config$pcaMapping$color)) NULL 
                            else .data[[analysis$config$pcaMapping$color]]),
                   shape = (if (is.null(analysis$config$pcaMapping$shape)) NULL 
                            else .data[[analysis$config$pcaMapping$shape]]),
                   text  = (if (is.null(analysis$config$pcaMapping$hover)) NULL 
                            else .data[[analysis$config$pcaMapping$hover]])),
               size = 5) +
    labs(color=analysis$config$pcaMapping$color,
         shape=analysis$config$pcaMapping$shape) +
    (if (is.null(analysis$config$pcaMapping$label)) NULL 
     else geom_text_repel(aes(label = .data[[analysis$config$pcaMapping$label]]),
                          size = 4, hjust = 0.5, vjust = -0.5, alpha=0.5)) +
    scale_x_continuous(expand = c(0.5,0)) +
    theme_bw() + 
    #ggtitle(paste0(analysis$config$analysis," PCA")) +
    (if (!is.null(analysis$config$pcaMapping$path)){
      geom_path(aes(linetype=.data[[analysis$config$pcaMapping$path]]))
    }) +
    theme(legend.key.width = unit(1.2, "cm")) +
    (if (!is.null(analysis$config$pcaMapping$ellipse)){
      stat_ellipse(aes(color=.data[[analysis$config$pcaMapping$ellipse]]), type="norm", level=0.67)
    })+
    theme(text = element_text(size=10)) #, arrow=arrow(ends="last", type="closed", length=unit(0.1, "inches")))
  
  
}

.pcaDataGKT <- function (object, intgroup = "condition", ntop = 500) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(pca$x, group = group, 
                  intgroup.df, name = colnames(object))
  attr(d, "percentVar") <- percentVar
  return(d)
}

.pcaPlotGKT <- function(object, intgroup = "condition", xpc = 1, ypc = 2, ntop = 500)
{
  d <- .pcaDataGKT(object, intgroup, ntop)
  percentVar <- round(100 * attr(d, "percentVar"))
  ggplot(d,aes_string(x = names(d)[xpc], y = names(d)[ypc])) + labs(x=paste0(names(d)[xpc],": ", percentVar[xpc], "% variance"), y=paste0(names(d)[ypc],": ", percentVar[ypc], "% variance"))
}
