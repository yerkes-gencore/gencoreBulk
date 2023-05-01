#' Generate volcano plot from DESeq2 results object
#'
#' Generate volcano plot from DESeq2 results object. Labels can be custom set,
#' or specified as the top N genes by adjusted p value. 
#'
#' @param result DESeqResults object
#' @param labels If integer, the top N genes (by adjusted p value) to label.
#'  If a character vector, labels those genes regardless of values.
#'  
#' @import EnhancedVolcano
#' @inheritParams EnhancedVolcano::EnhancedVolcano
#' 
#' @inheritDotParams EnhancedVolcano::EnhancedVolcano 
#' 
#' @inherit EnhancedVolcano::EnhancedVolcano return 
#' 
#' @export
generateVolcanoPlot <- function(result,
                                labels = 20,
                                FCcutoff = log2(1.3),
                                pCutoff = analysis$analysis_config$alpha,
                                title=NULL,
                                caption=NULL,
                                subtitle=NULL,
                                ...){
  volData <- result[!is.na(result$padj),]
  if (is.numeric(labels)) {
    volData <- volData[order(volData$padj),]
    labels <- rownames(volData[1:labels,])
  }
  volplot <- EnhancedVolcano(volData,
                             x = 'log2FoldChange',
                             y = 'padj',
                             lab = rownames(volData),
                             selectLab = labels,
                             drawConnectors = TRUE,
                             colConnectors = "lightgrey",
                             pCutoff = alpha,
                             FCcutoff = FCcutoff,
                             title = title,
                             caption = caption,
                             subtitle = subtitle,
                             ...)
  return(volplot)
}