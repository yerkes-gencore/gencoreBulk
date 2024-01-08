#' GSEA dotplot
#'
#' Make dotplot from fGSEA result, showing top n pathways
#'
#' @rdname gseaDotplot_joint
#'
#' @param gsea_results A table with multiple GSEA results from `fgsea::fgseaSimple()` 
#'  combined by `combine_GSEA_results()`
#' @param pathway_order Order of pathways to plot
#' @param x_order Order of comparisons on X axis
#' @param significance A vector of values to indicate significance with asterisks.
#'  Each subsequent value will add an extra asterisk. E.g. `c(0.05, 0.01)` will 
#'  give one asteriks to values below 0.05 and two asterisks to values below 0.01.
#'  If this is null, no asterisks will be plotted.
#'  @param p_val_col Column to use for significance values. Default 'pval'.
#'
#' @return A ggplot object
#' @export
#'
#' @importFrom rlang .data
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' combine_GSEA_results <- function(gsea_results,
#'                                  pathways){
#'  gsea_results <- lapply(gsea_results, function(x){x %>% filter(pathway %in% pathways)})
#'  gsea_results <- data.table::rbindlist(gsea_results, idcol='ID')
#' }
#'
#' pathways = c('REACTOME_SIGNALING_BY_GPCR', 'REACTOME_SIGNALING_BY_NTRKS')
#' gsea_results <- list('FAC vs Cont' = gsea_result, 'LiproxposFAC vs Cont' = gsea_result_2)
#' joint_GSEA_results <- combine_GSEA_results(gsea_results, pathways)
#' gseaDotplot_joint(joint_GSEA_results)
#' }
gseaDotplot_joint <- function(gsea_results,
                              pathway_order = NULL,
                              x_order = NULL,
                              significance = c(0.05, 0.01, 0.001),
                              p_val_col = 'pval'){
  if (!is.null(pathway_order)) {
    if (all(pathway_order %in% unique(gsea_results$pathway))){
      pathway_order <- order(factor(gsea_results$pathway, levels = pathway_order))
      gsea_results$pathway <- .wrap_underscore_strings_balance(gsea_results$pathway,36)
      ## reordering
      gsea_results <- gsea_results[pathway_order,]
      gsea_results$pathway <- factor(gsea_results$pathway, levels = unique(gsea_results$pathway))
    } else {
      warning('pathways specified in pathway_order not found, defaulting to arbitrary order')
    }
  } else {
    gsea_results$pathway <- .wrap_underscore_strings_balance(gsea_results$pathway,36)
  }
  
  if (!is.null(x_order)) {
    if (all(x_order %in% unique(gsea_results$ID))){
      gsea_results$ID <- factor(gsea_results$ID, levels = x_order)
    } else {
      warning('Specified x_order not all found in data, defaulting to arbitrary order')
    }
  }
  
  gsea_results$label <- NA
  caption <- ''
  if (!is.null(significance)) {
    if (is.numeric(significance)) {
      label <- '*'
      for (cutoff in significance) {
        gsea_results$label <- ifelse(gsea_results[[p_val_col]] < cutoff, label, gsea_results$label)
        caption <- paste(caption, label, '<', cutoff, ';', sep = ' ')
        label <- paste0(label, '*')
      }
      # gsea_results$label[is.numeric(gsea_results$label)] <- NA
    } else {
      stop('Significance argument should be a numeric vector')
    }
  } 
  ggplot(gsea_results, aes(x=.data$ID, y=.data$pathway, size=-log(.data[[p_val_col]]),
                           color=.data$NES, label = .data$label)) + 
    geom_point() +
    scale_color_gradient2(low="blue",
                          mid="white",
                          high="red",
                          midpoint=0,
                          breaks=c(-2,-1,0,1,2),
                          limits=c(min(gsea_results$NES,-1),
                                   max(gsea_results$NES,1))) +
    theme_classic(11) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_line(colour = "grey92"),#element_blank(),
          panel.grid.minor.y = element_line(colour = "grey92")) +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1)) +
    labs(x="Comparison",
         y="Gene set", 
         color = "Normalized\nenrichment\nscore",
         size="Nom p-val",
         title="GSEA pathway enrichments",
         caption = caption) +
    scale_radius(name="NOM p-val", 
                 range=c(1,8),
                 breaks=-log10(c(0.1,0.01,0.001,0.0001)),
                 labels=c(0.1,0.01,0.001,0.0001)) +
    geom_text(na.rm = TRUE, color = 'black', size = 3)
}