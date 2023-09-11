#' GSEA dotplot
#'
#' Make dotplot from fGSEA result, showing top n pathways
#'
#' @rdname gseaDotplot_joint
#'
#' @param gsea_results A table with multiple GSEA results from `fgsea::fgseaSimple()` 
#'  combined by `combine_GSEA_results()`
#' @param pathway_order Order of pathways to plot
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
                              pathway_order = NULL){
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
  
  
  
  ggplot(gsea_results) + 
    geom_point(aes(x=.data$ID, y=.data$pathway, size=-log(.data$pval), color=.data$NES)) +
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
         title="GSEA pathway enrichments") +
    scale_radius(name="NOM p-val", 
                 range=c(1,8),
                 breaks=-log10(c(0.1,0.01,0.001,0.0001)),
                 labels=c(0.1,0.01,0.001,0.0001))
}