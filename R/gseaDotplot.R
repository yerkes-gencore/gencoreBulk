#' GSEA dotplot
#'
#' Make dotplot from fGSEA result, showing top n pathways
#'
#' @rdname gseaDotplot
#' @param result A table with GSEA results from `fgsea::fgseaSimple()`
#' @param min_size Minimum size of gene set to plot
#' @param filter_source Character vector of pathway sources to filter for
#' @param signif_only If TRUE, only plot results with padj < `sig_cutoff`
#' @param top_n Show the top N results sorted by padj
#' @param sig_cutoff Threshold for significance, default to 0.05
#'
#' @return A ggplot object
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stringr str_replace_all str_wrap str_count
#' @importFrom utils head
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'   res <- runfgsea(result, gmt.file)
#'   gseaDotplot(res, filter_source = 'HALLMARK')
#' }
gseaDotplot_single <- function(result,
                        filter_source = NULL,
                        signif_only = FALSE,
                        top_n=20, 
                        min_size = 5,
                        sig_cutoff = 0.05){
  result <- result %>% 
    arrange(.data$pval, desc(.data$size)) %>% 
    mutate(perc = 100 * lengths(.data$leadingEdge) / .data$size) %>%
    mutate(name=paste0(.wrap_underscore_strings_balance(.data$pathway,36), "\nn=", .data$size)) %>%
    filter(.data$size >= min_size)
  if (!is.null(filter_source)){
    result <- result %>% 
      filter(.data$source %in% filter_source)
  }
  if (signif_only){
    result <- result %>%
      filter(.data$padj < sig_cutoff)
  }
  toppaths <- rbind(utils::head(result, n=top_n))
  
  toppaths$name <- factor(toppaths$name)
  dotplot <- ggplot(toppaths) + 
    geom_point(aes(x=.data$perc,
                   y=.data$name, 
                   size=-log10(.data$pval),
                   color=.data$NES)) +
    scale_color_gradient2(low="blue",
                          mid="white",
                          high="red",
                          midpoint=0,
                          breaks=c(-2,-1,0,1,2),
                          limits=c(min(toppaths$NES,-1),
                                   max(toppaths$NES,1))) +
    theme_classic(11) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_line(colour = "grey92"),
          panel.grid.minor.y = element_line(colour = "grey92")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x="% of genes in leading edge", y="Gene set",
         color = "Normalized\nenrichment\nscore", 
         size="Nom p-val", title="Top enriched pathways", 
         caption='GSEA p-value calculations are not continuous,\n
         there may be several or many pathways with the same p-value\n
         n = number of genes in pathway') + 
    scale_radius(name="NOM p-val",
                 range=c(1,8),
                 breaks=-log10(c(0.5,0.1,0.01,0.001)),
                 limits=c(0,3),
                 labels=c(0.5,0.1,0.01,0.001)) +
    scale_y_discrete(limits=toppaths$name)
  return(dotplot)
}