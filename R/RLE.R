plotRLE <- function(data, title){
  as_tibble(data, rownames = NA) %>% 
    pivot_longer(everything(), names_to = "Sample", values_to = "RLE") %>% 
    ggplot(aes(x=Sample,y=RLE)) +
    geom_hline(yintercept = 0, color = "red") + geom_violin(draw_quantiles = c(0.25,0.75), trim = TRUE, color = "lightgreen", alpha = 0.1) +
    geom_boxplot(alpha = 0) + 
    theme_bw() +
    theme(axis.text.x = element_text(vjust = 0.5,angle = 90),axis.title.x = element_blank(),aspect.ratio = 0.55) + 
    ggtitle(title) + 
    aes(fct_inorder(Sample))#+ scale_y_continuous(limits = c(-9,3.25), expand = c(0,0))
}