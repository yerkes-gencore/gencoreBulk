#' Parse STAR readsPerGene.out files
#' 
#' Returns the counts per gene and read count bins for a list of STAR output files
#' 
#' @param file.paths Named vector of file paths to STAR output files, where names
#'    correspond to sample names
#' @param library.type String of column to pull from in the output file. Options
#'  should be 'unstranded', 'sense', or 'antisense'
#'    
#' @returns List of 2 data frames, one containing read mapping bins, the other 
#'  read per gene counts
#'
#' @examples
#'    readfiles <- sapply(
#'      analysis$samplefileIDs,
#'      function(sid) {
#'        paste0(dir("my_STAR_output_dir",
#'                   pattern = sid, full.names = TRUE),
#'               "/", sid, 'readsPerGene.out.tab')})
#'               
#'    outs <- readCountFiles(readfiles, 'unstranded')
#' 
#' @export
parseReadPerGeneFiles <- function(file.paths, library.type = 'unstranded'){
  raw_out <- sapply(readfiles, 
                    function(x) {
                      read_tsv(x, 
                               col_names = c("gene_id",
                                             "unstranded_count",
                                             "sense_count",
                                             "antisense_count"),
                               col_types = c(gene_id = col_character(),
                                             unstranded_count = col_double(),
                                             sense_count = col_double(),
                                             antisense_count = col_double())) %>% 
                        select(gene_id, contains(library.type))},
                    simplify = FALSE,
                    USE.NAMES = TRUE)
  
  map_bins <- sapply(raw_out,
                     function(x) {x[c(1:4),][[paste0(library.type, '_count')]]},
                     USE.NAMES = TRUE)
  rownames(map_bins) <- unlist(raw_out[[1]][1:4,1])
  
  read_counts <- sapply(raw_out,
                        function(x) {x[-c(1:4),][[paste0(library.type, '_count')]]},
                        USE.NAMES = TRUE)
  rownames(read_counts) <- unlist(genes <- raw_out[[1]][-c(1:4),1])
  
  map_bins <- rbind(map_bins,"N_identified" = colSums2(read_counts))
  
  return(list(map_bins = map_bins, read_counts = read_counts))
}


#' Mapping bins plot
#'  
#' Plot values of mapping bins for each sample. Takes the output of `parseReadPerGeneFiles()`
#'
#' @param mapBins Data frame of 5 mapping bin categories for each sample
#' @param title   Optional title for the plot
#'  
#' @returns ggplot Plot
#' @export

mappingBinsPlot <- function(mapBins, title=''){
  data <- rownames_to_column(as_tibble(mapBins, rownames = NA),
                             var = "map_result") %>% 
    pivot_longer(!map_result, names_to = "SampleID", values_to = "count") 
  data$map_result <- unlist(lapply(data$map_result,
                                   function(x)toTitleCase(str_split(x, '_')[[1]][2])))
  data$map_result <- factor(data$map_result,
                            levels = c("Unmapped","Multimapping",
                                       "noFeature","Ambiguous","Identified"))
  
  ggplot(data, aes(x = SampleID, y = count, fill = map_result)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c('Unmapped' = 'grey',
                                 'Multimapping' = '#DB890E', 
                                 'noFeature' = '#C60A19', 
                                 'Ambiguous' = '#C4B31C', 
                                 'Identified' = '#209964' )) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    guides(fill=guide_legend(title="Map Result")) +
    aes(fct_inorder(SampleID)) + 
    labs(x="Sample", y='Number of reads') +
    ggtitle(title)
}


