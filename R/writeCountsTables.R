#'  Write data to csvs
#'
#'  Takes the raw counts matrix and an analysis object with yaml config 
#'  options saved to write data to named file outputs.
#'  
#' @param raw_counts matrix of raw count data
#' @param analysis object with config options saved from the template format
#' @param outdir folder to store outputs in, will be generated if it doesn't exist
#' @param generate_GSEA_app_tables write files for compatibility with GSEA desktop app
#' @param write_sample_table write sample metadata to a csv
#'
#' @returns NULL
#'
#' @examples
#'  writeCountTables(raw_counts, analysis, generate_GSEA_app_tables = TRUE)
#'  
#' @export


writeCountTables <- function(raw_counts,
                             analysis,
                             outdir = 'outputs',
                             generate_GSEA_app_tables = FALSE,
                             write_sample_table = TRUE) {
  if (!dir.exists(outdir)){ dir.create("outputs") }
  ## raw counts w/ gene symbols
  write_csv(as.data.frame(raw_counts) %>% 
              rownames_to_column(var = "gene_id"),
            file = paste0("outputs/raw_count_",
                          analysis$config$reference, "_",
                          analysis$config$analysis,".csv"),
            col_names = TRUE)
  if (write_sample_table) {
    write_csv(analysis$sampleTable,
              file = paste0("outputs/sample_table_",
                            analysis$config$analysis,".csv"))
  }
  if (generate_GSEA_app_tables) {
    analysis$assayRlogForGSEA <- assay(analysis$rld)
    analysis$assayRlogForGSEA <-
      analysis$assayRlogForGSEA[rowMeans(analysis$assayRlogForGSEA)>0,]
    ## GSEA cls file
    analysis$clsLinesGroup <-
      c(paste0(c(length(analysis$rldDrop$Group),
                 length(unique(analysis$rldDrop$Group)),1),
               collapse = " "),
        paste0(c("#",unique(as.vector(analysis$rldDrop$Group))),
               collapse = " "),
        paste0(analysis$rldDrop$Group, collapse = " "))
    write_lines(analysis$clsLinesGroup, 
                file = paste0("outputs/Group_",
                              analysis$config$analysis,".cls"))
    write_tsv(data.frame(
      Name = str_remove(rownames(analysis$assayRlogForGSEA), "[.].*"),
      Description = "na", analysis$assayRlogForGSEA),
      file = paste0("outputs/rlog_forGSEA_", 
                    analysis$config$analysis, ".txt"))
  }
}