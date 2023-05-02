#' Write fGSEA results to excel sheet
#' 
#' Take a list of fGSEA result tables returned from `runfgsea()`
#'  and write them to an excel worksheet, with one sheet for each result.
#'
#' @param results List of fGSEA results objects 
#' @param sheet_names List of names for sheets
#' @param output_name Name of file
#' @param outdir Location of file
#'
#' @return `NULL`
#' @export
#' 
#' @import openxlsx
#'
#' @examples
#' \dontrun{
#'  gsea_result_1 <- runfgsea(results$Group_FAC_vs_Cont,
#'   gmt.file,
#'   minSize = 1,
#'   breakdown_pathway_names = TRUE)
#'  
#'  gsea_result_2 <- runfgsea(results$Group_LiproxposFAC_vs_Cont,
#'   gmt.file,
#'   minSize = 1,
#'   breakdown_pathway_names = TRUE)
#' 
#'  gsea_results <- list('FAC vs Cont' = gsea_result_1, 'LiproxposFAC vs Cont' = gsea_result_2)
#'  writefGSEAResults(gsea_results)
#'}

writefGSEAResults <- function(results,
                              sheet_names = names(results),
                              output_name = paste0(analysis$analysis_config$analysis,
                                                   "_GSEA_results.xlsx"),
                              outdir = here::here('outputs')){
  outfile <- file.path(outdir, output_name)
  message(paste0('Writing results to ', outfile))
  wb <- openxlsx::createWorkbook(outfile)
  mapply(FUN=.addWorksheet_fGSEAres,
         result=results,
         sheet_name=sheet_names,
         MoreArgs = list(wb=wb)
  )
  saveWorkbook(wb, outfile, overwrite = TRUE)
}

.addWorksheet_fGSEAres <- function(wb,
                                   result,
                                   sheet_name){
  if (nchar(sheet_name)>31){
    sheet_name <- substr(sheet_name,1,31)
  }
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb,
                      sheet=sheet_name,
                      x=as.data.frame(result),
                      rowNames=FALSE)
  
}