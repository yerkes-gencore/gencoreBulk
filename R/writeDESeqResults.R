#' Write DESeqResults to excel sheet
#' 
#' Take a list of DESeqResults returned from `extractFilteredDESeqResults()`
#'  and write them to an excel worksheet, with one sheet for each result.
#'
#' @param results List of DESeqResult objects 
#' @param sheet_names List of names for sheets
#' @param output_name Name of file
#' @param outdir Location of file
#' @param padj_col Name of column with adjusted p values. For filtering and sorting
#' @param drop_NA Whether to exclude results with an adjusted p-value of `NA` 
#'
#' @return `NULL`
#' @export
#' 
#' @import openxlsx
#'
#' @examples
#' \dontrun{
#'  my_comparisons <- list('Group_FAC_vs_Cont',
#'  'Group_LiproxposFAC_vs_Cont',
#'  'Group_MTposFAC_vs_Cont')
#' 
#'  ## groups to include for greg's filter function, specific to each comparison
#'  ## should be levels of some group in colData(analysis$dds)
#'  filter_on_list <- list(c('FAC', 'Cont'),
#'                        c('LiproxposFAC', 'Cont'),
#'                        c('MTposFAC', 'Cont'))
#' 
#'  ## Get DESeq2 results for each comparison of interest
#'  ## intgroup specifies which metadata column your filter_on_list refers to
#'  results <- mapply(extractFilteredDESeqResults,
#'                    comparison = my_comparisons,
#'                    filt_groups = filter_on_list, 
#'                    intgroup = 'Group',
#'                    USE.NAMES = FALSE)
#'  names(results) <- lapply(results, function(x) x@metadata$contrast)
#'    writeDESeqResults(results)
#' }
writeDESeqResults <- function(results,
                              sheet_names = names(results),
                              output_name = "DEG_results.xlsx",
                              outdir = here::here('outputs'),
                              padj_col = 'padj',
                              drop_NA = TRUE){
  outfile <- file.path(outdir, output_name)
  message(paste0('Writing results to ', outfile))
  wb <- openxlsx::createWorkbook('ENPRC Gencore')
  mapply(FUN = .addWorksheet_DESeqres,
         result=results,
         sheet_name=sheet_names,
         padj_col = padj_col,
         MoreArgs = list(wb=wb,
                         drop_NA=drop_NA)
  )
  openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)
}

.addWorksheet_DESeqres <- function(wb,
                                   result,
                                   sheet_name,
                                   padj_col,
                                   drop_NA){
  if (nchar(sheet_name)>31){
    sheet_name <- substr(sheet_name,1,31)
  }
  if (drop_NA) {
    result <- result[!is.na(result[[padj_col]]),] 
  }
  result <- result[order(result[[padj_col]]),]
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb,
            sheet=sheet_name,
            x=as.data.frame(result),
            rowNames=TRUE)
}
