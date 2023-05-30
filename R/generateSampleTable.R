#' Generate minimal sampletable for use in QC template
#' 
#' A convenience function to start a sample table with the minimum
#'  requirements to run the QC template. You may need to manually adjust the 
#'  output to add additional metadata or trim sampleIDs. 
#'
#' @param config QC config with alignmentDir and STARreadSuffix specified
#' @param fileID_pattern Pattern to pull for FileID
#' @param sampleID_pattern Pattern to pull for SampleID
#' 
#'
#' @returns A dataframe with columns FileID and SampleID
#' @export
#'
#' @examples
#' \dontrun{
#' sampletable <- generateSampleTable(analysis$qc_config)
#' }
generateSampleTable <- function(config, 
                                fileID_pattern = '([^/])/.+',
                                sampleID_pattern = '^p\\d+.(s\\d+).+'){
  file_list <- dir(config$alignmentDir,
              pattern = config$STARreadSuffix,
              recursive = TRUE)
  fid <- gsub(pattern = fileID_pattern,
              replacement = '\\1', 
              x = file_list)
  sid <- gsub(pattern = sampleID_pattern,
              replacement = '\\1',
              x = file_list)
  return(data.frame(FileID = fid, SampleID = sid))
}