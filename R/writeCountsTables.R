#'  Write data to csvs
#'
#'  Takes the raw counts matrix and an analysis object with yaml config
#'  options saved to write data to named file outputs.
#'
#' @param analysis object with config options saved from the template format
#' @param outdir folder to store outputs in, will be generated if it doesn't exist
#' @param generate_GSEA_app_tables write files for compatibility with GSEA desktop app
#' @param write_sample_table write sample metadata to a csv
#'
#' @returns NULL
#'
#' @examples
#' \dontrun{
#' writeCountTables(raw_counts, analysis, generate_GSEA_app_tables = TRUE)
#' }
#' @importFrom here here
#' @importFrom readr write_csv write_tsv write_lines
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove
#'
#' @export


writeCountTables <- function(analysis,
                             outdir = here::here("outputs"),
                             generate_GSEA_app_tables = FALSE,
                             write_sample_table = TRUE) {
  if (!dir.exists(outdir)) {
    dir.create("outputs")
  }
  ## raw counts w/ gene symbols
  assays(analysis$dds)$counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    readr::write_csv(
      file = here::here(paste0(
        outdir, "/raw_count_",
        analysis$config$analysis, ".csv"
      )),
      col_names = TRUE
    )
  if (write_sample_table) {
    readr::write_csv(analysis$sampleTable,
      file = here::here(paste0(
        outdir, "/sample_table_",
        analysis$config$analysis, ".csv"
      ))
    )
  }
  if (generate_GSEA_app_tables) {
    assayRlogForGSEA <- assays(analysis$dds)$rld
    assayRlogForGSEA <-
      assayRlogForGSEA[rowMeans(assayRlogForGSEA) > 0, ]
    ## GSEA cls file
    analysis$clsLinesGroup <-
      c(
        paste0(
          c(
            length(assays(analysis$dds)$rld$Group),
            length(unique(assays(analysis$dds)$rld$Group)), 1
          ),
          collapse = " "
        ),
        paste0(c("#", unique(as.vector(assays(analysis$dds)$rld$Group))),
          collapse = " "
        ),
        paste0(assays(analysis$dds)$rld$Group, collapse = " ")
      )
    readr::write_lines(analysis$clsLinesGroup,
      file = here::here(paste0(
        outdir, "/Group_",
        analysis$config$analysis, ".cls"
      ))
    )
    readr::write_tsv(
      data.frame(
        Name = stringr::str_remove(rownames(assayRlogForGSEA), "[.].*"),
        Description = "na", assayRlogForGSEA
      ),
      file = here::here(paste0(
        outdir, "/rlog_forGSEA_",
        analysis$config$analysis, ".txt"
      ))
    )
  }
}
