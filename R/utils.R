#' @importFrom stringr str_replace_all str_wrap str_count

.wrap_underscore_strings_balance <- function(string, width) {
  stringr::str_replace_all(.str_wrap_balance(str_replace_all(string, "_", " "), width), " ", "_")
}

.str_wrap_balance <- function(string, width = 80, indent = 0, exdent = 0, USE.NAMES = FALSE) {
  out <- stringr::str_wrap(string, width, indent, exdent)
  vapply(out, function(string, width, indent, exdent) {
    wraps <- str_count(string, "\n")
    if (wraps > 0 && width > 1) {
      bwidth <- width
      repeat {
        bwidth <- bwidth - 1
        bstring <- stringr::str_wrap(string, bwidth, indent, exdent)
        bwraps <- stringr::str_count(bstring, "\n")
        if (bwraps > wraps || bwidth <= 1) break
        string <- bstring
      }
    }
    string
  }, character(1), width, indent, exdent, USE.NAMES = USE.NAMES)
}


#' Read an excel workbook with sheets to a list of tables
#'
#' Read an excel workbook with sheets to a list of tables
#'
#' @param filename File in format .xls or .xlsx
#' @param tibble BOOL, Return the data as a tibble instead of a data.frame
#' @param \dots Additional arguments passed to `openxlsx::read.xlsx()`
#'
#' @returns A list of data.frames (or tibbles)
#'
#' @import openxlsx
#' @export
#'
#' @examples
#' \dontrun{
#'   dge_results <- read_excel_allsheets(here('outputs/dge.xlsx'))
#' }
read_excel_allsheets <- function(filename,
                                 tibble = FALSE,
                                 ...) {
  # https://stackoverflow.com/a/12945838/15664425
  sheets <- openxlsx::getSheetNames(filename)
  x <- lapply(sheets, function(X) openxlsx::read.xlsx(filename, sheet = X, ...))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}