.wrap_underscore_strings_balance <- function(string, width) {
  stringr::str_replace_all(.str_wrap_balance(str_replace_all(string,"_"," "),width)," ","_")
}

.str_wrap_balance <- function(string, width = 80, indent = 0, exdent = 0,USE.NAMES = FALSE) {
  out <- stringr::str_wrap(string, width, indent, exdent)
  vapply(out, function(string,width,indent,exdent) {
    wraps <- str_count(string,"\n")
    if(wraps > 0 && width > 1) {
      bwidth <- width
      repeat {
        bwidth <- bwidth - 1
        bstring <- stringr::str_wrap(string, bwidth, indent, exdent)
        bwraps <- stringr::str_count(bstring,"\n")
        if(bwraps > wraps || bwidth <= 1) break
        string <- bstring
      }
    }
    string
  },character(1),width,indent,exdent,USE.NAMES = USE.NAMES)
}