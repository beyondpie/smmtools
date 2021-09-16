#' Read file by columns
#'
#' Ref: https://www.gastonsanchez.com/visually-enforced/how-to/2012/06/23/Read-file-by-columns/
#'
#' @param filenm str
#' @param cols vector[Integer]
#' @param sep str
#' @param ... others for read.csv
#' 
#' @return data.frame from read.csv
#' 
#' @export
read_columns <- function(filenm, cols = c(1), sep = ",", ...) {
  args <- c("-f", paste(cols, collapse = ","),
           "-d", paste0("'", sep, "'"), filenm)
  command <- "cut"
  command <- paste(c(shQuote(command), args), collapse = " ")
  tmp <- system(command, intern = TRUE)
  invisible(utils::read.table(pipe(command), sep = sep,...))
}
