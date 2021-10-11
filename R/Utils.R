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
  args <- c(
    "-f", paste(cols, collapse = ","),
    "-d", paste0("'", sep, "'"), filenm
  )
  command <- "cut"
  command <- paste(c(shQuote(command), args), collapse = " ")
  tmp <- system(command, intern = TRUE)
  invisible(utils::read.table(pipe(command), sep = sep, ...))
}

#' Load sparse matrix from mat file
#'
#' Each row of matfile: barcode[,genome_bin_index:counts]
#'
#' @param matfile str
#' @param ... for Matrix::sparseMatrix usage
#'
#' @return sparseMatrix
#'
#' @export
load_sparse_matrix <- function(matfile, ...) {
  lines <- data.table::fread(input = matfile, sep = "\n", header = FALSE)
  row_col_value <- lapply(seq(nrow(lines)), FUN = function(j) {
    l <- as.character(lines[j, 1])
    a <- unlist(strsplit(x = l, split = ","))
    barcode <- a[1]
    col_value <- t(vapply(a[-1], FUN = function(i) {
      invisible(as.numeric(unlist(strsplit(x = i, split = ":"))))
    }, FUN.VALUE = c(0, 0)))
    r <- cbind(j, col_value)
  })
  r <- as.data.frame(do.call(what = rbind, args = row_col_value))
  invisible(Matrix::sparseMatrix(i = r[, 1], j = r[, 2], x = r[, 3], ...))
}

#' Remove chr in GRSeqnames
#'
#' Ref: ArchR .convertGRSeqnames
#' @param gr GenomeRanges
#' @return gr2 GenomeRanges simplified
#' @export
simplifyGRSeqnames <- function(gr) {
  gr2 <- GenomicRanges::GRanges(
    seqnames = gsub("chr", "", GenomeInfoDb::seqnames(gr)),
    ranges = IRanges::ranges(gr), strand = GenomicRanges::strand(gr)
  )
  S4Vectors::mcols(gr2) <- S4Vectors::mcols(gr)
  return(gr2)
}

#' get subset of GR
#'
#' Ref: ArchR .subsetSeqnamesGR
#' @param gr GenomeRanges
#' @param names chromosome
#' @return gr2 subset of gr whose seqnames limited in names
#' @export
subsetSeqnamesGR <- function(gr, names) {
  gr <- gr[which(as.character(GenomeInfoDb::seqnames(gr)) %in% names), ]
  GenomeInfoDb::seqlevels(gr) <- as.character(unique(GenomeInfoDb::seqnames(gr)))
  return(gr)
}

#' Suppress message
#' Ref: ArchR .suppressAll
suppressAll <- function(expr = NULL) {
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}
