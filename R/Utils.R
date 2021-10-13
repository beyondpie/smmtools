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

#' tile ChromeSize
#' Ref: ArchR
#' @param chromSizes GenomeRanges
#' @param nChunk integer, default is 3
#' @return tileChromSizes GenomeRanges
#' @export
tileChrom <- function(chromSizes, nChunk = 3) {
  tileChromSizes <- unlist(GenomicRanges::tile(chromSizes, nChunk))
  S4Vectors::mcols(tileChromSizes)$chunkName <- paste0(
    GenomeInfoDb::seqnames(tileChromSizes), "#chunk",
    seq_along(tileChromSizes)
  )
  return(tileChromSizes)
}

#' Get fragments from a given chromsome in rawH5File.
#' Ref: ArchR .getFragsFromArrow
#' @return fragments IRanges
#' @export
getFragsOfAChromFromRawH5File <- function(rawH5File, chr="chr1", sampleName=NULL, nChunk=3) {
  o <- rhdf5::h5closeAll()
  nFrags <- 0
  fragList <- list()
  barcodeValueList <- list()
  barcodeList <- list()
  for(i in 1:nChunk) {
    barcodeValueList[[i]]  <- fastH5Read(file = rawH5File,
                             name = paste0("Fragments/", chr, "#chunk", i, "/BarcodeLength"),
                             method = "fast")
    nFrags <- sum(barcodeValueList[[i]]) + nFrags
    fragList[[i]] <- fastH5Read(file = rawH5File,
                                name = paste0("Fragments/", chr, "#chunk", i, "/Ranges"),
                                method = "fast")
    barcodeList[[i]] <- fastH5Read(file = rawH5File,
                                   name = paste0("Fragments/", chr, "#chunk", i, "/BarcodeValue"),
                                   method = "fast")
  }
  if(nFrags == 0) {
    output <- IRanges::IRanges(start = 1, end = 1)
    S4Vectors::mcols(output)$RG <- c("tmp")
    output <- output[-1, ]
    return(output)
  }
  frags <- do.call(what = rbind, args = fragList)
  output <- IRanges::IRanges(start = frags[,1], end = frags[,2])
  barcodes <- do.call("c", barcodeList)
  barcodeValue <- do.call("c", barcodeValueList)
  S4Vectors::mcols(output)$RG <- Rle(values = paste0(sampleName, "#", barcodes),
                                     lengths = barcodeValue)
  return(output)
}


#' Fast read H5 file
#' Ref: ArchR .h5read
#' @return results return by h5read
#' @export
fastH5Read <- function(
  file = NULL,
  name = NULL,
  method = "fast",
  index = NULL,
  start = NULL,
  block = NULL,
  count = NULL
  ){

  if(tolower(method) == "fast" & is.null(index) & is.null(start) & is.null(block) & is.null(count)){
    fid <- rhdf5::H5Fopen(file)
    dapl <- rhdf5::H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count)
  }
  o <- h5closeAll()
  return(res)
}
