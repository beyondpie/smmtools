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
#' @export
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
#' @param barcodes vector of strings, make sure it's not a data.frame with one column.
#' @return fragments IRanges
#' @importFrom rhdf5 h5ls h5read h5closeAll
#' @importFrom S4Vectors mcols mcols<- match
#' @export
getFragsOfAChrFromRawH5File <- function(rawH5File, chr="chr1", sampleName=NULL, barcodes = NULL) {
  o <- rhdf5::h5closeAll()
  nFrags <- 0
  fragList <- list()
  barcodeValueList <- list()
  barcodeList <- list()
  groups <- h5ls(file = rawH5File)
  groups <- groups[groups$group == "/Fragments" & groups$otype == "H5I_GROUP", "name"]
  groups <- groups[grep(pattern = paste0(chr, "#"), x = groups)]
  for(i in seq_along(groups)) {
    g <- groups[i]
    barcodeValueList[[i]]  <- h5read(file = rawH5File,
                             name = paste0("Fragments/", g, "/BarcodeLength"))
    h5closeAll()
    nFrags <- sum(barcodeValueList[[i]]) + nFrags
    fragList[[i]] <- h5read(file = rawH5File,
                            name = paste0("Fragments/", g, "/Ranges"))
    h5closeAll()

    barcodeList[[i]] <- h5read(file = rawH5File,
                               name = paste0("Fragments/", g, "/BarcodeValue"))
    h5closeAll()

  }
  if(nFrags == 0) {
    message(paste("No fragments for", chr, "of sample", sampleName))
    output <- IRanges::IRanges(start = 1, end = 1)
    mcols(output)$RG <- c("tmp")
    output <- output[-1, ]
    return(output)
  }
  frags <- do.call(what = rbind, args = fragList)
  output <- IRanges::IRanges(start = frags[,1], width = frags[,2])
  message(paste("Get", length(output), " fragments for", chr, "of sample", sampleName))
  barcodesFromH5 <- do.call("c", barcodeList)
  barcodeValue <- do.call("c", barcodeValueList)
  ## mcols(output)$RG <- S4Vectors::Rle(values = paste0(sampleName, "#", barcodesFromH5),
                                                ## lengths = barcodeValue)
  mcols(output)$RG <- S4Vectors::Rle(values = barcodesFromH5, lengths = barcodeValue)
  if(!is.null(barcodes)) {
    message(paste("Barcodes is not empty, and will select the fragments for these barcodes."))
    r <- output[BiocGenerics::which( match(mcols(output)$RG, barcodes, nomatch = 0) > 0 )]
    message(paste("After selection/filtering, we have", length(r),
                  "fragments left for", chr, "of sample", sampleName))
  } else {
    r <- output
  }
  h5closeAll()
  return(r)
}


#' Fast read H5 file
#' Ref: ArchR .h5read
#' @importFrom rhdf5 h5read h5closeAll H5Fopen H5Pcreate
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
    fid <- H5Fopen(file)
    dapl <- H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count)
  }
  o <- h5closeAll()
  return(res)
}

#' Remove sampleName in barcodes if they have
#' @param barcodes vector of strings
#' @param sep string, symbol used to separate sampleName and barcode, default is '#'
#' @export
removeSampleName <- function(barcodes, sep = "#") {
  if(is.data.frame(barcodes)) {
    message("Input Barcodes is data.frame, will treat the first column as barcodes.")
    barcodes <- barcodes[,1]
  }
  if(grepl(sep, barcodes[1], fixed = TRUE)) {
    r <- gsub(pattern = paste0(".+",sep), "", barcodes)
  } else {
    r <- barcodes
  }
  return(r)
}

#' @param bmat Matrix, cell by feature
#' @export
sampleBasedOnDepth <- function(bmat, n) {
  depths <- log(Matrix::rowSums(bmat) + 1, 10)
  dens <- stats::density(x = depths, bw = "nrd", adjust = 1)
  samplingProb <- 1 / (stats::approx(x = dens$x, y = dens$y, xout = depths)$y + .Machine$double.eps)
  idx <- sort(sample(x = seq_along(depths), size = min(n, nrow(bmat)), prob = samplingProb))
  return(idx)
}

#' @return matrix
#' @export
getNormOVE <- function(p1, p2){
  pp = tcrossprod(p1, p2);
	ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +
    matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
	ee = pp/(ss - pp)
	return(ee)	
}

#' @export
eig_decomp <- function(M, n_eigs) {
  sym <- Matrix::isSymmetric(M)
	n <- nrow(M)
	f <- function(x, A = NULL) {
    as.matrix(A %*% x)
  }
	wh <- if (sym) 'LA' else 'LM'
	#constraints: n >= ncv > nev
	ar <- igraph::arpack(f, extra = M, sym = sym, options = list(
		which = wh, n = n, ncv = min(n, 4*n_eigs), nev = n_eigs + 1))
	if (!sym) {
		ar$vectors <- Re(ar$vectors)
		ar$values  <- Re(ar$values)
	}
	if (length(dim(ar$vectors)) == 0L) {
		ar$vectors <- matrix(ar$vectors, ncol = 1L)
  }
	return(ar)
}
