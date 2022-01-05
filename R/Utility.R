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
readColumnsFromTextFile <- function(filenm, cols = c(1), sep = ",", ...) {
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
loadSparseMatrixFromMatTextFile <- function(matfile, ...) {
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
#' @param verbose bool, print logs when set is TRUE, default is TRUE.
#' @return fragments IRanges
#' @importFrom rhdf5 h5ls h5read 
#' @importFrom S4Vectors mcols mcols<- match
#' @export
getFragsOfAChrFromRawH5File <- function(rawH5File, chr="chr1", sampleName=NULL, barcodes = NULL,
                                        verbose = TRUE) {
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
    nFrags <- sum(barcodeValueList[[i]]) + nFrags
    fragList[[i]] <- h5read(file = rawH5File,
                            name = paste0("Fragments/", g, "/Ranges"))
    barcodeList[[i]] <- h5read(file = rawH5File,
                               name = paste0("Fragments/", g, "/BarcodeValue"))
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
  if (verbose) {
    message(paste("Get", length(output), " fragments for", chr, "of sample", sampleName))
  }
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

#' Sampling cells based on their sequencing depths without replace.
#'
#' @param bmat Matrix, cell by feature
#' @param n integer
#' @return vector of integers, index for the sampled cells
#' with length of min(n, ncells)
#' @export
sampleBasedOnDepth <- function(bmat, n) {
  depths <- log(Matrix::rowSums(bmat) + 1, 10)
  dens <- stats::density(x = depths, bw = "nrd", adjust = 1)
  samplingProb <- 1 / (stats::approx(x = dens$x, y = dens$y, xout = depths)$y + .Machine$double.eps)
  idx <- sort(sample(x = seq_along(depths), size = min(n, nrow(bmat)),
                     prob = samplingProb, replace = FALSE))
  return(idx)
}

#' NormOVE from SnapATAC
#' @return matrix
#' @export
getNormOVE <- function(p1, p2){
  pp = tcrossprod(p1, p2);
	ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +
    matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
	ee = pp/(ss - pp)
	return(ee)	
}

#' Get eigen decomposion
#' @param M matrix
#' @n_eigs integer, number of eigns
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

#' Get i from sparseMatrix
#' i index starts froms 1.
#' Ref: https://stackoverflow.com/questions/21099612/extract-i-and-j-from-a-sparse-matrix
#' @param sm sparseMatrix
#' @return vector of integer
#' @import Matrix
#' @export
getiFromSparseMatrix <- function(sm) {
  t <- as(sm, "TsparseMatrix")
  return(t@i+1)
}

#' Get j from sparseMatrix
#' j index starts froms 1.
#' Ref: https://stackoverflow.com/questions/21099612/extract-i-and-j-from-a-sparse-matrix
#' @param sm sparseMatrix
#' @return vector of integer
#' @import Matrix
#' @export
getjFromSparseMatrix <- function(sm) {
  t <- as(sm, "TsparseMatrix")
  return(t@j+1)
}
#' Get x/values from sparseMatrix
#' 
#' Ref: https://stackoverflow.com/questions/21099612/extract-i-and-j-from-a-sparse-matrix
#' @param sm sparseMatrix
#' @return numeric vector
#' @import Matrix
#' @export
getvalFromSparseMatrix <- function(sm) {
  t <- as(sm, "TsparseMatrix")
  return(t@x)
}

#' Transpose a sparseMatrix
#' A method for the generics function "t"
#' DEPRECATED: we can directly use Matrix::t.
#' @param sm sparseMatrix
#' @return sparseMatrix
t.sparseMatrix <- function(sm) {
  smi <- getiFromSparseMatrix(sm)
  smj <- getjFromSparseMatrix(sm)
  smx <- getvalFromSparseMatrix(sm)
  r <- Matrix::sparseMatrix(i = smj, j = smi, x = smx,
                            dims = c(ncol(sm), nrow(sm)), index1 = TRUE)
  return(r)
}

#' Get index in the subject vector for the query vector
#' @param query vector, same type as in subject
#' @param subject vector, same type as in query
#' @return integer vector, index in subject for the query vector
#' @export
getIndex <- function(query, subject) {
  t <- match(query, subject, nomatch = 0)
  r <- t[t>0]
  if(length(r) < 1){
    warning("Query has no element in subject. Empty vector will be returned.")
  }
  return(t[t > 0])
}

#' Plot correlations between adjacent two principle components.
#' Refer: SnapATAC plotDimReductPW
#' @param smat cell by feature
#' @param outpdf string, output of pdf file, default is NULL.
#' @param ndim integer, max dims to plot, default is 50
#' @param point_size double, default 0.5
#' @param point_color string, default "grey"
#' @param point_shape integer, default 19
#' @param point_alpha double, default 0.5
#' @param nsample integer, max number of random cells to use, default 10000
#' @param height double, height of pdf, default 7
#' @param width, double, width of pdf, default 7
#' @param ... parameters, used for plot.
#' @return value returned by graphics::par()
#' @export
plotDimReductPairwise <- function(smat, outpdf = NULL,ndim = 50, point_size = 0.5,
                                  point_color = "grey", point_shape = 19,
                                  point_alpha = 0.5, nsample = 10000,
                                  height = 7, width = 7, ...) {
  ## only less than 50 dims will be used.
  p <- min(ncol(smat), ndim, 50)
  if(!is.null(outpdf)) {
    outdir <- dirname(outpdf)
    if(!dir.exists(outdir)) {
      message(paste(outdir, "does not exist and will be created."))
      dir.create(outdir)
    }
    if(file.exists(outpdf)) {
      message(paste(outpdf, "exists and will be re-created."))
      file.remove(outpdf)
    }
    pdf(outpdf, width = width, height = height)
  }
  op <- par(mfrow = c(5, 5), oma = c(3, 3, 1,1) + 0.2,
            mar = c(0, 0, 1,1) + 0.2)
  pca_plot <- split(seq(p), ceiling(seq(p)/2))
  if( (p %% 2 ) == 1) {
    pca_plot <- pca_plot[1:length(pca_plot)-1]
  }
  for(x in pca_plot) {
    rows <- sort(sample(seq(nrow(smat)), size = min(nrow(smat), nsample)))
    plot(x = smat[rows, x[1]], y = smat[rows, x[2]],
         cex = point_size, col = scales::alpha(point_color, point_alpha),
         mtext(paste(paste("PCs", x[1]), x[2], sep = " vs "), side = 3),
         xlab = "", ylab = "", yaxt = "n", xaxt = "n"
         )       
         ## xlab = paste0("Dim", x[1]), ylab = paste0("Dim", x[2]), ...)
  }
  if(!is.null(outpdf)) {
    dev.off()
  }
  invisible(graphics::par(mfrow = c(1,1)))
}

#' Plot correlations between adjacent two principle components.
#' Refer: SnapATAC plotDimReductElbow
#' @param sdev numeric vector, standard deviation for PCs.
#' @param outpdf string, output of pdf file, default is NULL.
#' @param point_size double, default 1.5
#' @param point_color string, default "red"
#' @param point_shape integer, default 19
#' @param point_alpha double, default 1
#' @param height double, height of pdf, default 7
#' @param width, double, width of pdf, default 7
#' @param ... parameters, used for plot.
#' @return None
#' @export
plotDimReductElbow <- function(sdev, outpdf = NULL, point_size = 1.5,
                               point_shape = 19, point_color = "red",
                               point_alpha = 1, height= 7, width = 7, ...) {
  if(!is.null(outpdf)) {
    cleanOutfile(outpdf)
    pdf(outpdf, width = width, height = height)
  }
  plot(x = seq(length(sdev)), y = sdev, cex = point_size, pch = point_shape,
       col = scales::alpha(point_color, point_alpha), xlab = "Principle Components",
       ylab = "Standard Deviation of PCs.", ...)
  if(!is.null(outpdf)) {
    dev.off()
  }
}

#' Create outdir if not exist, and remove the old file if needed.
#' @param outf string
#' @param remove bool, default is TRUE
#' @return None
#' @export
cleanOutfile <- function(outf, remove = TRUE){
  outdir <- dirname(outf)
  if(!dir.exists(outdir)) {
    message(paste(outdir, "does not exist and will be created."))
    dir.create(outdir)
  }
  if(file.exists(outf) & remove) {
    message(paste(outf, "exists and will be removed."))
    file.remove(outf)
  }
}
