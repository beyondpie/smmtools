#' Generate Gene Matrix by directly mapping fragments onto genes.
#'
#' @importFrom rhdf5 h5createFile h5createDataset h5closeAll h5write
#' @import S4Vectors
#' @import GenomicRanges
#' @export
getGeneMatrix <- function(rawH5File, outdir, outfilenm,
                          genes = NULL, genenms = NULL,
                          barcodes = NULL,
                          genome = "mm10", sampleName = NULL,
                          excludeChr = c("chrM")) {
  tstart <- Sys.time()
  message(paste("Begin to run getGeneMatrix with genome:", genome))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(outdir, outfilenm)
  if (file.exists(outfile)) {
    message(paste(outfile, "exist, and remove it."))
    file.remove(outfile)
  }
  rhdf5::h5createFile(file = outfile)
  ## blacklist here, to be test if we need blacklist
  annotGenome <- getAnnotFromArchRData(tag = "genome", genome = genome)
  blacklist <- annotGenome$blacklist
  if (is.null(genes)){
    annotGenes <- getAnnotFromArchRData(tag = "gene", genome = genome)
    genes <- annotGenes$genes
  }
  if (is.null(genenms)) {
    genenms <- genes$symbol
  } else {
    elementMetadata(genes)$symbol <- genenms
  }
  
  chrnms <- as.character(genes@seqnames@values)
  chrnms <- chrnms[!(match(chrnms, excludeChr, nomatch = 0) > 0)]
  chrDFs <- lapply(chrnms, function(chr) {
    message(paste("current:", chr))
    g <- genes[genes@seqnames == chr]
    fragments <-
      getFragsOfAChrFromRawH5File(
        rawH5File = rawH5File,
        chr = chr,
        sampleName = sampleName,
        barcodes = barcodes
      )
    gf <- GRanges(seqnames = chr, ranges = fragments)
    black_ovs <- as.data.frame(GenomicRanges::findOverlaps(query = gf, subject = blacklist))
    if (nrow(black_ovs) > 1) {
      gf <- gf[-black_ovs$queryHits]
    }
    ovs <- as.data.frame(GenomicRanges::findOverlaps(query = gf, subject = g))
    if (nrow(ovs) < 1) {
      return(NULL)
    }
    df <- lapply(unique(ovs$subjectHits), function(i) {
      fid <- ovs$queryHits[ovs$subjectHits == i]
      f <- fragments[fid]
      frle <- as.data.frame(table(mcols(f)$RG))
      colnames(frle) <- c("values", "lengths")
      frle$values <- as.character(frle$values)
      rownames(frle) <- frle$values
      cols <- S4Vectors::match(frle$values, barcodes, nomatch = 0)
      cols <- cols[!(cols == 0)]
      values <- frle[barcodes[cols], "lengths"]
      r <- data.frame(
        i = rep(match(g[i]$symbol, genenms), length(cols)),
        j = cols,
        val = values
      )
      return(r)
    })
    df <- df[!sapply(df, is.null)]
    r <- do.call(rbind, df)
    return(r)
  })
  chrDFs <- chrDFs[!sapply(chrDFs, is.null)]
  wdf <- do.call(rbind, chrDFs)
  suppressALL(h5createDataset(
    file = outfile, dataset = "i", storage.mode = "integer",
    dims = c(nrow(wdf), 1), level = 0
  ))
  suppressALL(h5write(obj = wdf$i, file = outfile, name = "i"))

  suppressALL(h5createDataset(
    file = outfile, dataset = "j", storage.mode = "integer",
    dims = c(nrow(wdf), 1), level = 0
  ))
  suppressALL(h5write(obj = wdf$j, file = outfile, name = "j"))

  suppressALL(h5createDataset(
    file = outfile, dataset = "val", storage.mode = "integer",
    dims = c(nrow(wdf), 1), level = 0, fillValue = 0
  ))
  suppressALL(h5write(obj = wdf$val, file = outfile, name = "val"))

  suppressALL(h5write(obj = genenms, file = outfile, name = "gene"))
  suppressALL(h5write(obj = barcodes, file = outfile, name = "barcode"))
  h5closeAll()
  mat <- Matrix::sparseMatrix(
    i = wdf$i, j = wdf$j, x = wdf$val,
    dims = c(length(genes), length(barcodes))
  )
  rownames(mat) <- genenms
  colnames(mat) <- barcodes
  return(mat)
}

#' Load sparse matrix of gene matrix with both row and col names
#' @param geneMatrixH5File string, h5 file to read
#' @importFrom rhdf5 h5createFile h5closeAll 
#' @export
loadGeneMatrix <- function(geneMatrixH5File) {
  h5_barcode <- h5read(file = geneMatrixH5File, name = "barcode")
  h5_gene <- h5read(file = geneMatrixH5File, name = "gene")
  h5_i <- as.vector(h5read(file = geneMatrixH5File, name = "i"))
  h5_j <- as.vector(h5read(file = geneMatrixH5File, name = "j"))
  h5_val <- as.vector(h5read(file = geneMatrixH5File, name = "val"))

  h5_mat <- Matrix::sparseMatrix(i = h5_i, j = h5_j, x = h5_val,
                                 index1 = TRUE,
                                 dims = c(length(h5_gene), length(h5_barcode)))
  rownames(h5_mat) <- h5_gene
  colnames(h5_mat) <- h5_barcode
  h5closeAll()
  return(h5_mat)
}
