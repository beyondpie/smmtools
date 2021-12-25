#' Generate Gene Matrix by directly mapping fragments onto genes.
#'
#' FIXME: this function should be wrong at somewhere
#' since it's not consistante with SnapATAC
#'
#' @param rawH5File string
#' @param outfile string
#' @param barcodes vector of string, neccesary
#' @param genome string, dfeault is "mm10"
#' @param genes GenomicRanges, default is NULL, will load automatically using data under the package
#' Otherwise supported by the user.
#' @param genenms vector of string, names of genes, default is NULL
#' If we used the genes by ourselves, and genenms is NULL, we use genes$symbol to get the genenms.
#' This might be error: for example, if we use rtracklayer::import to generate the genes, then
#' they use genes$name to store gene names. So we need to check this mannually.
#' @param sampleName string, default is NULL
#' @param excludeChr vector of string, default is c("chrM")
#' @param compressLevel integer, used for rhdf5, default is 9.
#' @return sparseMatrix, gene by cell
#' @importFrom rhdf5 h5createFile h5createDataset
#' @import S4Vectors
#' @import GenomicRanges
#' @export
getGeneMatrix <- function(rawH5File, outfile,
                          barcodes,
                          genes = NULL, genenms = NULL,
                          genome = "mm10", sampleName = NULL,
                          excludeChr = c("chrM"),
                          compressLevel = 9) {
  tstart <- Sys.time()
  message(paste("Begin to run getGeneMatrix with genome:", genome))
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
    fragments <-
      getFragsOfAChrFromRawH5File(
        rawH5File = rawH5File,
        chr = chr,
        sampleName = sampleName,
        barcodes = barcodes
      )
    gf <- GRanges(seqnames = chr, ranges = fragments)
    black_ovs <- as.data.frame(GenomicRanges::findOverlaps(query = gf, subject = blacklist))
    if (nrow(black_ovs) >= 1) {
      gf <- gf[-black_ovs$queryHits]
    }
    g <- genes[genes@seqnames == chr]
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
    if(is.null(df)) {
      message(paste("No barcodes match for", chr, "from sample", sampleName))
      message("This might be error, please double check the barcodes. Now continue.")
    }
    r <- do.call(rbind, df)
    return(r)
  })
  
  chrDFs <- chrDFs[!sapply(chrDFs, is.null)]
  
  if (is.null(chrDFs)) {
    message(paste("WARNINIG: Get empty gmat! Will not write it into:", outfile))
    message(paste("This might because of no barcodes matched with the one you given."))
    return(NULL)
  }
  
  wdf <- do.call(rbind, chrDFs)
  
  message(paste("Finish getGeneMatrix, and now write it into:", outfile))
  outdir <- dirname(outfile)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (file.exists(outfile)) {
    message(paste(outfile, "exist, and remove it."))
    file.remove(outfile)
  }
  
  rhdf5::h5createFile(file = outfile)
  suppressAll(h5createDataset(
    file = outfile, dataset = "i", storage.mode = "integer",
    dims = c(nrow(wdf), 1), level = compressLevel
  ))
  suppressAll(h5write(obj = wdf$i, file = outfile, name = "i"))

  suppressAll(h5createDataset(
    file = outfile, dataset = "j", storage.mode = "integer",
    dims = c(nrow(wdf), 1), level = compressLevel
  ))
  suppressAll(h5write(obj = wdf$j, file = outfile, name = "j"))

  suppressAll(h5createDataset(
    file = outfile, dataset = "val", storage.mode = "integer",
    dims = c(nrow(wdf), 1), level = compressLevel, fillValue = 0
  ))
  suppressAll(h5write(obj = wdf$val, file = outfile, name = "val"))

  suppressAll(h5write(obj = genenms, file = outfile, name = "gene"))
  suppressAll(h5write(obj = barcodes, file = outfile, name = "barcode"))
  message(paste("Gmat has been written into:", outfile))
  mat <- Matrix::sparseMatrix(
    i = wdf$i, j = wdf$j, x = wdf$val,
    dims = c(length(genes), length(barcodes))
  )
  ## matrix allows the repeat rownames (data.frame not)
  rownames(mat) <- genenms
  colnames(mat) <- barcodes
  return(mat)
}

#' Load sparse matrix of gene matrix with both row and col names
#' @param geneMatrixH5File string, h5 file to read
#' @importFrom rhdf5 h5createFile
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
  rownames(h5_mat) <- as.vector(h5_gene)
  colnames(h5_mat) <- as.vector(h5_barcode)
  return(h5_mat)
}

#' Get Gmat from Snap object
#' @param snapFile string
#' @param geneBedFile string, default is NULL. Only used for mapping with blacklist
#' @param ibarcode vector of string, the barcodes we want to, we assume all of them are in
#' this snapFile. So use this carefully.
#' @param blacklistBedFile string, bed file of the blacklist, smmtools save one under data dir for mm10.
#' @param outfile string, default is NULL, otherwise please end of .h5, and the function will save gmat.
#' if the file already exists, the function will delete it.
#' @param compressLevel integer, default is 9, used for rhdf5.
#' @return sparseMatrix, cell by feature, with both row and colname and ordered by ibarcode if provided
#' @import Matrix
#' @export
getGmatFromSnap <- function(snapFile, geneBedFile = NULL,ibarcode = NULL, blacklistBedFile = NULL,
                            outfile = NULL, compressLevel = 9) {
  barcode <- rhdf5::h5read(file = snapFile, name = "/BD/name")
  ## /GM/name is saved as python binary string in SnapTools, but when loading into R, it's fine.
  ## i.e. there is b'' around a gene name.
  geneName <- as.character(rhdf5::h5read(file = snapFile, name = "/GM/name"))
  idx <- as.integer(rhdf5::h5read(file = snapFile, name = "/GM/idx"))
  idy <- as.integer(rhdf5::h5read(file = snapFile, name = "/GM/idy"))
  ## original count is saved as raw in R, use as.integer/as.numeric to convert it to number.
  count <- as.integer(rhdf5::h5read(file = snapFile, name = "/GM/count"))
  gmat <- Matrix::sparseMatrix(i = idx, j = idy, x = count, dims = c(length(barcode), length(geneName)))
  rownames(gmat) <- as.vector(as.character(barcode))
  colnames(gmat) <- as.vector(as.character(geneName))
  ## filter barcode
  if(is.null(ibarcode)) {
    gmat <- gmat[as.vector(ibarcode), ]
  }
  ## filter blacklist
  ## NOTE: temp remove this due to the duplicates happen in geneBedFile
  ## if( (!is.null(blacklistBedFile)) & (!is.null(geneBedFile)) ) {
  ##   blacklist <- read.table(file = blacklistBedFile, header = FALSE)
  ##   blackgr <- GenomicRanges::GRanges(seqnames = blacklist[,1],
  ##                                     ranges = IRanges::IRanges(start = blacklist[,2], end = blacklist[, 3]))
  ##   genelist <- read.table(file = geneBedFile, header = FALSE)
  ##   rownames(genelist) <- genelist[, 4]
  ##   ## select and reorder the gene list as geneName in snap
  ##   genelist <- genelist[geneName]
  ##   ggr <- GenomicRanges::GRanges(seqnames = genelist[,1],
  ##                                 ranges = IRanges::IRanges(start = genelist[,2], end = genelist[,3]))
  ##   black_ovs <- as.data.frame(GenomicRanges::findOverlaps(query = ggr, subject = blackgr))
  ##   if(nrow(black_ovs) >= 1) {
  ##     gmat <- gmat[,-sort(unique(black_ovs$queryHits))]
  ##   }
  ## }
  ## save gmat
  if(!is.null(outfile)) {
    outdir <- dirname(outfile)
    if (!dir.exists(outdir)) {
      message(paste(outdir, "not exists and create it."))
      dir.create(outdir)
    }
    if (file.exists(outfile)) {
      message(paste(outfile, "exists and remove it."))
      file.remove(outfile)
    }
    t <- as(gmat, "TsparseMatrix")
    rhdf5::h5createFile(outfile)
    suppressAll(rhdf5::h5createDataset(
      file = outfile, dataset = "i", storage.mode = "integer", dims = c(length(t@i),1), level = compressLevel))
    suppressAll(rhdf5::h5write(obj = t@i + 1, file = outfile, name = "i"))
    suppressAll(rhdf5::h5createDataset(
      file = outfile, dataset = "j", storage.mode = "integer", dims = c(length(t@i),1), level = compressLevel))
    suppressAll(rhdf5::h5write(obj = t@j + 1, file = outfile, name = "j"))
    suppressAll(rhdf5::h5createDataset(
      file = outfile, dataset = "val", storage.mode = "integer", dims = c(length(t@i),1), level = compressLevel))
    suppressAll(rhdf5::h5write(obj = t@x, file = outfile, name = "val"))
    suppressAll(rhdf5::h5write(obj = as.vector(colnames(gmat)), file = outfile, name = "gene"))
    suppressAll(rhdf5::h5write(obj = as.vector(rownames(gmat)), file = outfile, name = "barcode"))
    message(paste(outfile, " has saved the gmat."))
  }
  return(gmat)
}


#' Load SnapATAC Gmat from the h5 file.
#' 
#' @param gmatFile string
#' @return sparseMatrix, cell by gene with names
#' @export
loadSnapATACGmatFromFile <- function(gmatFile) {
  barcode <- rhdf5::h5read(file = gmatFile, name = "barcode")
  gene <- rhdf5::h5read(file = gmatFile, name = "gene")
  idx <- as.vector(rhdf5::h5read(file = gmatFile, name = "i"))
  idy <- as.vector(rhdf5::h5read(file = gmatFile, name = "j"))
  val <- as.vector(rhdf5::h5read(file = gmatFile, name = "val"))
  gmat <- Matrix::sparseMatrix(i = idx, j = idy, x = val, index1 = TRUE,
                               dims = c(length(barcode), length(gene)))
  rownames(gmat) <- as.vector(barcode)
  colnames(gmat) <- as.vector(gene)
  return(gmat)
}
