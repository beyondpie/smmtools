#' Generate matrix of couts for each tile per cell
#'
#' Ref: ArchR
#'
#' @param outfilenm string, end with .h5
#' @param barcodes vector of string, no need to append sampleName, will do that in this code block
#' @importFrom rhdf5 h5createFile h5createGroup h5createDataset h5closeAll h5write
#' @importFrom S4Vectors start end match DataFrame mcols Rle
#' @importFrom GenomicRanges seqnames split
#' @return  outfile string, where the matrix are saved 
#' @export
getTileMatrix <- function(rawH5File, outdir, outfilenm,
                          genome = "mm10",
                          sampleName = NULL, barcodes = NULL,
                          tileSize = 5000,
                          excludeChr = c("chrM")) {
  tstart <- Sys.time()
  message(paste("Begin to run getTileMatrix with tileSize: ", tileSize, "at", tstart))

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(outdir, outfilenm)
  if(file.exists(outfile)) {
    message(paste(outfile, "exist, and remove it."))
    file.remove(outfile)
  }
  rhdf5::h5createFile(file = outfile)
  ## h5createGroup(file = outfile, group = "TileMatrix")
  ## h5createGroup(file = outfile, group = "colnames")
  ## h5createGroup(file = outfile, group = "rownames")

  annotGenome <- getAnnotFromArchRData(tag = "genome", genome = genome)

  blacklist <- annotGenome$blacklist
  if (!is.null(blacklist)) {
    if (length(blacklist) > 0) {
      message(paste("Find blacklist for genome", genome))
      blacklist <- split(blacklist, seqnames(blacklist))
    }
  }else {
    message(paste("No blacklist for genome", genome))
  }
  
  chromSizes <- annotGenome$chromSizes
  chromLengths <- end(chromSizes)
  names(chromLengths) <- chromSizes@seqinfo@seqnames
  chromLengths <- chromLengths[!(match(names(chromLengths), excludeChr, nomatch = 0) > 0)]

  dfParams <- data.frame(
    seqnames = names(chromLengths),
    length = as.vector(chromLengths),
    tileSize = tileSize,
    binarize = FALSE,
    stringsAsFactors = FALSE
  )

  ## barcodes_with_sampleName <- barcodes
  ## if (!grepl(pattern = "#", barcodes[1], fixed = TRUE)) {
  ##   barcodes_with_sampleName <- paste0(sampleName, "#", barcodes)
  ## }

  featureDF <- lapply(seq_along(chromLengths), function(x) {
    DataFrame(seqnames = names(chromLengths)[x],
              idx = seq_len(trunc(chromLengths[x]) / tileSize + 1))
  }) %>% Reduce("rbind", .)
  featureDF$start <- (featureDF$idx - 1) * tileSize

  df <- data.frame(featureDF, stringsAsFactors = FALSE)
  supressAll(h5write(obj = df, file = outfile, name = "/FeatureDF"))

  ## record colnames for the matrix
  ## TODO: may have bugs when not all barcodes are in the raw h5.
  ## Should write this part later, and check if barcodes in the least.
  supressAll(h5write(obj = barcodes, file = outfile, name = "/Barcodes"))

  for (z in seq_along(chromLengths)) {
    chr <- names(chromLengths)[z]
    chrl <- chromLengths[z]
    message(paste(chr, "with length ", chrl))
    fragments <- getFragsOfAChrFromRawH5File(
      rawH5File = rawH5File, chr = chr,
      sampleName = sampleName,
      barcodes = barcodes
    )
    if(length(fragments) < 1) {
      message(paste(chr, "has no fragments left after filtering the barcodes."))
      next
    } else {
      message(paste(chr, "has", length(fragments), "left after filtering the barcodes."))
    }
    nTiles <- trunc(chrl / tileSize) + 1
    matchBarcodes <- match(mcols(fragments)$RG, barcodes)
    mat <- Matrix::sparseMatrix(
      i = as.integer(c(trunc(start(fragments) / tileSize),
                       trunc(end(fragments) / tileSize)) + 1),
      j = as.integer(as.vector(c(matchBarcodes, matchBarcodes))),
      x = rep(1, 2 * length(fragments)),
      dims = c(nTiles, length(barcodes))
    )
    ## remove blacklisted tiles
    if (!is.null(blacklist) & (length(blacklist) > 0)) {
      blacklistz <- blacklist[[chr]]
      if (length(blacklistz) > 0) {
        tile2 <- floor(tileSize / 2)
        blacklistIdx <- unique(
          trunc(start(unlist(GenomicRanges::slidingWindows(blacklistz, tile2, tile2))) / tileSize) + 1)
        blacklistIdx <- sort(blacklistIdx)
        idxToZero <- which(match((mat@i + 1), blacklistIdx, nomatch = 0) > 0)
        if (length(idxToZero) > 0) {
          mat@x[idxToZero] <- 0
          mat <- Matrix::drop0(mat)
        }
      }
    } ## end of remove blacked tiles
    ## save to file
    h5createGroup(file = outfile, group = chr)
    lengthI <- length(mat@i)
    supressAll(h5createDataset(
      file = outfile, dataset = paste0(chr, "/i"), storage.mode = "integer",
      dims = c(lengthI, 1), level = 0
    ))
    supressAll(h5write(obj = mat@i + 1, file = outfile, name = paste0(chr, "/i")))

    supressAll(h5createDataset(
      file = outfile, dataset = paste0(chr, "/x"), storage.mode = "double",
      dims = c(lengthI, 1), level = 0
    ))
    supressAll(h5write(obj = mat@x, file = outfile, name = paste0(chr, "/x")))

    ## #Convert Columns to Rle
    j <- Rle(findInterval(seq_along(mat@x) - 1, mat@p[-1]) + 1)
    ## effective ncols
    lengthRle <- length(j@lengths)

    supressAll(h5createDataset(
      file = outfile, dataset = paste0(chr, "/jLengths"), storage.mode = "integer",
      dims = c(lengthRle, 1), level = 0
    ))
    supressAll(h5write(obj = j@lengths, file = outfile, name = paste0(chr, "/jLengths")))

    supressAll(h5createDataset(
      file = outfile, dataset = paste0(chr, "/jValues"), storage.mode = "integer",
      dims = c(lengthRle, 1), level = 0
    ))
    ## corresoinds to cols with elements
    supressAll(h5write(obj = j@values, file = outfile, name = paste0(chr, "/jValues")))
  } ## end of for loop of chrs
  h5closeAll()
  return(outfile)
}

#' load TileMatrix into a sparseMatrix
#' Ref: ArchR
#'
#' @param barcodes vector of string, if not null, filter tilematrix by barcodes, default is NULL.
#' We assume barcodes are all in the tileMatrixFile. NOTE: No check for this.
#' @param binarize bool, default is FALSE, this is a simple binarization,
#' is not the version SnapATAC describes, should not be used.
#' TODO: remove this binarize param.
#' @return sparseMatrix, feature by cell, both row and colnames, ordered by barcodes if provides.
#' @export
loadTileMatrix <- function(tileMatrixFile, barcodes = NULL, binarize = FALSE) {
  colNames <- removeSampleName(barcodes = h5read(file = tileMatrixFile, name = "/Barcodes"))
  if(!is.null(barcodes)) {
    barcodes <- removeSampleName(barcodes = barcodes)
    idxCols <- which(colNames %in% barcodes)
  } else {
    idxCols <- seq_along(colNames)
  }
  featureDF <- h5read(file = tileMatrixFile, name = "/FeatureDF")
  chrs <- unique(featureDF$seqnames)
  matList <- lapply(seq_along(chrs), function(x) {
    chr <- chrs[x]
    idxRows <- featureDF[featureDF$seqnames == chr, "idx"]

    RleJ <- S4Vectors::Rle(values = h5read(file = tileMatrixFile, name = paste0("/", chr, "/jValues")),
                lengths = h5read(file = tileMatrixFile, name = paste0("/", chr, "/jLengths")))
    matchJ <- S4Vectors::match(RleJ, idxCols, nomatch = 0)
    idxJ <- BiocGenerics::which(matchJ > 0)
    j <- matchJ[idxJ]
    i <- h5read(file = tileMatrixFile, name = paste0("/", chr, "/i"))[idxJ]

    matchI <- match(i, idxRows, nomatch = 0)
    idxI <- which(matchI > 0)
    ## i <- i[idxI]
    j <- j[idxI]
    i <- matchI[idxI]
    
    if(binarize) {
      x <- rep(1, length(j))
    } else {
      x <- h5read(file = tileMatrixFile, name = paste0("/", chr, "/x"))[idxJ][idxI]
    }
    mat <- Matrix::sparseMatrix(
      i = as.integer(i),
      j = j,
      x = x,
      dims = c(length(idxRows), length(idxCols))
    )
    rownames(mat) <- paste0(chr,"_",featureDF[featureDF$seqnames == chr, "idx"])
    return(mat)
  })
  mat <- Reduce("rbind", matList)
  colnames(mat) <- colNames[idxCols]

  ## double check order
  rowNames <- paste0(featureDF$seqnames, "_",featureDF$idx)
  mat <- mat[rowNames, , drop = FALSE]
  if(!is.null(barcodes)) {
    ## one-dim array is not suitable for sparseMatrix
    ## use vector instead just in case
    mat <- mat[, as.vector(barcodes), drop = FALSE]
  }
  return(mat)
}
  
