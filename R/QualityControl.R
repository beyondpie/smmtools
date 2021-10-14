#' get barcodes-related fragments stats.
#'
#' Ref: ArchR
#' 
#' @param nChunk integer partition each chrom into nChunk pieces
#' @export
sumFragmentSingleThread <- function(tabixFile, chromsSizes, outdir,
                                    sampleName = NULL, barcodes = NULL,
                                    nChunk = 3) {
  rawH5File <- file.path(outdir, paste0(paste(sampleName, "tabix2H5", "nChunk", nChunk, sep = "_"), ".h5"))
  fragmentH5File <- file.path(outdir, paste0(sampleName, "sumFragment.h5"))
  o <- suppressAll(file.remove(fragmentH5File))
  o <- rhdf5::h5closeAll()
  o <- rhdf5::h5createFile(file = fragmentH5File)
  o <- rhdf5::h5createGroup(file = fragmentH5File, group = "Metadata")
  o <- rhdf5::h5createGroup(file = fragmentH5File, group = "Fragments")

  tileChromSizes <- tileChrom(chromSizes = chromSizes, nChunk = nChunk)
  ## Transform tab.gz to h5 file
  tabixToH5SingleThread(
    tabixFile = tabixFile, tileChromSizes = tileChromSizes,
    sampleName = sampleName, rawH5File = rawH5File,
    barcodes = barcodes
  )
  chunkName <- S4Vectors::mcols(x = tileChromSizes)$chunkName
  ## get nfrag per barcode
  dtNfragmentPerBarcode <- getNfragmentPerBarcode(chrRegions = chunkName, rawH5File = rawH5File)
  return(dtNfragmentPerBarcode)
}

getNfragmentPerBarcode <- function(chrRegions, rawH5File) {
  ## get nfrag per barcode
  dtList <- laaply(seq_along(chrRegions), function(x) {
    chrRegion <- chrRegions[x]
    ## fragmentRanges <- paste0("Fragments/", chrRegion, "/Ranges")
    barcodeLength <- rhdf5::h5read(
      file = rawH5File,
      name = paste0("Fragments/", chrRegion, "/BarcodeLength")
    )
    barcodeValue <- rhdf5::h5read(
      file = rawH5File,
      name = paste0("Fragments/", chrRegion, "/BarcodeValue")
    )
    dt <- NULL
    if ((length(barcodeValue) > 0) & (length(barcodeLength) > 0)) {
      dt <- data.table(
        values = barcodeValue,
        lengths = barcodeLength
      )
    }
  })
  names(dtList) <- chrRegions
  dt <- Reduce(f = "rbind", dtList)
  dt <- dt[, sum(lengths.V1), by = list(values.V1)]
  ## Order to reduce numebr of hyperslabs
  dt <- dt[order(V1, descreasing = TRUE)]
  return(dt)
}

#' Get TSS Enrichment scores
#'
#' Ref: ArchR
#' @return List two element: TSSE, TSSReads
#' @importFrom S4Vectors mcols split DataFrame queryHits subjectHits
#' @importFrom BiocGenerics strand match start end pmax
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges width ranges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom rhdf5 h5ls
#' @export
fastGetTSSEnrichmentMultiThreads <- function(TSS, barcodes,
                                             rawH5File, nChunkInRawH5File = 3,
                                             window = 101, norm = 100,
                                             flank = 2000, minNorm = 0.2, maxFragSize = NULL,
                                             sampleName = NULL,
                                             nthread = 2) {
  tstart <- Sys.time()
  message(paste("Get TSS Enrichment Scores starts at", tstart))
  TSS <- resize(x = TSS, width = 1, fix = "start")
  strand(x = TSS) <- "*"
  TSS <- unique(TSS)
  tssWindow <- resize(x = TSS, width = window, fix = "center")
  tssWindow$type <- "window"
  tssFlank <- c(
    GRanges(
      seqnames = seqnames(TSS),
      ranges = IRanges(
        start = end(TSS) + flank - norm + 1,
        end = end(TSS) + flank
      )
    ),
    GRanges(
      seqnames = seqnames(TSS),
      ranges = IRanges(
        start = start(TSS) - flank,
        end = start(TSS) - flank + norm - 1
      )
    )
  )
  tssFlank$type <- "flank"
  tssFeatures <- c(tssWindow, tssFlank)
  tssFeatureList <- split(x = tssFeature, f = seqnames(tssFeatures))
  ## get Available chrs
  groups <- h5ls(file = rawH5File)
  groups <- groups[groups$group == "/" & groups$otype == "H5I_GROUP", "name"]
  groups <- groups[grepl("Fragments", groups)]

  chrs <- gsub(
    pattern = "Fragments", replacement = "",
    unique(unlist(strsplit(x = groups, split = "#"))[1])
  )

  tssFeatureList <- split(x = tssFeature, f = seqnames(tssFeatures))
  tssFeatureList <- tssFeatureList[chrs]

  countDF <- mclapply(X = seq_along(tssFeatureList), FUN = function(i) {
    nWindow <- rep(0, length(barcodes))
    names(nWindow) <- barcodes
    ## clone in R
    nFlank <- nWindow
    feature <- tssFeatureList[[i]]
    fragments <- getFragsOfAChrFromRawH5File(
      rawH5File = rawH5File, chr = names(tssFeatureList)[i],
      sampleName = sampleName, nChunk = nChunkInRawH5File
    )
    if (length(fragments) == 0) {
      names(nWindow) <- NULL
      names(nFlank) <- NULL
      return(DataFrame(nWindow = nWindow, nFlank = nFlank))
    }
    if (length(fragments) > 0) {
      if (!is.null(maxFragSize)) {
        fragments <- fragments[width(fragments) <= maxFragSize]
      }
      mcols(fragments)$RG@values <- match(
        mcols(fragments)$RG@values,
        barcodes
      )
      mcols(feature)$RG@typeIdx <- match(mcols(feature)$type, c("window", "flank"))
      ## count each insertion
      for(y in seq_len(2)) {
        if(y == 1) {
          temp <- IRanges(start(fragments), width = 1)
        } else {
          temp <- IRanges(end(fragments), width = 1)
        }
        o <- findOverlaps(ranges(feature), temp)
        x <- as.vector(mcols(fragments)$RG[subjectHits(o)])
        x <-  (x >= 1) & (x <= length(barcodes))
        y <- mcols(feature)$typeIdx[queryHits(o)]
        y <- (y >= 1) & (y <= 2)
        mat <- y %*% t(x)
        nWindow <- nWindow + mat[1, ]
        nFlank <- nFlank + mat[2, ]
      }
    }
    names(nWindow) <- NULL
    names(nFlank) <- NULL
    return(DataFrame(nWindow = nWindow, nFlank = nFlank))
  }, mc.cores = nthread)
  cumDF <- countDF[[i]]
  for(i in 2:length(countDF)) {
    cumDF$nWindow <- cumDF$nWindow + countDF[[i]]$nWindow
    cumDF$nFlank <- cumDF$nFlank + countDF[[i]]$nFlank
  }
  cWn <- cumDF[,1] / window
  cFn <- cumDF[,2] / norm
  ## pmax: like Relu(vec - scalar) + scalar
  tssScores <- round(2 * cWn / pmax(cFn, minNorm), 3)
  names(tssScores) <- barcodes
  tend <- Sys.time()
  message(paste("Get TSS Enrichment Scores ends at", tend))
  return(list(TSSE = tssScores, TSSReads = cumDf[,1]))
}
