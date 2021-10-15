#' get barcodes-related fragments stats.
#'
#' Ref: ArchR
#' 
#' @param nChunk integer partition each chrom into nChunk pieces
#' @export
sumFragmentSingleThread <- function(tabixFile, chromSizes, outdir,
                                    sampleName = NULL, barcodes = NULL,
                                    nChunk = 3) {
  dir.create(outdir, showWarnings = TRUE, recursive = TRUE)
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
    sampleName = sampleName, outH5File = rawH5File,
    barcode = barcodes
  )
  chunkName <- S4Vectors::mcols(x = tileChromSizes)$chunkName
  ## get nfrag per barcode
  dtNfragmentPerBarcode <- getNfragmentPerBarcode(chrRegions = chunkName, rawH5File = rawH5File,
                                                  sampleName = sampleName)
  return(dtNfragmentPerBarcode)
}

#' get numebr of fragments per barcode
#' Ref: ArchR
#'
#' @param chrRegions array of characters, chrRegion names in rawH5File
#' @param rawH5File string 
#' @return tibble shape = [length(barcodes), 2], colnames = c("barcode", "count")
#' @importFrom dplyr %>% group_by summarise arrange desc
#' @export
getNfragmentPerBarcode <- function(chrRegions, rawH5File, sampleName = NULL) {
  ## get nfrag per barcode
  dtList <- lapply(seq_along(chrRegions), function(x) {
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
      dt <- data.table::data.table(
        values = paste(sampleName, barcodeValue, sep = "#"),
        lengths = barcodeLength
      )
    }
    invisible(dt)
  })
  names(dtList) <- chrRegions
  dt <- Reduce(f = "rbind", dtList)
  colnames(dt) <- c("barcode", "count")
  dt <- dt %>% group_by(barcode) %>% summarise(count = sum(count)) %>% arrange(desc(count))
  return(dt)
}

#' Get TSS Enrichment scores
#'
#' Ref: ArchR
#' @return List two element: TSSE, TSSReads
#' @importFrom S4Vectors mcols mcols<- split DataFrame queryHits subjectHits
#' @importFrom BiocGenerics strand<- match start end pmax 
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges width ranges resize
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
  ## get Available chrs
  groups <- h5ls(file = rawH5File)
  groups <- groups[groups$group == "/Fragments" & groups$otype == "H5I_GROUP", "name"]

  chrs <- unique(gsub(pattern = "#chunk\\d+",
                      replacement = "",
                      x = groups))
  
  tssFlank$type <- "flank"
  tssFeatures <- c(tssWindow, tssFlank)
  tssFeatureList <- split(x = tssFeatures, f = seqnames(tssFeatures))
  tssFeatureList <- tssFeatureList[chrs]

  countDF <- parallel::mclapply(X = seq_along(tssFeatureList), FUN = function(i) {
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
      mcols(feature)$typeIdx <- match(mcols(feature)$type, c("window", "flank"))
      ## count each insertion
      for(y in seq_len(2)) {
        if(y == 1) {
          temp <- IRanges(start(fragments), width = 1)
        } else {
          temp <- IRanges(end(fragments), width = 1)
        }
        o <- findOverlaps(ranges(feature), temp)
        ## CALL CPP
        mat <- tabulate2dCpp(
          x = as.vector(mcols(fragments)$RG[subjectHits(o)]),
          xmin = 1,
          xmax = length(barcodes),
          y = mcols(feature)$typeIdx[queryHits(o)],
          ymin = 1,
          ymax = 2
        )

        nWindow <- nWindow + mat[1, ]
        nFlank <- nFlank + mat[2, ]
      }
    }
    names(nWindow) <- NULL
    names(nFlank) <- NULL
    return(DataFrame(nWindow = nWindow, nFlank = nFlank))
  }, mc.cores = nthread)
  cumDF <- countDF[[1]]
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
  return(list(TSSE = tssScores, TSSReads = cumDF[,1]))
}
