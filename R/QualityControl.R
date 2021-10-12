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

#' @export
fastGetTSSEnrichmentMultiThreads <- function(TSS, barcodes,
                                             rawH5File,
                                             window = 101, norm = 100,
                                             flank = 2000, minNorm = 0.2, maxFragSize = NULL,
                                             nthread = 2) {
  tstart <- Sys.time()
  message(paste("Get TSS Encrichment Scores starts at", tstart))
  TSS <- resize(x = TSS, width = 1, fix = "start")
  BiocGenerics::strand(x = TSS) <- "*"
  TSS <- unique(TSS)
  tssWindow <- resize(x = TSS, width = window, fix = "center")
  tssWindow$type <- "window"
  tssFlank <- c(
    GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(TSS),
      ranges = IRanges::IRanges(start = BiocGenerics::end(TSS) + flank - norm + 1,
                                end = BiocGenerics::end(TSS) + flank)
    ),
    GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(TSS),
      ranges = IRanges::IRanges(
        start = BiocGenerics::start(TSS) - flank,
        end = BiocGenerics::start(TSS) - flank + norm - 1
      )
    )
  )
  tssFlank$type <- "flank"
  tssFeatures <- c(tssWindow, tssFlank)
  tssFeatureList <- S4Vectors::split(x = tssFeature, f = GenomeInfoDb::seqnames(tssFeatures))
  ## get Available chrs
  groups <- rhdf5::h5ls(file = rawH5File)
  groups <- groups[groups$group == "/" & groups$otype == "H5I_GROUP", "name"]
  groups <- groups[grepl("Fragments", groups)]
  
  chrs <- gsub(pattern = "Fragments", replacement = "",
               unique(unlist(strsplit(x = groups, split = "#"))[1]))
  
  tssFeatureList <- S4Vectors::split(x = tssFeature, f = GenomeInfoDb::seqnames(tssFeatures))
  tssFeatureList <- tssFeatureList[chrs]

  coutDF <- mclapply(X=seq_along(featureList), FUN = function(x) {
    nWindow <- rep(0, length(barcodes))
    names(nWindow) <- barcodes
    ## clone in R
    nFlank <- nWindow
    feature <- featureList[[x]]
    
  }, mc.cores = nthread)
}
