#' @param nChunk integer partition each chrom into nChunk pieces
#' @export
sumFragmentSingleThread <- function(tabixFile, chromsSizes, outdir,
                                    sampleName = NULL, barcodes = NULL,
                                    nChunk = 3) {
  outH5File = file.path(outdir, paste0(paste(sampleName,"tabix2H5", "nChunk", nChunk, sep="_"), ".h5"))
  tileChromSizes <- tileChrom(chromSizes = chromSizes, nChunk = nChunk)
  ## Transform tab.gz to h5 file
  tabixToH5SingleThread(tabixFile = tabixFile, tileChromSizes = tileChromSizes,
                        sampleName = sampleName, outH5File = outH5File
                        barcodes = barcodes)
  chunkName <- S4Vectors::mcols(x = tileChromSizes)$chunkName
  ## get nfrag per barcode
  dtList <- laaply(seq_along(chunkName), function(x) {
    chrRegion <- chunkName[x]
    ## fragmentRanges <- paste0("Fragments/", chrRegion, "/Ranges")
    barcodeLength <- rhdf5::h5read(file = outH5File, name = paste0("Fragments/", chrRegion, "/BarcodeLength"))
    barcodeValue <- rhdf5::h5read(file = outH5File, name = paste0("Fragments/", chrRegion, "/BarcodeValue"))
    dt <- NULL
    if ((length(barcodeValue) > 0) & (length(barcodeLength) > 0)) {
      dt <- data.table(values = barcodeValue,
                       lengths = barcodeLength)
    }
  })
  names(dtList) <- chunkName
  dt <- Reduce(f = "rbind", dtList)
  dt <- dt[, sum(lengths), by = list(values)]
  
}
