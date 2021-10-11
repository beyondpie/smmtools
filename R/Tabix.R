#' Read fragment information
#'
#' Ref: ArchR .tabixTotmp
#'
#' @param tabixFile string fragment file name, only tab/tab.gz file (.tabix should be in the same folder)
#' @param chromSizes GenomeRanges
#' @param sampleName string the sample corresponding to the fragment file
#' @param barcodes vector of string, optional, barcodes we want
#' @param outH5File string filename of the output h5file
#' @param nChunk integer partition each chrom into nChunk pieces
#' @return outputH5File string name of the output file
#' @export
tabixToh5SingleThread <- function(tabixFile, chromSizes = NULL,
                                  sampleName = NULL, barcodes = NULL,
                                  outH5File = tempfile(
                                    pattern = paste0("tmp-", sampleName, ".h5"),
                                    tmpdir = tempdir()
                                  ),
                                  nChunk = 3) {
  tstart <- Sys.time()
  message(paste("Tabix to h5 with single thread starts at", tstart))
  o <- rhdf5::h5closeAll()
  o <- rhdf5::h5createFile(file = outH5File)
  o <- rhdf5::h5createGroup(file = outH5File, group = "Fragments")
  o <- rhdf5::h5createGroup(file = outH5File, group = "Metadata")

  ## after unlist, tileChromSizes is GenomicRanges object
  tileChromSizes <- unlist(GenomicRanges::tile(chromSizes, nChunk))
  S4Vectors::mcols(tileChromSizes)$chunkName <- paste0(
    GenomeInfoDb::seqnames(tileChromSizes), "#chunk",
    seq_along(tileChromSizes)
  )

  for (x in seq_along(tileChromSizes)) {
    if (x %% 10 == 0) {
      sprintf("%s Reading TableFile %s Percent", sampleName, round(100 * x / length(tileChromSizes)), 3)
    }
    tileChromStringVector <- Rsamtools::scanTabix(file = tabixFile, param = tileChromSizes[x])[[1]]
    tmp <- read.table(textConnection(tileChromStringVector))
    dt <- data.table::data.table(start = as.integer(tmp$V2 + 1), end = tmp$V3, barcode = tmp$V4)
    message(paste0(
      sampleName, " Fragment-Chunk-(", x, " of ", length(tileChromSizes), ")-",
      nrow(dt)
    ))
    message(paste0(
      sampleName, " Barcodes-Chunk-(", x, " of ", length(tileChromSizes), ")-",
      unique(dt$barcode)
    ))
    # Care for Break Points from ArchR
    dt <- dt[dt$start >= S4Vectors::start(tileChromSizes[x]), ]

    if (!is.null(barcodes)) {
      dt <- dt[dt$barcode %in% barcodes, ]
    }
    # Order by barcodes
    data.table::setkey(dt, barcode)
    dt <- dt[order(barcode)]
    barcodeRle <- Rle(paste0(dt$barcode))
    fragmentRanges <- paste0("Fragments/", chrRegion, "/Ranges")
    barcodeLength <- paste0("Fragments/", chrRegion, "/BarcodeLenght")
    barcodeValue <- paste0("Fragments/", chrRegion, "/BarcodeValue")

    chrRegion <- S4Vectors::mcols(tileChromSizes)$chunkName[x]
    o <- rhdf5::h5createGroup(file = outH5File, group = paste0("Fragments/", chrRegion))
    o <- suppressAll(expr = rhdf5::h5createDataset(
      file = outH5File, dataset = fragmentRanges,
      storage.mode = "integer", dims = c(nrow(dt), 2),
      level = 0
    ))
    o <- suppressAll(expr = rhdf5::h5createDataset(
      file = outH5File, dataset = barcodeLength,
      storage.mode = "integer", dims = c(length(barcodeRle@lengths), 1), level = 0
    ))

    o <- suppressAll(expr = rhdf5::h5createDataset(
      file = outH5File,
      datast = barcodeValue,
      storage.mode = "character",
      dims = c(length(barcodeRle@lengths), 1), level = 0,
      size = max(nchar(barcoeRle@values)) + 1
    ))
    o <- rhdf5::h5write(obj = cbind(dt$start, dt$end-dt$start +1), file = outH5File, name = fragmentRanges)
    o <- rhdf5::h5write(obj = barcodeRle@lengths, file = outH5File, name = barcodeLength)
    o <- rhdf5::h5write(obj = barcodeRle@values, file = outH5File, name = barcodeValue)
  }
  estart <- Sys.time()
  message(paste("Tabix to h5 with single thread ends at", estart))
  return(outH5File)
}
