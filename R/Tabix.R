#' Read fragment information
#'
#' Ref: ArchR .tabixTotmp
#'
#' @param tabixFile string fragment file name, only tab/tab.gz file (.tabix should be in the same folder)
#' @param tileChromSizes GenomeRanges, tiled ChromSizes
#' @param sampleName string the sample corresponding to the fragment file
#' @param barcodes vector of string, optional, barcodes we want
#' @param outH5File string filename of the output h5file
#' @param nChunk integer partition each chrom into nChunk pieces
#' @return outputH5File string name of the output file
#' @export
tabixToH5SingleThread <- function(tabixFile, tileChromSizes,
                                  sampleName = NULL, barcodes = NULL,
                                  outH5File = tempfile(
                                    pattern = paste0("tmp-", sampleName, ".h5"),
                                    tmpdir = tempdir()
                                  )) {
  tstart <- Sys.time()
  message(paste("Tabix to h5 with single thread starts at", tstart), "...")
  o <- rhdf5::h5closeAll()
  o <- rhdf5::h5createFile(file = outH5File)
  o <- rhdf5::h5createGroup(file = outH5File, group = "Fragments")
  o <- rhdf5::h5createGroup(file = outH5File, group = "Metadata")

  for (x in seq_along(tileChromSizes)) {
    if (x %% 10 == 0) {
      message(paste("Reading", sampleName, "TableFile", round(100 * x / length(tileChromSizes), 3), "Percent."))
    }
    tileChromStringVector <- Rsamtools::scanTabix(file = tabixFile, param = tileChromSizes[x])[[1]]
    tmp <- utils::read.table(textConnection(tileChromStringVector))
    dt <- data.table::data.table(start = as.integer(tmp$V2 + 1), end = tmp$V3, barcode = tmp$V4)
    # Care for Break Points from ArchR
    dt <- dt[dt$start >= S4Vectors::start(tileChromSizes[x]), ]

    if (!is.null(barcodes)) {
      dt <- dt[dt$barcode %in% barcodes, ]
    }
    ## message(paste0(
    ##   sampleName, " Fragment-Chunk-(", x, " of ", length(tileChromSizes), ")-",
    ##   nrow(dt)
    ## ))
    ## message(paste0(
    ##   sampleName, " Barcodes-Chunk-(", x, " of ", length(tileChromSizes), ")-",
    ##   unique(dt$barcode)
    ## ))
    # Order by barcodes
    chrRegion <- S4Vectors::mcols(tileChromSizes)$chunkName[x]
    data.table::setkey(dt, barcode)
    dt <- dt[order(dt$barcode), ]
    barcodeRle <- S4Vectors::Rle(paste0(dt$barcode))
    fragmentRanges <- paste0("Fragments/", chrRegion, "/Ranges")
    barcodeLength <- paste0("Fragments/", chrRegion, "/BarcodeLength")
    barcodeValue <- paste0("Fragments/", chrRegion, "/BarcodeValue")

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
      dataset = barcodeValue,
      storage.mode = "character",
      dims = c(length(barcodeRle@lengths), 1), level = 0,
      size = max(nchar(barcodeRle@values)) + 1
    ))
    o <- rhdf5::h5write(obj = cbind(dt$start, dt$end-dt$start +1), file = outH5File, name = fragmentRanges)
    o <- rhdf5::h5write(obj = barcodeRle@lengths, file = outH5File, name = barcodeLength)
    o <- rhdf5::h5write(obj = barcodeRle@values, file = outH5File, name = barcodeValue)
  }
  estart <- Sys.time()
  message(paste("Tabix to h5 with single thread ends at", estart))
  return(outH5File)
}
