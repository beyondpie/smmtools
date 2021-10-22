#' Generate matrix of cell by genome bins for one sample
#'
#' Ref: ArchR
#'
#' @importFrom rhdf5 h5closeAll
#' @export
getTileMatrix <- function(rawH5File, outdir, outfilenm,
                          genome = "mm10",
                          sampleName = NULL, barcodes = NULL,
                          nChunk = 3, tileSize = 5000,
                          binarize = TRUE, excludeChr = c("chrM", "chrY"),
                          threads = 1) {
  tstart <- Sys.time()
  message(paste("Begin to run getTileMatrix  with tileSize: ", tileSize, "at", tstart))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(outdir, outfilenm)
  o <- h5closeAll()

  annotGenome <- getAnnotFromArchRData(tag = "genome", genome = genome)
  return()
}
