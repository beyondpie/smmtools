tabixToh5SingleThread <- function(tabixFile, sampleName = NULL, barcodes = NULL,
                                  outH5File = tempfile(pattern = paste0("tmp-", sampleName, ".h5"),
                                                       tmpdir = tempdir())
                                  nChunk = 3) {
  tstart <- Sys.time()
  message(paste("Tabix to h5 with single thread at", tstart))
  o <- rhdf5::h5closeAll()
  o <- rhdf5::h5createFile(file = outH5File)
  o <- rhdf5::h5createGroup(file = outH5File, group = "Fragments")
  o <- rhdf5::h5createGroup(file = outH5File, group = "Metadata")
  
  return(outH5File)
}
