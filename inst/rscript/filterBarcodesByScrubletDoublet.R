library(optparse)

option_list <- list(
  make_option(c("--doubletScoreFile"), type = "character"),
  make_option(c("--outBarcodeFile"), type = "character")
)
args <- parse_args(OptionParser(option_list = option_list))
doubletScores <- read.table(file = args$doubletScoreFile, header = TRUE,
                            sep = ",", comment.char = "")
out <- doubletScores$barcode[doubletScores$scrubletScore <= doubletScores$thresGMM]
invisible(write.table(x = out, file = args$outBarcodeFile,
                      quote = FALSE,
                      row.names = FALSE,
                      col.names = FALSE,
                      sep = ","))
