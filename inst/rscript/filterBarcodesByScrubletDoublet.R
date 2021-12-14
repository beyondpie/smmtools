library(optparse)

option_list <- list(
  make_option(c("--doubletScoreFile"), type = "character"),
  make_option(c("--outBarcodeFile"), type = "character"),
  make_option(c("--outPlotFile"), type = "character"),
  make_option(c("--sampleName"), type = "character")
)
args <- parse_args(OptionParser(option_list = option_list))
doubletScores <- read.table(file = args$doubletScoreFile, header = TRUE,
                            sep = ",", comment.char = "")

pdf(file = args$outPlotFile, width = 10, height = 8)
hist(doubletScores$scrubletScore, breaks = 50, xlab = "Scrublet Scores",
     main = paste("Histogram of Scrublet scores for sample", args$sampleName))
abline(v = doubletScores$thresGMM[1], col = "red", lwd=3, lty=2)
dev.off()

out <- doubletScores$barcode[doubletScores$scrubletScore <= doubletScores$thresGMM]
invisible(write.table(x = out, file = args$outBarcodeFile,
                      quote = FALSE,
                      row.names = FALSE,
                      col.names = FALSE,
                      sep = ","))
