library(optparse)

option_list <- list(
  make_option(c("--doubletScoreFile"), type = "character",
              default = "../out/scrublet/CEMBA210429_18B.csv"),
  make_option(c("--outPlotFile"), type = "character",
              default = "../out/plotScrubletScores/CEMBA210429_18B.pdf"),
  make_option(c("--sampleName"), type = "character", default = "CEMBA210429_18B")
)
args <- parse_args(OptionParser(option_list = option_list))
doubletScores <- read.table(file = args$doubletScoreFile, header = TRUE,
                            sep = ",", comment.char = "")
ncell <- nrow(doubletScores)
ratioGMM <- round(sum(doubletScores$scrubletScore > doubletScores$thresGMM) / ncell, 3)
ratioRaw <- round(sum(doubletScores$scrubletScore > doubletScores$thresScrublet) / ncell, 3)
ratioMixEM <- round(sum(doubletScores$scrubletScore > doubletScores$thresMixEM) / ncell, 3)

pdf(file = args$outPlotFile, width = 10, height = 8)
hist(doubletScores$scrubletScore, breaks = 50, xlab = "Scrublet Scores",
     main = paste("Histogram of Scrublet scores for sample", args$sampleName))
abline(v = doubletScores$thresGMM[1], col = "red", lwd=3, lty=2)
abline(v = doubletScores$thresScrublet[1], col = "blue", lwd = 3, lty = 2)
abline(v = doubletScores$thresMixEM[1], col = "green", lwd = 3, lty = 2)
legend("topright", legend = c("GMM", "Raw", "MixEM"),
       col = c("red", "blue", "green"), lty = c(2,2,2), cex = 0.8)
mtext(text = paste("Total cells", ncell,
                   "Doublet rate: GMM", paste0("(",round(doubletScores$thresGMM[1],2),")"),ratioGMM,
                   "; Raw", paste0("(",round(doubletScores$thresScrublet[1],2),")"), ratioRaw,
                   "; MixEM", paste0("(",round(doubletScores$thresMixEM[1],2),")"),ratioMixEM), side = 3)
dev.off()
