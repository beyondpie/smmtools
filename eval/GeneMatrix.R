## get gene matrix for doublet removement
library(rtracklayer)
library(rhdf5)
## remotes::install_github("beyondpie/smmtools")
## detach("package:smmtools", unload=TRUE)
library(smmtools)

## install.packages("plot3D")
## install.packages("doSNOW")
## remotes::install_github("r3fang/SnapATAC")
## library(SnapATAC)

## setwd("/mac/git-recipes/smmtools/eval")
rawH5File <- file.path("rawh5", "CEMBA201008_9F_raw.h5")
sumFragFile <- file.path("rawh5", "CEMBA201008_9F_smmtoolsQC.csv")
snapFile <- file.path("snap", "CEMBA201008_9F.snap")

outdir <- "gm"
outfilenm <- "CEMBA201008_9F.h5"
sampleName <- "CEMBA201008_9F"
genome <- "mm10"
excludeChr <- c("chrM", "chrY", "chrX", paste0("chr", 1:17))

## load geneup2k bed files
geneUp2K <- rtracklayer::import(file.path("snap", "gencode.vM16.geneUp2k.bed"), format = "bed")
genenms <- geneUp2K$name

## get barcodes filtered from smmtools
sumFrag <- read.table(file = sumFragFile, header = TRUE, sep = ",", comment.char = "")
index_from_nFrag <- which(sumFrag$nUniqFrag >= 1000)
index_from_qc <- which(sumFrag$TSSE >= 10)
index_join <- intersect(index_from_nFrag, index_from_qc)
barcodes <- sumFrag[index_join,"barcode"]

## get gmat from smmtools
gmat <- getGeneMatrix(rawH5File = rawH5File, outdir = outdir, outfilenm = outfilenm,
                      genes = geneUp2K, genenms = genenms,
                      genome = genome, excludeChr = excludeChr, barcodes = barcodes)

## test scrublet
gmat <- smmtools::loadGeneMatrix(geneMatrixH5File = file.path(outdir, outfilenm))
barcode <- colnames(gmat)
out <- smmtools::SnapATAC_runScrublet(
  mat = t(gmat), path_to_python = "/root/miniconda3/bin/python",
  expected_doublet_rate = 0.06,
  min_counts = 3,
  min_cells = 5,
  min_conv_pctl = 85,
  n_pc = 30)


invisible(write.table(
  x = data.frame(barcode = barcode,
                 scrubletScore = out[[4]],
                 thresGMM = rep(out[[1]], length(out[[4]]) )
  ),
  file = file.path(outdir, paste0(sampleName,"_scrublet.csv")),
  quote = FALSE,
  row.names = FALSE, col.names = TRUE,
  sep = ","
))

doubletScores <- read.table(file = file.path(outdir, paste0(sampleName,"_scrublet.csv")), header = TRUE,
                            sep = ",", comment.char = "")

pdf(file = file.path(outdir, paste0(sampleName,"_scrublet.pdf")), width = 10, height = 8)
hist(doubletScores$scrubletScore, breaks = 50, xlab = "Scrublet Scores",
     main = paste("Histogram of Scrublet scores for sample", sampleName))
dev.off()

out <- doubletScores$barcode[doubletScores$scrubletScore <= doubletScores$thresGMM]
invisible(write.table(x = out, file = file.path(outdir, paste0(sampleName,"_filterDoublet.csv")),
                      quote = FALSE,
                      row.names = FALSE,
                      col.names = FALSE,
                      sep = ","))




