## get gene matrix for doublet removement
library(smmtools)
detach("package:smmtools", unload=TRUE)
## library(smmtools)
library(SnapATAC)

setwd("/mac/git-recipes/smmtools/eval")
rawH5File <- file.path("rawh5", "CEMBA201008_9F_raw.h5")
sumFragFile <- file.path("rawh5", "CEMBA201008_9F_smmtoolsQC.csv")
snapFile <- file.path("snap", "CEMBA201008_9F.snap")

outdir <- "gm"
outfilenm <- "CEMBA201008_9F.h5"
sampleName <- "CEMBA201008_9F"
genome <- "mm10"
excludeChr <- c("chrM")

## get barcodes filtered from smmtools
sumFrag <- read.table(file = sumFragFile, header = TRUE, sep = ",", comment.char = "")
index_from_nFrag <- which(sumFrag$nUniqFrag >= 1000)
index_from_qc <- which(sumFrag$TSSE >= 10)
index_join <- intersect(index_from_nFrag, index_from_qc)
barcodes <- sumFrag[index_join]

gmat <- getGeneMatrix(rawH5File = rawH5File, outdir = outdir, outfilenm = outfilenm,
                      genome = genome, excludeChr = excludeChr, barcodes = barcodes)



