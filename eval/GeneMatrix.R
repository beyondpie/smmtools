## get gene matrix for doublet removement
library(rtracklayer)
## remotes::install_github("beyondpie/smmtools")
## detach("package:smmtools", unload=TRUE)
## library(smmtools)

## install.packages("plot3D")
## install.packages("doSNOW")
## remotes::install_github("r3fang/SnapATAC")
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

## load geneup2k bed files
geneUp2K <- rtracklayer::import(file.path("snap", "gencode.vM16.geneUp2k.bed"), format = "bed")
genenms <- geneUp2K$name

## get barcodes filtered from smmtools
sumFrag <- read.table(file = sumFragFile, header = TRUE, sep = ",", comment.char = "")
index_from_nFrag <- which(sumFrag$nUniqFrag >= 1000)
index_from_qc <- which(sumFrag$TSSE >= 10)
index_join <- intersect(index_from_nFrag, index_from_qc)
barcodes <- sumFrag[index_join,"barcode"]

## load snap
snap <- createSnap(file = snapFile, sample = sampleName, num.cores = 2)
snap <- snap[which(snap@barcode %in% barcodes), ]
snap <- SnapATAC::addGmatToSnap(obj = snap)

## get gmat from smmtools
gmat <- getGeneMatrix(rawH5File = rawH5File, outdir = outdir, outfilenm = outfilenm,
                      genes = geneUp2K, genenms = genenms,
                      genome = genome, excludeChr = excludeChr, barcodes = barcodes)



