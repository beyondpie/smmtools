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
sample <- "CEMBA210429_18B"
rawH5File <- file.path("rawh5", paste0(sample,"_raw.h5"))
sumFragFile <- file.path("rawh5", paste0(sample,"_smmtoolsQC.csv"))
snapFile <- file.path("snap", paste0(sample,".snap"))

outdir <- "gm"
outfilenm <- paste0(sample,".h5")
genome <- "mm10"

excludeChr <- c("chrM", "chrY", "chrX", paste0("chr", 1:17))
sampleName <- sample

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
genes <- geneUp2K

annotGenome <- getAnnotFromArchRData(tag = "genome", genome = genome)
blacklist <- annotGenome$blacklist
if (is.null(genes)){
  annotGenes <- getAnnotFromArchRData(tag = "gene", genome = genome)
  genes <- annotGenes$genes
}
if (is.null(genenms)) {
  genenms <- genes$symbol
} else {
  elementMetadata(genes)$symbol <- genenms
}

chrnms <- as.character(genes@seqnames@values)

chr <- "chr1"

g <- genes[genes@seqnames == chr]
fragments <-
  getFragsOfAChrFromRawH5File(
    rawH5File = rawH5File,
    chr = chr,
    sampleName = sampleName,
    barcodes = barcodes
  )
gf <- GRanges(seqnames = chr, ranges = fragments)
black_ovs <- as.data.frame(GenomicRanges::findOverlaps(query = gf, subject = blacklist))
if (nrow(black_ovs) >= 1) {
  gf <- gf[-black_ovs$queryHits]
}
ovs <- as.data.frame(GenomicRanges::findOverlaps(query = gf, subject = g))
if (nrow(ovs) < 1) {
  return(NULL)
}


## test original logic
i <- 1
fid <- ovs$queryHits[ovs$subjectHits == i]
f <- fragments[fid]
frle <- as.data.frame(table(mcols(f)$RG))
colnames(frle) <- c("values", "lengths")
frle$values <- as.character(frle$values)
rownames(frle) <- frle$values
cols <- S4Vectors::match(frle$values, barcodes, nomatch = 0)
cols <- cols[!(cols == 0)]
values <- frle[barcodes[cols], "lengths"]
r <- data.frame(
  i = rep(match(g[i]$symbol, genenms), length(cols)),
  j = cols,
  val = values
)

## generate gmat from rawh5
library(rhdf5)
library(S4Vectors)
library(GenomicRanges)

genes <- read.table(file.path("snap", "gencode.vM16.geneUp2k.bed"), header = FALSE)
geneGR <- GenomicRanges::GRanges(seqnames = genes[,1],
                                 ranges = IRanges::IRanges(start = genes[,2], end = genes[,3]))

gmat_smmtools <- getGeneMatrix(rawH5File = rawH5File, outfile = outfile,
                               genes = geneGR, genenms = genes[,4],
                               barcodes = barcodes, excludeChr = excludeChr)

gmat_smmtools <- t(as.matrix(gmat_smmtools))

## gmat_snapatac <- getGmatFromSnap(snapFile = file.path("snap/", paste0(sampleName, ".snap")),
##                         geneBedFile = file.path("../data/gencode.vM16.geneUp2k.bed"),
##                         ibarcode = barcodes, blacklistBedFile = "../data/mm10.blacklist.bed",
##                         outfile = file.path("gmat/", paste0(sampleName, ".h5")))
gmat_snapatac <- getGmatFromSnap(snapFile = file.path("snap/", paste0(sampleName, ".snap")),
                        ibarcode = barcodes,
                        outfile = file.path("gmat/", paste0(sampleName, ".h5")))
## compare with smmtools
colindex <- which(colSums(gmat_smmtools) > 1)
cols <- intersect(colnames(gmat_smmtools)[colindex], colnames(gmat_snapatac))
gmat_snapatac <- as.matrix(gmat_snapatac[rownames(gmat_smmtools), cols])
gmat_snapatac <- gmat_snapatac / rowSums(gmat_snapatac)
gmat_smmtools <- gmat_smmtools[, cols]
gmat_smmtools <- gmat_smmtools / rowSums(gmat_smmtools)


