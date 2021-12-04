## get gene matrix for doublet removement
library(smmtools)
library(rhdf5)
library(S4Vectors)
library(GenomicRanges)

setwd("/mac/git-recipes/smmtools/eval")
rawH5File <- file.path("rawh5", "4B_rep1_deep_tabix2H5_nChunk_3.h5")
barcodesFile <- file.path("rawh5", "4B_rep1_deep_barcodes_after_QC.txt")
barcodes <- read.table(file = barcodesFile, header = FALSE, comment.char = "")
barcodes <- removeSampleName(barcodes = barcodes, sep = "#")

outdir <- "gm"
outfilenm <- "4B_rep1.h5"

sampleName <- "4B_rep1_deep"
genome <- "mm10"
excludeChr <- c("chrM", "chrY")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outfile <- file.path(outdir, outfilenm)
if(file.exists(outfile)) {
  message(paste(outfile, "exist, and remove it."))
  file.remove(outfile)
}
rhdf5::h5createFile(file = outfile)
annotGenes <- getAnnotFromArchRData(tag = "gene", genome = genome)
## blacklist here, to be test if we need blacklist
annotGenome <- getAnnotFromArchRData(tag = "genome", genome = genome)
blacklist <- annotGenome$blacklist

genes <- annotGenes$genes
genenms <-genes$symbol
gnms <- as.character(genes@seqnames@values)
gnms <- gnms[!(match(gnms, excludeChr, nomatch = 0) > 0)]
a <- lapply(gnms, function(chr) {
  message(paste("current:",chr))
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
  if(nrow(black_ovs) > 1) {
    gf <- gf[-black_ovs$queryHits]
  }
  ovs <- as.data.frame(GenomicRanges::findOverlaps(query = gf, subject = g))
  if(nrow(ovs) < 1) {
    return(NULL)
  }
  ovs_g <- unique(ovs$subjectHits)
  rows <- lapply(ovs_g, function(i) {
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
    return(r)
  })
  rows <- rows[!sapply(rows, is.null)]
  r <- do.call(rbind, rows)
  return(r)
})
a <- a[!sapply(a, is.null)]
wdf <- do.call(rbind, a)
mat <- Matrix::sparseMatrix(i = wdf$i, j = wdf$j, x = wdf$val,
                            dims = c(length(genes), length(barcodes)))

rownames(mat) <- genenms
colnames(mat) <- barcodes

suppressALL(h5createDataset(file = outfile, dataset = "i", storage.mode = "integer", 
                            dims = c(nrow(wdf),1), level = 0))
suppressALL(h5write(obj = wdf$i, file = outfile, name = "i"))

suppressALL(h5createDataset(file = outfile, dataset = "j", storage.mode = "integer",
                            dims = c(nrow(wdf),1), level = 0))
suppressALL(h5write(obj = wdf$j, file = outfile, name = "j"))

suppressALL(h5createDataset(file = outfile, dataset = "val", storage.mode = "integer",
                            dims = c(nrow(wdf),1), level = 0, fillValue = 0))
suppressALL(h5write(obj = wdf$val, file = outfile, name = "val"))

suppressALL(h5write(obj = genenms, file = outfile, name = "gene"))
suppressALL(h5write(obj = barcodes, file = outfile, name = "barcode"))

h5closeAll()

## load h5
h5_barcode <- h5read(file = outfile, name = "barcode")
h5_gene <- h5read(file = outfile, name = "gene")
h5_i <- as.vector(h5read(file = outfile, name = "i"))
h5_j <- as.vector(h5read(file = outfile, name = "j"))
h5_val <- as.vector(h5read(file = outfile, name = "val"))

h5_mat <- Matrix::sparseMatrix(i = h5_i, j = h5_j, x = h5_val,
                               index1 = TRUE,
                               dims = c(length(h5_gene), length(h5_barcode)))
rownames(h5_mat) <- h5_gene
colnames(h5_mat) <- h5_barcode





