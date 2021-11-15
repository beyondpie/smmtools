library(ggplot2)
library(stringr)
library(data.table)
library(uwot)

## analysis QC: nfragment and TSSE

qc_cancers <- read.table(file = "cancer_tissue_qc.tsv", header = TRUE,
                         quote = "", comment.char = "")
qc_QY1287_Kai <- qc_cancers[qc_cancers$sample == "QY_1287",]

qc_QY1287_SnapATAC <- read.table(file = "out/QY_1287_QC.tsv", header = TRUE,
                                 quote = "", comment.char = "", sep = ",")
rownames(qc_QY1287_SnapATAC) <- qc_QY1287_SnapATAC$barcode
tsse_QY1287_smmtools <- read.table(file = "out/QY_1287_sumFragment.csv", header = TRUE,
                                   quote = "", comment.char = "", sep = ",")
rownames(tsse_QY1287_smmtools) <- tsse_QY1287_smmtools$barcode

fragment_bam210x <- data.frame(tsse_QY1287_smmtools$nUniqFrag,qc_QY1287_SnapATAC[rownames(tsse_QY1287_smmtools), "uniq"])
colnames(fragment_bam210x) <- c("smmtools", "SnapATAC")

plot(x = log10(fragment_bam210x$smmtools), y = log10(fragment_bam210x$SnapATAC))

## visualize the clustering result
