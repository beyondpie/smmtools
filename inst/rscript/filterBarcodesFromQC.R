library(optparse)

option_list <- list(
  make_option(c("--qcf"), type = "character"),
  make_option(c("--outf"), type = "character"),
  make_option(c("--TSSE"), type = "integer", default = 10),
  make_option(c("--nFrag"), type = "integer", default = 1000),
  make_option(c("--sampleName"), type = "character")
)

## TSSE Ref:
## https://www.encodeproject.org/atac-seq/

args <- parse_args(OptionParser(option_list = option_list))
sampleName <- args$sampleName

sumFrag <- read.table(file = args$qcf, header = TRUE, sep = ",", comment.char = "")

n <- nrow(sumFrag)
message(paste("From QC file of", sampleName, ": total barcodes", n, "."))

index_from_nFrag <- which(sumFrag$nUniqFrag >= args$nFrag)
message(paste("Based nFrag threshold", args$nFrag, length(index_from_nFrag), "are left."))

index_from_qc <- which(sumFrag$TSSE >= args$TSSE)
message(paste("Based TSSE threshold", args$TSSE, length(index_from_qc), "are left."))

index_join <- intersect(index_from_nFrag, index_from_qc)
message(paste("Finally", length(index_join), "barcodes are left."))

invisible(write.table(
  x = sumFrag$barcode[index_join], file = args$outf,
  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ","
))
