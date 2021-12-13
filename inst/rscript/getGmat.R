library(optparse)
## require rtracklayer package

option_list <- list(
  make_option(c("--rawH5File"), type = "character"),
  make_option(c("--outdir"), type = "character"),
  make_option(c("--outfnm"), type = "character"),
  make_option(c("--geneBedFile"), type = "character"),
  make_option(c("--excludeChrs"), type = "character", default = "chrM", help = "string seperated by comma"),
  make_option(c("--genome"), type = "character", default = "mm10")
  make_option(c("--barcodeFile"), type = "character", help = "no header, just one column of barcodes.")
)

args <- parse_args(OptionParser(option_list = option_list))

genes <- rtracklayer::import(args$geneBedFile, format = "bed")
genenms <- genes$name
barcodes <- as.vector(read.table(file = args$barcodeFile, header = FALSE, comment.char = ""))
excludeChr <- unlist(strsplit(args$excludeChr, split = ","))


invisible(smmtools::getGeneMatrix(rawH5File = args$rawH5File,
                        outdir = args$outdir,
                        outfilenm = args$outfnm,
                        genes = genes, genenms = genenms,
                        genome = args$genome,
                        excludeChr = excludeChr)

