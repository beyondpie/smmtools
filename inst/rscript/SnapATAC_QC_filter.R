suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("SnapATAC"))
library(tictoc)
library(stringr)

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required = TRUE, help = "input snap format file")
parser$add_argument("--fragment_num", default = 1000, help = "cutoff of fragment.num [default %(default)s]")
parser$add_argument("--umi_num", default = 500, help = "cutoff of UMI.num [default %(default)s]")
parser$add_argument("--mito_ratio", default = 1, help = "cutoff of mito.ratio [default %(default)s]")
parser$add_argument("--dup_ratio", default = 1, help = "cutoff of dup.ratio [default %(default)s]")
parser$add_argument("--umap_ratio", default = 0, help = "cutoff of umap.ratio [default %(default)s]")
parser$add_argument("--pair_ratio", default = 0, help = "cutoff of pair.ratio [default %(default)s]")
parser$add_argument("--bin_size", default = 5000, help = "binSize to use [default %(default)s]")
parser$add_argument("--black_list", required = TRUE, help = "black list file")
parser$add_argument("--pc_num", default = 50, help = "num of PCA components [default %(default)s]")
parser$add_argument("--cpu", default = 4, help = "# of cpus [default %(default)s]")
parser$add_argument("--ignoreChroms", default = "chrM|random")
parser$add_argument("--plotBarcodeRawSnapFile", required = TRUE)
parser$add_argument("--plotBarcodeQCFilteredSnapFile", required = TRUE)
parser$add_argument("--binCoverageFile", required = TRUE)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

snapF <- args$input
fragment_num <- as.numeric(args$fragment_num)
umi_num <- as.numeric(args$umi_num)
mito_ratio <- as.numeric(args$mito_ratio)
dup_ratio <- as.numeric(args$dup_ratio)
umap_ratio <- as.numeric(args$umap_ratio)
pair_ratio <- as.numeric(args$pair_ratio)
bin_size <- as.numeric(args$bin_size)
black_list <- args$black_list
pc_num <- as.numeric(args$pc_num)
cpus <- as.numeric(args$cpu)

pnames <- head(unlist(str_split(tail(unlist(str_split(snapF, "/")), 1), ".snap")), 1)
tic("readSnap")
x.sp <- createSnap(
  file = snapF,
  sample = pnames,
  num.cores = cpus
)
x.sp <- addBmatToSnap(
  x.sp,
  bin.size = bin_size,
  num.cores = cpus
)
toc()

## plot Barcode information
plotBarcode(x.sp,
  pdf.file.name = args$plotBarcodeFile, pdf.width = 7, pdf.height = 7,
  col = "grey", border = "grey", breaks = 50
)


tic("filterCells")
x.sp <- filterCells(x.sp,
  subset.names = c(
    "fragment.num",
    "UMI",
    "mito.ratio",
    "dup.ratio",
    "umap.ratio",
    "pair.ratio"
  ),
  low.thresholds = c(fragment_num, umi_num, 0, 0, umap_ratio, pair_ratio),
  high.thresholds = c(Inf, Inf, mito_ratio, dup_ratio, 1, 1)
)
toc()
plotBarcode(x.sp,
  pdf.file.name = args$plotBarcodeQCFilteredSnapFile, pdf.width = 7, pdf.height = 7,
  col = "grey", border = "grey", breaks = 50
)
## Output to stdout
summarySnap(x.sp)

## Matrix Binarization
tic("makeBinary")
x.sp <- makeBinary(x.sp, mat = "bmat")
toc()

## Feature Selection (filter bins according to black list and chrom)
tic("filterBins")
black_list <- read.table(black_list, sep = "\t")
library(GenomicRanges)
black_list.gr <- GRanges(
  black_list[, 1],
  IRanges(black_list[, 2], black_list[, 3])
)
idy1 <- queryHits(findOverlaps(x.sp@feature, black_list.gr))
idy2 <- grep(args$ignoreChroms, x.sp@feature)
idy <- unique(c(idy1, idy2))
x.sp <- x.sp[, -idy, mat = "bmat"]

## Filter Bins
x.sp <- filterBins(x.sp, low.threshold=-2, high.threshold=2, mat="bmat");
plotBinCoverage(
  x.sp,
  pdf.file.name = args$binCoverageFile,
  col = "grey",
  border = "grey",
  breaks = 10,
  xlim = c(-6, 6)
)

# filter empty rows
idx <- which(Matrix::rowSums(x.sp@bmat)!=0)
x.sp@metaData <- x.sp@metaData[idx,]
x.sp@barcode <- x.sp@barcode[idx]
x.sp@sample <- x.sp@sample[idx]
x.sp@file <- x.sp@file[idx]
x.sp@bmat <- x.sp@bmat[idx,]
toc()

## unmodified
outfname = paste(outF, ".pre.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".pre.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)




