library(optparse)
library(smmtools)

option_list <- list(
  make_option(c("--sampleNames"), type = "character",
              help = "sample names corresponding to tileMatrix files.",
              default = "QY_1287,QY_1288"),
  make_option(c("--tileMatrixDir"), type = "character",
              help = "directory storing the tileMatrix files.",
              default = "out"),
  make_option(c("--tileMatrixFNames"), type = "character",
              help = "file names seprated by comma.",
              default = "QY_1287_bmat.h5,QY_1288_bmat.h5"),
  make_option(c("--barcodeDir"), type = "character",
              help = "directory storing the barcode files.",
              default = "out"),
  make_option(c("--barcodeFNames"), type = "character",
              help = "barcode file names corrresponding tileMatrix files",
              default = "QY_1287_barcodes_after_doublet.csv,QY_1288_barcodes_after_doublet.csv"),
  make_option(c("--nLandmark"), type = "integer", default = 10000),
  make_option(c("--nSampleRegrn"),
    type = "integer", default = 2000,
    help = "nsamples for fitting regression to normalized jaccard matrix",
  ),
  make_option(c("--nPCDiffusionMaps"), type = "integer", default = 50,
              help = "n of principle components for diffusion map"),
  make_option(c("--nPCHarmony"), type = "integer", default = 30,
              help = "n of top components used for harmony"),
  make_option(c("--saveHarmonyFile"), type = "character", default = "./tmp.rds")
)

args <- parse_args(OptionParser(option_list = option_list))

## load matrix with filter barcodes
sampleNames <- unlist(strsplit(args$sampleNames, split = ",", fixed = TRUE))
tileMatrixFVec <- file.path(
  args$tileMatrixDir,
  unlist(strsplit(args$tileMatrixFNames, split = ",", fixed = TRUE))
)
barcodeFVec <- file.path(
  args$barcodeDir,
  unlist(strsplit(args$barcodeFNames, split = ",", fixed = TRUE))
)

if (length(tileMatrixFVec) != length(sampleNames)) {
  stop(paste(length(sampleNames), "sample names do not match", lenegth(tileMatrixFVec), "tileMatrix."))
}

if (length(tileMatrixFVec) != length(barcodeFVec)) {
  stop(paste(length(tileMatrixFVec), "tileMatirx vs", length(barcodeFVec), "barcode files."))
}

## cost large memory:
## two matrix roughly 4.4G
tileMatrixList <- lapply(seq_along(tileMatrixFVec), FUN = function(i) {
  barcodes <- read.csv(file = barcodeFVec[i], header = FALSE)[,1]
  message(paste("Loading tileMatrix:", tileMatrixFVec[i], length(barcodes),"nBarcodes."))
  tileMatrix <- loadTileMatrix(
    tileMatrixFile = tileMatrixFVec[i],
    barcodes = barcodes,
    binarize = FALSE
  )
  colnames(tileMatrix) <- paste0(sampleNames[i], "#", barcodes)
  return(tileMatrix)
})

## feature by cell
tileMatrix <- do.call(what = "cbind", args = tileMatrixList)
rm(tileMatrixList)
barcodes_with_sampleNames <- colnames(tileMatrix)
message(paste("After cbind the matrix,", length(barcodes_with_sampleNames), "barcodes are got in total."))

## binarize
## cell by matrix
bTileMatrix <- SnapATAC_BinarizeBmat(bmat = tileMatrix, z_threshold = 1.65, outlier = 1e-3)
message("Finish SnapATAC_BinarizeBmat.")
rm(tileMatrix)

## diffusion map

dimReductByDiffusionMaps <- SnapATAC_DiffusionMaps(
  bmatSnap = bTileMatrix$bmat,
  nLandmark = args$nLandmark,
  nPC = args$nPCDiffusionMaps,
  seed = 1
)
rm(bTileMatrix)
message(paste("Finish  diffusionmaps with", args$nLandmark, "landmarks", args$nPCDiffusionMaps, "dims."))
## run harmony
harmonyEmbed <- SnapATAC_runHarmony(
  mapSnapATAC = dimReductByDiffusionMaps,
  nPC = args$nPCHarmony,
  meta_data = gsub(pattern = "#.+", replacement = "", barcodes_with_sampleNames),
  vars_use = NULL,
  weightDimReduct = FALSE
)
rownames(harmonyEmbed) <- barcodes_with_sampleNames
message(paste("Finish SnapATAC_runHarmony", args$nPCHarmony, "dims are used."))
## save file
saveRDS(object = harmonyEmbed, file = args$saveHarmonyFile)


## run harmony locally on imac
## args$sampleNames <- "2C_rep1_deep,2C_rep2_deep,3C_rep1_deep,3C_rep2_deep,4B_rep1_deep,4B_rep2_deep,5D_rep1_deep,5D_rep2_deep"
## args$tileMatrixFNames <- "2C_rep1_deep_tileMatrix_5000.h5,2C_rep2_deep_tileMatrix_5000.h5,3C_rep1_deep_tileMatrix_5000.h5,3C_rep2_deep_tileMatrix_5000.h5,4B_rep1_deep_tileMatrix_5000.h5,4B_rep2_deep_tileMatrix_5000.h5,5D_rep1_deep_tileMatrix_5000.h5,5D_rep2_deep_tileMatrix_5000.h5"

## args$barcodeFNames <- "2C_rep1_deep_filteredBarcodesFromDoublet.csv,2C_rep2_deep_filteredBarcodesFromDoublet.csv,3C_rep1_deep_filteredBarcodesFromDoublet.csv,3C_rep2_deep_filteredBarcodesFromDoublet.csv,4B_rep1_deep_filteredBarcodesFromDoublet.csv,4B_rep2_deep_filteredBarcodesFromDoublet.csv,5D_rep1_deep_filteredBarcodesFromDoublet.csv,5D_rep2_deep_filteredBarcodesFromDoublet.csv"

