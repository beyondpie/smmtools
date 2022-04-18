#' Conver snap object to Seurat by gmat and non-linear embedding
#' provided by SnapATAC
#'
#' @param snap snap object defined in SnapATAC package
#' snap@gmat must be not empty; snap@smat@dmat should not be empty 
#' @param eigDims vector, used for choosing PCA components, default 1:50
#' @param assay characters, name used in Seurat object
#' @param pcaPrefix characters, default "SnapATAC_"
#' @return Seurat object
#' @import Seurat
#' @export
snapGmat2Seurat <- function(snap, eigDims = 1:50,
                            assay = "GeneScore",
                            pcaPrefix = "SnapATAC_") {
  # check snap@gmat
  # check snap@smat@dmat
  pcaUse <- snap@smat@dmat[, eigDims]
  metaData <- snap@metaData
  rownames(pcaUse) <- paste(metaData$sample, metaData$barcode, sep = ".")
  colnames(pcaUse) <- paste0(pcaPrefix, 1:ncol(pcaUse))
  rownames(metaData) <- paste(metaData$sample, metaData$barcode, sep = ".")
  gmatUse <- t(snap@gmat)
  colnames(gmatUse) <- paste(metaData$sample, metaData$barcode, sep = ".")
  snapSeurat <- CreateSeuratObject(counts = gmatUse, assay = assay)
  snapSeurat <- AddMetaData(snapSeurat, metadata = metaData)
  snapSeurat[["pca"]] <- new(Class = "DimReduc", cell.embeddings = pcaUse,
                             feature.loadings = matrix(0,0,0), feature.loadings.projected = matrix(0,0,0),
                             assay.used = assay, stdev = rep(1, ncol(pcaUse)),
                             key = pcaPrefix, jackstraw = new(Class = "JackStrawData"), misc = list())
  return(snapSeurat)
}
