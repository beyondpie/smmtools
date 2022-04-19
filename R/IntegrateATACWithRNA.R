#' Convert snap object to Seurat by gmat and non-linear embedding
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
                             feature.loadings = matrix(0,0,0),
                             feature.loadings.projected = matrix(0,0,0),
                             assay.used = assay, stdev = rep(1, ncol(pcaUse)),
                             key = pcaPrefix,
                             jackstraw = new(Class = "JackStrawData"), misc = list())
  return(snapSeurat)
}

#' Integration with single-cell RNA sequencing data
#'
#' @param snapSeurat Seurat object transformed from snap object
#' Assume meta.data has no column named ClusterName
#' @param rnaSeurat Seurat obejct, reference scRNA sequencing ddta
#' Assume meta.data has a column named ClusterName
#' @param outprefix characters, used for naming output file
#' @param outdir characters, if dir not exists with create it
#' @param snapAssay characters, snapSeurat assay name, default "GeneScore"
#' @return 
#' @import ggplot2 Seurat
#' @export
integrateWithScRNASeq <- function(snapSeurat,
                                  rnaSeurat,
                                  outprefix,
                                  outdir,
                                  eigDims = NULL,
                                  snapAssay = "GeneScore",
                                  preprocessSnap = TRUE,
                                  preprocessRNA = TRUE) {
  if (nrow(snapSeruat) < 1) {
    stop("No cells in snapSeurat.")
  }
  if (nrow(rnaSeurat) < 1 ) {
    stop("No cells in rnaSeurat.")
  }
  dir.create(path = outdir, showWarnings = FALSE)
  DefaultAssay(snapSeurat) <- snapAssay
  if(is.null(eigDims)) {
    eigDims <- 1:nrow(snapSeruat[["pca"]])
  }
  if (preprocessSnap) {
    snapSeurat <- NormalizeData(snapSeurat)
    snapSeurat <- ScaleData(snapSeurat, features = rownames(snapSeurat))
  }
  if (preprocessRNA) {
    rnaSeurat <- NormalizeData(rnaSeurat)
    rnaSeurat <- FindVariableFeatures(rnaSeurat)
    rnaSeurat <- ScaleData(rnaSeruat, features = rownames(rnaSeurat))
    rnaSeurat <- RunPCA(rnaSeurat, features = VariableFeatures(rnaSeurat))
  }
  anchors <- FindTransferAnchors(
    reference = rnaSeruat,
    query = snapSeurat, 
    features = VariableFeatures(rnaSeurat),
    reference.assay = "RNA",
    query.assay = snapAssay,
    reduction = "cca"
  )
  ## predict type
  transferLabel <- TransferData(
    anchorset = anchors,
    refdata = rnaSeurat,
    weight.reduction = snapSeurat[["pca"]],
    dims = eigDims
  )
  rownames(transferLabel) <- colnames(snapSeurat)
  t <- data.frame(
    predictId = transferLabel$predicted.id,
    predictMaxScore = apply(transferLabel[,-1], 1, max),
    row.names = colnames(snapSeurat)
  )
  snapSeurat <- AddMetaData(snapSeurat, metadata = t)
  ## impute gene expression
  refdata <- GetAssayData(
    object = rnaSeurat,
    assay = "RNA",
    slot = "data"
  )
  geneImpute <- TransferData(
    anchorset = anchors,
    refdata = refdata,
    weight.reduction = snapSeurat[["pca"]],
    dims = eigDims
  )
  snapSeurat[["RNA"]] <- CreateAssayObject(counts = geneImpute@data)
  rm(geneImpute)
  ## coEmbed
  coEmbed <- merge(x = snapSeurat, y = rnaSeurat)
  DefaultAssay(coEmbed) <- "RNA"
  coEmbed$tech <- with(coEmbed@meta.data, {
    ifelse(!is.na(ClusterName), "RNA", "ATAC")
  })
  ## save predict score histogram
  with_pdf(
    new = file.path(outdir,
                    paste(outprefix, "cellTypePredictScore.pdf", sep = ".")),
    code = {
      hist(snapSeurat$predictMaxScore,
        xlim = c(0, 1),
        xlab = "Cell Type Prediction Score",
        main = "Histogram of Cell Type Prediction Scores from scRNA-seq data")
      abline(v = 0.5, col = "red", lwd = 2, lty = 2)
    })

  ## output co-embedding tech labels
    p <- DimPlot(coEmbed, group.by = "tech") +
      ggtitle("Co embedding with technique label")
    ggsave(
      filename = file.path(outputDir,
        paste(prefix, "coEmbed.tech.umap.pdf", sep = ".")),
      plot = p)
  return(coEmbed)
}


