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
    feature.loadings = matrix(0, 0, 0),
    feature.loadings.projected = matrix(0, 0, 0),
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
#' Assume meta.data has a column for rnaType
#' @param eigDims dims used for snapSeurat pca,
#' if NULL (default), use all the dims.
#' @param snapAssay characters, snapSeurat assay name, default "GeneScore"
#' @param preprocessSnap bool, if need to normalize and scale snap, default TRUE
#' @param preprocessRNA bool, if need to normalize, scale rnaSeurat, default TRUE
#' @param rnaTypeColnm characters, colname for rnaSeurat's type, default "ClusterName"
#' @param reso numeric, clustering resolution parameter, default 0.5
#' @return Seurat object, merge snapSeurat and rnaSeurat
#' Add predict Id from rnaSeurat for snapSeurat on metaData,
#' Add imputed gene expression from rnaSeurat for snapSeurat with Assay: "RNA"
#' Clustering on the merged Seurat object as Co-Embedding Seurat object
#' @import ggplot2 Seurat
#' @export
integrateWithScRNASeq <- function(snapSeurat,
                                  rnaSeurat,
                                  eigDims = NULL,
                                  snapAssay = "GeneScore",
                                  preprocessSnap = TRUE,
                                  preprocessRNA = TRUE,
                                  rnaTypeColnm = "ClusterName",
                                  reso = 0.5) {
  if (nrow(snapSeurat) < 1) {
    stop("No cells in snapSeurat.")
  }
  if (nrow(rnaSeurat) < 1) {
    stop("No cells in rnaSeurat.")
  }
  DefaultAssay(snapSeurat) <- snapAssay
  if (is.null(eigDims)) {
    eigDims <- 1:nrow(snapSeurat[["pca"]])
  }
  if (preprocessSnap) {
    snapSeurat <- NormalizeData(snapSeurat)
    snapSeurat <- ScaleData(snapSeurat, features = rownames(snapSeurat))
  }
  if (preprocessRNA) {
    rnaSeurat <- NormalizeData(rnaSeurat)
    rnaSeurat <- FindVariableFeatures(rnaSeurat)
    rnaSeurat <- ScaleData(rnaSeurat, features = rownames(rnaSeurat))
    rnaSeurat <- RunPCA(rnaSeurat, features = VariableFeatures(rnaSeurat))
  }
  anchors <- FindTransferAnchors(
    reference = rnaSeurat,
    query = snapSeurat,
    features = VariableFeatures(rnaSeurat),
    reference.assay = "RNA",
    query.assay = snapAssay,
    reduction = "cca"
  )
  ## predict type
  transferLabel <- TransferData(
    anchorset = anchors,
    refdata = rnaSeurat@meta.data[, rnaTypeColnm],
    weight.reduction = snapSeurat[["pca"]],
    dims = eigDims
  )
  rownames(transferLabel) <- colnames(snapSeurat)
  t <- data.frame(
    predictId = transferLabel$predicted.id,
    predictMaxScore = apply(transferLabel[, -1], 1, max),
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
  coEmbed$tech <- ifelse(!is.na(coEmbed@meta.data[, rnaTypeColnm]), "RNA", "ATAC")
  coEmbed <- ScaleData(coEmbed,
    features = VariableFeatures(rnaSeurat),
    do.scale = FALSE)
  coEmbed <- RunPCA(coEmbed,
    features = VariableFeatures(rnaSeurat),
    verbose = FALSE)
  coEmbed <- RunUMAP(coEmbed, dims = eigDims)
  coEmbed <- FindNeighbors(coEmbed, dims = eigDims)
  coEmbed <- FindClusters(coEmbed, resolution = reso)
  return(coEmbed)
}

getPredictId <- function(ovlp, atacLabel = "MajorType", coEmbed) {
  nms <- rownames(ovlp)
  p1 <- vapply(nms, function(i) {
    colnames(ovlp)[which.max(ovlp[i, ])]
  }, "")
  p2 <- vapply(nms, function(i) {
    pred <- coEmbed$predictId[coEmbed[[atacLabel]] %in% i]
    stat <- table(pred)
    names(stat)[which.max(stat)]
  }, "")
  data.frame(fromOvlp = p1, fromCoEmbed = p2,
    row.names = nms)
}


# t1: table with 2 columns: coembed labels, raw labels
# t2: table with 2 columns: coembed labels, raw labels
cal_ovlpScore <- function(t1, t2) {
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x) {
    x / sum(x)
  })
  t2.pct <- apply(t2.table, 2, function(x) {
    x / sum(x)
  })
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(
    anno1 = as.character(),
    anno2 = as.character(), ovlpScore = as.numeric()
  )
  for (t1.label in t1.labels) {
    for (t2.label in t2.labels) {
      t1.pct.df <- data.frame(t1.pct[, t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[, t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- join(t1.pct.df, t2.pct.df, by = "ident", type = "full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1 = t1.label, anno2 = t2.label, ovlpScore = ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}

getOverlapMatrix <- function(metaF,
                             atacLabel = "MajorType", rnaLabel = "ClusterName") {
  ident2rna <- data.frame(idents = metaF$coembed.idents, rna_label = metaF[[rnaLabel]])
  ident2rna <- ident2rna[complete.cases(ident2rna), ]
  ident2atac <- data.frame(idents = metaF$coembed.idents, atac_label = metaF[[atacLabel]])
  ident2atac <- ident2atac[complete.cases(ident2atac), ]
  ovlp <- cal_ovlpScore(ident2rna, ident2atac)
  colnames(ovlp) <- c(rnaLabel, atacLabel, "ovlpScore")
  ovlp <- reshape2::dcast(ovlp, as.formula(paste(atacLabel, "~", rnaLabel, sep = " ")),
    value.var = "ovlpScore", fun.aggregate = identity, fill = 0.0)
  ovlp <- as.data.frame(ovlp)
  rnm <- ovlp[[atacLabel]]
  rownames(ovlp) <- rnm
  ## remove a name vector
  ovlp[[atacLabel]] <- NULL
  ovlp <- as.matrix(sapply(ovlp, as.numeric))
  rownames(ovlp) <- rnm
  return(ovlp)
}
