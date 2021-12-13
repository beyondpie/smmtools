#' 

#' @export
runUMAP <- function(mapSnapATAC, nPC = 20, seed = 10, weightDimReduct = FALSE) {
  dmat <- mapSnapATAC$dmat
  sdev <- mapSnapATAC$sdev
  nDim <- ncol(dmat)
  mat <- dmat[, 1:min(nPC, nDim)]
  if (weightDimReduct) {
    mat <- mat %*% diag(sdev[1:min(nPC, nDim)])
  }
  set.seed(seed = seed)
  message("Run UMAP")
  umap <- umap::umap(mat, quietly = TRUE)$layout
  colnames(umap) <- c("UMAP-1", "UMAP-2")
  return(umap)
}

#' Umap based on ArchR default parameters
#'
#' @export
SnapATAC_getArchRUmap <- function(mapSnapATAC, nPC = 30, seed = 10, weightDimReduct = FALSE,
                            umapNeighbors = 40, umapMindist = 0.4,
                            umapMetric = "cosine",
                            umapVerbose = FALSE,
                            umapFastSGD = TRUE,
                            umapNthread = 2
                            ) {
  dmat <- mapSnapATAC$dmat
  sdev <- mapSnapATAC$sdev
  nDim <- ncol(dmat)
  mat <- dmat[, 1:min(nPC, nDim)]
  if (weightDimReduct) {
    mat <- mat %*% diag(sdev[1:min(nPC, nDim)])
  }
  set.seed(seed = seed)
  uwotUmap <- uwot::umap(X = mat, n_neighbors = umapNeighbors,
                         metric = umapMetric, n_threads = umapNthread,
                         min_dist = umapMindist,
                         ret_model = FALSE, ret_nn = FALSE)
  return(uwotUmap)
}

#' Umap in SnapATAC.
#'
#' @export
SnapATAC_getDefaultUmap <- function(mapSnapATAC, nPC = 30, seed = 10){
}
