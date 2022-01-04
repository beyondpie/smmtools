#' Umap based on ArchR and Seuratdefault parameters
#'
#' Use uwot pacakge since it's more stable in R.
#' Uwot::umap default parameters are from ArchR.
#'
#' Ref: https://github.com/jlmelville/uwot
#' Ref: https://satijalab.org/seurat/reference/runumap
#' 
#' @param dmat matrix, cell by feature
#' @param sdev vector, scaling factors for dims, default is NULL
#' @param nPC integer, dims till the nPC-th dim will be used for UMAP, default is 30.
#' When nfeature in mat is less than nPC, nfeature will be used.
#' @param ncomp integer, defaul is 2
#' @param seed integer, default is 10
#' @param weightDimReduct bool, default is FALSE, together used with sdev.
#' @param umapNeighbors integer, default is 30
#' In umot::umap, default is 15.
#' @param umapMindist double, default is 0.3
#' In umot::umap, default is 0.01
#' @param umapMetric string, default is "cosine"
#'"euclidean" is the default in uwot::umap
#' @param umapVerbose bool, default is FALSE
#' @param umapFastSGD bool, default is TRUE
#' @param umapNthread integer, default is 1
#' @param a double, default is NULL
#' @param b double, default is NULL
#' NOTE: a = 1.8956, b = 0.8006 is the default UMAP
#' @param init string, default is "special"
#' NOTE: Use init = "spca" rather than init = "spectral" for better reproducibility
#' (although the latter is the default and preferred method for UMAP initialization).
#' @param retnn bool, default is FALSE
#' @param retmodel bool, default is FALSE
#' Set both retnn and retmodel as TRUE to save model
#' @return matrix, cell by ncomp
#' @export
runUmap <- function(dmat, sdev = NULL, nPC = 30,
                    ncomp = 2,
                    seed = 10, weightDimReduct = FALSE,
                    umapNeighbors = 30, umapMindist = 0.3,
                    umapMetric = "cosine",
                    umapVerbose = FALSE,
                    umapFastSGD = TRUE,
                    umapNthread = 1,
                    a = NULL,
                    b = NULL,
                    init = "spectral",
                    retnn = FALSE,
                    retmodel = FALSE) {
  nDim <- ncol(dmat)
  mat <- dmat[, 1:min(nPC, nDim)]
  if (weightDimReduct & (!is.null(sdev))) {
    message("Scale the dims based on sdev.")
    mat <- mat %*% diag(sdev[1:min(nPC, nDim)])
  }
  set.seed(seed = seed)
  uwotUmap <- uwot::umap(X = mat, n_neighbors = umapNeighbors,
                         metric = umapMetric, n_threads = umapNthread,
                         min_dist = umapMindist,
                         a = a, b = b, init = init,
                         ret_nn = retnn, ret_model = retmodel)
  return(uwotUmap)
}

