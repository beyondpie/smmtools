#' Perform batch effect using harmony
#'
#' @param mapSnapATAC list, output of SnapATAC_runDiffusionMaps, at least contains
#' - dmat matrix, cell by feature, dense matrix of embedding for the cells
#' - sdev numeric vector, standard deviation from runDiffusionMaps.
#' @param meta_data Either (1) Dataframe with variables to integrate or (2) vector with labels.
#' @param nPC integer, how many PCs used for Harmony, default is 30.
#' @param vars_use string, if meta_data is dataframe,
#' this defined which variable(s) to remove (character vector).
#' Default is NULL
#' @param weightDimReduct bool, if weighted the PCs by multiplying sdev, default is FALSE
#' @return Matrix, cell by updated principle components
#' @export
SnapATAC_runHarmony <- function(mapSnapATAC, meta_data,
                                nPC = 30, vars_use = NULL,
                                weightDimReduct = FALSE) {
  dmat <- mapSnapATAC$dmat
  sdev <- mapSnapATAC$sdev
  nCell <- nrow(dmat)
  nDim <- ncol(dmat)
  mat <- dmat[, 1:min(nPC, nDim)]
  if(weightDimReduct) {
    mat <- mat %*% diag(sdev[1:min(nPC, nDim)])
  }
  harmonyEmbed <- harmony::HarmonyMatrix(data_mat = mat, meta_data = meta_data,
                                         vars_use = vars_use, do_pca = FALSE)
  return(harmonyEmbed)
}
