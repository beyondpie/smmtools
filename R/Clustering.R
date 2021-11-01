#' Simplified runKNN in SnapATAC
#' @param nn_eps Error bound when performing nearest neighbor seach using RANN.
#' default of 0.0 implies exact nearest neighbor search
#' @export
runKNN <- function(smat, k = 20, nn_eps = 0.0){
  message(paste("Generate KNN with", k))
  ncell <- nrow(smat)
  if(ncell < k){
    message(paste("Ncell", ncell, "is smaller than K nearst neighbor", k))
    k <- ncell-1
    message("Set k as Ncell - 1.")
  }
  
  nnRanked <- RANN::nn2(data = smat, k = k, searchtype = "stardard",
                        eps = nn_eps)$nn.idx
  j <- as.numeric(t(nnRanked))
  i <- (seq_along(j)-1) %/% k + 1
  kmat <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(ncell, ncell))
  return(kmat)
}

#' @importFrom reticulate r_to_py py_to_r
#' @export
runLeiden <- function(kmat,
                      path_to_python = NULL,
                      reso = 1.0, seed = 10,
                      partitionType = "RB") {
  message(paste("Run Leiden for clustering with resolution", reso, "and partitionType", partitionType))
  if(!is.null(path_to_python)) {
    reticulate::use_python(path_to_python)
    message("Use the Python located in:", path_to_python, "\n")
  }
  setSessionTimeLimit(cpu = Inf, elapsed = Inf)
  ld <- reticulate::import(module = "smmuty", convert = FALSE)
  ldCluster <- as.factor(py_to_r(
    ld$leiden(knn = r_to_py(kmat), reso = reso, seed = seed, opt = partitionType)))
  
  message("Summary of clustering:")
  print(table(ldCluster))
  return(ldCluster)
}
