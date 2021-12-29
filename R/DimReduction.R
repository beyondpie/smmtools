#' DiffusionMaps method in SnapATAC
#'
#' @param bmatSnap sparse matrix, cell by feature
#' @param nLandmark integer, numebr of landmark cells, default is 10,000
#' @param n integer, used for fitting regression of random depth effect, default is 1000
#' @param outlier double, used for choosing threshold of normalized Jaccard similarity matrix,
#' default is 0.999
#' @param nPC integer, default is 30
#' @param seed integer, default is 1
#' @return list of two element
#' - dmat dense matrix, cell by principle components, keep the same order of cells
#' - sdev vector, length of principle components
#' @export
SnapATAC_DiffusionMaps <- function(bmatSnap, nLandmark = 10000,
                                   n = 1000, outlier = 0.999,
                                   nPC = 30, seed = 1) {
  set.seed(seed)
  nCell <- nrow(bmatSnap)

  idxLandmark <- sampleBasedOnDepth(bmat = bmatSnap, n = nLandmark)
  
  bmatLandmark <- bmatSnap[idxLandmark, ]
  mapLandmark <- SnapATAC_runDiffusionMaps(bmat = bmatLandmark,
                                           nPC = nPC, n = n, outlier = outlier)
  message("Finish getting landmark mappings.")
  rownames(mapLandmark$dmat) <- rownames(bmatSnap)[idxLandmark]
  if (nCell > nLandmark) {
    bmatQuery <- bmatSnap[-idxLandmark, ]
    mapQuery <- SnapATAC_runDiffusionMapsExtension(
      bmatLandmark = bmatLandmark,
      bmatQuery  = bmatQuery,
      mapLandmark = mapLandmark
    )
    message("Finish projecting query on landmark space.")
    rownames(mapQuery$dmat) <- rownames(bmatSnap)[-idxLandmark]
    mapMerge <- rbind(mapLandmark$dmat, mapQuery$dmat)
    mapReorder <- mapMerge[rownames(bmatSnap),]
  } else {
    mapReorder <- mapLandmark$dmat[rownames(bmatSnap),]
  }
  return(invisible(list(dmat = mapReorder, sdev = mapLandmark$sdev)))
}

#' Core of DiffusionMaps for landmark cells.
#' 
#' @param bmat sparse Matrix, cell by feature
#' @param nPC integer, default is 30
#' @param n integer, used for fitting regression of random depth effect, default is 1000.
#' @param outlier double, used for choosing threshold of normalized Jaccard similarity matrix,
#' default is 0.999
#' @return list of five elements
#' - dmat dense matrix, cell by features
#' - sdev vector, standard deviation
#' - cutoff numeric, cutoff of outiliers for the normalized Jaccard similarity matrix
#' - betas numeric vector of three elements, regression coefficients for reducing the random
#' depth effect in Jaccard similarity matrix.
#' - diagFactor Diagonal Matrix, inverse of square root of
#' the rowSums of normalized Jaccard similarity matrix.
#' @importFrom Matrix diag diag<-
#' @export
SnapATAC_runDiffusionMaps <- function(bmat, nPC = 30, n = 1000, outlier = 0.999) {
  message("Running SnapATAC Diffusion Maps...")
  if(any(Matrix::rowSums(bmat) == 0)) {
    stop("Some cells have no reads, please remove them firstly.")
  }
  message("Step1: calculate Jaccard similarity matrix.")
  A <- Matrix::tcrossprod(bmat, bmat)
  rowDepth <- Matrix::rowSums(A)
  J <- as.matrix(A / (2 * replicate(ncol(A), rowDepth) - A ))

  message("Step2: fit regression model to reduce cell depth effect.")
  idx <- sampleBasedOnDepth(bmat = bmat, n = 1000)
  subrmeans <- Matrix::rowMeans(bmat[idx, ])
  subNormOVE <- getNormOVE(p1 = subrmeans, p2 = subrmeans)
  jmat <- J[idx, idx]
  fitData <- data.frame(x = subNormOVE[upper.tri(subNormOVE)], y = jmat[upper.tri(jmat)])
  model <- stats::lm(formula = y ~ x + I(x^2), data = fitData)
  betas <- as.numeric(model$coefficients)
  message("Regression coefficients: ", paste(round(betas, 5), collapse = ", "))

  message("Step3: normalize Jaccard similarity matrix by reducing the fitted random depth.")
  rmean <- Matrix::rowMeans(bmat)
  allNormOVE <- getNormOVE(p1 = rmean, p2 = rmean)
  ## pred is a matrix
  pred <- betas[1]  + betas[2] * allNormOVE + betas[3] * (allNormOVE ** 2)
  normJ <- J / pred

  message("Step4: set cutoffs for the outliers of the normalized Jaccard similarity matrix.")
  cutoff <- quantile(normJ, outlier)
  normJ[normJ > cutoff] <- cutoff

  message("Step5: get eigin decomposition by igraph arpack.")
  diag(normJ) <- 0.0
  rsum <- Matrix::rowSums(normJ)
  diagFactor <- Matrix::Diagonal(x = 1 / sqrt(rsum))
  transition <- as.matrix(diagFactor %*% normJ %*% diagFactor)
  diag(transition) <- 0
	eig_transitions <- eig_decomp(transition, n_eigs = nPC)

  message(paste("Return the diffusion map result: dim from 2 to", nPC + 1))
  dmat <- eig_transitions$vectors[, 2:(nPC+1)]
  sdev <- eig_transitions$value[2:(nPC + 1)]
  return(invisible(list(dmat = dmat, sdev = sdev,
                        cutoff = cutoff, betas = betas,
                        diagFactor = diagFactor)))
}

#' Run DiffusionMaps for query cells baesd on the landmark results.
#' 
#' @param bmatLandmark matrix, cell by feature
#' @param bmatQuery matrix, cell by feature
#' @param mapLandmark results from difusionMaps for landmark cells.
#' @return list of two elements
#' - dmat matrix, nquery cell by feature
#' - sdev numeric vector, standard deviation for the query cells.
#' @export
SnapATAC_runDiffusionMapsExtension <- function(bmatLandmark, bmatQuery, mapLandmark) {
  if(ncol(bmatLandmark) != ncol(bmatQuery)) {
    stop("Features have different dims between Landmark and Query.")
  }
  message("Step1: calculate Jaccard similarity matrix.")
  A <- Matrix::tcrossprod(x = bmatQuery, y = bmatLandmark)
  rsumQuery <- Matrix::rowSums(bmatQuery)
  rsumLandmark <- Matrix::rowSums(bmatLandmark)
  J <- as.matrix(A / (replicate(ncol(A), rsumQuery) + t(replicate(nrow(A), rsumLandmark))-A))
  message("Step2: normalze Jaccard similarity matrix.")
  rmeanQuery <- Matrix::rowMeans(bmatQuery)
  rmeanLandmark <- Matrix::rowMeans(bmatLandmark)
  normOVE <- getNormOVE(p1 = rmeanQuery, p2 = rmeanLandmark)
  betas <- mapLandmark$betas
  pred <- betas[1] + betas[2] * normOVE + betas[3] * (normOVE ** 2)
  normJ <- J / pred
  message("Step3: set cutoffs for the outliers of the normalized jaccard similarity matrix.")
  cutoffLandmark <- mapLandmark$cutoff
  normJ[normJ > cutoffLandmark] <- cutoffLandmark
  message("Step4: project query to landmark map space.")
  diagFactorLandmark <- mapLandmark$diagFactor
  diagFactorQuery <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(normJ, na.rm = TRUE)))
  transition <-  diagFactorQuery %*% normJ %*% diagFactorLandmark
  dmatLandmark <- mapLandmark$dmat
  sdevLandmark <- mapLandmark$sdev
  dmatQuery <- as.matrix(t( t(as.matrix(transition) %*% dmatLandmark) / sdevLandmark ))
  sdevQuery <- sdevLandmark
  return(invisible(list(dmat = dmatQuery, sdev = sdevQuery)))
}
