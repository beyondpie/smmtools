#' @param bmat sparse matrix, feature by cell
#' @return dense matrix, cell by principle components, keep the same order of cells
#' @export dimension reduction in SnapATAC
SnapATAC_DiffusionMaps <- function(bmat, nLandmark = 10000, nPC = 30, seed = 1) {
  set.seed(1)
  ## change to cell by feature bmat used in SnapATAC
  bmatSnap <- Matrix::t(bmat)
  nCell <- nrow(bmatSnap)

  idxLandmark <- sampleBasedOnDepth(bmat = bmatSnap, n = nLandmark)
  
  bmatLandmark <- bmatSnap[idxLandmark, ]
  mapLandmark <- SnapATAC_runDiffusionMaps(bmat = bmatLandmark, nPC = nPC)
  rownames(mapLandmark) <- rownames(bmatSnap)[idxLandmark]
  if (nCell > nLandmark) {
    bmatQuery <- bmatSnap[-idxLandmark, ]
    mapQuery <- SnapATAC_runDiffusionMapsExtension(
      bmatLandmark = bmatLandmark,
      bmatQuery  = bmatQuery,
      mapLandmark = mapLandmark
    )
    rownames(mapQuery) <- rownames(bmatSnap)[-idxLandmark]
    mapMerge <- rbind(mapLandmark, mapQuery)
    mapReorder <- mapMerge[rownames(bmatSnap),]
    return(mapReorder)
  } else {
    return(mapLandmark[rownames(bmatSnap),])
  }
}

#' @param bmat sparse Matrix, cell by feature
#' @importFrom Matrix %*% diag diag<-
#' @export
SnapATAC_runDiffusionMaps <- function(bmat, nPC = 30, n = 1000, outlier = 0.999) {
  message("Running SnapATAC Diffusion Maps...")
  if(any(Matrix::rowSums(bmat) == 0)) {
    stop("Some cells have no reads, please remove them firstly.")
  }
  message("Step1: computing Jaccard similarity matrix.")
  A <- Matrix::tcrossprod(bmat, bmat)
  rowDepth <- Matrix::rowSums(A)
  J <- as.matrix(A / (2 * replicate(ncol(A), rowDepth) - A ))

  message("Step2: fitting regression model to reduce cell depth effect.")
  idx <- sampleBaseOnDepth(bmat = bmat, n = 1000)
  subNormOVE <- getNormOVE(avgDepths = Matrix::rowMeans(J[idx, idx]))
  fitData <- data.frame(x = subNormOVE[upper.tri(subNormOVE)], y = jmat[upper.tri(jmat)])
  model <- stats::lm(formula = y ~ x + I(x^2), data = fitData)
  betas <- as.numeric(model$coefficients)

  message("Step3: normalize Jaccard similarity matrix by reducing the fitted random depth.")
  rmean <- Matrix::rowMeans(J)
  allNormOVE <- getNormOVE(avgDepths = rmean)
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
	eig_transitions <- eig_decomp(transitions, n_eigs = nPC)

  message(paste("Return the diffusion map result: dim from 2 to", nPC + 1))
  dmat <- eig_transitions$vectors[, 2:(nPC+1)]
  sdev <- eig_transitions$value[2:(nPC + 1)]
  return(invisible(list(dmat = dmat, sdev = sdev)))
}
