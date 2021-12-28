#' Binarize Bmat according SnapATAC
#'
#' - Filter features if the covariances are smaller or larger than thresholds.
#' - Too high values will be set as zeros.
#'
#' @param bmat sparse matrix, feature by cells (default) or verse
#' In SnapATAC package, bmat is cell by feature. Be sure to claim the format
#' with the parameter cell_by_feature below.
#' @param cell_by_feature, bool, if bmat is cell by feature or not,
#' default is FALSE
#' @param mu double or NULL, mean of log10 count to get z score,
#' default is NULL, which will be estimated from the data.
#' @param sigma double or NULL, standard deviation of log10 count to get z score,
#' default is NULL, which will be estimated from the data.
#' @param z_threshold double, z score threshold, features highter than it or lower than -it,
#' will be removed, default is 1.65
#' @param outlier double, count larger than quantile of 1- outlier will be treated as zero,
#' default is 1e-3.
#' @return list of two elements
#' - bmat sparse matrix, cell by feature (consistent with SnapATAC)
#' - featureIndexKept integer vector, feature index kept after filtering
#' - mu double, mean for z_score
#' - sigma double, standard deviation for z_score
#' - z_threshold double, cutoff for z_threshold
#' - outlier double, cutoff for count
#' @export
SnapATAC_BinarizeBmat <- function(bmat, cell_by_feature = FALSE,
                                  mu = NULL, sigma = NULL,
                                  z_threshold = 1.65,
                                  outlier = 1e-3) {
  if (cell_by_feature) {
    bmat <- Matrix::t(bmat)
  }
  idx <- which(Matrix::rowSums(bmat) > 0)
  cov <- log(Matrix::rowSums(bmat)[idx] + 1, 10)
  if (is.null(mu)) {
    mu <- mean(cov)
  }
  if (is.null(sigma)) {
    sigma <- stats::sd(cov)
  }
  zcov <- (cov - mu) / sigma
  low <- -z_threshold
  up <- z_threshold
  idx2 <- which(zcov >= low & zcov <= up)
  idx <- idx[idx2]
  upd_bmat <- bmat[idx, ]

  count <- upd_bmat@x
  cutoff <- max(1, quantile(count, 1 - outlier))
  upd_bmat@x[upd_bmat@x > cutoff] <- 0
  upd_bmat@x[upd_bmat@x > 0] <- 1
  ## consistent with SnapATAC: cell by feature
  final_bmat <- Matrix::t(Matrix::drop0(upd_bmat))
  invisible(list(
    bmat = final_bmat,
    featureIndexKept = idx,
    mu = mu,
    sigma = sigma,
    z_threshold = z_threshold,
    outlier = outlier
  ))
}
