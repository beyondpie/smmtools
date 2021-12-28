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
#' @param nFragms data.frame, cols: barcode, nFragment
#' @return list of two elements
#' - bmat sparse matrix, cell by feature (consistent with SnapATAC)
#' - featureIndexKept integer vector, feature index kept after filtering
#' @export
SnapATAC_BinarizeBmat <- function(bmat, cell_by_feature = FALSE,
                                  z_threshold = 1.65,
                                  outlier = 1e-3) {
  if(cell_by_feature) {
    bmat <- t(bmat)
  }
  idx <- which(Matrix::rowSums(bmat) > 0)
	cov = log(Matrix::rowSums(bmat)[idx] + 1, 10)
	zcov = (cov - mean(cov)) / stats::sd(cov)
  low <- - z_threshold
  up <- z_threshold
	idx2 = which(zcov >= low & zcov <= up);
	idx = idx[idx2];
  upd_bmat <- bmat[idx, ]

  count <- upd_bmat@x;
  cutoff <- max(1, quantile(count, 1 - outlier))
  upd_bmat@x[upd_bmat@x > cutoff] <- 0
  upd_bmat@x[upd_bmat@x > 0] <- 1
  ## consistent with SnapATAC
  final_bmat <- Matrix::t(Matrix::drop0(upd_bmat))
  return(invisible(list(bmat = final_bmat, featureIndexKept = idx)))
}
