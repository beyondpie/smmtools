#' Doublet removement using Scrublet
#'
#' Need smmutil
#' 
#' @param bmat sparse matrix, feature by cell, transpose will be used when as input in scrublet 
#' @importFrom reticulate r_to_py py_ro_r import
#' @importFrom base nrow ncol
#' @return list, five fields
#'   1. threshold based on GMM for simulation scores [5]
#'   2. threshold estimated from Scrublets
#'   3. probs of real data based on GMM
#'   4. doublet scores of real data based on Scrublets
#'   5. doublet scores of simulation data based on Scrublets
#' @export
runScrublet <- function(bmat, path_to_python, expected_doublet_rate = 0.06,
                        min_counts = 3,
                        min_cells = 5L,
                        min_conv_pctl = 85,
                        n_pc = 30L) {
  cat("Epoch: loading python environments and packages ... \n", file = stderr())
  # load library and python env
  reticulate::use_python(path.to.python)
  cat("use the Python located in:", path.to.python, "\n")
  setSessionTimeLimit(cpu = Inf, elapsed = Inf)

  cat("Epoch: identify potential doublets ... \n")
  scr <- import(module = "smmuty", covert = FALSE)
  out <- src$detectDoublet(
    ## transpose since scrublet needs cell-by-feature matrix
    counts_matrix = r_to_py(t(bmat)),
    expected_doublet_rate = expected_doublet_rate,
    min_counts = min_counts,
    min_cells = as.integer(min_cells),
    min_gene_variability_pctl = min_conv_pctl,
    n_prin_comps = min(n_pc, nrow(bmat), ncol(bmat)))
  out <- py_to_r(out)
  ## summary
  cat("Epoch: summary doublets detection results ... \n")
  message(paste("Threshold from GMM:", round(out[[1]], 3)))
  message(paste("Threshold from Raw Scrublets", round(out[[2]], 3)))
  message(paste("Number of cells used:"), length(out[[3]]))
  message(paste("Doublet rate based on GMM:", round(length(which(out[[3]]>0.5)) / length(out[[3]]), 4)))
  return(out)
}



