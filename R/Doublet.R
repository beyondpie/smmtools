#' Doublet removement using Scrublet
#'
#' @param mat sparse matrix, cell by feature
#' @importFrom reticulate r_to_py py_to_r import
#' @return list, five fields
#'   1. threshold based on GMM for simulation scores [5]
#'   2. threshold estimated from Scrublets
#'   3. probs of real data based on GMM
#'   4. doublet scores of real data based on Scrublets
#'   5. doublet scores of simulation data based on Scrublets
#' @export
SnapATAC_runScrublet <- function(mat, path_to_python, expected_doublet_rate = 0.08,
                        min_counts = 3,
                        min_cells = 5L,
                        min_conv_pctl = 85,
                        n_pc = 30L) {
  message("Epoch: loading python environments and packages ... \n")
  # load library and python env
  reticulate::use_python(path_to_python)
  message("Use the Python located in:", path_to_python, "\n")
  setSessionTimeLimit(cpu = Inf, elapsed = Inf)

  message("Epoch: identify potential doublets ... \n")
  scr <- import(module = "smmuty", convert = FALSE)
  out <- scr$detectDoublet(
    ## scrublet needs cell-by-feature matrix
    counts_matrix = r_to_py(mat),
    expected_doublet_rate = expected_doublet_rate,
    min_counts = min_counts,
    min_cells = as.integer(min_cells),
    min_gene_variability_pctl = min_conv_pctl,
    n_prin_comps = as.integer(min(n_pc, nrow(mat), ncol(mat)))
  )
  out <- py_to_r(out)
  names(out) <- c("thresGMM", "thresScrublet", "probs", "scrubletScore", "simuScrubletScore")
  ## summary
  scores <- out$scrubletScore
  ncell <- length(out$scrubletScore)
  tgmm <- out$thresGMM
  tscrub <- out$thresScrublet
  message("Epoch: summary doublets detection results ... \n")
  message(paste("Threshold from GMM:", round(tgmm, 3)))
  message(paste("Threshold from Raw Scrublets", round(tscrub, 3)))
  message(paste("Number of cells in total:"), ncell)
  message(paste("Doublet rate based on GMM:", round(sum(scores > tgmm) / ncell, 4)))
  message(paste("Doublet rate based on RawScrublet:",
                round(sum(scores > tscrub) / ncell, 4)))
  return(out)
}

## getScrubletThresholdByEM <- function()
