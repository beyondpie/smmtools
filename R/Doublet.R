#' Doublet removement using Scrublet
#'
#' GMM is used by Kai Zhang in his Cell ATAC Atlas paper.
#' NormalMixEM is used by Yang Li in his CEMBA 1.0
#'
#' @param mat sparse matrix, cell by feature
#' @param path_to_python string
#' @param expected_doublet_rate double, default is 0.08
#' @param min_counts numeric, default is 3
#' @param min_cells integer, default is 5
#' @param min_conv_pctl numeric, default is 85
#' @param n_pc integer, default is 30
#' @importFrom reticulate r_to_py py_to_r import
#' @return list, five fields
#'   - thresGMM: threshold based on GMM for simulation scores
#'   - thresScrublet: threshold estimated from Scrublets
#'   - probs: probs of real data based on GMM
#'   - scrubletScore: doublet scores of real data based on Scrublets
#'   - simuScrubletScore: doublet scores of simulation data based on Scrublets
#'   - thresMixEM: threshold based normalmixEM for simulation scores.
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
  scores <- out$scrubletScore
  ncell <- length(out$scrubletScore)
  tgmm <- out$thresGMM
  tscrub <- out$thresScrublet
  tMixEM <- getScrubletThresholdByMixEM(scrubletSimScores = out$simuScrubletScore,
                                        defaultCutoff = tscrub)
  out$thresMixEM <- tMixEM
  ## summary
  message("Epoch: summary doublets detection results ... \n")
  message(paste("Threshold from GMM:", round(tgmm, 3)))
  message(paste("Threshold from Raw Scrublets", round(tscrub, 3)))
  message(paste("Number of cells in total:"), ncell)
  message(paste("Doublet rate based on GMM:", round(sum(scores > tgmm) / ncell, 4)))
  message(paste("Doublet rate based on RawScrublet:",round(sum(scores > tscrub) / ncell, 4)))
  message(paste("Doublet rate based on RawScrublet:", round(sum(scores > tMixEM) / ncell, 4)))
  ## get result
  return(out)
}

#' Another way to estimate the threshold for Scrublet
#' Method used by Yang Li in CEMBA 1.0.
#' @param scrubletSimScores numeric vector
#' @param defaultCutoff double, if mixEM failed, use this as threschold.
#' @param prob double, prob of a scrublet score treated as doublet, default is 0.5
#' @param lower double, used for uniroot
#' @return double 
#' @importFrom mixtools normalmixEM
#' @export
getScrubletThresholdByMixEM <- function(scrubletSimScores, defaultCutoff,
                                        prob = 0.5, lower = 0.2) {
  x <- scrubletSimScores[scrubletSimScores > 0]
  x <- x[is.finite(x)]
  model <- normalmixEM(x, k = 2)
  i <- which.min(model$mu)
  f <- function(x) {
    prob - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /
               (model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) +
                  model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
  }
  ## How about using try-catch?
  t <- try(stats::uniroot(f = f, prob = 0.5, lower = 0.2, upper = 1), silent = TRUE)
  if("try-error" %in% class(t)) {
    cutoff <- defaultCutoff
  } else {
    cutoff <- stats::uniroot(f = f, prob = 0.5, lower = 0.2, upper = 1)$root
  }
  return(cutoff)
}


