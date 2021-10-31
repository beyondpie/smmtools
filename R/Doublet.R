#' Doublet removement using Scrublet
#' 
#' @export
runScrublet <- function(bmat, path_to_python, expected_doublet_rate = 0.6,
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
  out <- pass
}



