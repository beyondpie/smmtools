library(optparse)

option_list <- list(
  make_option(c("--gmatFile"), type = "character"),
  make_option(c("--pathToPython"), type = "character"),
  make_option(c("--outScrubletFile"), type = "character", help = "seperated by comma"),
  make_option(c("--priorRate"), type = "double", default = 0.08),
  make_option(c("--minCount"), type = "double", default = 3),
  make_option(c("--minCells"), type = "integer", default = 5),
  make_option(c("--minConvPctl"), type = "double", default = 85),
  make_option(c("--nPC"), type = "integer", default = 30)
)

args <- parse_args(OptionParser(option_list = option_list))

message(paste0("Load gmat", args$gmatFile))
gmat <- smmtools::loadGeneMatrix(geneMatrixH5File = args$gmatFile)
barcode <- colnames(gmat)
message("Run Scrublet with GMM fitting ...")
mat <- t(as.matrix(gmat))
out <- smmtools::SnapATAC_runScrublet(
  mat = mat, path_to_python = args$pathToPython,
  expected_doublet_rate = args$priorRate,
  min_counts = args$minCount,
  min_cells = args$minCells,
  min_conv_pctl = args$minConvPctl,
  n_pc = args$nPC)

message("Save Scrublet results.")
invisible(write.table(
  x = data.frame(barcode = barcode,
                 scrubletScore = out[[4]],
                 thresGMM = rep(out[[1]], length(out[[4]]) ),
                 thresScrublet = rep(out[[2]], length(out[[4]]))
                 ),
  file = args$outScrubletFile,
  quote = FALSE,
  row.names = FALSE, col.names = TRUE,
  sep = ","
))
