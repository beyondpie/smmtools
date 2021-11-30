library(optparse)
library(smmtools)
option_list <- list(
  make_option(c("--tabixFile"),
    type = "character",
    help = "tab/tab.gz file from 10x genomics"
  ),
  make_option(c("--coverH5File"), type = "integer"),
  make_option(c("--genome"), type = "character"),
  make_option(c("--sampleName"), type = "character"),
  make_option(c("--rawH5File"), type = "character"),
  make_option(c("--sumFragFile", type = "character"))
)

args <- parse_args(OptionParser(option_list = option_list))
dir.create(path = dirname(args$sumFragFile), showWarnings = FALSE, recursive = TRUE)

invisible(sumFragmentSingleThread(
  genome = args$genome,
  sampleName = args$sampleName,
  rawH5File = args$rawH5File,
  sumFragFile = args$sumFragFile
))
