library(optparse)
library(smmtools)
option_list <- list(
  make_option(c("--tabixFile"), type = "character",
              help = "tab or tab.gz file from 10x genomics"),
  make_option(c("--rawH5File"), type = "character",
              help = "name of output file"),
  make_option(c("--genome"), type = "character", default = "mm10"),
  make_option(c("--sampleName"), type = "character"),
  make_option(c("--nChunk"), type = "integer", default = 3)
)

args <- parse_args(OptionParser(option_list = option_list))

dir.create(path = dirname(args$rawH5File), showWarnings = FALSE, recursive = TRUE)

annotGenome <- getAnnotFromArchRData(tag = "genome", genome = args$genome)
annotGene <- getAnnotFromArchRData(tag = "gene", genome = args$genome)
annotGene <- subsetGeneAnnoByGenomeAnno(
  geneAnnotation = annotGene,
  genomeAnnotation = annotGenome
)

tileChromSizes <- tileChrom(chromSizes = annotGenome$chromSizes, nChunk = args$nChunk)

invisible(tabixToH5SingleThread(
  tabixFile = args$tabixFile,
  tileChromSizes = tileChromSizes,
  sampleName = args$sampleName,
  outH5File = args$rawH5File,
  barcode = NULL
))


