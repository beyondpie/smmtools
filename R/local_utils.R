## use  package-local variables
## Ref: https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r/12605694
smmenv <- new.env(parent = emptyenv())
assign(
  x = geomeInfo,
  value = data.frame(
    genome = c("hg19", "hg38", "mm9", "mm10"),
    BSGenome = c(
      "BSgenome.Hsapiens.UCSC.hg19",
      "BSgenome.Hsapiens.UCSC.hg38",
      "BSgenome.Mmusculus.UCSC.mm10",
      "BSgenome.Mmusculus.UCSC.mm9"
    ),
    blasklist = c(
      "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz",
      "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
      "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
      "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm9-blacklist.bed.gz"
      ),
    TxDb = c(
      "TxDb.Hsapiens.UCSC.hg19.knownGene",
      "TxDb.Hsapiens.UCSC.hg38.knownGene",
      "TxDb.Mmusculus.UCSC.mm9.knownGene",
      "TxDb.Mmusculus.UCSC.mm10.knownGene"
    ),
    OrgDb = c(
      "org.Hs.eg.db",
      "org.Hs.eg.db",
      "org.Mm.eg.db",
      "org.Mm.eg.db"
    ),
    row.names = c("hg19", "hg38", "mm9", "mm10")
  ),
  envir = smmenv
)


#' Read file by columns
#'
#' Ref: https://www.gastonsanchez.com/visually-enforced/how-to/2012/06/23/Read-file-by-columns/
#'
#' @param filenm str
#' @param cols vector[Integer]
#' @param sep str
#' @param ... others for read.csv
#'
#' @return data.frame from read.csv
#'
#' @export
read_columns <- function(filenm, cols = c(1), sep = ",", ...) {
  args <- c(
    "-f", paste(cols, collapse = ","),
    "-d", paste0("'", sep, "'"), filenm
  )
  command <- "cut"
  command <- paste(c(shQuote(command), args), collapse = " ")
  tmp <- system(command, intern = TRUE)
  invisible(utils::read.table(pipe(command), sep = sep, ...))
}

#' Load sparse matrix from mat file
#'
#' Each row of matfile: barcode[,genome_bin_index:counts]
#'
#' @param matfile str
#' @param ... for Matrix::sparseMatrix usage
#'
#' @return sparseMatrix
#'
#' @export
load_sparse_matrix <- function(matfile, ... ) {
  lines <- data.table::fread(input = matfile, sep = "\n", header = FALSE)
  row_col_value <- lapply(seq(nrow(lines)), FUN = function(j) {
    l <- as.character(lines[j,1])
    a <- unlist(strsplit(x = l, split = ","))
    barcode <- a[1]
    col_value <- t(vapply(a[-1], FUN = function(i) {
      invisible(as.numeric(unlist(strsplit(x = i, split = ":"))))
    }, FUN.VALUE = c(0, 0)))
    r <- cbind(j, col_value)
  })
  r <- as.data.frame(do.call(what = rbind, args = row_col_value))
  invisible(Matrix::sparseMatrix(i = r[, 1], j = r[, 2], x = r[, 3],...))
}

## TODO: support mtx in R instead of the sparse matrix we defined before.
## This is more general, can be directly used in Python.


#' Install genome.
#'
#' Install the genome from BiocManager if the corresponding BSgenome does not exists.
#' Ref: ArchR addArchRGenome
#' 
#' @param genomes str
#' @return BSgenomeName str
#'
#' @export
installGenome <- function(genome) {
  lowergenome <- tolower(genome)
  if(lowergenome %ni% smmenv$genomeInfo$genome){
    stop(paste0(genome, " not in ", paste0(smmenv$genomeInfo$genome, collapse = ",")))
  }
  BSgenomeNm <- smmenv$genoemInfo[lowergenome, "BSGenome"]
  if(!requireNamespace(BSgenomeNm, quitely = TRUE)) {
    BiocManager::install(BSgenomeNM)
  }
  return(BSgenomeNM)
}

#' get black list
#'
#' Return GRanges. Empty GRanges will be return if no name matches the genome or network issusue.
#'
#' Ref: ArchR .getBlacklist
#'
#' @param genome str
#' @return GRanges
#'
#' @export

getBlacklist <- function(genome) {
  lowergenome <- tolower(genome)
  if (lowergenoem %in% smmenv$genomeInfo$genome) {
    bl <- tryCatch({
      blacklist <- rtracklayer::import.bed(smmenv$genoemInfo[lowergenome, "blacklist"])
    }, error = function(x) {
      message("Find the blacklist, but failed at download. Return Empty GRanges.")
      GenomicRanges::GRanges()
    })
    return(bl)
  }
  message(paste0(genome, " not in ", paste0(smmenv$genomeInfo$genome, collapse = ",")))
  return(GenomicRanges::GRanges())
}


#' get TxDb object, which records gene/transcript coordinates


