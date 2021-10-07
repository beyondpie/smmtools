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

#' Install genome/TxDb/OrgDb
#'
#' Install the genome from BiocManager if the corresponding BSgenome does not exists.
#' Ref: ArchR addArchRGenome
#' 
#' @param genome str
#' @param dbnm str
#' @return name str
#' @export
installGenomeRelatedDatabase <- function(genome, dbnm) {
  lowergenome <- tolower(genome)
  lowerdbnm <- tolower(dbnm)
  if(!(lowergenome %in% smmenv$genomeInfo$genome)){
    stop(paste0(genome, " not in ", paste0(smmenv$genomeInfo$genome, collapse = ",")))
  }

  if(!(lowerdbnm %in% colnames(smmenv$genoemInfo))) {
    stop(paste0(lowerdbnm, " not in ", paste0(colnames(smmenv$genoemInfo), collapse = ",")))
  }
  fulldb <- smmenv$genomeInfo[lowergenoem, lowerdbnm]
  if(!requireNamespace(fulldb, quitely = TRUE)) {
    BiocManager::install(fulldb, update = FALSE)
  }
  return(fulldb)
  
}

#' get black list
#'
#' Return GRanges. Empty GRanges will be return if no name matches the genome or network issusue.
#' Ref: ArchR .getBlacklist
#'
#' @param genome str
#' @return GRanges
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

#' filter ChrGR
#'
#' Ref: ArchR filterChrGR
#'
#' @param gr GRanges
#' @param remove charactervector extra seqlevels to be removed
#' @param underscore bool remove seqlevels containing an underscore
#' @param standard bool keep all the standard genomes
#' @return gr GRanges
#' @export
filterChrGR <- function(gr = NULL, remove = NULL, understcore = TRUE){
  if(is.null(gr)) {
    return(gr)
  }
  if(standard){
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
  }
  chrNames <- GenomeInfoDb::seqlevels(gr)
  chrRemove <- c()
  #first we remove all chr with an underscore
  if(underscore){
    chrRemove <- c(chrRemove, which(grepl("_", chrNames)))
  }
  #next we remove all chr specified in remove
  if(!is.null(remove)){
    chrRemove <- c(chrRemove, which(chrNames %in% remove))
  }
  
  if(length(chrRemove) > 0){
    chrKeep <- chrNames[-chrRemove]
  }else{
    chrKeep <- chrNames
  }
  #this function restores seqlevels
  GenomeInfoDb::seqlevels(gr, pruning.mode="coarse") <- chrKeep
  return(gr)
}

#' get TSS from TxDb
#'
#' Ref: ArchR createGeneAnnotation
#'
#' @param TxDb TxDbobjet  library(TxDb.Mmusculus.UCSC.mm10.knownGene); TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' @return TSS GenomicRanges
#' @export
#' @examples
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' TSS <- getTSS(TxDb)
getTSS <- function(TxDb) {
  return(BiocGenerics::unique(IRanges::resize(x = GenomicFeatures::transcripts(TxDb), width = 1, fix = "start")))
}
