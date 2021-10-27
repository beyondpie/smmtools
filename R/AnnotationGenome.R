#' @useDynLib smmtools
#' @importFrom Rcpp sourceCpp
NULL

## use  package-local variables
## Ref: https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r/12605694
smmenv <- new.env(parent = emptyenv())
assign(
  x = "genomeInfo",
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
    ArchRgeneAnnotation = c(
      "ArchRgeneAnnoHg19",
      "ArchRgeneAnnoHg38",
      "ArchRgeneAnnoMm9",
      "ArchRgeneAnnoMm10"
    ),
    ArchRgenomeAnnotation = c(
      "ArchRgenomeAnnoHg19",
      "ArchRgenomeAnnoHg38",
      "ArchRgenomeAnnoMm9",
      "ArchRgenomeAnnoMm10"
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
  if (lowergenome %in% smmenv$genomeInfo$genome) {
    bl <- tryCatch({
      blacklist <- rtracklayer::import.bed(smmenv$genomeInfo[lowergenome, "blacklist"])
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
#' @param TxDb TxDbobjet
#' @return TSS GenomicRanges
#' @export
#' @examples
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' TSS <- getTSSFromTxDb(TxDb)
getTSSFromTxDb <- function(TxDb) {
  return(BiocGenerics::unique(IRanges::resize(x = GenomicFeatures::transcripts(TxDb), width = 1, fix = "start")))
}

#' get gene or genome Annotation from ArchR data.
#'
#' Ref: ArchR getArchRGenome
#'
#' @param tag string, gene or genome
#' @param genome string genomenm
#' @return annotation SimpleGenomicRangesList
#' @export
#' @examples
#' getAnnotFromArchRData(tag = "gene", genome = "mm10")
getAnnotFromArchRData <- function(tag, genome) {
  lowertag <- tolower(tag)
  if (!(lowertag %in% c("gene", "genome"))) {
    stop(paste0(lowertag, " not in ", paste0(c("gene", "genome"), collapse = ",")))
  }
  lowergenome <- tolower(genome)
  if(!(lowergenome %in% smmenv$genomeInfo$genome)){
    stop(paste0(genome, " not in ", paste0(smmenv$genomeInfo$genome, collapse = ",")))
  }
  ArchRtag <- "ArchRgeneAnnotation"
  if(tag == "genome") {
    ArchRtag <- "ArchRgenomeAnnotation"
  }
  AnnoName <- smmenv$genomeInfo[genome, ArchRtag]
  eval(parse(text = paste0("data(", AnnoName, ', package = "smmtools")')))
  return(eval(parse(text = gsub("ArchR", "", AnnoName))))
}


#' subset GeneAnnotation by GenomeAnnotation
#'
#' Ref: ArchR .validGeneAnnoByGenomeAnno
#' @param geneAnnotation SimpleGenomicRangesList for genes
#' @param genomeAnnotation SimpleGenomicRangesList for genomic
#' @return geneAnnotation SimpleGenomicRangesList limited by genomicAnnotation
#' @export
subsetGeneAnnoByGenomeAnno <- function(geneAnnotation, genomeAnnotation) {
  genomeChrs <- unique(GenomeInfoDb::seqnames(genomeAnnotation$chromSize))
  geneChrs <- unique(GenomeInfoDb::seqnames(geneAnnotation$genes))
  if(!all(geneChrs %in% genomeChrs)) {
    geneNotIn <- geneChrs[!(geneChrs %in% genomeChrs)]
    message("Found Gene Seqnames not in GenomeAnnotation chromSizes, Removing: ", paste0(geneNotIn, collapse = ","))
    geneAnnotation$genes <- subsetSeqnamesGR(gr = geneAnnotation$genes, names = genomeChrs)
  }
  exonChrs <- unique(GenomeInfoDb::seqnames(geneAnnotation$genes))
  if(!all(exonChrs %in% genomeChrs)) {
    exonNotIn <- exonChrs[!(exonChrs %in% genomeChrs)]
    message("Found Exon Seqnames not in GenomeAnnotation chromSizes, Removing: ", paste0(exonNotIn, collapse = ","))
    geneAnnotation$exons <- subsetSeqnamesGR(gr = geneAnnotation$exons, names = genomeChrs)
  }
  TSSChrs <- unique(GenomeInfoDb::seqnames(geneAnnotation$genes))
  if(!all(TSSChrs %in% genomeChrs)) {
    TSSNotIn <- TSSChrs[!(TSSChrs %in% genomeChrs)]
    message("Found TSS Seqnames not in GenomeAnnotation chromSizes, Removing: ", paste0(TSSNotIn, collapse = ","))
    geneAnnotation$TSS <- subsetSeqnamesGR(gr = geneAnnotation$TSS, names = genomeChrs)
  }
  return(geneAnnotation)
}

