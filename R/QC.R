#' Remove chr in GRSeqnames
#' 
#' Ref: ArchR .convertGRSeqnames
#' @param gr GenomeRanges
#' @return gr2 GenomeRanges simplified
#' @export
simplifyGRSeqnames <- function(gr) {
  gr2 <- GenomicRanges::GRanges(seqnames = gsub("chr", "", GenomeInfoDb::seqnames(gr)),
                                ranges = IRanges::ranges(gr), strand = GenomicRanges::strand(gr))
  S4Vectors::mcols(gr2) <- S4Vectors::mcols(gr)
  return(gr2)
}


getFragmentInfo <- function(fragmentFile) {
  data.table::fread(file = fragmentFile, )
}
