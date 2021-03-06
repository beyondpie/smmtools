% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TileMatrix.R
\name{getBmatFromSnap}
\alias{getBmatFromSnap}
\title{Get Bmat from Snap object.}
\usage{
getBmatFromSnap(
  snapFile,
  binSize = 5000,
  chrs = c(paste0("chr", 1:19), "chrX", "chrY"),
  ibarcode = NULL,
  blacklistBedFile = NULL,
  outfile = NULL,
  compressLevel = 9
)
}
\arguments{
\item{chrs}{vector of string, the chroms we want. Since there might exist chroms like
chrM, chrY_JH584303_random, chrUn_GL456390, so here we use explicitly chrs we want to
filter chroms. The default is specific for mm10.}

\item{ibarcode}{vector of string, the barcodes we want to, we assume all of them are in
this snapFile. So use this carefully.}

\item{blacklistBedFile}{string, bed file of the blacklist, smmtools save one under data dir for mm10.}
}
\value{
sparseMatrix, cell by feature, with both row and colname and ordered by ibarcode if provided
}
\description{
Note: this function is validated with SnapATAC::addBmatFromSnap
To include all the chromes in snap file,
chrs <- unique(unlist(lapply(strsplit(snap@feature$name, ":"), function(x){x\link{1}})))
where snap is a snap object created by SnapATAC.
}
