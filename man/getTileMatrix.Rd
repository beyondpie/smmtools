% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TileMatrix.R
\name{getTileMatrix}
\alias{getTileMatrix}
\title{Generate matrix of couts for each tile per cell}
\usage{
getTileMatrix(
  rawH5File,
  outfile,
  genome = "mm10",
  sampleName = NULL,
  barcodes = NULL,
  tileSize = 5000,
  excludeChr = c("chrM")
)
}
\arguments{
\item{rawH5File}{string}

\item{outfile}{string}

\item{genome}{string, default is "mm10"}

\item{sampleName}{string, default is NULL}

\item{barcodes}{vector of string, no need to append sampleName, will do that in this code block}

\item{tileSize}{integer, default is 5000,}

\item{excludeChr}{vector of string, default is c("chrM")}
}
\value{
string, where the matrix are saved
}
\description{
Ref: ArchR
}
