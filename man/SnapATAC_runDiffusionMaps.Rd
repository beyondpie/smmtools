% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DimReduction.R
\name{SnapATAC_runDiffusionMaps}
\alias{SnapATAC_runDiffusionMaps}
\title{Core of DiffusionMaps for landmark cells.}
\usage{
SnapATAC_runDiffusionMaps(bmat, nPC = 30, n = 1000, outlier = 0.999)
}
\arguments{
\item{bmat}{sparse Matrix, cell by feature, after binarization}

\item{nPC}{integer, default is 30}

\item{n}{integer, used for fitting regression of random depth effect, default is 1000.}

\item{outlier}{double, used for choosing threshold of normalized Jaccard similarity matrix,
default is 0.999}
}
\value{
list of five elements
\itemize{
\item dmat dense matrix, cell by features
\item sdev vector, standard deviation
\item cutoff numeric, cutoff of outiliers for the normalized Jaccard similarity matrix
\item betas numeric vector of three elements, regression coefficients for reducing the random
depth effect in Jaccard similarity matrix.
\item diagFactor Diagonal Matrix, inverse of square root of
the rowSums of normalized Jaccard similarity matrix.
}
}
\description{
Core of DiffusionMaps for landmark cells.
}
