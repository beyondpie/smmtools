
<!-- README.md is generated from README.Rmd. Please edit that file -->

# smmtools

<!-- badges: start -->
<!-- badges: end -->

The goal of **smmtools** is to analyze 10x multi-omic ATAC sequencing
data with SnapATAC pipelines. This package is mainly based on ArchR,
SnapATAC and SnapTools. The significant difference is that we
decentralize the data representation, and make smmtools as a pool of
different functional utilities. So that users can combine different
functions and complete the tasks. Also it’s quite easy to hack smmtools
by either adding new functions or change some functions.

# Installation

## Install python package *smmuty*

``` shell
## install the local package smmuty with utilities like scrublet, lenden algorithm support.
pip install smmuty
```

## Install R package of smmtools

``` r
# install.packages("devtools")
devtools::install_github("beyondpie/smmtools")
```
