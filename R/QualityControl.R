#' Transform original tabix File to H5File, then get barcodes-related fragments stats.
#'
#' - Ref: ArchR
#' - Method is run on single core, but the speed should be OK.
#' - No filtering barcodes here unless we provide barcode array ourselves.
#' - Side Effect:
#'   - output rawH5File
#'   - output sumFragment File (descreasingly ordered by nUniqFrag)
#'
#' @param tabixFile string, 10x cell ranger result tab/tab.gz file name
#' @param outfilenm string, output file name, if NULL, a default one will be genereated.
#' @param coverH5File bool, whether to cover the rawH5File if exists
#' @param genome string, namely mm9, mm10, hg19 or hg38
#' @param sampleName string, name for the tabix file
#' @param barcodes vector of string, limite fragment on specific barcodes
#' @param nChunk integer, partition each chrom into nChunk pieces
#' @param tsseWindow integer
#' @param tsseNorm interger
#' @param tsseFlank integer
#' @param tsseMinNorm float
#' @param tsseMaxFragSize integer
#' @return data.frame with cols: barcode, nUniqFrag, TSSE, TSSReads
#' @export
sumFragmentSingleThread <- function(rawH5File,
                                    sumFragFile,
                                    genome,
                                    sampleName = NULL,
                                    barcodes = NULL,
                                    nChunk = 3,
                                    tsseWindow = 101,
                                    tsseNorm = 100,
                                    tsseFlank = 2000,
                                    tsseMinNorm = 0.2,
                                    tsseMaxFragSize = NULL) {
  annotGenome <- getAnnotFromArchRData(tag = "genome", genome = genome)
  annotGene <- getAnnotFromArchRData(tag = "gene", genome = genome)
  annotGene <- subsetGeneAnnoByGenomeAnno(
    geneAnnotation = annotGene,
    genomeAnnotation = annotGenome
  )
  tileChromSizes <- tileChrom(chromSizes = annotGenome$chromSizes, nChunk = nChunk)
  chunkName <- S4Vectors::mcols(x = tileChromSizes)$chunkName
  ## get nfrag per barcode
  nFragPerBarcode <- getNfragmentPerBarcode(
    chrRegions = chunkName, rawH5File = rawH5File,
    sampleName = sampleName
  )
  TSSEnrich <- fastGetTSSEnrichmentSingleThread(
    TSS = annotGene$TSS, barcodes = nFragPerBarcode$barcode,
    rawH5File = rawH5File,
    window = tsseWindow, norm = tsseNorm, flank = tsseFlank,
    minNorm = tsseMinNorm, maxFragSize = tsseMaxFragSize,
    sampleName = sampleName
  )
  sumFrag <- as.data.frame(merge(x = nFragPerBarcode, y = TSSEnrich, by = "barcode", all = TRUE))
  sumFrag <- sumFrag[order(sumFrag$nUniqFrag, decreasing = TRUE), ]
  invisible(write.table(
    x = sumFrag,
    file = sumFragFile,
    quote = FALSE,
    sep = ",", row.names = FALSE, col.names = TRUE
  ))
  return(sumFrag)
}

#' get numebr of fragments per barcode
#' Ref: ArchR
#'
#' @param chrRegions array of characters, chrRegion names in rawH5File
#' @param rawH5File string
#' @return tibble shape = [length(barcodes), 2], colnames = c("barcode", "count")
#' @importFrom dplyr %>% group_by summarise arrange desc
#' @export
getNfragmentPerBarcode <- function(chrRegions, rawH5File, sampleName = NULL) {
  ## get nfrag per barcode
  dtList <- lapply(seq_along(chrRegions), function(x) {
    chrRegion <- chrRegions[x]
    ## fragmentRanges <- paste0("Fragments/", chrRegion, "/Ranges")
    barcodeLength <- rhdf5::h5read(
      file = rawH5File,
      name = paste0("Fragments/", chrRegion, "/BarcodeLength")
    )
    barcodeValue <- rhdf5::h5read(
      file = rawH5File,
      name = paste0("Fragments/", chrRegion, "/BarcodeValue")
    )
    dt <- NULL
    if ((length(barcodeValue) > 0) & (length(barcodeLength) > 0)) {
      dt <- data.table::data.table(
        ## values = paste(sampleName, barcodeValue, sep = "#"),
        values = barcodeValue,
        lengths = barcodeLength
      )
    }
    invisible(dt)
  })
  names(dtList) <- chrRegions
  dt <- Reduce(f = "rbind", dtList)
  colnames(dt) <- c("barcode", "nUniqFrag")
  dt <- dt %>%
    group_by(barcode) %>%
    summarise(nUniqFrag = sum(nUniqFrag)) %>%
    arrange(desc(nUniqFrag))
  return(dt)
}

#' Get TSS Enrichment scores
#'
#' Ref: ArchR
#' @return List two element: TSSE, TSSReads
#' @importFrom S4Vectors mcols mcols<- split DataFrame queryHits subjectHits
#' @importFrom BiocGenerics strand<- match start end pmax
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges width ranges resize
#' @importFrom GenomeInfoDb seqnames
#' @importFrom rhdf5 h5ls
#' @return dataframe with cols, barcode,TSSE,TSSRead
#' @export
fastGetTSSEnrichmentSingleThread <- function(TSS, barcodes,
                                             rawH5File,
                                             window = 101, norm = 100,
                                             flank = 2000, minNorm = 0.2, maxFragSize = NULL,
                                             sampleName = NULL) {
  tstart <- Sys.time()
  message(paste("Get TSS Enrichment Scores starts at", tstart))
  TSS <- resize(x = TSS, width = 1, fix = "start")
  strand(x = TSS) <- "*"
  TSS <- unique(TSS)
  tssWindow <- resize(x = TSS, width = window, fix = "center")
  tssWindow$type <- "window"
  tssFlank <- c(
    GRanges(
      seqnames = seqnames(TSS),
      ranges = IRanges(
        start = end(TSS) + flank - norm + 1,
        end = end(TSS) + flank
      )
    ),
    GRanges(
      seqnames = seqnames(TSS),
      ranges = IRanges(
        start = start(TSS) - flank,
        end = start(TSS) - flank + norm - 1
      )
    )
  )
  ## get Available chrs
  groups <- h5ls(file = rawH5File)
  groups <- groups[groups$group == "/Fragments" & groups$otype == "H5I_GROUP", "name"]

  chrs <- unique(gsub(
    pattern = "#chunk\\d+",
    replacement = "",
    x = groups
  ))

  tssFlank$type <- "flank"
  tssFeatures <- c(tssWindow, tssFlank)
  tssFeatureList <- split(x = tssFeatures, f = seqnames(tssFeatures))
  tssFeatureList <- tssFeatureList[chrs]

  countDF <- lapply(X = seq_along(tssFeatureList), FUN = function(i) {
    nWindow <- rep(0, length(barcodes))
    names(nWindow) <- barcodes
    nFlank <- rep(0, length(barcodes))
    names(nFlank) <- barcodes
    fragments <- getFragsOfAChrFromRawH5File(
      rawH5File = rawH5File, chr = names(tssFeatureList)[i],
      sampleName = sampleName, barcodes = NULL
    )
    if (length(fragments) == 0) {
      names(nWindow) <- NULL
      names(nFlank) <- NULL
      return(DataFrame(nWindow = nWindow, nFlank = nFlank))
    }
    if (length(fragments) > 0) {
      if (!is.null(maxFragSize)) {
        fragments <- fragments[width(fragments) <= maxFragSize]
      }
      mcols(fragments)$RG@values <- match(
        mcols(fragments)$RG@values,
        barcodes
      )
      feature <- tssFeatureList[[i]]
      mcols(feature)$typeIdx <- match(mcols(feature)$type, c("window", "flank"))
      ## count each insertion
      for (y in seq_len(2)) {
        if (y == 1) {
          temp <- IRanges(start(fragments), width = 1)
        } else {
          temp <- IRanges(end(fragments), width = 1)
        }
        o <- findOverlaps(ranges(feature), temp)
        ## CALL CPP
        mat <- tabulate2dCpp(
          x = as.vector(mcols(fragments)$RG[subjectHits(o)]),
          xmin = 1,
          xmax = length(barcodes),
          y = mcols(feature)$typeIdx[queryHits(o)],
          ymin = 1,
          ymax = 2
        )

        nWindow <- nWindow + mat[1, ]
        nFlank <- nFlank + mat[2, ]
      }
    }
    names(nWindow) <- NULL
    names(nFlank) <- NULL
    return(DataFrame(nWindow = nWindow, nFlank = nFlank))
  })
  cumDF <- countDF[[1]]
  for (i in 2:length(countDF)) {
    cumDF$nWindow <- cumDF$nWindow + countDF[[i]]$nWindow
    cumDF$nFlank <- cumDF$nFlank + countDF[[i]]$nFlank
  }
  cWn <- cumDF[, 1] / window
  cFn <- cumDF[, 2] / norm
  ## pmax: like Relu(vec - scalar) + scalar
  tssScores <- round(2 * cWn / pmax(cFn, minNorm), 3)
  names(tssScores) <- barcodes
  tend <- Sys.time()
  message(paste("Get TSS Enrichment Scores ends at", tend))
  return(data.frame(barcode = barcodes, TSSE = tssScores, TSSRead = cumDF[, 1]))
}

#' @importFrom S4Vectors mcols mcols<- split DataFrame queryHits subjectHits
#' @importFrom BiocGenerics strand<- match start end pmax
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges width ranges resize
#' @importFrom GenomeInfoDb seqnames
#' @importFrom rhdf5 h5ls
#' @importFrom dplyr %>% group_by summarise arrange desc
#' @export
getFragSize <- function(rawH5File, pal = NULL, barcodes = NULL, sampleName = NULL, maxSize = 750) {
  groups <- h5ls(file = rawH5File)
  groups <- groups[groups$group == "/Fragments" & groups$otype == "H5I_GROUP", "name"]
  chrs <- unique(gsub(
    pattern = "#chunk\\d+",
    replacement = "",
    x = groups
  ))

  for (i in seq_along(chrs)) {
    if (i == 1) {
      fsi <- getFragsOfAChrFromRawH5File(
        rawH5File = rawH5File, chr = chrs[i],
        sampleName = sampleName, barcodes = barcodes
      ) %>%
        width() %>%
        tabulate(nbins = maxSize)
    } else {
      fsi <- fsi + getFragsOfAChrFromRawH5File(
        rawH5File = rawH5File, chr = chrs[i],
        sampleName = sampleName, barcodes = barcodes
      ) %>%
        width() %>%
        tabulate(nbins = maxSize)
    }
  }
  df <- data.frame(fragmentSize = seq_along(fsi), fragmentPercent = round(100 * fsi / sum(fsi), 4))
  if (is.null(pal)) {
    pal <- paletteDiscrete(values = c(1))
  }
  p <- ggplot(df, aes(fragmentSize, fragmentPercent)) +
    geom_line(size = 1) +
    theme_ArchR() +
    xlab("ATAC-seq Fragment Size (bp)") +
    ylab("Percentage of Fragments") +
    scale_color_manual(values = pal) +
    scale_y_continuous(limits = c(0, max(df$fragmentPercent) * 1.05), expand = c(0, 0)) +
    scale_x_continuous(limits = c(min(df$fragmentSize), max(df$fragmentSize)), expand = c(0, 0))
  return(list(df = df, p = p))
}
