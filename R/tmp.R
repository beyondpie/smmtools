.tabixToTmp <- function(
  tabixFile = NULL, 
  sampleName = NULL,
  validBC = NULL,
  tmpFile = .tempfile(pattern = paste0("tmp-",sampleName,"-arrow"), fileext=".arrow"),
  chromSizes = NULL, 
  nChunk = 3,
  gsubExpression = NULL, 
  verbose = TRUE,
  prefix = "",
  tstart = NULL,
  threads = 1,
  logFile = NULL
  ){

  .requirePackage("Rsamtools", source = "bioc")

  #######################################################################################################
  # We will dump a chunked genome into an Hdf5 file in a memory efficient manner!
  #######################################################################################################

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  tstart2 <- Sys.time()
  
  .logDiffTime(sprintf("%s Tabix Bed To Temporary File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  o <- h5closeAll()
  o <- h5createFile(tmpFile)
  o <- h5createGroup(tmpFile, paste0("Fragments"))
  o <- h5createGroup(tmpFile, paste0("Metadata"))
  o <- h5write(obj = "Arrow", file = tmpFile, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = tmpFile, name = "ArchRVersion")
  o <- h5write(obj = "tmp", file = tmpFile, name = "Metadata/Sample")

  tileChromSizes <- unlist(GenomicRanges::tile(chromSizes, nChunk))
  mcols(tileChromSizes)$chunkName <- paste0(seqnames(tileChromSizes),"#chunk",seq_along(tileChromSizes))

  .logThis(tileChromSizes, name = "tileChromSizes", logFile = logFile)

  readTiledChrom <- .safelapply(seq_along(tileChromSizes), function(x){

    tryCatch({

      errorCheck <- 0

      if(threads == 1){
        if(x %% 10 == 0){
          .logDiffTime(sprintf("%s Reading TabixFile %s Percent", prefix, round(100*x/length(tileChromSizes)),3), t1 = tstart, 
            verbose = verbose,  logFile = logFile)
        }
      }else{
        if(x %% (2 * threads + 1) == 0){
          .logDiffTime(sprintf("%s Reading TabixFile %s Percent", prefix, round(100*x/length(tileChromSizes)),3), t1 = tstart, 
                    verbose = verbose, logFile = logFile)
        }
      }

      dt <- tryCatch({
        Rsamtools::scanTabix(tabixFile, param = tileChromSizes[x])[[1]] %>%
          textConnection %>% 
          {tryCatch(read.table(.), error = function(e) NULL)} %>% 
          {data.table(V2=.$V2 + 1, V3=.$V3, V4=.$V4)}
      }, error = function(f){
          NULL
      })

      if(is.null(dt)){

        dt <- tryCatch({
          Rsamtools::scanTabix(tabixFile, param = .convertGRSeqnames(tileChromSizes[x], method = "remove"))[[1]] %>%
            textConnection %>% 
            {tryCatch(read.table(.), error = function(e) NULL)} %>% 
            {data.table(V2=.$V2 + 1, V3=.$V3, V4=.$V4)}
        }, error = function(f){
            NULL
        })

        if(!is.null(dt)){
          .logMessage(msg = paste0(prefix, " found fragments when removed chromosome prefix : ", paste0(tileChromSizes[x])), logFile = logFile)          
        }

      }

      if(x == 1){
        .logThis(dt, name = paste0(prefix, " .tabixToTmp Fragments-Chunk-(",x," of ",length(tileChromSizes),")-", tileChromSizes[x]), logFile = logFile)
        .logThis(unique(dt$V4), name = paste0(prefix, " .tabixToTmp Barcodes-Chunk-(",x," of ",length(tileChromSizes),")-", tileChromSizes[x]), logFile = logFile)
      }

      if(is.null(dt)){
        return(list(tmpChrFile = NULL, errorCheck = errorCheck))
      }

      #Care for Break Points
      dt <- dt[dt$V2 >= start(tileChromSizes[x]),]

      if(!is.null(gsubExpression)){
        dt$V4 <- gsub(gsubExpression, "", dt$V4)
      }

      #Check for valid barcodes
      if(!is.null(validBC)){
        dt <- dt[dt$V4 %in% validBC, ]
      }

      if(all(!is.null(dt), nrow(dt) > 0)){

        errorCheck <- errorCheck + 1

        if(threads == 1){

          #Order by bc
          setkey(dt, V4)
          dt <- dt[order(V4)]
          RG <- Rle(paste0(dt$V4))

          chrTmp <- mcols(tileChromSizes)$chunkName[x]
          chrPos <- paste0("Fragments/",chrTmp,"/Ranges")
          chrRGLengths <- paste0("Fragments/",chrTmp,"/RGLengths")
          chrRGValues <- paste0("Fragments/",chrTmp,"/RGValues")
          lengthRG <- length(RG@lengths)
          o <- h5createGroup(tmpFile, paste0("Fragments/",chrTmp))
          o <- .suppressAll(h5createDataset(tmpFile, chrPos, storage.mode = "integer", dims = c(nrow(dt), 2), level = 0))
          o <- .suppressAll(h5createDataset(tmpFile, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
          o <- .suppressAll(h5createDataset(tmpFile, chrRGValues, storage.mode = "character", 
            dims = c(lengthRG, 1), level = 0, size = max(nchar(RG@values)) + 1))
          o <- h5write(obj = cbind(dt$V2,dt$V3 - dt$V2 + 1), file = tmpFile, name = chrPos)
          o <- h5write(obj = RG@lengths, file = tmpFile, name = chrRGLengths)
          o <- h5write(obj = RG@values, file = tmpFile, name = chrRGValues)

          rm(dt, RG)
          gc()

          return(list(tmpChrFile = NULL, errorCheck = errorCheck))

        }else{

          chrTmp <- mcols(tileChromSizes)$chunkName[x]

          #Temporary File
          tmpChrFile <- paste0(gsub(".arrow", "", tmpFile), ".", chrTmp, ".arrow")
          if(file.exists(tmpChrFile)){
            file.remove(tmpChrFile)
          }

          o <- h5createFile(tmpChrFile)

          #Order by bc
          setkey(dt, V4)
          dt <- dt[order(V4)]
          RG <- Rle(paste0(dt$V4))

          chrPos <- paste0(chrTmp, "._.Ranges")
          chrRGLengths <- paste0(chrTmp, "._.RGLengths")
          chrRGValues <- paste0(chrTmp, "._.RGValues")
          lengthRG <- length(RG@lengths)

          o <- .suppressAll(h5createDataset(tmpChrFile, chrPos, storage.mode = "integer", dims = c(nrow(dt), 2), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGLengths, storage.mode = "integer", dims = c(lengthRG, 1), level = 0))
          o <- .suppressAll(h5createDataset(tmpChrFile, chrRGValues, storage.mode = "character", 
            dims = c(lengthRG, 1), level = 0, size = max(nchar(RG@values)) + 1))
          
          o <- h5write(obj = cbind(dt$V2,dt$V3 - dt$V2 + 1), file = tmpChrFile, name = chrPos)
          o <- h5write(obj = RG@lengths, file = tmpChrFile, name = chrRGLengths)
          o <- h5write(obj = RG@values, file = tmpChrFile, name = chrRGValues)

          rm(dt, RG)
          gc()

          return(list(tmpChrFile = tmpChrFile, errorCheck = errorCheck))

        }

      }

    }, error = function(e){

      errorList <- list(
        x = x,
        tileChromSizes = tileChromSizes[x],
        fragments = if(exists("dt", inherits = FALSE)) dt else "Error with Fragments!"
      )

      .logError(e, fn = ".tabixToTmp", info = prefix, errorList = errorList, logFile = logFile)

    })


  }, threads = threads)

  if(threads > 1){

    #Parallel Linkage Hdf5
    .logDiffTime(sprintf("%s Parallel Hdf5 Linkage Temporary File", prefix), t1 = tstart, verbose = FALSE, logFile = logFile)

    file.remove(tmpFile)
    o <- h5closeAll()
    fid <- H5Fcreate(tmpFile)

    o <- h5createGroup(fid, paste0("Fragments"))
    o <- h5createGroup(fid, paste0("Metadata"))
    o <- h5write(obj = "Arrow", file = fid, name = "Class")
    o <- h5write(obj = paste0(packageVersion("ArchR")), file = fid, name = "ArchRVersion")
    o <- h5write(obj = "tmp", file = fid, name = "Metadata/Sample")

    tmpChrFiles <- lapply(readTiledChrom, function(x) x$tmpChrFile) %>% unlist

    for(i in seq_along(tmpChrFiles)){

      tmpChrFilei <- tmpChrFiles[i]
      splitNames <- stringr::str_split(h5ls(tmpChrFilei)$name, "._.", simplify=TRUE)
      chunkName <- splitNames[,1]
      group <- splitNames[,2]

      o <- h5createGroup(fid, paste0("Fragments/", chunkName[1]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[1],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[1]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[2],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[2]))

      H5Lcreate_external(target_file_name = tmpChrFilei, 
                         target_obj_name = h5ls(tmpChrFilei)$name[3],
                         link_loc = fid,
                         link_name = paste0("Fragments/", chunkName[1], "/", group[3]))

    }

    H5Fclose(fid)

  }

  errorCheck <- sum(lapply(readTiledChrom, function(x) x$errorCheck) %>% unlist)

  if(errorCheck == 0){
    if(!is.null(validBC)){
      .logStop("No fragments found! Possible error with validBarcodes!", logFile = logFile)
    }else{
      .logStop("No fragments found!", logFile = logFile)
    }
  }

  .logDiffTime(sprintf("%s Successful creation of Temporary File", prefix), t1 = tstart, verbose = verbose, logFile = logFile)

  return(tmpFile)

}


.fastFragmentInfo <- function(
  ArrowFile = NULL,
  cellNames = .availableCells(ArrowFile),
  nucLength = 147,
  prefix = NULL,
  logFile = NULL
  ){

  .logDiffTime(paste0(prefix, " Computing fragment size info!"), t1 = tstart, verbose = FALSE, logFile = logFile)

  #Info to get
  matNuc <- matrix(0, nrow = length(cellNames), ncol = 3)
  nFrags <- rep(0, length(cellNames))
  fragDist <- rep(0, 1000)

  chrArrow <- .availableChr(ArrowFile)
  for(x in seq_along(chrArrow)){

    countFrags <- tryCatch({

      #Read in frags
      fragx <- .getFragsFromArrow(ArrowFile = ArrowFile, chr = chrArrow[x], out = "IRanges", cellNames = cellNames)

      if(length(fragx) > 0){

        mcols(fragx)$RG@values <- S4Vectors::match(mcols(fragx)$RG@values, cellNames)
        nFragsx <- S4Vectors:::tabulate(mcols(fragx)$RG, nbins = length(cellNames))

        #Get Distributions
        fragDistx <- tabulate(width(fragx), nbins = 1000)
        w <- trunc(width(fragx)/nucLength) + 1
        w[w > 3] <- 3

        #Get Nuc Info
        matNucx <- tabulate2dCpp(
          w, xmin = 1, xmax = 3, 
          as.integer(mcols(fragx)$RG), ymin = 1, ymax = length(cellNames)
        )   

        list(nFrags = nFragsx, matNuc = matNucx, fragDist = fragDistx)

      }else{

        list(nFrags = NULL, matNuc = NULL, fragDist = NULL)

      }


    }, error = function(e){

      errorList <- list(
        x = x,
        chr = chrArrow[x],
        fragments = if(exists("fragx", inherits = FALSE)) fragx else "Error with Fragments!",
        nFrags = if(exists("nFragsx", inherits = FALSE)) nFragsx else "Error with Counting Fragment Barcodes!",
        fragDist = if(exists("fragDistx", inherits = FALSE)) fragDistx else "Error with Counting Fragment Width Numbers!",
        matNucx = if(exists("matNucx", inherits = FALSE)) matNucx else "Error with Counting Fragment Nucleosome Spanning Numbers!"
      )

      .logError(e, fn = ".fastFragmentInfo", info = prefix, errorList = errorList, logFile = logFile)

    })

    if(!is.null(countFrags[[1]])){
      nFrags <- nFrags + countFrags[[1]]
    }
    
    if(!is.null(countFrags[[2]])){
      matNuc <- matNuc + countFrags[[2]]
    }

    if(!is.null(countFrags[[3]])){
      fragDist <- fragDist + countFrags[[3]]
    }

  }

  df <- DataFrame(matNuc)
  colnames(df) <- c("nMonoFrags", "nDiFrags", "nMultiFrags")
  df$cellNames <- cellNames
  df$nFrags <- nFrags
  df <- df[,c("cellNames","nFrags","nMonoFrags", "nDiFrags", "nMultiFrags")]

  out <- list(dfSummary = df, fragDistribution  = fragDist)

  .logDiffTime(paste0(prefix, " Finished computing fragment size info!"), t1 = tstart, verbose = FALSE, logFile = logFile)

  return(out)

}

.fastTSSEnrichment <- function(
  TSS = NULL, 
  ArrowFile = NULL, 
  cellNames = NULL, 
  window = 101, 
  norm = 100, 
  flank = 2000, 
  minNorm = 0.2, #Handles low cell reads inflated TSS values
  maxFragSize = NULL,
  threads = 1,
  prefix = NULL,
  logFile = NULL
  ){

  tstart <- Sys.time()

  #Validate
  ArrowFile <- .validArrow(ArrowFile)
  TSS <- .validGRanges(TSS)

  if(is.null(cellNames)){
    cellNames <- .availableCells(ArrowFile)
  }

  #Create Window and Flank
  TSS <- resize(TSS, 1, fix = "start")
  strand(TSS) <- "*"
  TSS <- unique(TSS)
  tssWindow <- resize(TSS, window, "center")
  tssWindow$type <- "window"
  tssFlank <- c(
    #Positive Flank
    GRanges(seqnames(TSS), IRanges(end(TSS) + flank - norm + 1, end(TSS) + flank)),
    #Negative Flank
    GRanges(seqnames(TSS), IRanges(start(TSS) - flank, start(TSS) - flank + norm - 1))
  )
  tssFlank$type <- "flank"
  tssFeatures <- c(tssWindow, tssFlank)
  #.logThis(tssFeatures, paste0(prefix, " tssFeatures"), logFile = logFile)

  #Counting
  countList <- .fastTSSCounts(
    feature = tssFeatures, 
    ArrowFile = ArrowFile, 
    cellNames = cellNames, 
    maxFragSize = maxFragSize,
    threads = threads, 
    prefix = prefix,
    logFile = logFile
  )

  #Normalize per BP
  cWn <- countList$nWindow / window
  cFn <- countList$nFlank / norm

  #Compute scores
  tssScores <- 2 * cWn / (pmax(cFn, minNorm))
  names(tssScores) <- cellNames
  tssScores <- round(tssScores, 3)

  .logDiffTime(paste0(prefix, " Computed TSS Scores!"), t1 = tstart, verbose = FALSE, logFile = logFile)

  return(list(tssScores=tssScores, tssReads=cWn * window))

}

.fastTSSCounts <- function(
  feature = NULL, 
  ArrowFile = NULL, 
  cellNames = NULL,
  maxFragSize = NULL, 
  threads = 1,
  prefix = NULL,
  logFile = NULL
  ){
  
  tstart1 <- Sys.time()
  featureList <- split(feature, seqnames(feature))
  chrArrow <- .availableChr(ArrowFile)
  featureList <- featureList[chrArrow]
  if(length(featureList)==0){
    stop("Error No Overlap in Chromosomes and TSS Chromosomes!")
  }

  #Count
  countDF <- .safelapply(seq_along(featureList), function(x){

    tryCatch({

      #Count Vector
      nWindow <- rep(0, length(cellNames))
      names(nWindow) <- cellNames

      nFlank <- rep(0, length(cellNames))
      names(nFlank) <- cellNames

      ###############################################################################
      # Get Fragments
      ###############################################################################
      featurex <- featureList[[x]]
      fragments <- .getFragsFromArrow(
        ArrowFile = ArrowFile, 
        chr = names(featureList)[x], 
        out = "IRanges", 
        cellNames = cellNames
      )

      if(length(fragments) > 0){
        if(!is.null(maxFragSize)){
          fragments <- fragments[width(fragments) <= maxFragSize]
        }
      }

      if(length(fragments) > 0){

        mcols(fragments)$RG@values <- match(mcols(fragments)$RG@values, cellNames)
        mcols(featurex)$typeIdx <- match(mcols(featurex)$type, c("window", "flank"))
        
        ###############################################################################
        # Count Each Insertion
        ###############################################################################
        for(y in seq_len(2)){
            
          if(y==1){
            temp <- IRanges(start(fragments), width=1)
          }else if(y==2){
            temp <- IRanges(end(fragments), width=1)
          }
          stopifnot(length(temp) == length(fragments))

          o <- findOverlaps(ranges(featurex), temp)
          remove(temp)
          gc()
          
          mat <- tabulate2dCpp(
                x = as.vector(mcols(fragments)$RG[subjectHits(o)]),
                xmin = 1,
                xmax = length(cellNames),
                y = mcols(featurex)$typeIdx[queryHits(o)],
                ymin = 1,
                ymax = 2
            )

          #Add To
          nWindow <- nWindow + mat[1, ]
          nFlank <- nFlank + mat[2, ]
          rm(o, mat)

        }
        
        rm(fragments)
        gc()

      }

      names(nWindow) <- NULL
      names(nFlank) <- NULL

      DataFrame(nWindow = nWindow, nFlank = nFlank)

    }, error = function(e){

      errorList <- list(
        x = x,
        chr = names(featureList)[x],
        fragments = if(exists("fragments", inherits = FALSE)) fragments else "Error with Fragments!",
        features = if(exists("featurex", inherits = FALSE)) featurex else "Error with Features!" 
      )

      .logError(e, fn = ".fastTSSCounts", info = prefix, errorList = errorList, logFile = logFile)

    })


  }, threads = threads)

  for(i in seq_along(countDF)){
    if(i == 1){
      cDF <- countDF[[i]]
    }else{
      cDF$nWindow <- cDF$nWindow + countDF[[i]]$nWindow
      cDF$nFlank <- cDF$nFlank + countDF[[i]]$nFlank
    }
  }

  rm(countDF)
  gc()

  out <- list(nWindow = cDF[,1] , nFlank = cDF[, 2])

  return(out)

}

.fastFeatureCounts <- function(
  feature = NULL, 
  ArrowFile = NULL, 
  cellNames = NULL, 
  threads = 1,
  prefix = NULL,
  logFile = NULL
  ){

  tstart1 <- Sys.time()
  featureNames <- unique(mcols(feature)$FeatureName)
  featureList <- split(feature, seqnames(feature))
  chrArrow <- .availableChr(ArrowFile)
  featureList <- featureList[chrArrow]
  if(length(featureList)==0){
    .logStop("Error No Overlap in Chromosomes and Feature Chromosomes!", logFile = logFile)
  }

  #Count
  countMat <- .safelapply(seq_along(featureList), function(x){

    tryCatch({

      m <- matrix(0, nrow = length(featureNames), ncol = length(cellNames))
      rownames(m) <- featureNames

      ###############################################################################
      # Get Fragments
      ###############################################################################
      featurex <- featureList[[x]]
      fragments <- .getFragsFromArrow(
        ArrowFile = ArrowFile, 
        chr = names(featureList)[x], 
        out = "IRanges", 
        cellNames = cellNames
      )

      if(length(fragments) > 0){

        mcols(fragments)$RG@values <- match(mcols(fragments)$RG@values, cellNames)
        mcols(featurex)$typeIdx <- match(mcols(featurex)$FeatureName, featureNames)
        
        ###############################################################################
        # Count Each Insertion
        ###############################################################################
        for(y in seq_len(2)){
            
          if(y==1){
            temp <- IRanges(start(fragments), width=1)
          }else if(y==2){
            temp <- IRanges(end(fragments), width=1)
          }
          stopifnot(length(temp) == length(fragments))

          o <- findOverlaps(ranges(featurex), temp)
          remove(temp)
          gc()
          
          m <- m + tabulate2dCpp(
                x = as.vector(mcols(fragments)$RG[subjectHits(o)]),
                xmin = 1,
                xmax = length(cellNames),
                y = mcols(featurex)$typeIdx[queryHits(o)],
                ymin = 1,
                ymax = nrow(m)
            )

          rm(o)

        }
        
        rm(fragments)
        gc()

      }

      m

    }, error = function(e){

      errorList <- list(
        x = x,
        chr = names(featureList)[x],
        fragments = if(exists("fragments", inherits = FALSE)) fragments else "Error with Fragments!",
        features = if(exists("featurex", inherits = FALSE)) featurex else "Error with Features!" 
      )

      .logError(e, fn = ".fastFeatureCounts", info = prefix, errorList = errorList, logFile = logFile)

    })

  }, threads = threads)

  countMat <- Reduce("+", countMat)
  colnames(countMat) <- cellNames

  countMat

}

.isTabix <- function(file = NULL){
  tryCatch({
    TabixFile(file)
    TRUE
  }, error = function(x){
    tryCatch({
      if(getArchRVerbose()) message("Attempting to index ", file," as tabix..")
      indexTabix(file, format = "bed")
      TRUE
    }, error = function(y){
      FALSE
    })
  })
}
