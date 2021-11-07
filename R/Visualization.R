#' @export
runUMAP <- function(mapSnapATAC, nPC = 20, seed = 10, weightDimReduct = FALSE) {
  dmat <- mapSnapATAC$dmat
  sdev <- mapSnapATAC$sdev
  nDim <- ncol(dmat)
  mat <- dmat[, 1:min(nPC, nDim)]
  if (weightDimReduct) {
    mat <- mat %*% diag(sdev[1:min(nPC, nDim)])
  }
  set.seed(seed = seed)
  message("Run UMAP")
  umap <- umap::umap(mat, quietly = TRUE)$layout
  colnames(umap) <- c("UMAP-1", "UMAP-2")
  return(umap)
}

#' Create a color panel
#' Ref: From SnapATAC
#' @param num.color Number of diffenrent colors to return
#' @export
createColorPanel <- function(num.color) {
  colPanel <- c(
    "grey", "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
    "#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744",
    "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
    "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
    "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
  )
  if (num.color > length(colPanel)) {
    colPanel <- c(colPanel, colVector(num.color - length(colPanel)))
  } else {
    colPanel <- colPanel[1:num.color]
  }
  return(colPanel)
}


#' Create a color panel
#' Ref: From SnapATAC
#' @param num.color Number of diffenrent colors to return
#' @param type Type of colors to generate
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @export
colVector <- function(num.color = 60, type = c("qual", "div", "seq")) {
  type <- match.arg(type)
  set.seed(10)
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == type, ]
  set.seed(10)
  col_vector <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  return(col_vector)
}

#' Ref: From SnapATAC
findCentrod <- function(x, y){
	x.ls = split(data.frame(x),y);
	centroid.ls = lapply(split(data.frame(x),y), function(xx) apply(xx, 2, median))
	centroid.df = data.frame(do.call(rbind, centroid.ls))
	centroid.df$Y = names(centroid.ls);
	
	return(centroid.df);		
}

#' @importFrom grDevices xy.coords 
#' @importFrom graphics strwidth strheight 
textHalo <- function(
	x, y=NULL, 
	labels, 
	col='white', 
	bg='black', 
	r=0.1,
	... 
){

	theta= seq(0, 2*pi, length.out=50);
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')

    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
}


#' Ref: From SnapATAC
#' @export
plot2DEmbed <- function(dmat,
                        point.size = 1,
                        point.shape = 19,
                        point.alpha = 0.8,
                        clusterLabel = NULL,
                        text.add = TRUE,
                        text.size = 1,
                        text.color = "black",
                        text.halo.add = TRUE,
                        text.halo.color = "white",
                        text.halo.width = 0.2,
                        legend.add = FALSE,
                        legend.pos = c("bottomleft", "bottom", "left", "topleft", "top", "topright", "right", "center"),
                        legend.text.size = 1,
                        legend.text.color = "black",
                        down.sample = 10000,
                        pdf.file.name = NULL,
                        pdf.width = 7,
                        pdf.height = 7,
                        myxlims = NULL,
                        myylims = NULL,
                        ...) {
  ncell <- nrow(dmat)
  down.sample <- as.integer(down.sample)
  data.use <- as.data.frame(dmat)
  if (is.null(myxlims)) {
    xlims <- c(-max(abs(data.use[, 1])) * 1.05, max(abs(data.use[, 1])) * 1.2)
  } else {
    xlims <- myxlims
  }
  if (is.null(myylims)) {
    ylims <- c(-max(abs(data.use[, 2])) * 1.05, max(abs(data.use[, 2])) * 1.2)
  } else {
    ylims <- myylims
  }

  cluster <- clusterLabel
  if (((x <- length(cluster)) == 0L) | (is.null(clusterLabel))) {
    warning("cluster does not exist, text.add is ignored")
    text.add <- FALSE
  }

  if (length(cluster) != 0L) {
    data.use$col <- factor(cluster)
  } else {
    data.use$col <- factor(1)
  }

  if (!is.null(pdf.file.name)) {
    if (file.exists(pdf.file.name)) {
      warning("pdf.file already exists and remove it.")
      file.remove(pdf.file.name)
    } else {
      if (!file.create(pdf.file.name)) {
        stop("cannot create pdf.file, not a directory")
      }
      file.remove(pdf.file.name)
    }
    pdf(pdf.file.name, width = pdf.width, height = pdf.height)
  }

  legend.pos <- match.arg(legend.pos)
  down.sample <- min(down.sample, ncell)
  idx.ds <- sort(sample(seq(ncell), down.sample))
  data.use <- data.use[idx.ds, , drop = FALSE]

  colPanel <- createColorPanel(length(unique(data.use$col)))
  graphics::plot(
    data.use[, c(1, 2)],
    cex = point.size,
    pch = point.shape,
    col = scales::alpha(colPanel[factor(data.use$col)], point.alpha),
    bty = "l",
    font.lab = 2,
    col.axis = "darkgrey",
    xlim = xlims,
    ylim = ylims,
    ...
  )
  graphics::box(bty = "l", lwd = 2)
  if (text.add) {
    xx <- findCentrod(data.use[, c(1, 2)], data.use$col)
    textHalo(
      x = xx[, 1], y = xx[, 2], labels = xx[, 3], col = text.color, bg = text.halo.color,
      r = text.halo.width, cex = text.size
    )
  }

  if (legend.add) {
    legend.ncol <- as.integer(length(levels(data.use$col)) / 25) + 1
    graphics::legend(
      "topright",
      legend = levels(factor(data.use$col)),
      col = colPanel,
      pch = point.shape,
      pt.cex = 1,
      bty = "n",
      cex = legend.text.size,
      text.col = legend.text.color,
      horiz = FALSE,
      ncol = legend.ncol
    )
  }
  if (!is.null(pdf.file.name)) {
    dev.off()
  }
}
