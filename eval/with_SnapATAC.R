library(ggplot2)
library(stringr)
library(scales)
library(data.table)
library(uwot)
library(RColorBrewer)


# analysis QC: nfragment and TSSE
qc_cancers <- read.table(file = "cancer_tissue_qc.tsv", header = TRUE,
                         quote = "", comment.char = "")

qc_QY1287_Taiji <- qc_cancers[qc_cancers$sample == "QY_1287",]
rownames(qc_QY1287_Taiji) <- qc_QY1287_Taiji$barcode

qc_QY1287_SnapATAC <- read.table(file = "out/QY_1287_QC.tsv", header = TRUE,
                                 quote = "", comment.char = "", sep = ",")
rownames(qc_QY1287_SnapATAC) <- qc_QY1287_SnapATAC$barcode
tsse_QY1287_smmtools <- read.table(file = "out/QY_1287_sumFragment.csv", header = TRUE,
                                   quote = "", comment.char = "", sep = ",")
rownames(tsse_QY1287_smmtools) <- tsse_QY1287_smmtools$barcode

fragment_bam210x <- data.frame(tsse_QY1287_smmtools$nUniqFrag,
                               qc_QY1287_SnapATAC[rownames(tsse_QY1287_smmtools), "uniq"])
rownames(fragment_bam210x) <- rownames(tsse_QY1287_smmtools)
colnames(fragment_bam210x) <- c("smmtools", "SnapATAC")

## check fragment from smmtools, and fragment from SnapTools.
## this helps us to know if tabix2h5 file is OK.
plot(x = log10(fragment_bam210x$smmtools), y = log10(fragment_bam210x$SnapATAC))

## check SnapTools and Taiji' number of fragments.
fragment_snap_and_Taiji <- data.frame(qc_QY1287_Taiji$log10frag,
                                      log10(qc_QY1287_SnapATAC[rownames(qc_QY1287_Taiji), "uniq"]))
colnames(fragment_snap_and_Taiji) <- c("Taiji", "SnapATAC")
plot(x = fragment_snap_and_Taiji$Taiji, fragment_snap_and_Taiji$SnapATAC)

par(mfrow = c(1,3))
plot(x = log10(fragment_bam210x[rownames(qc_QY1287_Taiji), "smmtools"]),
     y = log10(fragment_bam210x[rownames(qc_QY1287_Taiji), "SnapATAC"]))
plot(x = fragment_snap_and_Taiji$Taiji, fragment_snap_and_Taiji$SnapATAC)
plot(x = log10(fragment_bam210x[rownames(qc_QY1287_Taiji), "smmtools"]),
     y = fragment_snap_and_Taiji$Taiji)

## check TSSE scores
tsse_smmtools_and_Taiji <- data.frame(tsse_QY1287_smmtools[rownames(qc_QY1287_Taiji), "TSSE"],
                                      qc_QY1287_Taiji$TSSe)
colnames(tsse_smmtools_and_Taiji) <- c("smmtools", "Taiji")
plot(x = tsse_smmtools_and_Taiji$smmtools,
     y = tsse_smmtools_and_Taiji$Taiji, pch = 16,
      col = alpha("blue", 0.4), xlim = c(2, 40), ylim = c(2,40))
abline(coef = c(0,1), col = "red")

## visualize the clustering result
harmonyMatrix <- readRDS(file = "out/harmony_matrix.rds")
## one sample QY1288
harmonyMatrix <- harmonyMatrix[grepl(pattern = "QY_1288", x = rownames(harmonyMatrix)), ]

RB_cluster <- as.data.frame(fread(file = "out/RB_cluster.csv", sep = ","))
rownames(RB_cluster) <- RB_cluster$barcode
RB_cluster <- RB_cluster[rownames(harmonyMatrix),]

umap <- uwot::umap(harmonyMatrix, init = "spectral")
umap <- as.data.frame(umap, row.names = rownames(harmonyMatrix))
colnames(umap) <- c("UMAP1", "UMAP2")
umap$cluster <- as.factor(RB_cluster$cluster)
umap$sample <- as.factor(gsub("#.*", "", rownames(umap)))

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

colVector <- function(num.color = 60, type = c("qual", "div", "seq")) {
  type <- match.arg(type)
  set.seed(10)
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == type, ]
  set.seed(10)
  col_vector <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  return(col_vector)
}

cols <- colVector(length(unique(umap$cluster)))
p <- ggplot(data = umap, aes(UMAP1, UMAP2, colour = cluster)) +
  geom_point() +
  scale_fill_manual(values = cols,aesthetics = c("colour", "fill"))

b <- ggplot(data = umap, aes(UMAP1, UMAP2, colour = sample)) +
  geom_point()

