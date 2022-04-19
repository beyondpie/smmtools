test_that("snapGmat2Seurat works", {
  snap <- readRDS(test_path("testdata", "testSnap.rds"))
  seurat <- snapGmat2Seurat(snap = snap)
  expect_s4_class(seurat, "Seurat")
})

test_that("integrateWithScRNASeq works", {
  snap <- readRDS(test_path("testdata", "testSnap.rds"))
  snapSeurat <- snapGmat2Seurat(snap = snap)
  mbSeurat <- readRDS(test_path("testdata", "testMouseBrainOrgSeurat.rds"))
  coEmbed <- integrateWithScRNASeq(snapSeurat = snapSeurat,
                                   rnaSeurat = mbSeurat,
                                   eigDims = 1:20,
                                   snapAssay = "GeneScore",
                                   preprocessSnap = T,
                                   preprocessRNA = T,
                                   reso = 0.5)
  expect_equal(ncol(coEmbed), ncol(snapSeurat) + ncol(mbSeurat))
})
