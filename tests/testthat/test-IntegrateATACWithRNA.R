test_that("snapGmat2Seurat works", {
  snap <- readRDS(test_path("testdata", "testSnap.rds"))
  seurat <- snapGmat2Seurat(snap = snap)
  expect_s4_class(seurat, "Seurat")
})
