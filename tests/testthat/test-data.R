test_that("the TNBC panel has the documented form", {
  expect_equal(dim(tnbc_metabric), c(299L, 1000L))
  expect_false(anyNA(tnbc_metabric))
  expect_null(rownames(tnbc_metabric))
  expect_false(anyDuplicated(colnames(tnbc_metabric)) > 0)
  expect_true(is.numeric(tnbc_metabric))
})
