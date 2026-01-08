test_that("pam_fixed respects fixed medoids", {
  d <- as.matrix(dist(matrix(runif(20), ncol = 2)))
  res <- pam_fixed(d, k = 2, fixed_ids = 1)
  expect_true(1 %in% res$medoids)
})

