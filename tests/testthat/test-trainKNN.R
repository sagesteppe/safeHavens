test_that("trainKNN partitions data and fits a knn model", {

  skip_if_not_installed("caret")
  set.seed(123)

  x <- data.frame(
    ID = factor(rep(c("A", "B"), each = 10)),
    v1 = rnorm(20),
    v2 = rnorm(20)
  )

  res <- suppressWarnings(trainKNN(x, split_prop = 0.7))

  expect_type(res, "list")
  expect_named(res, c("fit.knn", "confusionMatrix"))
  expect_s3_class(res$fit.knn, "train")
  expect_s3_class(res$confusionMatrix, "confusionMatrix")
})
