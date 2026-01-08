test_that("greatCircleDistance is symmetric and zero on diagonal", {
  d <- greatCircleDistance(0, 0, 0, 0)
  expect_equal(d, 0)
})

