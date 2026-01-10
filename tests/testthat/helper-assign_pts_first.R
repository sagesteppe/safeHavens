make_test_data <- function() {
  x <- data.frame(
    g1 = c(1, 2, 3, 4),
    g2 = c(4, 3, 2, 1),
    g3 = c(2, 2, 2, 2)
  )
  props <- c(g1 = 2, g2 = 1, g3 = 1)
  nf_pct <- c(g1 = 0.1, g2 = 0.9, g3 = 0.5)
  list(x = x, props = props, nf_pct = nf_pct)
}

test_that("assign_pts_frst returns expected structure", {
  d <- make_test_data()
  res <- assign_pts_frst(d$x, d$props, d$nf_pct)

  expect_s3_class(res, "data.frame")
  expect_true(all(c("Assignment", "ID") %in% names(res)))
  expect_equal(res$ID, seq_len(nrow(d$x)))
})

test_that("number of assignments does not exceed props", {
  d <- make_test_data()
  res <- assign_pts_frst(d$x, d$props, d$nf_pct)

  assigned <- table(na.omit(res$Assignment))

  for (nm in names(assigned)) {
    expect_lte(assigned[[nm]], d$props[[nm]])
  }
})

test_that("points are not assigned more than once", {
  d <- make_test_data()
  res <- assign_pts_frst(d$x, d$props, d$nf_pct)

  expect_false(any(duplicated(res$ID[!is.na(res$Assignment)])))
})

test_that("props larger than available points do not error", {
  x <- data.frame(a = c(1, 2), b = c(2, 1))
  props <- c(a = 10, b = 10)
  nf_pct <- c(a = 0.9, b = 0.1)

  res <- assign_pts_frst(x, props, nf_pct)

  expect_equal(sum(!is.na(res$Assignment)), nrow(x))
})

test_that("zero props results in no assignments for that grid", {
  x <- data.frame(a = c(1, 2), b = c(2, 1))
  props <- c(a = 0, b = 2)
  nf_pct <- c(a = 0.9, b = 0.1)

  res <- assign_pts_frst(x, props, nf_pct)

  expect_false(any(res$Assignment == "a", na.rm = TRUE))
})

test_that("loop continues safely after all points assigned", {
  x <- data.frame(a = c(1, 2), b = c(2, 1))
  props <- c(a = 2, b = 2)
  nf_pct <- c(a = 0.9, b = 0.1)

  expect_silent(assign_pts_frst(x, props, nf_pct))
})

