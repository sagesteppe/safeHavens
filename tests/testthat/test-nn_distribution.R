make_pres_sf <- function(n = 5) {
  sf::st_as_sf(
    data.frame(
      id = seq_len(n),
      x = seq_len(n),
      y = seq_len(n)
    ),
    coords = c("x", "y"),
    crs = 5070
  )
}

test_that("nn_distribution returns distances for each index", {

  pres <- make_pres_sf(5)
  folds <- list(c(1, 3, 5))

  res <- nn_distribution(folds, pres)

  expect_type(res, "double")
  expect_length(res, length(unlist(folds)))
})

test_that("nn_distribution handles multiple folds", {

  pres <- make_pres_sf(6)
  folds <- list(c(1, 2), c(4, 6))

  res <- nn_distribution(folds, pres)

  expect_length(res, 4)
})
