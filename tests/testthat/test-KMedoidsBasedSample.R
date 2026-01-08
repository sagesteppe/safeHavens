library(testthat)

test_that("KMedoidsBasedSample returns expected structure", {

  # Create 5 points
  n_sites <- 5
  df <- data.frame(
    site_id = 1:n_sites,
    lat = runif(n_sites, 0, 10),
    lon = runif(n_sites, 0, 10),
    required = c(TRUE, rep(FALSE, n_sites-1)),
    coord_uncertainty = 0
  )

  dist_mat <- sapply(1:nrow(df), function(i) {
    greatCircleDistance(df$lat[i], df$lon[i], df$lat, df$lon)
  })

  input_data <- list(distances = dist_mat, sites = df)

  res <- KMedoidsBasedSample(input_data, n = 3, n_bootstrap = 10, n_restarts = 1, n_local_search_iter = 5, verbose = FALSE)

  expect_true(all(c("input_data", "selected_sites", "stability_score", "stability", "settings") %in% names(res)))
  expect_equal(nrow(res$input_data), n_sites)
  expect_true(all(res$input_data$selected %in% c(TRUE, FALSE)))
  expect_true(is.numeric(res$stability_score))
})

test_that("KMedoidsBasedSample respects required sites", {

  df <- data.frame(
    site_id = 1:4,
    lat = c(0, 1, 2, 3),
    lon = c(0, 1, 2, 3),
    required = c(TRUE, FALSE, FALSE, FALSE),
    coord_uncertainty = 0
  )

  dist_mat <- sapply(1:nrow(df), function(i) {
    greatCircleDistance(df$lat[i], df$lon[i], df$lat, df$lon)
  })

  input_data <- list(distances = dist_mat, sites = df)

  res <- KMedoidsBasedSample(input_data, n = 2, n_bootstrap = 5, n_restarts = 1, n_local_search_iter = 3, verbose = FALSE)

  # required site should always be in the solution
  expect_true(all(df$site_id[df$required] %in% res$selected_sites))
})
