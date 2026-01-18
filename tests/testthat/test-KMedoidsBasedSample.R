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

  dist_mat <- sapply(seq_len(nrow(df)), function(i) {
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

  dist_mat <- sapply(seq_len(nrow(df)), function(i) {
    greatCircleDistance(df$lat[i], df$lon[i], df$lat, df$lon)
  })

  input_data <- list(distances = dist_mat, sites = df)

  res <- KMedoidsBasedSample(input_data, n = 2, n_bootstrap = 5, n_restarts = 1, n_local_search_iter = 3, verbose = FALSE)

  # required site should always be in the solution
  expect_true(all(df$site_id[df$required] %in% res$selected_sites))
})

test_that("required column is created if missing", {

  df <- data.frame(
    site_id = 1:5,
    lat = runif(5),
    lon = runif(5),
    coord_uncertainty = 0
  )

  dist_mat <- sapply(seq_len(nrow(df)), function(i) {
    greatCircleDistance(df$lat[i], df$lon[i], df$lat, df$lon)
  })

  res <- KMedoidsBasedSample(
    input_data = list(distances = dist_mat, sites = df),
    n = 2,
    n_bootstrap = 1,
    n_restarts = 1,
    verbose = FALSE
  )

  expect_true("required" %in% names(res$input_data))
})

test_that("dropout disabled works", {

  df <- data.frame(
    site_id = 1:6,
    lat = runif(6),
    lon = runif(6),
    required = FALSE,
    coord_uncertainty = 0
  )

  dist_mat <- sapply(seq_len(nrow(df)), function(i) {
    greatCircleDistance(df$lat[i], df$lon[i], df$lat, df$lon)
  })

  res <- KMedoidsBasedSample(
    list(distances = dist_mat, sites = df),
    n = 3,
    dropout_prob = 0,
    n_bootstrap = 2,
    n_restarts = 1,
    verbose = FALSE
  )

  expect_length(res$selected_sites, 3)
})


test_that("environmental distance skips jitter", {

  df <- data.frame(
    site_id = 1:5,
    lat = runif(5),
    lon = runif(5),
    required = FALSE,
    coord_uncertainty = 50000
  )

  dist_mat <- matrix(runif(25), 5)

  res <- KMedoidsBasedSample(
    list(distances = dist_mat, sites = df),
    n = 2,
    n_bootstrap = 1,
    n_restarts = 1,
    distance_type = "environmental",
    verbose = FALSE
  )

  expect_length(res$selected_sites, 2)
})

test_that("pam_fixed respects fixed medoids", {

  set.seed(1)
  dist_mat <- matrix(runif(25), 5)
  diag(dist_mat) <- 0

  res <- pam_fixed(dist_mat, k = 2, fixed_ids = 1)

  expect_true(1 %in% res$medoids)
  expect_length(res$medoids, 2)
})


test_that("greatCircleDistance returns zero for identical points", {
  d <- greatCircleDistance(0, 0, 0, 0)
  expect_equal(d, 0)
})


test_that("jitter_coords handles zero and NA uncertainty", {

  j <- jitter_coords(
    lat = c(0, 1),
    lon = c(0, 1),
    uncertainty_m = c(NA, 0)
  )

  expect_equal(j$jittered_lat, c(0, 1))
  expect_equal(j$jittered_lon, c(0, 1))
})
