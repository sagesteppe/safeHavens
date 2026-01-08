test_that("bootstrap iteration returns valid solution size", {
  d <- as.matrix(dist(matrix(runif(40), ncol = 2)))
  df <- data.frame(lat = runif(20), lon = runif(20))
  
  res <- run_bootstrap_iteration(
    distances = d,
    sites_df = df,
    n = 3,
    seeds = integer(0),
    available_sites = 1:20,
    uncertain_idx = integer(0),
    n_local_search_iter = 10,
    n_restarts = 1,
    distance_type = "geographic"
  )
  
  expect_length(res$solution, 3)
})

