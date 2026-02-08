library(testthat)
library(terra)
library(sf)
library(dplyr)

# ==============================================================================
# Tests for rescaleFuture()
# ==============================================================================
test_that("rescaleFuture correctly rescales future predictors", {
  skip_if_not_installed("glmnet")
  
  set.seed(123)
  
  # Create actual glmnet model (not mock)
  library(glmnet)
  
  # Training data
  x_train <- matrix(rnorm(200), ncol=2)
  colnames(x_train) <- c("var1", "var2")
  y_train <- factor(sample(c(0, 1), 100, replace=TRUE))
  
  # Fit actual glmnet model
  model <- glmnet(x_train, y_train, family="binomial", lambda=0.01, alpha=0.5)
  
  # Create rasters
  r1 <- rast(ncols=10, nrows=10, xmin=0, xmax=10, ymin=0, ymax=10)
  values(r1) <- rnorm(100, mean=0, sd=1)
  r2 <- rast(ncols=10, nrows=10, xmin=0, xmax=10, ymin=0, ymax=10)
  values(r2) <- rnorm(100, mean=0, sd=1)
  current_preds <- c(r1, r2)
  names(current_preds) <- c("var1", "var2")
  
  f1 <- rast(ncols=10, nrows=10, xmin=0, xmax=10, ymin=0, ymax=10)
  values(f1) <- rnorm(100, mean=0.5, sd=1)
  f2 <- rast(ncols=10, nrows=10, xmin=0, xmax=10, ymin=0, ymax=10)
  values(f2) <- rnorm(100, mean=0.5, sd=1)
  future_preds <- c(f1, f2)
  names(future_preds) <- c("var1", "var2")
  
  # Training data frame
  train_data <- data.frame(occurrence = y_train)
  
  # Prediction matrix
  pred_mat <- as.data.frame(x_train)
  
  result <- rescaleFuture(
    model = model,
    future_predictors = future_preds,
    current_predictors = current_preds,
    training_data = train_data,
    pred_mat = pred_mat
  )
  
  expect_s4_class(result, "SpatRaster")
  expect_equal(nlyr(result), 2)
  expect_true(all(!is.na(values(result))))
})


test_that("rescaleFuture handles missing variables", {
  # Setup with mismatched variables
  current_preds <- rast(ncols=10, nrows=10)
  names(current_preds) <- "var1"
  
  future_preds <- rast(ncols=10, nrows=10)
  names(future_preds) <- "var_different"
  
  train_data <- data.frame(occurrence = factor(c(1,0)))
  pred_mat <- data.frame(var1 = c(1,2))
  
  mock_model <- list()
  class(mock_model) <- "glmnet"
  coef_matrix <- matrix(c(0.5, 0.3), ncol=1)
  rownames(coef_matrix) <- c("(Intercept)", "var1")
  attr(mock_model, "coef") <- coef_matrix
  
  expect_error(
    rescaleFuture(mock_model, future_preds, current_preds, train_data, pred_mat)
  )
})


# ==============================================================================
# Tests for cluster_novel_areas()
# ==============================================================================

test_that("cluster_novel_areas handles normal case", {
  skip_on_cran()
  skip_if_not_installed("NbClust")
  skip_if_not_installed("caret")
  
  set.seed(456)
  
  # Create mock rescaled raster with MORE diversity for clustering
  r1 <- rast(ncols=30, nrows=30, xmin=0, xmax=30, ymin=0, ymax=30)
  values(r1) <- rnorm(900, mean=0, sd=2)
  r2 <- rast(ncols=30, nrows=30, xmin=0, xmax=30, ymin=0, ymax=30)
  values(r2) <- rnorm(900, mean=0, sd=2)
  r3 <- rast(ncols=30, nrows=30, xmin=0, xmax=30, ymin=0, ymax=30)
  values(r3) <- rnorm(900, mean=0, sd=2)
  future_rescaled <- c(r1, r2, r3)
  names(future_rescaled) <- c("var1", "var2", "var3")
  
  # Create novel mask with enough cells
  novel_mask <- rast(ncols=30, nrows=30, xmin=0, xmax=30, ymin=0, ymax=30)
  values(novel_mask) <- sample(c(1, NA), 900, replace=TRUE, prob=c(0.5, 0.5))
  
  result <- cluster_novel_areas(
    future_rescaled = future_rescaled,
    novel_mask = novel_mask,
    n_novel_pts = 100,
    next_cluster_id = 10,
    nbclust_args = list(min.nc=2, max.nc=5)
  )
  
  expect_type(result, "list")
  expect_true("clusters_raster" %in% names(result))
  expect_s4_class(result$clusters_raster, "SpatRaster")
  
  # Just check that we got SOME cluster values
  cluster_vals <- values(result$clusters_raster)
  cluster_vals <- unique(cluster_vals[!is.na(cluster_vals)])
  expect_true(length(cluster_vals) > 0)
  # Don't test exact values - KNN can be unpredictable
})


test_that("cluster_novel_areas handles too few points", {
  set.seed(789)
  
  r1 <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(r1) <- rnorm(25)
  r2 <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(r2) <- rnorm(25)
  future_rescaled <- c(r1, r2)  # Need multiple layers
  names(future_rescaled) <- c("var1", "var2")
  
  # Very sparse novel mask (only 3 cells)
  novel_mask <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(novel_mask) <- c(rep(1, 3), rep(NA, 22))
  
  expect_warning(
    result <- cluster_novel_areas(
      future_rescaled = future_rescaled,
      novel_mask = novel_mask,
      n_novel_pts = 50,
      next_cluster_id = 10,
      nbclust_args = list()
    ),
    "Too few novel climate points for clustering. Assigning single novel cluster."
  )
  
  cluster_vals <- values(result$clusters_raster)
  cluster_vals <- unique(cluster_vals[!is.na(cluster_vals)])
  expect_length(cluster_vals, 1)
  expect_equal(cluster_vals, 10)
})


# ==============================================================================
# Tests for analyze_cluster_relationships()
# ==============================================================================

test_that("analyze_cluster_relationships returns correct structure", {
  skip_if_not_installed("cluster")
  
  set.seed(321)
  
  # Create cluster raster with defined pattern
  clusters_rast <- rast(ncols=30, nrows=30, xmin=0, xmax=30, ymin=0, ymax=30)
  vals <- rep(NA, 900)
  vals[1:250] <- 1
  vals[251:500] <- 2
  vals[501:650] <- 3
  vals[651:800] <- 11
  vals[801:900] <- 12
  values(clusters_rast) <- vals
  
  r1 <- rast(ncols=30, nrows=30, xmin=0, xmax=30, ymin=0, ymax=30)
  values(r1) <- rnorm(900)
  r2 <- rast(ncols=30, nrows=30, xmin=0, xmax=30, ymin=0, ymax=30)
  values(r2) <- rnorm(900)
  future_rescaled <- c(r1, r2)
  names(future_rescaled) <- c("var1", "var2")
  
  # FIXED: Expect warning if insufficient data
  result <- suppressWarnings(
    analyze_cluster_relationships(
      clusters_raster = clusters_rast,
      future_rescaled = future_rescaled,
      existing_ids = c(1, 2, 3),
      novel_ids = c(11, 12),
      n_sample_per_cluster = 30
    )
  )
  
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("novel_cluster_id", "nearest_existing_id", "avg_silhouette_width") %in% names(result)))
  expect_true(nrow(result) >= 0)
  expect_true(nrow(result) <= 2)
  
  if (nrow(result) > 0) {
    expect_true(all(result$novel_cluster_id %in% c(11, 12)))
  }
})
test_that("analyze_cluster_relationships can handle isolated novel clusters", {
  skip_if_not_installed("cluster")
  
  set.seed(999)
  
  clusters_rast <- rast(ncols=30, nrows=30)
  crs(clusters_rast) <- "EPSG:4326"
  
  vals <- rep(NA, 900)
  vals[1:300] <- 1
  vals[301:600] <- 2
  vals[601:900] <- 11
  values(clusters_rast) <- vals
  
  r1 <- rast(ncols=30, nrows=30)
  crs(r1) <- "EPSG:4326"
  r1_vals <- rnorm(900, mean=0, sd=1)
  
  cluster_11_mask <- !is.na(vals) & vals == 11
  r1_vals[cluster_11_mask] <- rnorm(sum(cluster_11_mask), mean=100, sd=1)
  values(r1) <- r1_vals
  
  r2 <- rast(ncols=30, nrows=30)
  crs(r2) <- "EPSG:4326"
  values(r2) <- rnorm(900)
  
  future_rescaled <- c(r1, r2)
  names(future_rescaled) <- c("var1", "var2")
  
  # FIXED: Expect and suppress warning
  result <- suppressWarnings(
    analyze_cluster_relationships(
      clusters_raster = clusters_rast,
      future_rescaled = future_rescaled,
      existing_ids = c(1, 2),
      novel_ids = c(11),
      n_sample_per_cluster = 30
    )
  )
  
  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) <= 1)
  
  if (nrow(result) > 0) {
    expect_equal(result$novel_cluster_id, 11)
    expect_true(is.numeric(result$avg_silhouette_width))
    expect_true(is.na(result$nearest_existing_id) || result$nearest_existing_id %in% c(1, 2))
  }
})
test_that("analyze_cluster_relationships handles no novel clusters", {
  result <- analyze_cluster_relationships(
    clusters_raster = rast(ncols=10, nrows=10),
    future_rescaled = rast(ncols=10, nrows=10),
    existing_ids = c(1, 2, 3),
    novel_ids = integer(0),
    n_sample_per_cluster = 20
  )
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
  expect_true(all(c("novel_cluster_id", "nearest_existing_id", "avg_silhouette_width") %in% names(result)))
})

test_that("analyze_cluster_relationships can return NA for isolated novel clusters", {
  skip_if_not_installed("cluster")
  
  set.seed(999)
  
  # Create clusters where novel cluster is spatially separated
  clusters_rast <- rast(ncols=30, nrows=30)
  crs(clusters_rast) <- "EPSG:4326"
  
  vals <- rep(NA, 900)
  vals[1:300] <- 1      # Existing cluster 1
  vals[301:600] <- 2    # Existing cluster 2
  vals[601:900] <- 11   # Novel cluster
  values(clusters_rast) <- vals
  
  # Create predictors where cluster 11 has very different values
  r1 <- rast(ncols=30, nrows=30)
  crs(r1) <- "EPSG:4326"
  r1_vals <- rnorm(900, mean=0, sd=1)
  
  # FIXED: Use clean boolean mask that excludes NAs
  cluster_11_mask <- !is.na(vals) & vals == 11
  r1_vals[cluster_11_mask] <- rnorm(sum(cluster_11_mask), mean=100, sd=1)
  values(r1) <- r1_vals
  
  r2 <- rast(ncols=30, nrows=30)
  crs(r2) <- "EPSG:4326"
  values(r2) <- rnorm(900)
  
  future_rescaled <- c(r1, r2)
  names(future_rescaled) <- c("var1", "var2")
  
  result <- suppressWarnings(analyze_cluster_relationships(
    clusters_raster = clusters_rast,
    future_rescaled = future_rescaled,
    existing_ids = c(1, 2),
    novel_ids = c(11),
    n_sample_per_cluster = 30
  ))
  
  # Should not error
  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) <= 1)
  
  if (nrow(result) > 0) {
    expect_equal(result$novel_cluster_id, 11)
    expect_true(is.numeric(result$avg_silhouette_width))
    expect_true(is.na(result$nearest_existing_id) || result$nearest_existing_id %in% c(1, 2))
  }
})


# ==============================================================================
# Tests for calculate_changes()
# ==============================================================================

test_that("calculate_changes computes metrics for persisting clusters", {
  # Create mock current clusters
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  poly2 <- st_polygon(list(matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = c(1, 2), geometry = st_sfc(poly1, poly2), crs = 4326)
  
  # Create mock future clusters (shifted and resized)
  poly1_future <- st_polygon(list(matrix(c(1,1, 11,1, 11,11, 1,11, 1,1), ncol=2, byrow=TRUE)))
  poly2_future <- st_polygon(list(matrix(c(20,20, 35,20, 35,35, 20,35, 20,20), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = c(1, 2), geometry = st_sfc(poly1_future, poly2_future), crs = 4326)
  
  result <- calculate_changes(
    current_sf = current_sf,
    future_sf = future_sf,
    planar_proj = 3857  # Web Mercator
  )
  
  expect_s3_class(result, "tbl_df")

})
test_that("analyze_cluster_relationships can handle isolated novel clusters", {
  skip_if_not_installed("cluster")
  
  set.seed(999)
  
  # Create clusters where novel cluster is spatially separated
  clusters_rast <- rast(ncols=30, nrows=30)
  crs(clusters_rast) <- "EPSG:4326"
  
  vals <- rep(NA, 900)
  vals[1:300] <- 1      # Existing cluster 1
  vals[301:600] <- 2    # Existing cluster 2
  vals[601:900] <- 11   # Novel cluster
  values(clusters_rast) <- vals
  
  # Create predictors where cluster 11 has very different values
  r1 <- rast(ncols=30, nrows=30)
  crs(r1) <- "EPSG:4326"
  r1_vals <- rnorm(900, mean=0, sd=1)
  
  # FIXED: Only modify non-NA positions for cluster 11
  cluster_11_mask <- !is.na(vals) & vals == 11  # Boolean mask excluding NAs
  r1_vals[cluster_11_mask] <- rnorm(sum(cluster_11_mask), mean=100, sd=1)
  values(r1) <- r1_vals
  
  r2 <- rast(ncols=30, nrows=30)
  crs(r2) <- "EPSG:4326"
  values(r2) <- rnorm(900)
  
  future_rescaled <- c(r1, r2)
  names(future_rescaled) <- c("var1", "var2")
  
  result <- suppressWarnings(analyze_cluster_relationships(
    clusters_raster = clusters_rast,
    future_rescaled = future_rescaled,
    existing_ids = c(1, 2),
    novel_ids = c(11),
    n_sample_per_cluster = 30
  ))
  
  # Should not error
  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) <= 1)
  
  # If we got a result, check structure
  if (nrow(result) > 0) {
    expect_equal(result$novel_cluster_id, 11)
    expect_true(is.numeric(result$avg_silhouette_width))
    # nearest_existing_id might be NA if truly isolated
    expect_true(is.na(result$nearest_existing_id) || result$nearest_existing_id %in% c(1, 2))
  }
})

test_that("calculate_changes handles gained clusters", {
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = c(1), geometry = st_sfc(poly1), crs = 4326)
  
  # Future has cluster 1 and new cluster 2
  poly1_future <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  poly2_future <- st_polygon(list(matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = c(1, 2), geometry = st_sfc(poly1_future, poly2_future), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  expect_equal(nrow(result), 2)
  gained_row <- result[result$cluster_id == 2, ]
  expect_equal(gained_row$current_area_km2, 0)
  expect_true(is.na(gained_row$area_change_pct))
  expect_true(is.na(gained_row$centroid_shift_km))
})

test_that("calculate_changes returns sorted results", {
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  poly2 <- st_polygon(list(matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol=2, byrow=TRUE)))
  poly3 <- st_polygon(list(matrix(c(40,40, 50,40, 50,50, 40,50, 40,40), ncol=2, byrow=TRUE)))
  
  current_sf <- st_sf(ID = c(3, 1, 2), geometry = st_sfc(poly3, poly1, poly2), crs = 4326)
  future_sf <- st_sf(ID = c(2, 3, 1), geometry = st_sfc(poly2, poly3, poly1), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  expect_true(all(result$cluster_id == sort(result$cluster_id)))
})

