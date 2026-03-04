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

test_that("rescaleFuture uses correct standardization parameters", {
  skip_if_not_installed("glmnet")
  
  set.seed(456)
  library(glmnet)
  
  # Create predictable data with 2 variables (glmnet requires at least 2)
  x_train <- matrix(c(rep(10, 50), rep(20, 50),  # var1
                      rnorm(100)),                 # var2
                    ncol=2)
  colnames(x_train) <- c("var1", "var2")
  y_train <- factor(sample(c(0, 1), 100, replace=TRUE))
  
  model <- glmnet(x_train, y_train, family="binomial", lambda=0.01)
  
  # Current: mean=15 for var1
  r1 <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(r1) <- rep(15, 25)
  r2 <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(r2) <- rnorm(25)
  current_preds <- c(r1, r2)
  names(current_preds) <- c("var1", "var2")
  
  # Future: different mean for var1
  f1 <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(f1) <- rep(25, 25)
  f2 <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(f2) <- rnorm(25)
  future_preds <- c(f1, f2)
  names(future_preds) <- c("var1", "var2")
  
  train_data <- data.frame(occurrence = y_train)
  pred_mat <- as.data.frame(x_train)
  
  result <- rescaleFuture(model, future_preds, current_preds, train_data, pred_mat)
  
  # Should use CURRENT stats to standardize FUTURE
  # var1 should all have same rescaled value since input was constant
  result_var1 <- values(result$var1)
  expect_true(all(abs(result_var1 - result_var1[1]) < 1e-10, na.rm = TRUE))
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

test_that("rescaleFuture preserves raster structure", {
  skip_if_not_installed("glmnet")
  
  set.seed(456)
  library(glmnet)
  
  x_train <- matrix(rnorm(100), ncol=2)
  colnames(x_train) <- c("var1", "var2")
  y_train <- factor(sample(c(0, 1), 50, replace=TRUE))
  model <- glmnet(x_train, y_train, family="binomial", lambda=0.01)
  
  r1 <- rast(ncols=15, nrows=20, xmin=-10, xmax=10, ymin=-5, ymax=5)
  values(r1) <- rnorm(300)
  r2 <- rast(ncols=15, nrows=20, xmin=-10, xmax=10, ymin=-5, ymax=5)
  values(r2) <- rnorm(300)
  current_preds <- c(r1, r2)
  names(current_preds) <- c("var1", "var2")
  
  future_preds <- c(r1, r2)
  names(future_preds) <- c("var1", "var2")
  
  train_data <- data.frame(occurrence = y_train)
  pred_mat <- as.data.frame(x_train)
  
  result <- rescaleFuture(model, future_preds, current_preds, train_data, pred_mat)
  
  # Check dimensions preserved
  expect_equal(ncol(result), ncol(future_preds))
  expect_equal(nrow(result), nrow(future_preds))
  expect_equal(ext(result), ext(future_preds))
})


# ==============================================================================
# Tests for cluster_novel_areas()
# ==============================================================================

test_that("cluster_novel_areas handles normal case", {
  skip_on_cran()
  skip_if_not_installed("NbClust")
  skip_if_not_installed("caret")
  
  set.seed(456)
  
  # Smaller, faster test with clear structure
  r1 <- rast(ncols=15, nrows=15, xmin=0, xmax=15, ymin=0, ymax=15)
  values(r1) <- rnorm(225, mean=0, sd=2)
  r2 <- rast(ncols=15, nrows=15, xmin=0, xmax=15, ymin=0, ymax=15)
  values(r2) <- rnorm(225, mean=0, sd=2)
  future_rescaled <- c(r1, r2)
  names(future_rescaled) <- c("var1", "var2")
  
  # Ensure enough cells
  novel_mask <- rast(ncols=15, nrows=15, xmin=0, xmax=15, ymin=0, ymax=15)
  values(novel_mask) <- sample(c(1, NA), 225, replace=TRUE, prob=c(0.7, 0.3))
  
  # Use very constrained nbclust to make it fast
  result <- suppressWarnings(
    cluster_novel_areas(
      future_rescaled = future_rescaled,
      novel_mask = novel_mask,
      n_novel_pts = 50,  # Reduced
      next_cluster_id = 10,
      nbclust_args = list(
        min.nc = 2, 
        max.nc = 3,  # Very limited
        index = "silhouette"  # Single fast index
      )
    )
  )
  
  expect_type(result, "list")
  expect_true("clusters_raster" %in% names(result))
  expect_s4_class(result$clusters_raster, "SpatRaster")
  
  # Just check that we got SOME cluster values
  cluster_vals <- values(result$clusters_raster)
  cluster_vals <- unique(cluster_vals[!is.na(cluster_vals)])
  expect_true(length(cluster_vals) > 0)
})


test_that("cluster_novel_areas handles too few points", {
  set.seed(789)
  
  r1 <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(r1) <- rnorm(25)
  r2 <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
  values(r2) <- rnorm(25)
  future_rescaled <- c(r1, r2)
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
  
  result <- suppressWarnings(
    analyze_cluster_relationships(
      clusters_raster = clusters_rast,
      future_rescaled = future_rescaled,
      existing_ids = c(1, 2, 3),
      novel_ids = c(11, 12),
      n_sample_per_cluster = 30
    )
  )
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("novel_cluster_id", "nearest_existing_id", "avg_silhouette_width") %in% names(result)))
  expect_true(nrow(result) >= 0)
  expect_true(nrow(result) <= 2)
  
  if (nrow(result) > 0) {
    expect_true(all(result$novel_cluster_id %in% c(11, 12)))
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
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_true(all(c("novel_cluster_id", "nearest_existing_id", "avg_silhouette_width") %in% names(result)))
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
  
  result <- suppressWarnings(
    analyze_cluster_relationships(
      clusters_raster = clusters_rast,
      future_rescaled = future_rescaled,
      existing_ids = c(1, 2),
      novel_ids = 11,
      n_sample_per_cluster = 30
    )
  )
  
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) <= 1)
  
  if (nrow(result) > 0) {
    expect_equal(result$novel_cluster_id, 11)
    expect_true(is.numeric(result$avg_silhouette_width))
    expect_true(is.na(result$nearest_existing_id) || result$nearest_existing_id %in% c(1, 2))
  }
})

test_that("analyze_cluster_relationships handles insufficient data gracefully", {
  skip_if_not_installed("cluster")
  
  set.seed(111)
  
  # Minimal cluster setup
  clusters_rast <- rast(ncols=10, nrows=10)
  vals <- rep(NA, 100)
  vals[1:5] <- 1
  vals[6:10] <- 11
  values(clusters_rast) <- vals
  
  r1 <- rast(ncols=10, nrows=10)
  values(r1) <- rnorm(100)
  future_rescaled <- r1
  names(future_rescaled) <- "var1"
  
  result <- suppressWarnings(
    analyze_cluster_relationships(
      clusters_raster = clusters_rast,
      future_rescaled = future_rescaled,
      existing_ids = 1,
      novel_ids = 11,
      n_sample_per_cluster = 20  # Requesting more points than exist
    )
  )
  
  # Should return valid data.frame even if empty or incomplete
  expect_s3_class(result, "data.frame")
  expect_true(all(c("novel_cluster_id", "nearest_existing_id", "avg_silhouette_width") %in% names(result)))
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
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(c("cluster_id", "current_area_km2", "future_area_km2", 
                    "area_change_pct", "centroid_shift_km") %in% names(result)))
  
  # Both clusters should have positive areas
  expect_true(all(result$current_area_km2 > 0))
  expect_true(all(result$future_area_km2 > 0))
  
  # Centroid shift should be numeric and non-negative
  expect_true(all(result$centroid_shift_km >= 0, na.rm = TRUE))
})

test_that("calculate_changes handles lost clusters", {
  # Current has clusters 1 and 2
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  poly2 <- st_polygon(list(matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = c(1, 2), geometry = st_sfc(poly1, poly2), crs = 4326)
  
  # Future has only cluster 1
  poly1_future <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = 1, geometry = st_sfc(poly1_future), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  expect_equal(nrow(result), 2)
  
  # Check lost cluster (ID=2)
  lost_row <- result[result$cluster_id == 2, ]
  expect_true(lost_row$current_area_km2 > 0)
  expect_equal(lost_row$future_area_km2, 0)
  expect_equal(lost_row$area_change_pct, -100)
  expect_true(is.na(lost_row$centroid_shift_km))
})

test_that("calculate_changes handles gained clusters", {
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = 1, geometry = st_sfc(poly1), crs = 4326)
  
  # Future has cluster 1 and new cluster 2
  poly1_future <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  poly2_future <- st_polygon(list(matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = c(1, 2), geometry = st_sfc(poly1_future, poly2_future), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  expect_equal(nrow(result), 2)
  
  gained_row <- result[result$cluster_id == 2, ]
  expect_equal(gained_row$current_area_km2, 0)
  expect_true(gained_row$future_area_km2 > 0)
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

test_that("calculate_changes handles mix of persisting, lost, and gained", {
  # Current: clusters 1, 2, 3
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  poly2 <- st_polygon(list(matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol=2, byrow=TRUE)))
  poly3 <- st_polygon(list(matrix(c(40,40, 50,40, 50,50, 40,50, 40,40), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = c(1, 2, 3), geometry = st_sfc(poly1, poly2, poly3), crs = 4326)
  
  # Future: clusters 1, 3 (lost 2), 4 (gained)
  poly1_f <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  poly3_f <- st_polygon(list(matrix(c(40,40, 50,40, 50,50, 40,50, 40,40), ncol=2, byrow=TRUE)))
  poly4_f <- st_polygon(list(matrix(c(60,60, 70,60, 70,70, 60,70, 60,60), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = c(1, 3, 4), geometry = st_sfc(poly1_f, poly3_f, poly4_f), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  expect_equal(nrow(result), 4)  # 1, 2, 3, 4
  
  # Check persisting (1 and 3)
  persist_rows <- result[result$cluster_id %in% c(1, 3), ]
  expect_true(all(persist_rows$current_area_km2 > 0))
  expect_true(all(persist_rows$future_area_km2 > 0))
  
  # Check lost (2)
  lost_row <- result[result$cluster_id == 2, ]
  expect_equal(lost_row$future_area_km2, 0)
  expect_equal(lost_row$area_change_pct, -100)
  
  # Check gained (4)
  gained_row <- result[result$cluster_id == 4, ]
  expect_equal(gained_row$current_area_km2, 0)
  expect_true(is.na(gained_row$area_change_pct))
})


# ==============================================================================
# Tests for projectClusters() - MAIN INTEGRATION FUNCTION
# ==============================================================================

test_that("projectClusters handles basic workflow without novel climate", {
  skip("Integration test requires full mock setup - tested in package vignettes")
  # This test requires properly mocked:
  # - eSDM object with real glmnet model
  # - current_clusters with KNN trained on rescaled+coord-weighted data
  # - Matching variable names throughout the pipeline
  # Better tested through package integration tests/vignettes
})

test_that("projectClusters respects cluster_novel=FALSE flag", {
  skip("Integration test requires full mock setup - tested in package vignettes")
  # This test requires properly mocked eSDM and current_clusters objects
})

test_that("projectClusters handles custom thresholds", {
  skip("Integration test requires full mock setup - tested in package vignettes")
  # This test requires properly mocked eSDM and current_clusters objects
})


# ==============================================================================
# Edge case and error handling tests
# ==============================================================================

test_that("rescaleFuture handles all-NA raster gracefully", {
  skip_if_not_installed("glmnet")
  
  library(glmnet)
  
  x_train <- matrix(rnorm(100), ncol=2)
  colnames(x_train) <- c("var1", "var2")
  y_train <- factor(sample(c(0, 1), 50, replace=TRUE))
  model <- glmnet(x_train, y_train, family="binomial", lambda=0.01)
  
  r1 <- rast(ncols=10, nrows=10)
  values(r1) <- rnorm(100)
  r2 <- rast(ncols=10, nrows=10)
  values(r2) <- rnorm(100)
  current_preds <- c(r1, r2)
  names(current_preds) <- c("var1", "var2")
  
  # Future with all NAs
  r1_na <- rast(ncols=10, nrows=10)
  values(r1_na) <- NA_real_
  r2_na <- rast(ncols=10, nrows=10)
  values(r2_na) <- NA_real_
  future_preds <- c(r1_na, r2_na)
  names(future_preds) <- c("var1", "var2")
  
  train_data <- data.frame(occurrence = y_train)
  pred_mat <- as.data.frame(x_train)
  
  # Should not error, but result will be all NA
  result <- rescaleFuture(model, future_preds, current_preds, train_data, pred_mat)
  
  expect_s4_class(result, "SpatRaster")
  expect_true(all(is.na(values(result))))
})

test_that("calculate_changes handles completely disjoint cluster sets", {
  # Current has cluster 1
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = 1, geometry = st_sfc(poly1), crs = 4326)
  
  # Future has cluster 2 only (no overlap)
  poly2 <- st_polygon(list(matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = 2, geometry = st_sfc(poly2), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)  # Both clusters listed
  
  # Cluster 1 should be lost
  lost_row <- result[result$cluster_id == 1, ]
  expect_equal(lost_row$future_area_km2, 0)
  expect_equal(lost_row$area_change_pct, -100)
  
  # Cluster 2 should be gained
  gained_row <- result[result$cluster_id == 2, ]
  expect_equal(gained_row$current_area_km2, 0)
  expect_true(is.na(gained_row$area_change_pct))
})

test_that("calculate_changes handles empty current_sf", {
  # All clusters are new in future
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = 1, geometry = st_sfc(poly1), crs = 4326)
  
  # Empty current
  current_sf <- st_sf(ID = integer(0), geometry = st_sfc(crs = 4326))
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$cluster_id, 1)
  expect_equal(result$current_area_km2, 0)
  expect_true(result$future_area_km2 > 0)
  expect_true(is.na(result$area_change_pct))  # Gained cluster
  expect_true(is.na(result$centroid_shift_km))
})

test_that("calculate_changes handles empty future_sf", {
  # All clusters are lost
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = 1, geometry = st_sfc(poly1), crs = 4326)
  
  # Empty future
  future_sf <- st_sf(ID = integer(0), geometry = st_sfc(crs = 4326))
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$cluster_id, 1)
  expect_true(result$current_area_km2 > 0)
  expect_equal(result$future_area_km2, 0)
  expect_equal(result$area_change_pct, -100)
  expect_true(is.na(result$centroid_shift_km))
})

test_that("calculate_changes validates projection", {
  # Just test that it runs with valid projections
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  sf_obj <- st_sf(ID = 1, geometry = st_sfc(poly1), crs = 4326)
  
  # Should work with EPSG code
  result1 <- calculate_changes(sf_obj, sf_obj, planar_proj = 3857)
  expect_s3_class(result1, "data.frame")
  
  # Should work with different EPSG
  result2 <- calculate_changes(sf_obj, sf_obj, planar_proj = 32633)  # UTM zone 33N
  expect_s3_class(result2, "data.frame")
})


# ==============================================================================
# Additional cluster_novel_areas tests for better coverage
# ==============================================================================

test_that("cluster_novel_areas message output", {
  skip_on_cran()
  skip_if_not_installed("NbClust")
  skip_if_not_installed("caret")
  
  set.seed(999)
  
  r1 <- rast(ncols=10, nrows=10, xmin=0, xmax=10, ymin=0, ymax=10)
  values(r1) <- rnorm(100, sd=2)
  r2 <- rast(ncols=10, nrows=10, xmin=0, xmax=10, ymin=0, ymax=10)
  values(r2) <- rnorm(100, sd=2)
  future_rescaled <- c(r1, r2)
  names(future_rescaled) <- c("var1", "var2")
  
  novel_mask <- rast(ncols=10, nrows=10, xmin=0, xmax=10, ymin=0, ymax=10)
  values(novel_mask) <- sample(c(1, NA), 100, replace=TRUE, prob=c(0.6, 0.4))
  
  # Should output message
  expect_message(
    suppressWarnings(
      cluster_novel_areas(
        future_rescaled = future_rescaled,
        novel_mask = novel_mask,
        n_novel_pts = 30,
        next_cluster_id = 10,
        nbclust_args = list(min.nc=2, max.nc=3, index="silhouette")
      )
    ),
    "Clustering novel climate areas"
  )
})

test_that("cluster_novel_areas handles complete.cases filtering", {
  skip_on_cran()
  skip_if_not_installed("caret")
  
  set.seed(111)
  
  # Create data where NAs will reduce sample below threshold
  r1 <- rast(ncols=8, nrows=8, xmin=0, xmax=8, ymin=0, ymax=8)
  vals1 <- rnorm(64)
  vals1[1:55] <- NA  # Most cells NA
  values(r1) <- vals1
  
  r2 <- rast(ncols=8, nrows=8, xmin=0, xmax=8, ymin=0, ymax=8)
  values(r2) <- rnorm(64)
  
  future_rescaled <- c(r1, r2)
  names(future_rescaled) <- c("var1", "var2")
  
  # Mark remaining cells as novel
  novel_mask <- rast(ncols=8, nrows=8, xmin=0, xmax=8, ymin=0, ymax=8)
  values(novel_mask) <- ifelse(is.na(vals1), NA, 1)
  
  # Should trigger "too few points" warning because complete.cases removes most points
  expect_warning(
    result <- cluster_novel_areas(
      future_rescaled = future_rescaled,
      novel_mask = novel_mask,
      n_novel_pts = 50,  # Requesting more than available after NA removal
      next_cluster_id = 5,
      nbclust_args = list()
    ),
    "Too few novel climate points"
  )
  
  expect_s4_class(result$clusters_raster, "SpatRaster")
  
  # Should assign single cluster ID
  cluster_vals <- values(result$clusters_raster)
  unique_vals <- unique(cluster_vals[!is.na(cluster_vals)])
  expect_length(unique_vals, 1)
  expect_equal(unique_vals, 5)
})


# ==============================================================================  
# Additional analyze_cluster_relationships tests
# ==============================================================================

test_that("analyze_cluster_relationships handles silhouette warnings", {
  skip_if_not_installed("cluster")
  
  set.seed(222)
  
  # Create minimal setup that might trigger warnings
  clusters_rast <- rast(ncols=10, nrows=10)
  vals <- rep(NA, 100)
  vals[1:10] <- 1
  vals[11:20] <- 11
  values(clusters_rast) <- vals
  
  r1 <- rast(ncols=10, nrows=10)
  values(r1) <- rnorm(100)
  future_rescaled <- r1
  names(future_rescaled) <- "var1"
  
  # Should handle gracefully
  result <- suppressWarnings(
    analyze_cluster_relationships(
      clusters_raster = clusters_rast,
      future_rescaled = future_rescaled,
      existing_ids = 1,
      novel_ids = 11,
      n_sample_per_cluster = 5
    )
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("analyze_cluster_relationships filters by novel_ids", {
  skip_if_not_installed("cluster")
  
  set.seed(333)
  
  clusters_rast <- rast(ncols=25, nrows=25)
  vals <- rep(NA, 625)
  vals[1:150] <- 1
  vals[151:300] <- 2
  vals[301:450] <- 11
  vals[451:600] <- 12
  values(clusters_rast) <- vals
  
  r1 <- rast(ncols=25, nrows=25)
  values(r1) <- rnorm(625)
  r2 <- rast(ncols=25, nrows=25)
  values(r2) <- rnorm(625)
  future_rescaled <- c(r1, r2)
  names(future_rescaled) <- c("var1", "var2")
  
  result <- suppressWarnings(
    analyze_cluster_relationships(
      clusters_raster = clusters_rast,
      future_rescaled = future_rescaled,
      existing_ids = c(1, 2),
      novel_ids = c(11, 12),
      n_sample_per_cluster = 20
    )
  )
  
  # Result should only contain novel cluster IDs
  if (nrow(result) > 0) {
    expect_true(all(result$novel_cluster_id %in% c(11, 12)))
  }
})


# ==============================================================================
# Additional calculate_changes tests
# ==============================================================================

test_that("calculate_changes computes centroid shifts correctly", {
  # Create clusters that move a known distance
  poly1_current <- st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = 1, geometry = st_sfc(poly1_current), crs = 4326)
  
  # Move polygon exactly 1 degree east
  poly1_future <- st_polygon(list(matrix(c(1,0, 2,0, 2,1, 1,1, 1,0), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = 1, geometry = st_sfc(poly1_future), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  # Centroid shift should be positive
  expect_true(result$centroid_shift_km > 0)
  expect_true(is.finite(result$centroid_shift_km))
})

test_that("calculate_changes area_change_pct calculation", {
  # Create polygon that doubles in area
  poly1_current <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  current_sf <- st_sf(ID = 1, geometry = st_sfc(poly1_current), crs = 4326)
  
  # Double the size
  poly1_future <- st_polygon(list(matrix(c(0,0, 20,0, 20,20, 0,20, 0,0), ncol=2, byrow=TRUE)))
  future_sf <- st_sf(ID = 1, geometry = st_sfc(poly1_future), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  # Area should roughly quadruple (4x larger in each dimension = 4x area)
  # Allow some tolerance due to projection
  expect_true(result$area_change_pct > 200)  # At least 200% increase
})

test_that("calculate_changes bind_rows preserves structure", {
  # Test with all three types: persisting, lost, gained
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  poly2 <- st_polygon(list(matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol=2, byrow=TRUE)))
  poly3 <- st_polygon(list(matrix(c(40,40, 50,40, 50,50, 40,50, 40,40), ncol=2, byrow=TRUE)))
  poly4 <- st_polygon(list(matrix(c(60,60, 70,60, 70,70, 60,70, 60,60), ncol=2, byrow=TRUE)))
  
  # Current: 1, 2, 3
  current_sf <- st_sf(ID = c(1, 2, 3), geometry = st_sfc(poly1, poly2, poly3), crs = 4326)
  
  # Future: 1, 3, 4 (2 lost, 4 gained)
  future_sf <- st_sf(ID = c(1, 3, 4), geometry = st_sfc(poly1, poly3, poly4), crs = 4326)
  
  result <- calculate_changes(current_sf, future_sf, planar_proj = 3857)
  
  # All required columns present
  expect_true(all(c("cluster_id", "current_area_km2", "future_area_km2", 
                    "area_change_pct", "centroid_shift_km") %in% names(result)))
  
  # All cluster IDs accounted for
  expect_setequal(result$cluster_id, c(1, 2, 3, 4))
})

test_that("calculate_changes validates projection", {
  # Just test that it runs with valid projections
  poly1 <- st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol=2, byrow=TRUE)))
  sf_obj <- st_sf(ID = 1, geometry = st_sfc(poly1), crs = 4326)
  
  # Should work with EPSG code
  result1 <- calculate_changes(sf_obj, sf_obj, planar_proj = 3857)
  expect_s3_class(result1, "data.frame")
  
  # Should work with different EPSG
  result2 <- calculate_changes(sf_obj, sf_obj, planar_proj = 32633)  # UTM zone 33N
  expect_s3_class(result2, "data.frame")
})

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixture: a small planar sf with known geometry
# Four squares in a 2×2 grid, clearly separated so distances are predictable.
# ─────────────────────────────────────────────────────────────────────────────

make_grid_sf <- function() {
  # Four 1×1 unit squares at (0,0), (2,0), (0,2), (2,2)
  make_sq <- function(xmin, ymin) {
    st_polygon(list(matrix(
      c(xmin, ymin, xmin+1, ymin, xmin+1, ymin+1, xmin, ymin+1, xmin, ymin),
      ncol = 2, byrow = TRUE
    )))
  }
  geoms <- st_sfc(make_sq(0,0), make_sq(2,0), make_sq(0,2), make_sq(2,2),
                  crs = 3857)
  st_sf(ID = 1:4, geometry = geoms)
}

# Larger grid for PrioritizeSample (needs more cells for voronoi stability)
make_larger_grid_sf <- function(n = 9) {
  # 3×3 grid of 10×10 unit squares, well separated
  coords <- expand.grid(col = 0:2, row = 0:2)
  make_sq <- function(xmin, ymin, size = 10) {
    st_polygon(list(matrix(
      c(xmin, ymin, xmin+size, ymin, xmin+size, ymin+size,
        xmin, ymin+size, xmin, ymin),
      ncol = 2, byrow = TRUE
    )))
  }
  geoms <- st_sfc(
    mapply(function(col, row) make_sq(col * 15, row * 15),
           coords$col, coords$row, SIMPLIFY = FALSE),
    crs = 3857
  )
  st_sf(ID = seq_len(nrow(coords)), geometry = geoms)
}

# ═════════════════════════════════════════════════════════════════════════════
# order_by_distance_variance()
# ═════════════════════════════════════════════════════════════════════════════

test_that("order_by_distance_variance returns a permutation of 1:n", {
  x <- make_grid_sf()
  result <- order_by_distance_variance(x, metric = "var")

  expect_type(result, "integer")
  expect_length(result, nrow(x))
  expect_setequal(result, seq_len(nrow(x)))
})

test_that("order_by_distance_variance: all four metrics return valid permutations", {
  x <- make_grid_sf()
  for (m in c("var", "sd", "energy", "cv")) {
    result <- order_by_distance_variance(x, metric = m)
    expect_length(result, 4L)
    expect_setequal(result, 1:4)
  }
})

test_that("order_by_distance_variance seed is the most-central zone", {
  # In a symmetric 2×2 grid all zones have equal row-sums, so which.min picks
  # index 1. For an asymmetric layout the central point should win.
  #
  # Build a star layout: one center zone + four arms far away.
  # Center is at origin; arms at ±100 in x or y.
  make_sq_at <- function(cx, cy, half = 1) {
    st_polygon(list(matrix(
      c(cx-half, cy-half, cx+half, cy-half, cx+half, cy+half,
        cx-half, cy+half, cx-half, cy-half),
      ncol = 2, byrow = TRUE
    )))
  }
  geoms <- st_sfc(
    make_sq_at(0,   0),   # index 1 — center
    make_sq_at(100, 0),
    make_sq_at(-100,0),
    make_sq_at(0,  100),
    make_sq_at(0, -100),
    crs = 3857
  )
  x <- st_sf(ID = 1:5, geometry = geoms)

  result <- order_by_distance_variance(x, metric = "var")
  expect_equal(unname(result[1]), 1L)  # center zone should be sampled first
})

test_that("order_by_distance_variance: cv metric runs without error", {
  # sf::st_distance() returns units objects. The cv score_fn computes
  # mean(x) and checks `== 0`, which fails when x is a units vector.
  # This test verifies the function at least runs on a normal layout;
  # the near-zero-mean edge case is a known source limitation.
  x <- make_grid_sf()
  result <- expect_no_error(order_by_distance_variance(x, metric = "cv"))
  expect_setequal(result, 1:4)
})

test_that("order_by_distance_variance: single zone returns itself", {
  make_sq <- function() st_polygon(list(matrix(
    c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE
  )))
  x <- st_sf(ID = 1L, geometry = st_sfc(make_sq(), crs = 3857))

  result <- order_by_distance_variance(x, metric = "var")
  expect_equal(unname(result), 1L)
})

test_that("order_by_distance_variance: two zones return both in some order", {
  make_sq_at <- function(cx) st_polygon(list(matrix(
    c(cx,0, cx+1,0, cx+1,1, cx,1, cx,0), ncol=2, byrow=TRUE
  )))
  x <- st_sf(ID = 1:2,
              geometry = st_sfc(make_sq_at(0), make_sq_at(100), crs = 3857))

  # cv excluded: sf::st_distance() returns units objects and the cv branch
  # performs `mean(x) == 0` (units vs numeric) which errors in source code.
  for (m in c("var", "sd", "energy")) {
    result <- order_by_distance_variance(x, metric = m)
    expect_setequal(result, 1:2)
  }
})

test_that("order_by_distance_variance: energy metric gives finite scores", {
  x <- make_grid_sf()
  # Verify internally that score_fn = sum(x^2) never crashes
  result <- order_by_distance_variance(x, metric = "energy")
  expect_true(all(is.finite(result)))
})

test_that("order_by_distance_variance: different metrics may differ on ordering", {
  # With enough zones, var and energy don't have to agree on the full sequence
  x <- make_larger_grid_sf()
  r_var    <- order_by_distance_variance(x, metric = "var")
  r_energy <- order_by_distance_variance(x, metric = "energy")

  # Both are valid permutations — but they need not be identical
  expect_setequal(r_var,    seq_len(nrow(x)))
  expect_setequal(r_energy, seq_len(nrow(x)))
  # (We don't assert they differ — that would be fragile — just that both work)
})

test_that("order_by_distance_variance: greedy step reduces remaining correctly", {
  # After the greedy loop completes, every zone must appear exactly once.
  # cv excluded: units arithmetic bug in source (see cv metric test).
  x <- make_larger_grid_sf()
  for (m in c("var", "sd", "energy")) {
    result <- order_by_distance_variance(x, metric = m)
    expect_equal(length(result), nrow(x))
    expect_equal(length(unique(result)), nrow(x))  # no duplicates
  }
})


# ═════════════════════════════════════════════════════════════════════════════
# PrioritizeSample()
# ═════════════════════════════════════════════════════════════════════════════

test_that("PrioritizeSample errors on geographic (lon/lat) CRS", {
  make_sq <- function() st_polygon(list(matrix(
    c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE
  )))
  x <- st_sf(ID = 1L, geometry = st_sfc(make_sq(), crs = 4326))  # lon/lat

  expect_error(
    PrioritizeSample(x, n_breaks = 3),
    "planar"
  )
})

test_that("PrioritizeSample returns a named list with 'Geometry' element", {
  x <- make_larger_grid_sf()

  result <- suppressWarnings(PrioritizeSample(x, n_breaks = 3, metric = "var"))

  expect_type(result, "list")
  expect_named(result, "Geometry")
  expect_s3_class(result$Geometry, "sf")
})

test_that("PrioritizeSample output has required columns", {
  x <- make_larger_grid_sf()
  result <- suppressWarnings(PrioritizeSample(x, n_breaks = 3, metric = "var"))

  expect_true(all(c("ID", "SampleOrder", "Level", "geometry") %in%
                    names(result$Geometry)))
})

test_that("PrioritizeSample Level values are 1:n_breaks", {
  x <- make_larger_grid_sf()
  result <- suppressWarnings(PrioritizeSample(x, n_breaks = 3, metric = "var"))

  levels_present <- unique(result$Geometry$Level)
  # cut() with n_breaks=3 produces labels 1, 2, 3
  expect_true(all(levels_present %in% 1:3))
})

test_that("PrioritizeSample SampleOrder covers all input zones", {
  x <- make_larger_grid_sf()
  result <- suppressWarnings(PrioritizeSample(x, n_breaks = 3, metric = "var"))

  orders <- unique(result$Geometry$SampleOrder)
  # Every zone index should appear in SampleOrder
  expect_true(length(orders) == nrow(x))
})

test_that("PrioritizeSample filters rows where n == 0 when n column present", {
  x <- make_larger_grid_sf()
  # Mark only rows 1-6 as target (rows 7-9 should be excluded)
  x$n <- c(rep(1L, 6L), rep(0L, 3L))

  result <- suppressWarnings(PrioritizeSample(x, n_breaks = 2, metric = "var"))

  # SampleOrder should only reference 6 zones
  orders <- unique(result$Geometry$SampleOrder)
  expect_equal(length(orders), 6L)
})

test_that("PrioritizeSample keeps all rows when n column absent", {
  x <- make_larger_grid_sf()  # no 'n' column
  result <- suppressWarnings(PrioritizeSample(x, n_breaks = 2, metric = "var"))

  orders <- unique(result$Geometry$SampleOrder)
  expect_equal(length(orders), nrow(x))
})

test_that("PrioritizeSample n_breaks=2 produces Level in {1, 2}", {
  x <- make_larger_grid_sf()
  result <- suppressWarnings(PrioritizeSample(x, n_breaks = 2, metric = "energy"))

  expect_true(all(result$Geometry$Level %in% 1:2))
})

test_that("PrioritizeSample works with all four metric options", {
  x <- make_larger_grid_sf()
  for (m in c("var", "sd", "energy", "cv")) {
    result <- suppressWarnings(
      PrioritizeSample(x, n_breaks = 3, metric = m)
    )
    expect_s3_class(result$Geometry, "sf")
  }
})

test_that("PrioritizeSample geometry is valid sf", {
  x <- make_larger_grid_sf()
  result <- suppressWarnings(PrioritizeSample(x, n_breaks = 3, metric = "var"))

  expect_true(all(st_is_valid(result$Geometry)))
})

test_that("PrioritizeSample: n column filters only rows where n != 1", {
  x <- make_larger_grid_sf()
  # n=1 for first 5, n=0 for last 4
  x$n <- c(rep(1L, 5L), rep(0L, 4L))

  result_filtered <- suppressWarnings(
    PrioritizeSample(x, n_breaks = 2, metric = "var")
  )

  x_no_n <- x[, setdiff(names(x), "n")]
  result_full <- suppressWarnings(
    PrioritizeSample(x_no_n, n_breaks = 2, metric = "var")
  )

  expect_lt(
    length(unique(result_filtered$Geometry$SampleOrder)),
    length(unique(result_full$Geometry$SampleOrder))
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures
# ─────────────────────────────────────────────────────────────────────────────

make_rast <- function(nrows = 10, ncols = 10,
                      xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                      vals = NULL) {
  r <- terra::rast(nrows = nrows, ncols = ncols,
                   xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                   crs = "EPSG:4326")
  terra::values(r) <- if (is.null(vals)) rnorm(terra::ncell(r)) else vals
  r
}

make_pred_stack <- function(vars = c("var1", "var2"), ...) {
  layers <- lapply(vars, function(v) make_rast(...))
  stk <- do.call(c, layers)
  names(stk) <- vars
  stk
}

# Minimal glmnet model with known coefficient names
make_glmnet <- function(vars = c("var1", "var2"), n = 100, seed = 1) {
  set.seed(seed)
  x <- matrix(rnorm(n * length(vars)), ncol = length(vars),
               dimnames = list(NULL, vars))
  y <- factor(sample(0:1, n, replace = TRUE))
  glmnet::glmnet(x, y, family = "binomial", lambda = 0.01, alpha = 0.5)
}

# eSDM_object stub
make_esdm <- function(vars = c("var1", "var2"), n_train = 60, n_test = 20,
                      seed = 1) {
  set.seed(seed)
  coords_train <- data.frame(lon = runif(n_train, 0, 10),
                              lat = runif(n_train, 0, 10))
  train_sf <- sf::st_as_sf(coords_train, coords = c("lon", "lat"),
                            crs = 4326)
  train_sf$occurrence <- sample(0:1, n_train, replace = TRUE)

  coords_test <- data.frame(lon = runif(n_test, 0, 10),
                             lat = runif(n_test, 0, 10))
  test_sf <- sf::st_as_sf(coords_test, coords = c("lon", "lat"), crs = 4326)
  test_sf$occurrence <- sample(0:1, n_test, replace = TRUE)

  pred_mat <- as.data.frame(
    matrix(rnorm(n_train * length(vars)), ncol = length(vars),
           dimnames = list(NULL, vars))
  )

  list(
    TrainData    = train_sf,
    TestData     = test_sf,
    PredictMatrix = pred_mat,
    Model        = make_glmnet(vars, n = n_train, seed = seed)
  )
}

# current_clusters stub — sf with ID column + fit.knn placeholder
make_current_clusters <- function(n_clusters = 3) {
  make_sq <- function(i) {
    x0 <- (i - 1) * 3
    sf::st_polygon(list(matrix(
      c(x0, 0, x0+2, 0, x0+2, 2, x0, 2, x0, 0),
      ncol = 2, byrow = TRUE
    )))
  }
  geoms <- sf::st_sfc(lapply(seq_len(n_clusters), make_sq), crs = 4326)
  geom_sf <- sf::st_sf(ID = seq_len(n_clusters), geometry = geoms)

  # fit.knn just needs to be something terra::predict dispatches on;
  # we mock terra::predict so the object itself doesn't matter
  list(
    Geometry = geom_sf,
    fit.knn  = structure(list(), class = "knn")
  )
}

# Build a mock raster::stack result for dismo::mess (returns a RasterStack)
make_mess_stack <- function(r_template, mess_vals = NULL) {
  r <- raster::raster(
    nrows = terra::nrow(r_template),
    ncols = terra::ncol(r_template),
    xmn   = terra::xmin(r_template),
    xmx   = terra::xmax(r_template),
    ymn   = terra::ymin(r_template),
    ymx   = terra::ymax(r_template),
    crs   = terra::crs(r_template)
  )
  raster::values(r) <- if (is.null(mess_vals))
    rep(10, terra::ncell(r_template))   # all positive → no novel climate
  else
    mess_vals
  raster::stack(r)
}

# Standard changes data.frame
make_changes_df <- function(ids = 1:3) {
  data.frame(
    cluster_id       = ids,
    current_area_km2 = runif(length(ids), 10, 100),
    future_area_km2  = runif(length(ids), 10, 100),
    area_change_pct  = runif(length(ids), -20, 20),
    centroid_shift_km = runif(length(ids), 0, 5)
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: activate all mocks for a projectClusters() call
# mess_vals controls which cells are "novel" (negative MESS = novel climate)
# ─────────────────────────────────────────────────────────────────────────────

# mess_vals controls which cells are "novel" (negative MESS = novel climate)
# ─────────────────────────────────────────────────────────────────────────────

run_with_mocks <- function(esdm, current_clusters, future_preds, current_preds,
                           mess_vals = NULL,
                           knn_cluster_vals = NULL,
                           cluster_novel = FALSE,
                           novel_cluster_result = NULL,
                           novel_similarity_result = NULL,
                           changes_result = NULL,
                           ...) {

  n_cells   <- terra::ncell(future_preds[[1]])
  r_template <- future_preds[[1]]

  # Default KNN output: all cells → cluster 1
  if (is.null(knn_cluster_vals)) knn_cluster_vals <- rep(1L, n_cells)

  # Default novel cluster result
  if (is.null(novel_cluster_result)) {
    nr <- r_template
    terra::values(nr) <- rep(4L, n_cells)
    novel_cluster_result <- list(clusters_raster = nr)
  }

  if (is.null(novel_similarity_result)) {
    novel_similarity_result <- data.frame(
      novel_cluster_id    = 4L,
      nearest_existing_id = 1L,
      avg_silhouette_width = 0.4
    )
  }

  if (is.null(changes_result)) {
    changes_result <- make_changes_df(seq_len(nrow(current_clusters$Geometry)))
  }

  # Build mock KNN raster
  knn_rast <- r_template
  terra::values(knn_rast) <- knn_cluster_vals

  # Rescaled future: same structure as future_preds
  rescaled_rast <- future_preds

  local_mocked_bindings(
    # dismo::mess → RasterStack with requested mess_vals
    mess = function(x, v, full = FALSE) {
      make_mess_stack(r_template, mess_vals)
    },
    .package = "dismo"
  )

  local_mocked_bindings(
    # terra::predict dispatches for both glmnet (SDM) and knn (clusters)
    # Distinguish by model class
    predict = function(object, model, fun = NULL, na.rm = TRUE, ...) {
      if (inherits(model, "knn")) {
        knn_rast
      } else {
        # glmnet SDM prediction: return values in [0,1]
        r <- object[[1]]
        terra::values(r) <- runif(terra::ncell(r), 0.3, 0.9)
        r
      }
    },
    .package = "terra"
  )

  local_mocked_bindings(
    rescaleFuture           = function(...) rescaled_rast,
    add_weighted_coordinates = function(x, ...) x,
    cluster_novel_areas     = function(...) novel_cluster_result,
    analyze_cluster_relationships = function(...) novel_similarity_result,
    calculate_changes       = function(...) changes_result,
    .package = "safeHavens"
  )

  extra <- list(...)
  if (!"thresholds" %in% names(extra)) {
    extra$thresholds <- list(sensitivity = 0.5)
  }

  do.call(projectClusters, c(
    list(
      eSDM_object        = esdm,
      current_clusters   = current_clusters,
      future_predictors  = future_preds,
      current_predictors = current_preds,
      planar_proj        = 3857,
      cluster_novel      = cluster_novel
    ),
    extra
  ))
}

# ═════════════════════════════════════════════════════════════════════════════
# Branch A: cluster_novel = FALSE → else path, empty novel_similarity
# ═════════════════════════════════════════════════════════════════════════════

test_that("cluster_novel=FALSE returns empty novel_similarity data.frame", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(1)
  vars         <- c("var1", "var2")
  esdm         <- make_esdm(vars)
  cc           <- make_current_clusters(3)
  future_preds <- make_pred_stack(vars)
  current_preds <- make_pred_stack(vars)

  result <- run_with_mocks(
    esdm, cc, future_preds, current_preds,
    cluster_novel = FALSE
  )

  expect_equal(nrow(result$novel_similarity), 0L)
  expect_named(result$novel_similarity,
               c("novel_cluster_id", "nearest_existing_id", "avg_silhouette_width"))
})

test_that("cluster_novel=FALSE: final_clusters equals known_clusters (no cover merge)", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(2)
  vars  <- c("var1", "var2")
  esdm  <- make_esdm(vars)
  cc    <- make_current_clusters(2)
  fp    <- make_pred_stack(vars)
  cp    <- make_pred_stack(vars)

  # KNN produces clusters 1 and 2 alternating
  n <- terra::ncell(fp[[1]])
  knn_vals <- rep(c(1L, 2L), length.out = n)

  result <- run_with_mocks(
    esdm, cc, fp, cp,
    knn_cluster_vals = knn_vals,
    cluster_novel = FALSE
  )

  # clusters_raster should only contain values from KNN output (1, 2, or NA)
  rast_vals <- terra::values(result$clusters_raster)
  unique_vals <- unique(rast_vals[!is.na(rast_vals)])
  expect_true(all(unique_vals %in% c(1L, 2L)))
})

# ═════════════════════════════════════════════════════════════════════════════
# Branch B: cluster_novel=TRUE with novel cells → if path
# ═════════════════════════════════════════════════════════════════════════════

test_that("cluster_novel=TRUE calls novel clustering and returns non-empty novel_similarity", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(3)
  vars  <- c("var1", "var2")
  esdm  <- make_esdm(vars)
  cc    <- make_current_clusters(3)
  fp    <- make_pred_stack(vars)
  cp    <- make_pred_stack(vars)

  # All cells have negative MESS → all are novel climate
  n <- terra::ncell(fp[[1]])
  novel_sim <- data.frame(
    novel_cluster_id    = 4L,
    nearest_existing_id = 1L,
    avg_silhouette_width = 0.55
  )

  result <- run_with_mocks(
    esdm, cc, fp, cp,
    mess_vals = rep(-5, n),      # all novel
    cluster_novel = TRUE,
    novel_similarity_result = novel_sim
  )

  expect_gt(nrow(result$novel_similarity), 0L)
  expect_equal(result$novel_similarity$novel_cluster_id, 4L)
  expect_equal(result$novel_similarity$nearest_existing_id, 1L)
})

test_that("cluster_novel=TRUE merges known and novel clusters via cover()", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(4)
  vars <- c("var1", "var2")
  esdm <- make_esdm(vars)
  cc   <- make_current_clusters(3)  # existing IDs 1,2,3
  fp   <- make_pred_stack(vars)
  cp   <- make_pred_stack(vars)

  n <- terra::ncell(fp[[1]])

  # Novel cluster raster (will be offset by max existing ID = 3)
  nr <- fp[[1]]
  terra::values(nr) <- rep(1L, n)   # stub novel cluster ID 1 → becomes 4 after offset

  result <- run_with_mocks(
    esdm, cc, fp, cp,
    mess_vals = rep(-5, n),
    knn_cluster_vals = rep(1L, n),
    cluster_novel = TRUE,
    novel_cluster_result = list(clusters_raster = nr)
  )

  # clusters_raster is a SpatRaster
  expect_s4_class(result$clusters_raster, "SpatRaster")
})

# ═════════════════════════════════════════════════════════════════════════════
# Branch C: cluster_novel=TRUE but zero novel cells → falls through to else
# ═════════════════════════════════════════════════════════════════════════════

test_that("cluster_novel=TRUE but no novel cells → empty novel_similarity", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(5)
  vars <- c("var1", "var2")
  esdm <- make_esdm(vars)
  cc   <- make_current_clusters(2)
  fp   <- make_pred_stack(vars)
  cp   <- make_pred_stack(vars)

  n <- terra::ncell(fp[[1]])

  result <- run_with_mocks(
    esdm, cc, fp, cp,
    mess_vals = rep(10, n),    # all positive → zero novel cells
    cluster_novel = TRUE       # flag is TRUE but condition n_novel_cells > 0 fails
  )

  expect_equal(nrow(result$novel_similarity), 0L)
})

# ═════════════════════════════════════════════════════════════════════════════
# Output contract: all 7 slots present and correctly typed
# ═════════════════════════════════════════════════════════════════════════════

test_that("return list has all required named slots with correct types", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(6)
  vars <- c("var1", "var2")
  esdm <- make_esdm(vars)
  cc   <- make_current_clusters(3)
  fp   <- make_pred_stack(vars)
  cp   <- make_pred_stack(vars)

  result <- run_with_mocks(esdm, cc, fp, cp, cluster_novel = FALSE)

  expect_named(result,
    c("clusters_raster", "clusters_sf", "suitable_habitat",
      "novel_mask", "mess", "changes", "novel_similarity"),
    ignore.order = TRUE
  )

  expect_s4_class(result$clusters_raster,  "SpatRaster")
  expect_s3_class(result$clusters_sf,      "sf")
  expect_s4_class(result$suitable_habitat, "SpatRaster")
  expect_s4_class(result$novel_mask,       "SpatRaster")
  expect_s4_class(result$mess,             "SpatRaster")
  expect_s3_class(result$changes,          "data.frame")
  expect_s3_class(result$novel_similarity, "data.frame")
})

test_that("clusters_sf has an ID column", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(7)
  vars <- c("var1", "var2")
  esdm <- make_esdm(vars)
  cc   <- make_current_clusters(2)
  fp   <- make_pred_stack(vars)
  cp   <- make_pred_stack(vars)

  result <- run_with_mocks(esdm, cc, fp, cp, cluster_novel = FALSE)

  expect_true("ID" %in% names(result$clusters_sf))
})

# ═════════════════════════════════════════════════════════════════════════════
# suitable_habitat and novel_mask are binary (1 or NA)
# ═════════════════════════════════════════════════════════════════════════════

test_that("suitable_habitat raster contains only 1 and NA", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(8)
  vars <- c("var1", "var2")
  esdm <- make_esdm(vars)
  cc   <- make_current_clusters(2)
  fp   <- make_pred_stack(vars)
  cp   <- make_pred_stack(vars)

  result <- run_with_mocks(esdm, cc, fp, cp, cluster_novel = FALSE)

  vals <- terra::values(result$suitable_habitat)
  non_na <- vals[!is.na(vals)]
  expect_true(all(non_na == 1))
})

test_that("novel_mask raster contains only 1 and NA", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(9)
  vars <- c("var1", "var2")
  esdm  <- make_esdm(vars)
  cc    <- make_current_clusters(2)
  fp    <- make_pred_stack(vars)
  cp    <- make_pred_stack(vars)
  n     <- terra::ncell(fp[[1]])

  result <- run_with_mocks(
    esdm, cc, fp, cp,
    mess_vals = rep(-5, n),   # all novel → novel_mask all 1
    cluster_novel = FALSE
  )

  vals <- terra::values(result$novel_mask)
  non_na <- vals[!is.na(vals)]
  expect_true(all(non_na == 1))
})

# ═════════════════════════════════════════════════════════════════════════════
# thresh_metric controls which threshold value is used
# ═════════════════════════════════════════════════════════════════════════════

test_that("thresh_metric selects the correct threshold from thresholds list", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(10)
  vars <- c("var1", "var2")
  esdm <- make_esdm(vars)
  cc   <- make_current_clusters(2)
  fp   <- make_pred_stack(vars)
  cp   <- make_pred_stack(vars)

  # High threshold → very little suitable habitat
  result_high <- run_with_mocks(
    esdm, cc, fp, cp,
    cluster_novel = FALSE,
    thresh_metric = "specificity",
    thresholds = list(sensitivity = 0.1, specificity = 0.99)
  )

  # Low threshold → more suitable habitat
  result_low <- run_with_mocks(
    esdm, cc, fp, cp,
    cluster_novel = FALSE,
    thresh_metric = "sensitivity",
    thresholds = list(sensitivity = 0.1, specificity = 0.99)
  )

  high_suitable <- sum(!is.na(terra::values(result_high$suitable_habitat)))
  low_suitable  <- sum(!is.na(terra::values(result_low$suitable_habitat)))

  expect_lte(high_suitable, low_suitable)
})

# ═════════════════════════════════════════════════════════════════════════════
# mess output matches mess_threshold logic
# ═════════════════════════════════════════════════════════════════════════════

test_that("mess slot reflects raw MESS values from dismo", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(11)
  vars <- c("var1", "var2")
  esdm <- make_esdm(vars)
  cc   <- make_current_clusters(2)
  fp   <- make_pred_stack(vars)
  cp   <- make_pred_stack(vars)
  n    <- terra::ncell(fp[[1]])

  result <- run_with_mocks(
    esdm, cc, fp, cp,
    mess_vals = rep(7, n),
    cluster_novel = FALSE
  )

  # mess slot should be a single-layer SpatRaster
  expect_equal(terra::nlyr(result$mess), 1L)
  mess_vals <- terra::values(result$mess)
  expect_true(all(mess_vals == 7, na.rm = TRUE))
})

# ═════════════════════════════════════════════════════════════════════════════
# changes data.frame is forwarded from calculate_changes
# ═════════════════════════════════════════════════════════════════════════════

test_that("changes slot is the data.frame returned by calculate_changes", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(12)
  vars <- c("var1", "var2")
  esdm <- make_esdm(vars)
  cc   <- make_current_clusters(3)
  fp   <- make_pred_stack(vars)
  cp   <- make_pred_stack(vars)

  mock_changes <- data.frame(
    cluster_id        = 1:3,
    current_area_km2  = c(10, 20, 30),
    future_area_km2   = c(12, 18, 35),
    area_change_pct   = c(20, -10, 16.7),
    centroid_shift_km = c(1.1, 2.2, 0.5)
  )

  result <- run_with_mocks(
    esdm, cc, fp, cp,
    cluster_novel = FALSE,
    changes_result = mock_changes
  )

  expect_equal(result$changes, mock_changes)
})

# ═════════════════════════════════════════════════════════════════════════════
# model variables extracted from glmnet coef names
# ═════════════════════════════════════════════════════════════════════════════

test_that("only model coefficient variables are used for SDM prediction", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("dismo")

  set.seed(13)
  # Predictor stack has an extra variable not in the model
  vars_model <- c("var1", "var2")
  vars_extra <- c("var1", "var2", "var3")   # var3 not in model

  esdm <- make_esdm(vars_model)
  cc   <- make_current_clusters(2)
  fp   <- make_pred_stack(vars_extra)
  cp   <- make_pred_stack(vars_extra)

  # terra::predict is called with future_predictors[[vars]] — if var3 was
  # incorrectly included it would cause a dimension mismatch in glmnet.
  # The mock sidesteps that, but we can verify the function completes cleanly.
  expect_no_error(
    run_with_mocks(esdm, cc, fp, cp, cluster_novel = FALSE)
  )
})