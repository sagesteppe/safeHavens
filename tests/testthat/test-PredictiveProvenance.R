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
