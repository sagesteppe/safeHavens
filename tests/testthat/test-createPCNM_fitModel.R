# Test suite for createPCNM_fitModel helper functions
# Uses dismo package example data (Bradypus occurrence data)

library(testthat)
library(terra)
library(sf)
library(dplyr)

# Setup: Load example data used across tests
setup_pcnm_test_data <- function() {
  # Load Bradypus occurrence data
  x <- read.csv(file.path(system.file(package = "dismo"), 'ex', 'bradypus.csv'))
  x <- x[, c('lon', 'lat')]
  x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)
  
  # Load predictor rasters
  files <- list.files(
    path = file.path(system.file(package = "dismo"), 'ex'), 
    pattern = 'grd', 
    full.names = TRUE
  )
  predictors <- terra::rast(files)
  
  # Add occurrence column and extract predictors
  set.seed(1)
  x$occurrence <- factor(
    sample(c(0, 1), nrow(x), replace = TRUE, prob = c(0.4, 0.6))
  )

  x <- terra::extract(predictors, x, bind = TRUE) |>
    sf::st_as_sf()
  
  # Take subset for faster testing
 # x_subset <- x[1:30, ]
  
  # Planar projection for distance calculations
  planar_proj <- 3857  # Web Mercator
  
  list(
    occurrences = x,
    predictors = predictors,
    planar_proj = planar_proj
  )
}

# ==============================================================================
# Tests for calculate_distance_matrix()
# ==============================================================================

test_that("calculate_distance_matrix returns numeric matrix", {
  data <- setup_pcnm_test_data()
  
  result <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  # Should be numeric matrix
  expect_true(is.matrix(result))
  expect_type(result, "double")
})

test_that("calculate_distance_matrix returns square matrix", {
  data <- setup_pcnm_test_data()
  
  result <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  # Should be square matrix (n x n)
  n_points <- nrow(data$occurrences)
  expect_equal(nrow(result), n_points)
  expect_equal(ncol(result), n_points)
})

test_that("calculate_distance_matrix has zero diagonal", {
  data <- setup_pcnm_test_data()
  
  result <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  # Diagonal should be zero (distance from point to itself)
  diagonal <- diag(result)
  expect_true(all(diagonal == 0))
})

test_that("calculate_distance_matrix is symmetric", {
  data <- setup_pcnm_test_data()
  
  result <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  # Distance matrix should be symmetric (check values, ignore dimnames)
  expect_equal(as.numeric(result), as.numeric(t(result)))
})

test_that("calculate_distance_matrix contains positive distances", {
  data <- setup_pcnm_test_data()
  
  result <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  # All distances should be non-negative
  expect_true(all(result >= 0))
  
  # Off-diagonal elements should be positive
  # (assuming points are not duplicates)
  off_diag <- result[upper.tri(result)]
  expect_true(all(off_diag > 0))
})

test_that("calculate_distance_matrix respects planar projection", {
  data <- setup_pcnm_test_data()
  
  # Use different projections
  result_3857 <- calculate_distance_matrix(data$occurrences, 3857)
  result_5070 <- calculate_distance_matrix(data$occurrences, 5070)
  
  # Results should differ due to different projections
  # but same structure
  expect_equal(dim(result_3857), dim(result_5070))
  expect_false(all(result_3857 == result_5070))
})

# ==============================================================================
# Tests for create_pcnm_vectors()
# ==============================================================================

test_that("create_pcnm_vectors returns data frame with correct dimensions", {
  skip_if_not_installed("vegan")
  
  data <- setup_pcnm_test_data()
  dis <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  nv = 10

  result <- create_pcnm_vectors(dis, n_vectors = nv)
  
  # Should be a data frame
  expect_s3_class(result, "data.frame")
  
  # Should have at most n_vectors columns (vegan may return fewer)
  expect_true(ncol(result) <= nv)
  expect_true(ncol(result) > 0)
  
  # Should have same number of rows as points
  expect_equal(nrow(result), nrow(data$occurrences))
})

test_that("create_pcnm_vectors respects n_vectors parameter", {
  skip_if_not_installed("vegan")
  
  data <- setup_pcnm_test_data()
  dis <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  for (n in c(3, 6, 9)) {
    result <- create_pcnm_vectors(dis, n_vectors = n)
    expect_equal(ncol(result), n)
  }
})

test_that("create_pcnm_vectors produces orthogonal eigenvectors", {
  skip_if_not_installed("vegan")
  
  data <- setup_pcnm_test_data()
  dis <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  result <- create_pcnm_vectors(dis, n_vectors = 5)
  
  # Eigenvectors should be approximately orthogonal
  # Calculate correlation matrix
  cor_mat <- cor(result)
  
  # Off-diagonal elements should be near zero
  diag(cor_mat) <- 0
  expect_true(max(abs(cor_mat)) < 0.1)
})

test_that("create_pcnm_vectors produces numeric values", {
  skip_if_not_installed("vegan")
  
  data <- setup_pcnm_test_data()
  dis <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  result <- create_pcnm_vectors(dis, n_vectors = 2)

  # All values should be numeric
  expect_true(all(sapply(result, is.numeric)))
  
  # No missing values
  expect_true(all(!is.na(result)))
})

# ==============================================================================
# Tests for select_pcnm_features()
# ==============================================================================

test_that("select_pcnm_features returns character vector", {
  skip_if_not_installed("vegan")
  skip_if_not_installed("caret")
  #skip("Feature selection is slow - run manually")
  
  data <- setup_pcnm_test_data()
  dis <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  pcnm_df <- create_pcnm_vectors(dis, n_vectors = 20)
  
  # Create simple CV structure
  ctrl <- caret::rfeControl(
    method = "cv",
    number = 3,
    functions = caret::lrFuncs,
    verbose = FALSE
  )
  
  # Create simple CV indices
  cv_indices <- list(
    indx_train = caret::createFolds(data$occurrences$occurrence, k = 3)
  )
  
  result <- select_pcnm_features(
    pcnm_df,
    data$occurrences$occurrence,
    cv_indices
  )
  
  # Should select between 1 and 5 features
  expect_true(length(result) >= 1)
  expect_true(length(result) <= 5)
  
  # All selected features should be in original data
  expect_true(all(result %in% names(pcnm_df)))
})

# ==============================================================================
# Tests for combine_predictors()
# ==============================================================================

test_that("combine_predictors merges environmental and spatial data", {
  # Create mock data
  env_preds <- data.frame(
    bio1 = rnorm(30),
    bio2 = rnorm(30),
    bio3 = rnorm(30)
  )
  
  pcnm_df <- data.frame(
    PCNM1 = rnorm(30),
    PCNM2 = rnorm(30),
    PCNM3 = rnorm(30)
  )
  
  selected <- c("PCNM1", "PCNM3")
  
  result <- combine_predictors(env_preds, pcnm_df, selected)
  
  # Should be data frame
  expect_s3_class(result, "data.frame")
  
  # Should have env + selected PCNM columns
  expect_equal(ncol(result), ncol(env_preds) + length(selected))
  
  # Should have all environmental columns
  expect_true(all(names(env_preds) %in% names(result)))
  
  # Should have selected PCNM columns
  expect_true(all(selected %in% names(result)))
})

test_that("combine_predictors handles single PCNM variable", {
  env_preds <- data.frame(
    bio1 = rnorm(30),
    bio2 = rnorm(30)
  )
  
  pcnm_df <- data.frame(
    PCNM1 = rnorm(30),
    PCNM2 = rnorm(30)
  )
  
  selected <- "PCNM1"  # Single variable
  
  result <- combine_predictors(env_preds, pcnm_df, selected)
  
  # Should have correct dimensions
  expect_equal(ncol(result), 3)  # 2 env + 1 PCNM
  
  # Last column should be named correctly
  expect_equal(names(result)[3], "PCNM1")
})

test_that("combine_predictors preserves row count", {
  n_rows <- 25
  
  env_preds <- data.frame(
    bio1 = rnorm(n_rows),
    bio2 = rnorm(n_rows)
  )
  
  pcnm_df <- data.frame(
    PCNM1 = rnorm(n_rows),
    PCNM2 = rnorm(n_rows)
  )
  
  result <- combine_predictors(env_preds, pcnm_df, c("PCNM1", "PCNM2"))
  
  expect_equal(nrow(result), n_rows)
})

# ==============================================================================
# Tests for fit_combined_model()
# ==============================================================================

test_that("fit_combined_model returns list with required components", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")
  #skip("Model fitting is slow - run manually")
  
  # Create mock data
  combined_preds <- data.frame(
    bio1 = rnorm(50),
    bio2 = rnorm(50),
    PCNM1 = rnorm(50)
  )
  
  occurrence <- factor(sample(0:1, 50, replace = TRUE))
  
  cv_indices <- list(
    indx_train = caret::createFolds(occurrence, k = 3)
  )
  
  result <- fit_combined_model(combined_preds, occurrence, cv_indices)
  
  # Should return list
  expect_type(result, "list")
  expect_named(result, c("cv_model", "glmnet_model"))
  
  # Check model classes
  expect_s3_class(result$cv_model, "train")
  expect_s3_class(result$glmnet_model, "glmnet")
})

# ==============================================================================
# Tests for interpolate_pcnm_to_raster()
# ==============================================================================

test_that("interpolate_pcnm_to_raster returns SpatRaster", {
  skip_if_not_installed("fields")
  
  data <- setup_pcnm_test_data()
  
  # Create simple PCNM vector
  pcnm_vector <- rnorm(nrow(data$occurrences))
  coords <- sf::st_coordinates(data$occurrences)
  template <- data$predictors[[1]]
  
  result <- interpolate_pcnm_to_raster(pcnm_vector, coords, template)
  
  # Should return SpatRaster
  expect_s4_class(result, "SpatRaster")
})

test_that("interpolate_pcnm_to_raster matches template extent", {
  skip_if_not_installed("fields")
  
  data <- setup_pcnm_test_data()
  
  pcnm_vector <- rnorm(nrow(data$occurrences))
  coords <- sf::st_coordinates(data$occurrences)
  template <- data$predictors[[1]]
  
  result <- interpolate_pcnm_to_raster(pcnm_vector, coords, template)
  
  # Extent should match template
  expect_equal(terra::ext(result), terra::ext(template))
  
  # Resolution should match template
  expect_equal(terra::res(result), terra::res(template))
})

test_that("interpolate_pcnm_to_raster produces smooth surface", {
  skip_if_not_installed("fields")
  
  data <- setup_pcnm_test_data()
  
  pcnm_vector <- rnorm(nrow(data$occurrences))
  coords <- sf::st_coordinates(data$occurrences)
  template <- data$predictors[[1]]
  
  result <- interpolate_pcnm_to_raster(pcnm_vector, coords, template)
  
  # Should have values (not all NA)
  vals <- terra::values(result, na.rm = TRUE)
  expect_true(length(vals) > 0)
  
  # Values should be continuous/smooth (check range is reasonable)
  expect_true(is.finite(min(vals)))
  expect_true(is.finite(max(vals)))
})

# ==============================================================================
# Tests for create_pcnm_rasters()
# ==============================================================================

test_that("create_pcnm_rasters returns SpatRaster with correct layers", {
  skip_if_not_installed("fields")
  
  data <- setup_pcnm_test_data()
  
  # Create simple PCNM data
  pcnm_df <- data.frame(
    PCNM1 = rnorm(nrow(data$occurrences)),
    PCNM2 = rnorm(nrow(data$occurrences)),
    PCNM3 = rnorm(nrow(data$occurrences))
  )
  
  selected <- c("PCNM1", "PCNM3")
  
  result <- create_pcnm_rasters(
    pcnm_df,
    selected,
    data$occurrences,
    data$predictors[[1]]
  )
  
  # Should return SpatRaster
  expect_s4_class(result, "SpatRaster")
  
  # Should have correct number of layers
  expect_equal(terra::nlyr(result), length(selected))
  
  # Layer names should match selected
  expect_equal(names(result), selected)
})

test_that("create_pcnm_rasters handles single PCNM variable", {
  skip_if_not_installed("fields")
  
  data <- setup_pcnm_test_data()
  
  pcnm_df <- data.frame(
    PCNM1 = rnorm(nrow(data$occurrences)),
    PCNM2 = rnorm(nrow(data$occurrences))
  )
  
  selected <- "PCNM1"  # Single variable
  
  result <- create_pcnm_rasters(
    pcnm_df,
    selected,
    data$occurrences,
    data$predictors[[1]]
  )
  
  # Should return SpatRaster with one layer
  expect_s4_class(result, "SpatRaster")
  expect_equal(terra::nlyr(result), 1)
  expect_equal(names(result), selected)
})

test_that("create_pcnm_rasters matches template properties", {
  skip_if_not_installed("fields")
  
  data <- setup_pcnm_test_data()
  
  pcnm_df <- data.frame(
    PCNM1 = rnorm(nrow(data$occurrences)),
    PCNM2 = rnorm(nrow(data$occurrences))
  )
  
  selected <- c("PCNM1", "PCNM2")
  template <- data$predictors[[1]]
  
  result <- create_pcnm_rasters(
    pcnm_df,
    selected,
    data$occurrences,
    template
  )
  
  # Extent should match
  expect_equal(terra::ext(result), terra::ext(template))
  
  # Resolution should match
  expect_equal(terra::res(result), terra::res(template))
  
  # CRS should match
  expect_equal(terra::crs(result), terra::crs(template))
})

# ==============================================================================
# Integration test for createPCNM_fitModel
# ==============================================================================

test_that("createPCNM_fitModel full workflow", {
  skip_if_not_installed("vegan")
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("fields")
  skip_if_not_installed("CAST")
  #skip("Full integration test - slow, run manually")
  
  data <- setup_pcnm_test_data()
  
  # Prepare environmental predictors
  env_preds <- sf::st_drop_geometry(
    data$occurrences[, names(data$predictors)]
  )
  
  # Create CV structure
  cv_indices <- CAST::knndm(data$occurrences, data$predictors, k = 3)
  
  ctrl <- caret::rfeControl(
    method = "LGOCV",
    repeats = 3,
    number = 5,
    functions = caret::lrFuncs,
    index = cv_indices$indx_train,
    verbose = FALSE
  )
  
  result <- suppressWarnings(
    createPCNM_fitModel(
      x = data$occurrences,
      planar_proj = data$planar_proj,
      ctrl = ctrl,
      indices_knndm = cv_indices,
      sub = env_preds,
      predictors = data$predictors
    )
  )
  
  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("mod", "pred_mat", "cv_model", "pcnm"))
  
  # Check model objects
  expect_s3_class(result$mod, "glmnet")
  expect_s3_class(result$cv_model, "train")
  expect_s3_class(result$pred_mat, "data.frame")
  expect_s4_class(result$pcnm, "SpatRaster")
  
  # pred_mat should combine env and PCNM
  expect_true(ncol(result$pred_mat) > ncol(env_preds))
  
  # PCNM rasters should have same extent as predictors
  expect_equal(terra::ext(result$pcnm), terra::ext(data$predictors))
})

# ==============================================================================
# Property-based tests
# ==============================================================================

test_that("distance matrix triangle inequality holds", {
  data <- setup_pcnm_test_data()
  
  result <- calculate_distance_matrix(data$occurrences[1:10, ], data$planar_proj)
  
  # Check triangle inequality: d(i,k) <= d(i,j) + d(j,k)
  # Sample a few triplets to test
  for (i in 1:3) {
    for (j in 1:3) {
      for (k in 1:3) {
        if (i != j && j != k && i != k) {
          expect_true(result[i, k] <= result[i, j] + result[j, k] + 1e-10)
        }
      }
    }
  }
})

test_that("PCNM vectors capture spatial structure", {
  skip_if_not_installed("vegan")
  
  data <- setup_pcnm_test_data()
  dis <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  pcnm_vectors <- create_pcnm_vectors(dis, n_vectors = 5)
  
  # First PCNM vector should capture broad-scale spatial pattern
  # It should have spatial autocorrelation
  coords <- sf::st_coordinates(data$occurrences)
  
  # Correlation between PCNM1 and latitude should exist
  cor_lat <- cor(pcnm_vectors[, 1], coords[, 2])
  
  # Should have some correlation (positive or negative)
  expect_true(abs(cor_lat) > 0.1)
})

test_that("combined predictors maintain data integrity", {
  n <- 30
  env_preds <- data.frame(
    bio1 = rnorm(n),
    bio2 = rnorm(n)
  )
  
  pcnm_df <- data.frame(
    PCNM1 = rnorm(n),
    PCNM2 = rnorm(n),
    PCNM3 = rnorm(n)
  )
  
  result <- combine_predictors(env_preds, pcnm_df, c("PCNM1", "PCNM3"))
  
  # Original environmental data should be unchanged
  expect_equal(result[, 1:2], env_preds)
  
  # Selected PCNM data should match
  expect_equal(result$PCNM1, pcnm_df$PCNM1)
  expect_equal(result$PCNM3, pcnm_df$PCNM3)
})

test_that("interpolated rasters preserve spatial trends", {
  skip_if_not_installed("fields")
  
  data <- setup_pcnm_test_data()
  
  # Create PCNM vector with known spatial pattern (latitude gradient)
  coords <- sf::st_coordinates(data$occurrences)
  pcnm_vector <- coords[, 2]  # Use latitude as PCNM
  
  result <- interpolate_pcnm_to_raster(
    pcnm_vector,
    coords,
    data$predictors[[1]]
  )
  
  # Extract values at original points
  extracted <- terra::extract(result, data$occurrences, ID = FALSE)
  
  # Should be correlated with original values
  cor_val <- cor(extracted[, 1], pcnm_vector, use = "complete.obs")
  
  # High correlation expected for smooth interpolation
  expect_true(cor_val > 0.8)
})

# ==============================
# Additional tests for coverage
# ==============================
# ==============================================================================
# Tests for create_pcnm_vectors
# ==============================================================================

test_that("create_pcnm_vectors fails on truly empty matrix", {
  skip_if_not_installed("vegan")
  
  empty_dis <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(create_pcnm_vectors(empty_dis), "Cannot compute PCNM: distance matrix is empty")
})

test_that("create_pcnm_vectors returns correct number of columns", {
  skip_if_not_installed("vegan")
  
  data <- setup_pcnm_test_data()
  dis <- calculate_distance_matrix(data$occurrences, data$planar_proj)
  
  # Request fewer than available
  res <- create_pcnm_vectors(dis, n_vectors = 5)
  expect_lte(ncol(res), 5)
  expect_gte(ncol(res), 1)
  
  # Request more than available (should just return all)
  res2 <- create_pcnm_vectors(dis, n_vectors = 50)
  expect_lte(ncol(res2), nrow(data$occurrences) - 1) # pcnm returns <= n-1 vectors
})


# ==============================================================================
# Tests for combine_predictors & fit_combined_model with single predictor
# ==============================================================================

test_that("combine_predictors handles single PCNM column", {
  env <- data.frame(var1 = rnorm(5))
  pcnm_num <- data.frame(PCNM1 = rnorm(5)) # keep as data.frame
  combined <- combine_predictors(env, pcnm_num, "PCNM1")
  
  expect_s3_class(combined, "data.frame")
  expect_equal(ncol(combined), 2)
  expect_equal(names(combined)[2], "PCNM1")
})

test_that("fit_combined_model works with single predictor", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")
  
  preds <- data.frame(x1 = rnorm(20))
  y <- factor(sample(0:1, 20, replace = TRUE))
  cv <- list(indx_train = list(1:15))
  
  # Should run without error
  res <- suppressWarnings(fit_combined_model(preds, y, cv))
  expect_type(res, "list")
  expect_named(res, c("cv_model", "glmnet_model"))
})

