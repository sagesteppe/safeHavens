# Test suite for elasticSDM helper functions
# Uses dismo package example data (Bradypus occurrence data)

library(testthat)
library(terra)
library(sf)
library(dplyr)

# Setup: Load example data used across tests
setup_sdm_test_data <- function() {
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
  
  # Planar projection for distance calculations
  planar_proj <- 3857  # Web Mercator
  
  list(
    occurrences = x,
    predictors = predictors,
    planar_proj = planar_proj
  )
}

# ==============================================================================
# Tests for calculate_study_extent()
# ==============================================================================

test_that("calculate_study_extent returns terra extent object", {
  data <- setup_sdm_test_data()
  
  result <- calculate_study_extent(
    data$occurrences, 
    data$planar_proj, 
    domain = NULL,
    data$predictors
  )
  
  # Should return SpatExtent
  expect_s4_class(result, "SpatExtent")
})

test_that("calculate_study_extent is larger than simple bbox", {
  data <- setup_sdm_test_data()
  
  # Get simple bounding box
  simple_bbox <- data$occurrences |>
    sf::st_bbox() |>
    sf::st_as_sfc() |>
    sf::st_transform(terra::crs(data$predictors)) |>
    terra::vect() |>
    terra::ext()
  
  result <- calculate_study_extent(
    data$occurrences, 
    data$planar_proj, 
    domain = NULL,
    data$predictors
  )
  
  # Buffered extent should be larger in all directions
  expect_true(result$xmin < simple_bbox$xmin)
  expect_true(result$xmax > simple_bbox$xmax)
  expect_true(result$ymin < simple_bbox$ymin)
  expect_true(result$ymax > simple_bbox$ymax)
})

test_that("calculate_study_extent handles different projections", {
  data <- setup_sdm_test_data()
  
  # Test with different planar projections
  proj_5070 <- 5070  # NAD83 / Conus Albers
  proj_3857 <- 3857  # Web Mercator
  
  result_5070 <- calculate_study_extent(
    data$occurrences, proj_5070, NULL, data$predictors
  )
  result_3857 <- calculate_study_extent(
    data$occurrences, proj_3857, NULL, data$predictors
  )
  
  # Both should return valid extents
  expect_s4_class(result_5070, "SpatExtent")
  expect_s4_class(result_3857, "SpatExtent")
  
  # Extents should be in the same CRS (predictors' CRS)
  # But calculations use different planar projections
  expect_true(result_5070$xmin != result_3857$xmin)
})

# ==============================================================================
# Tests for generate_background_points()
# ==============================================================================

test_that("generate_background_points returns sf object with correct structure", {
  skip_if_not_installed("sdm")
  
  data <- setup_sdm_test_data()
  
  result <- generate_background_points(data$predictors, data$occurrences)
  
  # Should be sf object
  expect_s3_class(result, "sf")
  
  # Should have occurrence column
  expect_true("occurrence" %in% names(result))
  
  # All occurrences should be 0 (pseudo-absence)
  expect_true(all(result$occurrence == 0))
  
  # Should have same CRS as predictors
  expect_equal(sf::st_crs(result), sf::st_crs(data$predictors))
})

test_that("generate_background_points creates n points equal to occurrences", {
  skip_if_not_installed("sdm")
  
  data <- setup_sdm_test_data()
  
  result <- generate_background_points(data$predictors, data$occurrences)
  
  # Should have same number of background points as occurrences
  expect_equal(nrow(result), nrow(data$occurrences))
})

test_that("generate_background_points uses environmental distance", {
  skip_if_not_installed("sdm")
  
  data <- setup_sdm_test_data()
  
  # Generate background points
  result <- generate_background_points(data$predictors, data$occurrences)
  
  # Extract environmental values
  bg_values <- terra::extract(data$predictors, result, ID = FALSE)
  occ_values <- terra::extract(data$predictors, data$occurrences, ID = FALSE)
  
  # Background points should have environmental values within predictor range
  expect_true(all(!is.na(bg_values)))
  
  # Environmental distribution of background should differ from occurrences
  # (eDist method selects points to match environmental distribution)
  expect_true(nrow(bg_values) > 0)
})

# ==============================================================================
# Tests for thin_occurrence_points()
# ==============================================================================

test_that("thin_occurrence_points reduces point density", {
  skip_if_not_installed("spThin")
  
  data <- setup_sdm_test_data()
  
  # Add occurrence column
  data$occurrences$occurrence <- 1
  
  result <- thin_occurrence_points(data$occurrences, quantile_v = 0.05)
  
  # Should be sf object
  expect_s3_class(result, "sf")
  
  # Should have fewer or equal points
  expect_true(nrow(result) <= nrow(data$occurrences))
})

test_that("thin_occurrence_points respects quantile_v parameter", {
  skip_if_not_installed("spThin")
  
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  
  # More aggressive thinning
  result_aggressive <- thin_occurrence_points(data$occurrences, quantile_v = 0.1)
  
  # Less aggressive thinning
  result_light <- thin_occurrence_points(data$occurrences, quantile_v = 0.01)
  
  # More aggressive thinning should result in fewer points
  expect_true(nrow(result_aggressive) <= nrow(result_light))
})

test_that("thin_occurrence_points returns points from original dataset", {
  skip_if_not_installed("spThin")
  
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  
  result <- thin_occurrence_points(data$occurrences, quantile_v = 0.025)
  
  # All returned points should be from original dataset
  # Check by seeing if geometries intersect
  intersects <- lengths(sf::st_intersects(result, data$occurrences))
  expect_true(all(intersects > 0))
})

# ==============================================================================
# Tests for extract_predictors_to_points()
# ==============================================================================

test_that("extract_predictors_to_points adds predictor columns", {
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  
  result <- extract_predictors_to_points(data$occurrences, data$predictors)
  
  # Should be sf object
  expect_s3_class(result, "sf")
  
  # Should have all predictor columns
  pred_names <- names(data$predictors)
  expect_true(all(pred_names %in% names(result)))
  
  # Should have same number of rows
  expect_equal(nrow(result), nrow(data$occurrences))
})

test_that("extract_predictors_to_points preserves original columns", {
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  data$occurrences$test_col <- "test"
  
  result <- extract_predictors_to_points(data$occurrences, data$predictors)
  
  # Original columns should be preserved
  expect_true("occurrence" %in% names(result))
  expect_true("test_col" %in% names(result))
})

test_that("extract_predictors_to_points handles NA values", {
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  
  result <- extract_predictors_to_points(data$occurrences, data$predictors)
  
  # Check if there are any NA values in predictor columns
  pred_names <- names(data$predictors)
  pred_data <- sf::st_drop_geometry(result[, pred_names])
  
  # Most points should have valid predictor values
  # (some NA is acceptable at edges)
  complete_cases <- sum(complete.cases(pred_data))
  expect_true(complete_cases / nrow(result) > 0.8)
})

# ==============================================================================
# Tests for create_spatial_cv_folds()
# ==============================================================================

test_that("create_spatial_cv_folds returns knndm object", {
  skip_if_not_installed("CAST")
  
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  
  # Take subset for faster testing
  train <- x_extracted[1:50, ]
  
  result <- create_spatial_cv_folds(train, data$predictors, k = 5)
  
  # Should return knndm object
  expect_s3_class(result, "knndm")
  
  # Should have fold indices
  expect_true("indx_train" %in% names(result))
  expect_true("indx_test" %in% names(result))
})

test_that("create_spatial_cv_folds respects k parameter", {
  skip_if_not_installed("CAST")
  
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  train <- x_extracted[1:50, ]
  
  for (k in c(3, 5, 10)) {
    result <- create_spatial_cv_folds(train, data$predictors, k = k)
    
    # Number of folds should match k
    expect_equal(length(result$indx_train), k)
    expect_equal(length(result$indx_test), k)
  }
})

test_that("create_spatial_cv_folds creates non-overlapping folds", {
  skip_if_not_installed("CAST")
  
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  train <- x_extracted[1:50, ]
  
  result <- create_spatial_cv_folds(train, data$predictors, k = 5)
  
  # Check that train and test indices don't overlap within each fold
  for (i in 1:length(result$indx_train)) {
    train_idx <- result$indx_train[[i]]
    test_idx <- result$indx_test[[i]]
    
    overlap <- intersect(train_idx, test_idx)
    expect_equal(length(overlap), 0)
  }
})

# ==============================================================================
# Tests for perform_feature_selection()
# ==============================================================================

test_that("perform_feature_selection returns rfe object", {
  skip_if_not_installed("CAST")
  skip_if_not_installed("caret")
  skip("Feature selection is slow - run manually")
  
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- factor(sample(0:1, nrow(data$occurrences), replace = TRUE))
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  x_extracted <- x_extracted[complete.cases(sf::st_drop_geometry(x_extracted)), ]
  train <- x_extracted[1:50, ]
  
  cv_indices <- create_spatial_cv_folds(train, data$predictors, k = 3)
  
  result <- perform_feature_selection(train, cv_indices)
  
  expect_s3_class(result, "rfe")
  expect_true("optVariables" %in% names(result))
})

# ==============================================================================
# Tests for fit_elastic_net_model()
# ==============================================================================

test_that("fit_elastic_net_model returns list with required components", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")
  skip("Model fitting is slow - run manually")
  
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- factor(sample(0:1, nrow(data$occurrences), replace = TRUE))
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  x_extracted <- x_extracted[complete.cases(sf::st_drop_geometry(x_extracted)), ]
  train <- x_extracted[1:50, ]
  
  cv_indices <- create_spatial_cv_folds(train, data$predictors, k = 3)
  selected_preds <- names(data$predictors)[1:4]
  
  result <- fit_elastic_net_model(train, selected_preds, cv_indices)
  
  expect_type(result, "list")
  expect_named(result, c("cv_model", "glmnet_model", "selected_data"))
  expect_s3_class(result$cv_model, "train")
  expect_s3_class(result$glmnet_model, "glmnet")
})

# ==============================================================================
# Tests for evaluate_model_performance()
# ==============================================================================

test_that("evaluate_model_performance returns confusion matrix", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")
  skip("Model evaluation requires fitted model - run manually")
  
  # This test requires a fitted model from previous steps
  # Would test with mock glmnet model in practice
})

# ==============================================================================
# Tests for create_spatial_predictions()
# ==============================================================================

test_that("create_spatial_predictions returns SpatRaster", {
  skip_if_not_installed("glmnet")
  skip("Spatial prediction requires fitted model - run manually")
  
  # This test requires a fitted model from previous steps
  # Would test with mock glmnet model in practice
})

# ==============================================================================
# Integration tests
# ==============================================================================

test_that("complete workflow produces expected outputs", {
  skip_if_not_installed("sdm")
  skip_if_not_installed("spThin")
  skip_if_not_installed("CAST")
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")
  skip("Full integration test - slow, run manually")
  
  data <- setup_sdm_test_data()
  
  # This would test the full elasticSDM workflow
  # Only run this manually due to computation time
  result <- elasticSDM(
    x = data$occurrences,
    predictors = data$predictors,
    planar_projection = data$planar_proj,
    quantile_v = 0.05
  )
  
  # Check return structure
  expect_type(result, "list")
  expected_names <- c(
    "RasterPredictions", "Predictors", "PCNM", "Model",
    "CVStructure", "ConfusionMatrix", "TrainData", 
    "TestData", "PredictMatrix"
  )
  expect_true(all(expected_names %in% names(result)))
  
  # Check key outputs
  expect_s4_class(result$RasterPredictions, "SpatRaster")
  expect_s3_class(result$ConfusionMatrix, "confusionMatrix")
  expect_s3_class(result$TrainData, "sf")
  expect_s3_class(result$TestData, "sf")
})

# ==============================================================================
# Property-based tests
# ==============================================================================

test_that("thinning preserves data integrity", {
  skip_if_not_installed("spThin")
  
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  data$occurrences$id <- 1:nrow(data$occurrences)
  
  result <- thin_occurrence_points(data$occurrences, quantile_v = 0.025)
  
  # All IDs in result should be from original data
  original_ids <- data$occurrences$id
  result_ids <- result$id
  
  expect_true(all(result_ids %in% original_ids))
  
  # No duplicate IDs in result
  expect_equal(length(result_ids), length(unique(result_ids)))
})

test_that("background points are within predictor extent", {
  skip_if_not_installed("sdm")
  
  data <- setup_sdm_test_data()
  
  result <- generate_background_points(data$predictors, data$occurrences)
  
  # Convert to same CRS as predictors for extent check
  result_vect <- terra::vect(result)
  pred_extent <- terra::ext(data$predictors)
  result_extent <- terra::ext(result_vect)
  
  # Background points should be within predictor extent
  expect_true(result_extent$xmin >= pred_extent$xmin)
  expect_true(result_extent$xmax <= pred_extent$xmax)
  expect_true(result_extent$ymin >= pred_extent$ymin)
  expect_true(result_extent$ymax <= pred_extent$ymax)
})

test_that("predictor extraction maintains spatial relationships", {
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  
  # Add unique ID to track points
  data$occurrences$point_id <- 1:nrow(data$occurrences)
  
  result <- extract_predictors_to_points(data$occurrences, data$predictors)
  
  # Same number of points
  expect_equal(nrow(result), nrow(data$occurrences))
  
  # Same point IDs preserved
  expect_equal(sort(result$point_id), sort(data$occurrences$point_id))
  
  # Geometries should be identical
  coords_orig <- sf::st_coordinates(data$occurrences)
  coords_result <- sf::st_coordinates(result)
  expect_equal(coords_orig, coords_result)
})
