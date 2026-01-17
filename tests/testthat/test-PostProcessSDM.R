# tests/testthat/test-PostProcessSDM.R

# Helper to set up test data from extdata
setup_sdModel_test_data <- function() {
  rds_path <- system.file("extdata", "sdModel.rds", package = "safeHavens")
  if (rds_path == "") {
    skip("sdModel.rds not found in package extdata")
  }
  
  sdModel <- readRDS(rds_path)
  sdModel$RasterPredictions <- terra::unwrap(sdModel$RasterPredictions)
  sdModel$Predictors <- terra::unwrap(sdModel$Predictors)
  
  list(
    sdModel = sdModel,
    rast_cont = sdModel$RasterPredictions,
    predictors = sdModel$Predictors,
    train = sdModel$TrainData,
    test = sdModel$TestData,
    planar_proj = if(!is.null(sdModel$planar_proj)) sdModel$planar_proj else "EPSG:5070"
  )
}

# Test helper function: nn_distribution ----

test_that("nn_distribution computes nearest neighbor distances correctly", {
  skip_if_not_installed("sf")
  
  # Minimal presence data
  pres <- sf::st_as_sf(
    data.frame(x = c(1, 2, 3, 4, 5), y = c(1, 2, 3, 4, 5), occurrence = 1),
    coords = c("x", "y"), crs = 4326
  )
  
  # Mock cv fold indices
  folds <- list(1:3, 4:5)
  
  dists <- lapply(folds, nn_distribution, y = pres)
  
  # Check that each element is numeric and has expected length
  expect_true(all(sapply(dists, is.numeric)))
  expect_equal(length(dists[[1]]), 3)
  expect_equal(length(dists[[2]]), 2)
  
  # Distances should be non-negative
  expect_true(all(unlist(dists) >= 0))
})

test_that("nn_distribution returns units object from sf", {
  skip_if_not_installed("sf")
  
  pres <- sf::st_as_sf(
    data.frame(x = 1:4, y = 1:4, occurrence = 1),
    coords = c("x", "y"), crs = 4326
  )
  
  result <- nn_distribution(list(1:2), pres)
  
  # sf::st_distance returns units object
  expect_s3_class(result, "units")
})

# Test PostProcessSDM main function ----

test_that("PostProcessSDM runs with real data from extdata", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  result <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    thresh_metric = 'sensitivity',
    quant_amt = 0.25,
    planar_projection = data$planar_proj
  )
  
  expect_type(result, "list")
  expect_named(result, c("FinalRasters", "Threshold"))
})

test_that("PostProcessSDM returns correct raster structure", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  result <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    planar_projection = data$planar_proj
  )
  
  # Check FinalRasters
  expect_s4_class(result$FinalRasters, "SpatRaster")
  expect_equal(terra::nlyr(result$FinalRasters), 4)
  expect_equal(
    names(result$FinalRasters), 
    c('Predictions', 'Threshold', 'Clipped', 'Supplemented')
  )
})

test_that("PostProcessSDM threshold statistics are returned", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  result <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    planar_projection = data$planar_proj
  )
  
  # Check Threshold statistics
  expect_s3_class(result$Threshold, "data.frame")
  expect_true("sensitivity" %in% names(result$Threshold))
})

test_that("PostProcessSDM works with different threshold metrics", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  metrics <- c("kappa", "spec_sens", "no_omission")
  
  for (metric in metrics) {
    result <- PostProcessSDM(
      rast_cont = data$rast_cont,
      test = data$test,
      train = data$train,
      thresh_metric = metric,
      planar_projection = data$planar_proj
    )
    
    expect_s4_class(result$FinalRasters, "SpatRaster")
    expect_named(result, c("FinalRasters", "Threshold"))
  }
})

test_that("PostProcessSDM binary raster has only 1 and NA values", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  result <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    planar_projection = data$planar_proj
  )
  
  # Check Threshold layer (binary)
  binary_vals <- unique(terra::values(result$FinalRasters[['Threshold']]))
  binary_vals <- binary_vals[!is.na(binary_vals)]
  
  expect_true(all(binary_vals %in% c(1)))
  expect_length(binary_vals, 1)
})

test_that("PostProcessSDM clipped raster is subset of binary", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  result <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    planar_projection = data$planar_proj
  )
  
  # Clipped should have fewer or equal cells than Threshold
  n_threshold <- sum(!is.na(terra::values(result$FinalRasters[['Threshold']])))
  n_clipped <- sum(!is.na(terra::values(result$FinalRasters[['Clipped']])))
  
  expect_true(n_clipped <= n_threshold)
})

test_that("PostProcessSDM supplemented adds cells back", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  result <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    planar_projection = data$planar_proj
  )
  
  # Supplemented should have more or equal cells than Clipped
  n_clipped <- sum(!is.na(terra::values(result$FinalRasters[['Clipped']])))
  n_supplemented <- sum(!is.na(terra::values(result$FinalRasters[['Supplemented']])))
  
  expect_true(n_supplemented >= n_clipped)
})

test_that("PostProcessSDM respects quant_amt parameter", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  result1 <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    quant_amt = 0.1,
    planar_projection = data$planar_proj
  )
  
  result2 <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    quant_amt = 0.5,
    planar_projection = data$planar_proj
  )
  
  # Different quantiles should produce different clipping extents
  n_clipped1 <- sum(!is.na(terra::values(result1$FinalRasters[['Clipped']])))
  n_clipped2 <- sum(!is.na(terra::values(result2$FinalRasters[['Clipped']])))
  
  # Higher quantile = larger buffer = more cells (generally)
  expect_true(n_clipped1 != n_clipped2)
})

test_that("PostProcessSDM handles edge case with no outside points", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  skip("Edge case - may not apply to all datasets")
  
  data <- setup_sdModel_test_data()
  
  # This test verifies the function handles the case where all presence
  # points fall within suitable habitat after thresholding
  result <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    thresh_metric = 'no_omission',  # Very permissive threshold
    planar_projection = data$planar_proj
  )
  
  expect_s4_class(result$FinalRasters, "SpatRaster")
})

test_that("PostProcessSDM preserves CRS throughout processing", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  data <- setup_sdModel_test_data()
  
  result <- PostProcessSDM(
    rast_cont = data$rast_cont,
    test = data$test,
    train = data$train,
    planar_projection = data$planar_proj
  )
  
  # All output layers should have the same CRS as input
  original_crs <- terra::crs(data$rast_cont)
  
  for (lyr in names(result$FinalRasters)) {
    expect_equal(
      terra::crs(result$FinalRasters[[lyr]]), 
      original_crs,
      label = paste("CRS for layer", lyr)
    )
  }
})

# Minimal synthetic data tests ----

test_that("PostProcessSDM works with minimal synthetic data", {
  skip_if_not_installed("sf")
  skip_if_not_installed("dismo")
  skip_if_not_installed("CAST")
  
  # Create minimal raster
  r <- terra::rast(ncols = 20, nrows = 20, xmin = -120, xmax = -119,
                   ymin = 45, ymax = 46, crs = "EPSG:4326")
  terra::values(r) <- runif(400, 0, 1)
  names(r) <- "s0"
  
  # Create minimal test/train data
  set.seed(123)
  coords <- data.frame(
    x = runif(30, -120, -119),
    y = runif(30, 45, 46)
  )
  
  train <- coords[1:20, ] |>
    dplyr::mutate(occurrence = rep(c(1, 0), each = 10))
  
  test <- coords[21:30, ] |>
    dplyr::mutate(occurrence = rep(c(1, 0), each = 5))

  train <- sf::st_as_sf(train, coords = c("x","y"), crs = 4326)
  test  <- sf::st_as_sf(test, coords = c("x","y"), crs = 4326)
  
  result <- PostProcessSDM(
    rast_cont = r,
    test = test,
    train = train,
    planar_projection = "EPSG:32610"
  )
  
  expect_type(result, "list")
  expect_named(result, c("FinalRasters", "Threshold"))
  expect_equal(terra::nlyr(result$FinalRasters), 4)
})