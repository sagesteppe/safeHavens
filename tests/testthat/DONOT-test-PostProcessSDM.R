library(terra)
library(sf)
library(dplyr)

setup_sdModel_test_data <- function() {
  # Path to extdata inside the installed package
  rds_path <- system.file("extdata", "sdModel.rds", package = "safeHavens")
  
  if (rds_path == "") {
    stop("sdModel.rds not found in package extdata")
  }
  
  sdModel <- readRDS(rds_path)
  
  # Optionally extract predictors & other components if you need them in tests
  list(
    sdModel = sdModel,
    predictors = sdModel$Predictors,  # adjust based on slot names
    occurrences = dplyr::bind_rows(sdModel$TrainData, sdModel$TestData),
    planar_proj = sdModel$planar_proj  # if stored
  )
}

# ==============================================================================
# Helper: nn_distribution()
# ==============================================================================

test_that("nn_distribution computes nearest neighbor distances correctly", {
  # Minimal presence data
  pres <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, occurrence = 1),
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

# ==============================================================================
# PostProcessSDM minimal test using sdModel fixture
# ==============================================================================

test_that("PostProcessSDM handles small SDM object", {
  data <- setup_sdModel_test_data()
  sdModel <- data$sdModel
  predictors <- data$predictors
  occurrences <- data$occurrences
  planar_proj <- data$planar_proj

  rast_cont <- sdModel$RasterPredictions
  test_data <- sdModel$TestData
  train_data <- sdModel$TrainData
  
  result <- PostProcessSDM(
    rast_cont = rast_cont,
    test = test_data,
    train = train_data,
    planar_projection = planar_proj
  )
  
  # Structure checks
  expect_type(result, "list")
  expect_named(result, c("FinalRasters", "Threshold"))
  
  fr <- result$FinalRasters
  expect_s4_class(fr, "SpatRaster")
  expect_equal(terra::nlyr(fr), 4)
  expect_equal(names(fr), c("Predictions", "Threshold", "Clipped", "Supplemented"))
  
  # Threshold table
  expect_s3_class(result$Threshold, "data.frame")
  expect_true("sensitivity" %in% names(result$Threshold))
  
  # Raster values reasonable
  vals <- terra::values(fr[[1]], na.rm = TRUE)
  expect_true(all(vals >= 0 & vals <= 1))
})

# ==============================================================================
# PostProcessSDM edge cases
# ==============================================================================

test_that("PostProcessSDM works with planar projection as EPSG numeric", {
  data <- setup_sdModel_test_data()
  sdModel <- data$sdModel
  predictors <- data$predictors
  occurrences <- data$occurrences
  planar_proj <- data$planar_proj

  rast_cont <- sdModel$RasterPredictions
  test_data <- sdModel$TestData
  train_data <- sdModel$TrainData

  
  # Use UTM zone for planar projection
  planar_proj <- 32611
  
  result <- expect_silent(
    PostProcessSDM(
      rast_cont = rast_cont,
      test = test_data,
      train = test_data,
      planar_projection = planar_proj
    )
  )
  
  fr <- result$FinalRasters
  expect_s4_class(fr, "SpatRaster")
  expect_equal(terra::nlyr(fr), 4)
})

