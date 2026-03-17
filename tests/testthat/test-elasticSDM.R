# Test suite for elasticSDM helper functions
# Uses dismo package example data (Bradypus occurrence data)

library(terra)
library(sf)
library(dplyr)

# Setup: Load example data used across tests
setup_sdm_test_data <- function() {
  x <- read.csv(file.path(system.file(package = "dismo"), 'ex', 'bradypus.csv'))
  x <- x[, c('lon', 'lat')]
  x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)

  files <- list.files(
    path = file.path(system.file(package = "dismo"), 'ex'),
    pattern = 'grd',
    full.names = TRUE
  )
  predictors <- terra::rast(files)

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

  expect_s4_class(result, "SpatExtent")
})

test_that("calculate_study_extent is larger than simple bbox", {
  data <- setup_sdm_test_data()

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

  expect_true(result$xmin < simple_bbox$xmin)
  expect_true(result$xmax > simple_bbox$xmax)
  expect_true(result$ymin < simple_bbox$ymin)
  expect_true(result$ymax > simple_bbox$ymax)
})

test_that("calculate_study_extent handles different projections", {
  data <- setup_sdm_test_data()

  result_5070 <- calculate_study_extent(
    data$occurrences, 5070, NULL, data$predictors
  )
  result_3857 <- calculate_study_extent(
    data$occurrences, 3857, NULL, data$predictors
  )

  expect_s4_class(result_5070, "SpatExtent")
  expect_s4_class(result_3857, "SpatExtent")
  expect_true(result_5070$xmin != result_3857$xmin)
})

# ==============================================================================
# Tests for generate_background_points()
# ==============================================================================

test_that("generate_background_points returns sf object with correct structure", {
  skip_if_not_installed("sdm")

  data <- setup_sdm_test_data()
  result <- generate_background_points(data$predictors, data$occurrences)

  expect_s3_class(result, "sf")
  expect_true("occurrence" %in% names(result))
  expect_true(all(result$occurrence == 0))
  expect_equal(sf::st_crs(result), sf::st_crs(data$predictors))
})

test_that("generate_background_points creates n points equal to occurrences", {
  skip_if_not_installed("sdm")

  data <- setup_sdm_test_data()
  result <- generate_background_points(data$predictors, data$occurrences)

  expect_equal(nrow(result), nrow(data$occurrences))
})

test_that("generate_background_points respects fact parameter", {
  skip_if_not_installed("sdm")

  data <- setup_sdm_test_data()

  result_1x <- generate_background_points(data$predictors, data$occurrences, fact = 1.0)
  result_2x <- generate_background_points(data$predictors, data$occurrences, fact = 2.0)

  expect_equal(nrow(result_1x), nrow(data$occurrences))
  expect_equal(nrow(result_2x), nrow(data$occurrences) * 2)
})

test_that("generate_background_points resample=TRUE branch runs and adds extra points", {
  skip_if_not_installed("sdm")

  data <- setup_sdm_test_data()

  # resample = TRUE takes the density-based resampling branch
  result <- generate_background_points(
    data$predictors,
    data$occurrences,
    fact = 1.0,
    resample = TRUE
  )

  expect_s3_class(result, "sf")
  expect_true("occurrence" %in% names(result))
  expect_true(all(result$occurrence == 0))
  # Should still produce ~n points (sdm_bg_pts_n + dens_bg_pts_n = n_requested)
  expect_true(nrow(result) > 0)
})

test_that("generate_background_points resample=TRUE points are within predictor extent", {
  skip_if_not_installed("sdm")

  data <- setup_sdm_test_data()

  result <- generate_background_points(
    data$predictors,
    data$occurrences,
    resample = TRUE
  )

  result_vect <- terra::vect(result)
  pred_extent <- terra::ext(data$predictors)
  result_extent <- terra::ext(result_vect)

  expect_true(result_extent$xmin >= pred_extent$xmin)
  expect_true(result_extent$xmax <= pred_extent$xmax)
  expect_true(result_extent$ymin >= pred_extent$ymin)
  expect_true(result_extent$ymax <= pred_extent$ymax)
})

test_that("background points are within predictor extent", {
  skip_if_not_installed("sdm")

  data <- setup_sdm_test_data()
  result <- generate_background_points(data$predictors, data$occurrences)

  result_vect <- terra::vect(result)
  pred_extent <- terra::ext(data$predictors)
  result_extent <- terra::ext(result_vect)

  expect_true(result_extent$xmin >= pred_extent$xmin)
  expect_true(result_extent$xmax <= pred_extent$xmax)
  expect_true(result_extent$ymin >= pred_extent$ymin)
  expect_true(result_extent$ymax <= pred_extent$ymax)
})

# ==============================================================================
# Tests for thin_occurrence_points()
# ==============================================================================

test_that("thin_occurrence_points reduces point density", {
  skip_if_not_installed("spThin")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1

  result <- thin_occurrence_points(data$occurrences, quantile_v = 0.05)

  expect_s3_class(result, "sf")
  expect_true(nrow(result) <= nrow(data$occurrences))
})

test_that("thin_occurrence_points respects quantile_v parameter", {
  skip_if_not_installed("spThin")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1

  result_aggressive <- thin_occurrence_points(data$occurrences, quantile_v = 0.1)
  result_light      <- thin_occurrence_points(data$occurrences, quantile_v = 0.01)

  expect_true(nrow(result_aggressive) <= nrow(result_light))
})

test_that("thin_occurrence_points returns points from original dataset", {
  skip_if_not_installed("spThin")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1

  result <- thin_occurrence_points(data$occurrences, quantile_v = 0.025)

  intersects <- lengths(sf::st_intersects(result, data$occurrences))
  expect_true(all(intersects > 0))
})

test_that("thinning preserves data integrity", {
  skip_if_not_installed("spThin")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  data$occurrences$id <- seq_len(nrow(data$occurrences))

  result <- thin_occurrence_points(data$occurrences, quantile_v = 0.025)

  expect_true(all(result$id %in% data$occurrences$id))
  expect_equal(length(result$id), length(unique(result$id)))
})

# ==============================================================================
# Tests for extract_predictors_to_points()
# ==============================================================================

test_that("extract_predictors_to_points adds predictor columns", {
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1

  result <- extract_predictors_to_points(data$occurrences, data$predictors)

  expect_s3_class(result, "sf")
  expect_true(all(names(data$predictors) %in% names(result)))
  expect_equal(nrow(result), nrow(data$occurrences))
})

test_that("extract_predictors_to_points preserves original columns", {
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  data$occurrences$test_col <- "test"

  result <- extract_predictors_to_points(data$occurrences, data$predictors)

  expect_true("occurrence" %in% names(result))
  expect_true("test_col" %in% names(result))
})

test_that("extract_predictors_to_points handles NA values gracefully", {
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1

  result <- extract_predictors_to_points(data$occurrences, data$predictors)

  pred_names <- names(data$predictors)
  pred_data  <- sf::st_drop_geometry(result[, pred_names])

  complete_cases <- sum(complete.cases(pred_data))
  expect_true(complete_cases / nrow(result) > 0.8)
})

test_that("predictor extraction maintains spatial relationships", {
  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  data$occurrences$point_id <- seq_len(nrow(data$occurrences))

  result <- extract_predictors_to_points(data$occurrences, data$predictors)

  expect_equal(nrow(result), nrow(data$occurrences))
  expect_equal(sort(result$point_id), sort(data$occurrences$point_id))
  expect_equal(sf::st_coordinates(data$occurrences), sf::st_coordinates(result))
})

# ==============================================================================
# Tests for create_spatial_cv_folds()
# ==============================================================================

test_that("create_spatial_cv_folds returns knndm object", {
  skip_if_not_installed("CAST")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- 1
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  train <- x_extracted[1:50, ]

  result <- create_spatial_cv_folds(train, data$predictors, k = 5)

  expect_s3_class(result, "knndm")
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

  for (i in seq_along(result$indx_train)) {
    overlap <- intersect(result$indx_train[[i]], result$indx_test[[i]])
    expect_equal(length(overlap), 0)
  }
})

# ==============================================================================
# Tests for perform_feature_selection()
# ==============================================================================

test_that("perform_feature_selection returns rfe object", {
  skip_if_not_installed("CAST")
  skip_if_not_installed("caret")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- factor(sample(0:1, nrow(data$occurrences), replace = TRUE))
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  x_extracted <- x_extracted[complete.cases(sf::st_drop_geometry(x_extracted)), ]
  train <- x_extracted[1:50, ]

  cv_indices <- create_spatial_cv_folds(train, data$predictors, k = 3)
  result <- suppressWarnings(perform_feature_selection(train, cv_indices))

  expect_s3_class(result, "rfe")
  expect_true("optVariables" %in% names(result))
})

# ==============================================================================
# Tests for fit_elastic_net_model()
# ==============================================================================

test_that("fit_elastic_net_model returns list with required components", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- factor(sample(0:1, nrow(data$occurrences), replace = TRUE))
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  x_extracted <- x_extracted[complete.cases(sf::st_drop_geometry(x_extracted)), ]
  train <- x_extracted[1:50, ]

  cv_indices      <- create_spatial_cv_folds(train, data$predictors, k = 3)
  selected_preds  <- names(data$predictors)[1:4]

  result <- fit_elastic_net_model(train, selected_preds, cv_indices)

  expect_type(result, "list")
  expect_named(result, c("cv_model", "glmnet_model", "selected_data"))
  expect_s3_class(result$cv_model, "train")
  expect_s3_class(result$glmnet_model, "glmnet")
})

test_that("fit_elastic_net_model glmnet uses bestTune lambda and alpha", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- factor(sample(0:1, nrow(data$occurrences), replace = TRUE))
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  x_extracted <- x_extracted[complete.cases(sf::st_drop_geometry(x_extracted)), ]
  train <- x_extracted[1:50, ]

  cv_indices     <- create_spatial_cv_folds(train, data$predictors, k = 3)
  selected_preds <- names(data$predictors)[1:4]

  result <- fit_elastic_net_model(train, selected_preds, cv_indices)

  # lambda stored on glmnet object should match what caret selected
  expect_equal(result$glmnet_model$lambda, result$cv_model$bestTune$lambda)
  # selected_data should only contain the predictor columns, no geometry/occurrence
  expect_false("occurrence" %in% names(result$selected_data))
  expect_false("geometry" %in% names(result$selected_data))
  expect_true(all(selected_preds %in% names(result$selected_data)))
})

# ==============================================================================
# Tests for evaluate_model_performance()
# ==============================================================================

test_that("evaluate_model_performance returns a confusionMatrix", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- factor(sample(0:1, nrow(data$occurrences), replace = TRUE))
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  x_extracted <- x_extracted[complete.cases(sf::st_drop_geometry(x_extracted)), ]

  index <- unlist(caret::createDataPartition(x_extracted$occurrence, p = 0.8))
  train <- x_extracted[index, ]
  test  <- x_extracted[-index, ]

  selected_preds <- names(data$predictors)[1:4]
  cv_indices     <- create_spatial_cv_folds(train, data$predictors, k = 3)
  model_results  <- fit_elastic_net_model(train, selected_preds, cv_indices)

  cm <- evaluate_model_performance(
    model_results$glmnet_model,
    test,
    data$predictors,
    selected_preds
  )

  expect_s3_class(cm, "confusionMatrix")
  expect_true("table" %in% names(cm))
  expect_true("overall" %in% names(cm))
})

# ==============================================================================
# Tests for create_spatial_predictions()
# ==============================================================================

test_that("create_spatial_predictions returns a SpatRaster", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- factor(sample(0:1, nrow(data$occurrences), replace = TRUE))
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  x_extracted <- x_extracted[complete.cases(sf::st_drop_geometry(x_extracted)), ]
  train <- x_extracted[1:50, ]

  selected_preds <- names(data$predictors)[1:4]
  cv_indices     <- create_spatial_cv_folds(train, data$predictors, k = 3)
  model_results  <- fit_elastic_net_model(train, selected_preds, cv_indices)

  rast <- create_spatial_predictions(
    model_results$glmnet_model,
    data$predictors,
    selected_preds
  )

  expect_s4_class(rast, "SpatRaster")
  # Values should be probabilities bounded 0-1
  vals <- terra::values(rast, na.rm = TRUE)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("create_spatial_predictions matches predictor extent and resolution", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")

  data <- setup_sdm_test_data()
  data$occurrences$occurrence <- factor(sample(0:1, nrow(data$occurrences), replace = TRUE))
  x_extracted <- extract_predictors_to_points(data$occurrences, data$predictors)
  x_extracted <- x_extracted[complete.cases(sf::st_drop_geometry(x_extracted)), ]
  train <- x_extracted[1:50, ]

  selected_preds <- names(data$predictors)[1:4]
  cv_indices     <- create_spatial_cv_folds(train, data$predictors, k = 3)
  model_results  <- fit_elastic_net_model(train, selected_preds, cv_indices)

  rast <- create_spatial_predictions(
    model_results$glmnet_model,
    data$predictors,
    selected_preds
  )

  expect_equal(terra::ext(rast),   terra::ext(data$predictors))
  expect_equal(terra::res(rast),   terra::res(data$predictors))
  expect_equal(terra::nrow(rast),  terra::nrow(data$predictors))
  expect_equal(terra::ncol(rast),  terra::ncol(data$predictors))
})

# ==============================================================================
# Integration tests
# ==============================================================================

test_that("complete workflow produces expected outputs including AOA", {
  skip_if_not_installed("sdm")
  skip_if_not_installed("spThin")
  skip_if_not_installed("CAST")
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")

  data <- setup_sdm_test_data()

  result <- suppressWarnings(
    elasticSDM(
      x = data$occurrences,
      predictors = data$predictors,
      planar_projection = data$planar_proj,
      quantile_v = 0.05
    )
  )

  expect_type(result, "list")

  # All previously expected names still present
  expected_names <- c(
    "RasterPredictions", "Predictors", "PCNM", "Model",
    "CVStructure", "ConfusionMatrix", "TrainData",
    "TestData", "PredictMatrix",
    "AOA", "AOA_Diagnostics"           # new AOA outputs
  )
  expect_true(all(expected_names %in% names(result)))

  # Core spatial outputs
  expect_s4_class(result$RasterPredictions, "SpatRaster")
  expect_s3_class(result$ConfusionMatrix, "confusionMatrix")
  expect_s3_class(result$TrainData, "sf")
  expect_s3_class(result$TestData, "sf")

  # AOA surface: should be a 3-layer SpatRaster named AOA, DI, LPD
  expect_s4_class(result$AOA, "SpatRaster")
  expect_equal(terra::nlyr(result$AOA), 3)
  expect_equal(names(result$AOA), c("AOA", "DI", "LPD"))

  # AOA diagnostics
  expect_type(result$AOA_Diagnostics, "list")
  expect_true("threshold" %in% names(result$AOA_Diagnostics))
  expect_true("AOA_coverage" %in% names(result$AOA_Diagnostics))
  expect_true(result$AOA_Diagnostics$AOA_coverage >= 0 &
                result$AOA_Diagnostics$AOA_coverage <= 1)
})

test_that("complete workflow with PCNM=FALSE produces expected outputs", {
  skip_if_not_installed("sdm")
  skip_if_not_installed("spThin")
  skip_if_not_installed("CAST")
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")

  data <- setup_sdm_test_data()

  result <- suppressWarnings(
    elasticSDM(
      x = data$occurrences,
      predictors = data$predictors,
      planar_projection = data$planar_proj,
      quantile_v = 0.05,
      PCNM = FALSE
    )
  )

  expect_type(result, "list")
  expect_null(result$PCNM)
  expect_s4_class(result$RasterPredictions, "SpatRaster")
  expect_s4_class(result$AOA, "SpatRaster")
})