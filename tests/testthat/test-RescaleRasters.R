library(testthat)
library(glmnet)
library(terra)

test_that("RescaleRasters rescales predictors and returns valid coefficients", {

  set.seed(42)

  # -------------------------------
  # Mock training data
  # -------------------------------
  n <- 30
  training_data <- data.frame(
    occurrence = sample(0:1, n, replace = TRUE),
    var1 = rnorm(n),
    var2 = rnorm(n)
  )

  x <- as.matrix(training_data[, c("var1", "var2")])
  y <- training_data$occurrence

  # -------------------------------
  # Fit glmnet model
  # -------------------------------
  model <- glmnet::glmnet(
    x, y,
    alpha = 0.5,
    family = "binomial"
  )

  # -------------------------------
  # Raster predictors
  # -------------------------------
  r1 <- terra::rast(nrows = 5, ncols = 5, vals = rnorm(25))
  r2 <- terra::rast(nrows = 5, ncols = 5, vals = rnorm(25))
  predictors <- c(r1, r2)
  names(predictors) <- c("var1", "var2")

  # -------------------------------
  # Prediction matrix
  # -------------------------------
  pred_mat <- x

  # -------------------------------
  # Run function
  # -------------------------------

  out <- RescaleRasters(
    model = model,
    predictors = predictors,
    training_data = training_data,
    pred_mat = pred_mat
  )

  # -------------------------------
  # Structural tests
  # -------------------------------
  expect_type(out, "list")
  expect_named(out, c("RescaledPredictors", "BetaCoefficients"))

  expect_s4_class(out$RescaledPredictors, "SpatRaster")
  expect_true(is.data.frame(out$BetaCoefficients))

  # -------------------------------
  # Coefficient table invariants
  # -------------------------------
  coef_tab <- out$BetaCoefficients

  expect_true(all(c(
    "Variable",
    "Coefficient",
    "BetaCoefficient"
  ) %in% names(coef_tab)))

  # Intercept removed
  expect_false("(Intercept)" %in% coef_tab$Variable)

  # Correct number of rows
  expect_equal(nrow(coef_tab), 2)

  # Beta coefficients are finite (can be zero!)
  expect_true(all(is.finite(coef_tab$BetaCoefficient)))

  # -------------------------------
  # Raster invariants
  # -------------------------------
  rescaled <- out$RescaledPredictors

  expect_equal(
    sort(names(rescaled)),
    sort(coef_tab$Variable)
  )

  # Each layer should be numeric and finite
  for (i in seq_len(terra::nlyr(rescaled))) {
    vals <- terra::values(rescaled[[i]])
    expect_true(all(is.finite(vals[!is.na(vals)])))
  }


  # Raster layer names and coefficient variables must match exactly
  expect_setequal(
    names(rescaled),
    coef_tab$Variable
  )

  # Names must be unique
  expect_equal(
    length(names(rescaled)),
    length(unique(names(rescaled)))
  )

})
