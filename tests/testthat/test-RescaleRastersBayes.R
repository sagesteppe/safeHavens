test_that("RescaleRasters_bayes rescales predictors correctly (mean)", {

  skip_if_not_installed("terra")

  fake_coef <- data.frame(
    Variable   = c("var1", "var2"),
    Estimate   = c(0.5, -1),
    Est.Error  = c(0.1, 0.2),
    Q2.5       = c(0.3, -1.4),
    Q97.5      = c(0.7, -0.6),
    BetaWeight = c(0.5, -1)
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r1 <- terra::rast(nrows = 2, ncols = 2, vals = c(1,2,3,4))
  r2 <- terra::rast(nrows = 2, ncols = 2, vals = c(2,4,6,8))
  predictors <- c(r1, r2)
  names(predictors) <- c("var1", "var2")

  pred_mat <- data.frame(
    var1 = c(1,2,3,4),
    var2 = c(2,4,6,8)
  )

  out <- RescaleRasters_bayes(
    model = NULL,
    predictors = predictors,
    training_data = NULL,
    pred_mat = pred_mat
  )

  expect_true(inherits(out$RescaledPredictors, "SpatRaster"))
  expect_equal(names(out$RescaledPredictors), c("var1", "var2"))

  # manual check for first value
  mu1 <- mean(pred_mat$var1)
  sd1 <- sqrt(mean((pred_mat$var1 - mu1)^2))
  expected <- ((1 - mu1)/sd1) * abs(0.5)

  expect_equal(
    terra::values(out$RescaledPredictors[[1]])[1],
    expected,
    tolerance = 1e-6
  )
})

test_that("Zero SD variable produces zero raster and warns", {

  fake_coef <- data.frame(
    Variable   = "var1",
    Estimate   = 1,
    Est.Error  = 0.1,
    Q2.5       = 0.5,
    Q97.5      = 1.5,
    BetaWeight = 1
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 2, vals = 5)
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(5,5,5,5))

  expect_warning(
    out <- RescaleRasters_bayes(NULL, r, NULL, pred_mat),
    "zero SD"
  )

  expect_true(all(terra::values(out$RescaledPredictors) == 0))
})

test_that("Errors if raster missing model predictor", {

  fake_coef <- data.frame(
    Variable   = "var1",
    Estimate   = 1,
    Est.Error  = 0.1,
    Q2.5       = 0.5,
    Q97.5      = 1.5,
    BetaWeight = 1
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 2, vals = 1)
  names(r) <- "not_var1"

  pred_mat <- data.frame(var1 = c(1,2,3,4))

  expect_error(
    RescaleRasters_bayes(NULL, r, NULL, pred_mat),
    "not in `predictors` raster"
  )
})

test_that("beta_summary switches correctly", {

  local_mocked_bindings(
    extract_posterior_betas = function(model, beta_summary) {
      data.frame(
        Variable   = "var1",
        Estimate   = if (beta_summary == "Q2.5") 10 else 2,
        Est.Error  = 0.1,
        Q2.5       = 10,
        Q97.5      = 20,
        BetaWeight = if (beta_summary == "Q2.5") 10 else 2
      )
    }
  )

  r <- terra::rast(nrows = 2, ncols = 1, vals = c(1, 2))
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(1, 2))

  out <- RescaleRasters_bayes(
    NULL, r, NULL, pred_mat,
    beta_summary = "Q2.5"
  )

  expect_equal(out$BetaCoefficients$Q2.5, 10)
})
test_that("Uncertainty layer is added when requested", {

  fake_coef <- data.frame(
    Variable   = "var1",
    Estimate   = 1,
    Est.Error  = 0.5,
    Q2.5       = 0.2,
    Q97.5      = 1.5,
    BetaWeight = 1
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 2, vals = c(1,2,3,4))
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(1,2,3,4))

  out <- RescaleRasters_bayes(
    NULL,
    r,
    NULL,
    pred_mat,
    include_uncertainty = TRUE
  )

  expect_true("coef_uncertainty" %in%
                names(out$RescaledPredictors))

  expect_true(inherits(out$UncertaintyLayer, "SpatRaster"))
})

test_that("Intercept and GP terms are removed", {

  fake_coef <- data.frame(
    Variable   = c("Intercept", "sgp(x)", "var1"),
    Estimate   = c(0, 0, 1),
    Est.Error  = c(0, 0, 0.1),
    Q2.5       = c(0, 0, 0.5),
    Q97.5      = c(0, 0, 1.5),
    BetaWeight = c(0, 0, 1)
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 1, vals = c(1,2))
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(1,2))

  out <- RescaleRasters_bayes(NULL, r, NULL, pred_mat)

  expect_equal(names(out$RescaledPredictors), "var1")
})

test_that("Errors if no environmental predictors found", {

  fake_coef <- data.frame(
    Variable   = "Intercept",
    Estimate   = 1,
    Est.Error  = 0.1,
    Q2.5       = 0.5,
    Q97.5      = 1.5,
    BetaWeight = 1
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 1, vals = c(1,2))
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(1,2))

  expect_error(
    RescaleRasters_bayes(NULL, r, NULL, pred_mat),
    "No environmental predictor betas found"
  )
})

test_that("Uncertainty scaling handles zero span", {

  fake_coef <- data.frame(
    Variable   = "var1",
    Estimate   = 1,
    Est.Error  = 0,  # zero uncertainty
    Q2.5       = 0.5,
    Q97.5      = 1.5,
    BetaWeight = 1
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 1, vals = c(1,2))
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(1,2))

  out <- RescaleRasters_bayes(
    NULL, r, NULL, pred_mat,
    include_uncertainty = TRUE
  )

  expect_true(inherits(out$UncertaintyLayer, "SpatRaster"))
})

test_that("gp_x and gp_y are ignored", {

  fake_coef <- data.frame(
    Variable   = c("var1", "gp_x"),
    Estimate   = c(1, 999),
    Est.Error  = c(0.1, 0.1),
    Q2.5       = c(0.5, 0.5),
    Q97.5      = c(1.5, 1.5),
    BetaWeight = c(1, 999)
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 1, vals = c(1,2))
  names(r) <- "var1"

  pred_mat <- data.frame(
    var1 = c(1,2),
    gp_x = c(10,20)
  )

  out <- RescaleRasters_bayes(NULL, r, NULL, pred_mat)

  expect_equal(names(out$RescaledPredictors), "var1")
})

test_that("Errors listing multiple missing raster predictors", {

  fake_coef <- data.frame(
    Variable   = c("var1", "var2"),
    Estimate   = c(1, 2),
    Est.Error  = c(0.1, 0.1),
    Q2.5       = c(0.5, 1.5),
    Q97.5      = c(1.5, 2.5),
    BetaWeight = c(1, 2)
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 1, vals = c(1,2))
  names(r) <- "var1"   # var2 missing

  pred_mat <- data.frame(
    var1 = c(1,2),
    var2 = c(3,4)
  )

  expect_error(
    RescaleRasters_bayes(NULL, r, NULL, pred_mat),
    "var2"
  )
})

test_that("BetaCoefficients returns expected columns only", {

  fake_coef <- data.frame(
    Variable   = "var1",
    Estimate   = 1,
    Est.Error  = 0.1,
    Q2.5       = 0.5,
    Q97.5      = 1.5,
    BetaWeight = 999
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 1, vals = c(1,2))
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(1,2))

  out <- RescaleRasters_bayes(NULL, r, NULL, pred_mat)

  expect_equal(
    colnames(out$BetaCoefficients),
    c("Variable", "Estimate", "Est.Error", "Q2.5", "Q97.5")
  )
})

test_that("Predictor order is preserved", {

  fake_coef <- data.frame(
    Variable   = c("varB", "varA"),
    Estimate   = c(2, 1),
    Est.Error  = c(0.1, 0.1),
    Q2.5       = c(1.5, 0.5),
    Q97.5      = c(2.5, 1.5),
    BetaWeight = c(2, 1)
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r1 <- terra::rast(nrows = 2, ncols = 1, vals = c(1,2))
  r2 <- terra::rast(nrows = 2, ncols = 1, vals = c(3,4))
  predictors <- c(r2, r1)
  names(predictors) <- c("varB", "varA")

  pred_mat <- data.frame(
    varB = c(1,2),
    varA = c(3,4)
  )

  out <- RescaleRasters_bayes(NULL, predictors, NULL, pred_mat)

  expect_equal(names(out$RescaledPredictors), c("varB", "varA"))
})

test_that("Return structure is stable", {

  fake_coef <- data.frame(
    Variable   = "var1",
    Estimate   = 1,
    Est.Error  = 0.1,
    Q2.5       = 0.5,
    Q97.5      = 1.5,
    BetaWeight = 1
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 2, ncols = 1, vals = c(1,2))
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(1,2))

  out <- RescaleRasters_bayes(NULL, r, NULL, pred_mat)

  expect_named(
    out,
    c("RescaledPredictors", "BetaCoefficients", "UncertaintyLayer")
  )
})

test_that("uncertainty_wt scales uncertainty layer", {

  fake_coef <- data.frame(
    Variable   = "var1",
    Estimate   = 1,
    Est.Error  = 0.5,
    Q2.5       = 0.5,
    Q97.5      = 1.5,
    BetaWeight = 1
  )

  local_mocked_bindings(
    extract_posterior_betas = function(...) fake_coef
  )

  r <- terra::rast(nrows = 3, ncols = 1, vals = c(1,2,4))
  names(r) <- "var1"

  pred_mat <- data.frame(var1 = c(1,2,4))

  out1 <- RescaleRasters_bayes(
    NULL, r, NULL, pred_mat,
    include_uncertainty = TRUE,
    uncertainty_wt = 1
  )

  out2 <- RescaleRasters_bayes(
    NULL, r, NULL, pred_mat,
    include_uncertainty = TRUE,
    uncertainty_wt = 3
  )

  max1 <- terra::global(out1$UncertaintyLayer, "max", na.rm = TRUE)[1,1]
  max2 <- terra::global(out2$UncertaintyLayer, "max", na.rm = TRUE)[1,1]

  expect_gt(max2, max1)
})
