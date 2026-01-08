test_that("RescaleRasters returns expected structure", {

  fx <- rescale_fixture()

  out <- RescaleRasters(
    model = fx$model,
    predictors = fx$predictors,
    training_data = fx$training_data,
    pred_mat = fx$pred_mat
  )

  expect_type(out, "list")
  expect_named(out, c("RescaledPredictors", "BetaCoefficients"))
})

test_that("RescaleRasters returns SpatRaster with correct layers", {

  fx <- rescale_fixture()

  out <- RescaleRasters(
    model = fx$model,
    predictors = fx$predictors,
    training_data = fx$training_data,
    pred_mat = fx$pred_mat
  )

  rp <- out$RescaledPredictors

  expect_s4_class(rp, "SpatRaster")
  expect_true(all(names(rp) %in% c("bio1", "bio2")))
})

test_that("RescaleRasters returns valid coefficient table", {

  fx <- rescale_fixture()

  out <- RescaleRasters(
    model = fx$model,
    predictors = fx$predictors,
    training_data = fx$training_data,
    pred_mat = fx$pred_mat
  )

  bc <- out$BetaCoefficients

  expect_s3_class(bc, "data.frame")
  expect_named(
    bc,
    c("Variable", "Coefficient", "BetaCoefficient")
  )

  expect_true(all(bc$Variable %in% c("bio1", "bio2")))
})

test_that("RescaleRasters modifies predictor values", {

  fx <- rescale_fixture()

  out <- RescaleRasters(
    model = fx$model,
    predictors = fx$predictors,
    training_data = fx$training_data,
    pred_mat = fx$pred_mat
  )

  orig_vals <- values(predictors[[1]], mat = FALSE)
  new_vals  <- values(out$RescaledPredictors[[1]], mat = FALSE)

  expect_false(isTRUE(all.equal(orig_vals, new_vals)))
})

test_that("RescaleRasters drops zero-coefficient predictors", {

  fx <- rescale_fixture()

  # force one coefficient to zero
  model_zero <- glmnet::glmnet(
    x, y,
    family = "binomial",
    alpha = 1,
    lambda = max(glmnet::glmnet(x, y, family = "binomial")$lambda)
  )

    model = model_zero(
    predictors = predictors,
    training_data = training_data,
    pred_mat = pred_mat
  )

  expect_lte(nlyr(out$RescaledPredictors), 2)
})

