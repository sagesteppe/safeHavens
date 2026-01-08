test_that("writeSDMresults writes expected outputs when all inputs supplied", {

  skip_if_not_installed("terra")

  td <- tempfile()
  dir.create(td)

  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  taxon <- "TestTaxon"

  # ---- mock inputs ----
  cv_model <- list(a = 1, b = 2)

  pcnm <- terra::rast(
    ncols = 5, nrows = 5,
    xmin = 0, xmax = 5,
    ymin = 0, ymax = 5
  )
  terra::values(pcnm)  <- 1

  model <- list(model = "glmnet")

  cm <- list(
    table = matrix(c(10, 2, 3, 15), nrow = 2),
    byClass = c(Sensitivity = 0.8, Specificity = 0.9)
  )

  coef_tab <- data.frame(
    term = c("x1", "x2"),
    estimate = c(0.3, -0.1)
  )

  f_rasts <- terra::rast(
    ncols = 5, nrows = 5,
    xmin = 0, xmax = 5,
    ymin = 0, ymax = 5
  )
  terra::values(f_rasts)  <- 1

  thresh <- matrix(c(0.4, 0.6), nrow = 1)
  colnames(thresh) <- c("spec_sens", "maxKappa")

  # ---- run ----
  expect_no_error(
    writeSDMresults(
      path = td,
      taxon = taxon,
      cv_model = cv_model,
      pcnm = pcnm,
      model = model,
      cm = cm,
      coef_tab = coef_tab,
      f_rasts = f_rasts,
      thresh = thresh
    )
  )

  # ---- expectations ----
  expect_true(dir.exists(file.path(td, "Fitting")))
  expect_true(dir.exists(file.path(td, "PCNM")))
  expect_true(dir.exists(file.path(td, "Model")))
  expect_true(dir.exists(file.path(td, "Evaluation")))
  expect_true(dir.exists(file.path(td, "Threshold")))
  expect_true(dir.exists(file.path(td, "Raster")))

  expect_true(file.exists(file.path(td, "Fitting",  "TestTaxon-Fitting.rds")))
  expect_true(file.exists(file.path(td, "PCNM",     "TestTaxon-PCNM.tif")))
  expect_true(file.exists(file.path(td, "Model",    "TestTaxon-Model.rds")))
  expect_true(file.exists(file.path(td, "Model",    "TestTaxon-Coefficients.csv")))
  expect_true(file.exists(file.path(td, "Evaluation", "TestTaxon-CMatrix.csv")))
  expect_true(file.exists(file.path(td, "Evaluation", "TestTaxon-Metrics.csv")))
  expect_true(file.exists(file.path(td, "Threshold",  "TestTaxon-Threshold.csv")))
  expect_true(file.exists(file.path(td, "Raster",     "TestTaxon.tif")))
})


test_that("writeSDMresults runs when optional inputs are missing", {

  td <- withr::local_tempdir()

  expect_no_error(
    writeSDMresults(
      path = td,
      taxon = "TestTaxon"
    )
  )

  # nothing should be created
  expect_false(dir.exists(file.path(td, "Model")))
  expect_false(dir.exists(file.path(td, "Raster")))
})



