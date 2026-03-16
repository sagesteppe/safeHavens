library(testthat)

# =============================================================================
# Shared fixtures
# =============================================================================

make_predictors <- function(n_vars = 3, nrows = 10, ncols = 10) {
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = nrows, ncols = ncols, crs = "EPSG:4326"
  )
  layers <- lapply(seq_len(n_vars), function(i) {
    p <- r; terra::values(p) <- rnorm(terra::ncell(r)); p
  })
  stk <- do.call(c, layers)
  names(stk) <- paste0("bio", seq_len(n_vars))
  stk
}

make_points <- function(n = 80, with_occurrence = TRUE) {
  set.seed(42)
  coords <- data.frame(
    lon = runif(n, -119, -111),
    lat = runif(n, 31, 39)
  )
  sf_pts <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  if (with_occurrence) sf_pts$occurrence <- sample(0:1, n, replace = TRUE)
  sf_pts
}

make_mock_brmsfit <- function() {
  structure(list(), class = "brmsfit")
}

make_mock_loo <- function() {
  structure(
    list(
      estimates = matrix(c(-50, 2), nrow = 1,
                         dimnames = list("elpd_loo", c("Estimate", "SE"))),
      pointwise = matrix(0, nrow = 10, ncol = 1)
    ),
    class = c("psis_loo", "loo")
  )
}

make_mock_aoa <- function(r_template) {
  aoa <- r_template; terra::values(aoa) <- sample(0:1, terra::ncell(r_template), TRUE)
  di  <- r_template; terra::values(di)  <- runif(terra::ncell(r_template))
  lpd <- r_template; terra::values(lpd) <- runif(terra::ncell(r_template), 0, 100)
  structure(
    list(AOA = aoa, DI = di, LPD = lpd, parameters = list(threshold = 0.5)),
    class = "aoa"
  )
}

stub_epred <- function(object, newdata, ndraws = 200, ...) {
  n      <- nrow(newdata)
  ndraws <- max(ndraws %||% 200L, 1L)
  matrix(runif(ndraws * n, 0.3, 0.7), nrow = ndraws)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Standard brms + CAST mocks for a two-predictor (bio1, bio2) run.
# Wraps expr with local_mocked_bindings; caller supplies the expression.
with_standard_mocks <- function(expr, fixef_names = c("Intercept", "bio1", "bio2"),
                                fixef_vals = c(-0.1, 0.5, 0.3)) {
  fixef_mat <- cbind(
    Estimate    = fixef_vals,
    `Est.Error` = rep(0.1, length(fixef_vals))
  )
  rownames(fixef_mat) <- fixef_names

  local_mocked_bindings(
    brm               = function(...) make_mock_brmsfit(),
    loo               = function(...) make_mock_loo(),
    posterior_epred   = stub_epred,
    posterior_summary = function(...) matrix(0, 2, 4),
    rhat              = function(...) c(a = 1.001),
    neff_ratio        = function(...) c(a = 0.6),
    as_draws_df       = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params       = function(...) data.frame(Parameter = rep("divergent__", 200),
                                                 Value     = rep(0, 200)),
    fixef             = function(object, summary = TRUE, ...) fixef_mat,
    .package = "brms"
  )
  local_mocked_bindings(
    aoa    = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm  = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    .package = "CAST"
  )
  expr
}


# =============================================================================
# 1. Pure-logic helpers (no brms / Stan needed)
# =============================================================================

test_that("run_pca_on_raster handles edge cases", {
  r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  pred1 <- r; terra::values(pred1) <- 1:25
  pred2 <- r; terra::values(pred2) <- 25:1
  predictors <- c(pred1, pred2)
  names(predictors) <- c("a", "b")

  result <- run_pca_on_raster(predictors, pca_axes = 5)

  expect_equal(length(result$pc_names), 2)   # capped at 2 available
  expect_s4_class(result$raster, "SpatRaster")
  expect_s3_class(result$pca_model, "prcomp")
})

test_that("build_priors creates valid prior objects", {
  hs_prior <- build_priors(
    prior_type = "horseshoe", prior_scale = 1,
    pred_names = c("var1", "var2", "var3"), n_obs = 100, p0 = 2, gp_prior = NULL
  )
  expect_s3_class(hs_prior, "brmsprior")
  expect_true(any(grepl("horseshoe", hs_prior$prior)))

  norm_prior <- build_priors(
    prior_type = "normal", prior_scale = 2,
    pred_names = c("var1", "var2"), n_obs = 50, p0 = NULL, gp_prior = NULL
  )
  expect_s3_class(norm_prior, "brmsprior")
  expect_true(any(grepl("normal", norm_prior$prior)))

  student_prior <- build_priors(
    prior_type = "student", prior_scale = 1.5,
    pred_names = "var1", n_obs = 75, p0 = NULL, gp_prior = NULL
  )
  expect_s3_class(student_prior, "brmsprior")
  expect_true(any(grepl("student_t", student_prior$prior)))
})

test_that("add_planar_coords attaches gp_x and gp_y in km", {
  pts    <- data.frame(x = c(-120, -115), y = c(35, 38))
  sf_pts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = 4326)

  result <- add_planar_coords(sf_pts, 3857, scale_km = TRUE)

  expect_true("gp_x" %in% names(result))
  expect_true("gp_y" %in% names(result))
  expect_true(all(abs(result$gp_x) < 20000))   # km, not metres
})

test_that("add_coord_weights_to_points attaches weighted coordinate columns", {
  sample_pts <- terra::vect(
    cbind(x = c(-120, -115, -110), y = c(35, 37, 39)),
    crs = "EPSG:4326"
  )
  pt_env <- data.frame(PC1 = c(0.5, 1.2, -0.3), PC2 = c(0.1, 0.8, 0.4))

  result <- add_coord_weights_to_points(
    sample_pts, pt_env, coord_wt = 2.0, env_vars = c("PC1", "PC2")
  )

  expect_true("coord_x_w" %in% names(result))
  expect_true("coord_y_w" %in% names(result))
  expect_equal(nrow(result), 3)
})

test_that("prediction raster coordinate construction produces km-scale values", {
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 5, ncols = 5, crs = "EPSG:4326"
  )
  terra::values(r) <- 1:25

  coords_df   <- data.frame(
    x = rep(seq(-120, -110, length.out = 5), 5),
    y = rep(seq(30, 40, length.out = 5), each = 5)
  )
  coords_vect <- terra::vect(coords_df, geom = c("x", "y"), crs = "EPSG:4326")
  coords_proj <- terra::project(coords_vect, "EPSG:3857")
  coords_km   <- terra::crds(coords_proj) / 1000

  expect_equal(nrow(coords_km), 25)
  expect_equal(ncol(coords_km), 2)
  expect_true(all(abs(coords_km) < 50000))
})

test_that("AOA coverage statistic is a valid probability", {
  mock_aoa <- terra::rast(nrows = 10, ncols = 10)
  terra::values(mock_aoa) <- sample(0:1, 100, replace = TRUE, prob = c(0.3, 0.7))

  aoa_coverage <- sum(terra::values(mock_aoa), na.rm = TRUE) /
                  sum(!is.na(terra::values(mock_aoa)))

  expect_true(aoa_coverage >= 0 && aoa_coverage <= 1)
  expect_type(aoa_coverage, "double")
})

test_that("feature selection min_ffs_var guard resets to 2 when >= nlyr", {
  pred_names  <- c("PC1", "PC2", "PC3")
  min_ffs_var <- 10
  nlyr        <- 3
  if (min_ffs_var > nlyr) min_ffs_var <- 2
  expect_equal(min_ffs_var, 2)
})

test_that("compute_top3_clusters replicates rank1 when only one cluster seen", {
  n_pts   <- 10
  n_draws <- 20
  dc      <- matrix(5, nrow = n_pts, ncol = n_draws)

  result <- compute_top3_clusters(dc, n_pts, n_draws)

  expect_equal(result$labels_matrix[1, 1], 5)
  expect_equal(result$labels_matrix[1, 2], 5)
  expect_equal(result$labels_matrix[1, 3], 5)
  expect_equal(result$pct_matrix[1, 1], 100)
  expect_equal(result$pct_matrix[1, 2], 100)
  expect_equal(result$pct_matrix[1, 3], 100)
})

test_that("compute_top3_clusters handles mixed assignments", {
  n_pts   <- 5
  n_draws <- 100
  dc      <- matrix(NA_integer_, nrow = n_pts, ncol = n_draws)
  dc[1, ] <- sample(c(1, 2, 3), n_draws, replace = TRUE, prob = c(0.7, 0.2, 0.1))
  dc[2, ] <- rep(5, n_draws)

  result <- compute_top3_clusters(dc, n_pts, n_draws)

  expect_equal(length(unique(result$labels_matrix[1, ])), 3)
  expect_true(all(result$labels_matrix[2, ] == 5))
  expect_equal(result$pct_matrix[2, 1], 100)
})

test_that("compute_aoa_bayes correctly filters smooth and intercept terms", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  set.seed(60)
  n       <- 30
  coords  <- data.frame(lon = runif(n, -120, -115), lat = runif(n, 35, 40))
  train_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  train_sf$occurrence <- sample(0:1, n, replace = TRUE)
  train_sf$bio1 <- rnorm(n)
  train_sf$gp_x <- runif(n)
  train_sf$gp_y <- runif(n)

  r     <- terra::rast(xmin = -120, xmax = -115, ymin = 35, ymax = 40,
                       nrows = 4, ncols = 4, crs = "EPSG:4326")
  bio1r <- r; terra::values(bio1r) <- rnorm(16)
  predictors <- c(bio1r); names(predictors) <- "bio1"

  mock_fit <- structure(list(), class = "brmsfit")
  rn       <- c("Intercept", "bio1", "s(gp_x,gp_y)_1", "sgp_x_2", "sgp_y_3")
  mock_fixef <- cbind(
    Estimate    = c(-0.2, 0.9, 0.01, 0.02, -0.01),
    `Est.Error` = rep(0.1, 5)
  )
  rownames(mock_fixef) <- rn

  captured_vars <- NULL
  mock_aoa_result <- list(AOA = r, DI = r, LPD = r, parameters = list(threshold = 0.5))
  terra::values(mock_aoa_result$AOA) <- sample(0:1, 16, replace = TRUE)
  terra::values(mock_aoa_result$DI)  <- runif(16)
  terra::values(mock_aoa_result$LPD) <- runif(16)
  class(mock_aoa_result) <- "aoa"

  local_mocked_bindings(
    fixef = function(object, summary = TRUE, ...) mock_fixef,
    .package = "brms"
  )
  local_mocked_bindings(
    aoa = function(...) {
      captured_vars <<- list(...)$variables
      mock_aoa_result
    },
    .package = "CAST"
  )

  compute_aoa_bayes(
    model              = mock_fit,
    train              = train_sf,
    predictors         = predictors,
    cv_folds           = list(indx_train = list(1:24), indx_test = list(25:30)),
    use_posterior_mean = TRUE,
    LPD                = TRUE
  )

  expect_equal(captured_vars, "bio1")
})

test_that("compute_aoa_bayes errors when no env vars match raster", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  n       <- 20
  coords  <- data.frame(lon = runif(n, -120, -115), lat = runif(n, 35, 40))
  train_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  train_sf$occurrence <- sample(0:1, n, replace = TRUE)
  train_sf$gp_x <- runif(n)
  train_sf$gp_y <- runif(n)

  r     <- terra::rast(xmin = -120, xmax = -115, ymin = 35, ymax = 40,
                       nrows = 3, ncols = 3, crs = "EPSG:4326")
  bio1r <- r; terra::values(bio1r) <- rnorm(9)
  predictors <- c(bio1r); names(predictors) <- "bio1"

  mock_fit  <- structure(list(), class = "brmsfit")
  rn        <- c("Intercept", "sgp_x_1")
  mock_fixef <- cbind(Estimate = c(-0.1, 0.01))
  rownames(mock_fixef) <- rn

  local_mocked_bindings(
    fixef = function(object, summary = TRUE, ...) mock_fixef,
    .package = "brms"
  )

  expect_error(
    compute_aoa_bayes(
      model      = mock_fit,
      train      = train_sf,
      predictors = predictors,
      cv_folds   = list(indx_train = list(1:16), indx_test = list(17:20)),
      use_posterior_mean = TRUE
    ),
    "No environmental predictors matched"
  )
})

test_that("compute_aoa_bayes uses equal weights when use_posterior_mean = FALSE", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  set.seed(70)
  n       <- 30
  coords  <- data.frame(lon = runif(n, -120, -115), lat = runif(n, 35, 40))
  train_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  train_sf$occurrence <- sample(0:1, n, replace = TRUE)
  train_sf$PC1 <- rnorm(n); train_sf$PC2 <- rnorm(n)
  train_sf$gp_x <- runif(n); train_sf$gp_y <- runif(n)

  r     <- terra::rast(xmin = -120, xmax = -115, ymin = 35, ymax = 40,
                       nrows = 3, ncols = 3, crs = "EPSG:4326")
  PC1r  <- r; terra::values(PC1r) <- rnorm(9)
  PC2r  <- r; terra::values(PC2r) <- rnorm(9)
  predictors <- c(PC1r, PC2r); names(predictors) <- c("PC1", "PC2")

  mock_fit  <- structure(list(), class = "brmsfit")
  mock_fixef <- cbind(
    Estimate    = c(-0.2, 0.9, -0.4),
    `Est.Error` = c(0.1, 0.1, 0.1)
  )
  rownames(mock_fixef) <- c("Intercept", "PC1", "PC2")

  captured_weight <- NULL
  mock_aoa_result <- list(AOA = r, DI = r, LPD = r, parameters = list(threshold = 0.5))
  terra::values(mock_aoa_result$AOA) <- sample(0:1, 9, replace = TRUE)
  terra::values(mock_aoa_result$DI)  <- runif(9)
  terra::values(mock_aoa_result$LPD) <- runif(9)
  class(mock_aoa_result) <- "aoa"

  local_mocked_bindings(
    fixef = function(object, summary = TRUE, ...) mock_fixef,
    .package = "brms"
  )
  local_mocked_bindings(
    aoa = function(...) {
      captured_weight <<- list(...)$weight
      mock_aoa_result
    },
    .package = "CAST"
  )

  compute_aoa_bayes(
    model              = mock_fit,
    train              = train_sf,
    predictors         = predictors,
    cv_folds           = list(indx_train = list(1:24), indx_test = list(25:30)),
    use_posterior_mean = FALSE,
    LPD                = TRUE
  )

  expect_true(!is.null(captured_weight))
  expect_equal(as.numeric(captured_weight), c(1, 1))
})


# =============================================================================
# 2. check_convergence (mocked brms internals)
# =============================================================================

test_that("check_convergence returns expected structure", {
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")

  mock_fit <- structure(list(), class = "brmsfit")
  set.seed(1)
  n_draws     <- 2000
  mock_draws  <- posterior::as_draws_df(
    list(b_Intercept = rnorm(n_draws, 0, 0.1), b_PC1 = rnorm(n_draws, 0.5, 0.05))
  )
  mock_nuts   <- data.frame(
    Parameter = rep(c("divergent__", "other__"), each = n_draws),
    Value     = c(rep(0, n_draws), rep(1, n_draws))
  )

  local_mocked_bindings(
    posterior_summary = function(...) matrix(0, nrow = 2, ncol = 4,
                                              dimnames = list(c("b_Intercept", "b_PC1"), NULL)),
    rhat              = function(...) c(b_Intercept = 1.002, b_PC1 = 1.003),
    neff_ratio        = function(...) c(b_Intercept = 0.6, b_PC1 = 0.55),
    as_draws_df       = function(...) mock_draws,
    nuts_params       = function(...) mock_nuts,
    .package = "brms"
  )

  result <- check_convergence(mock_fit)

  expect_type(result, "list")
  expect_named(result, c("max_Rhat", "min_BulkESS", "min_TailESS", "n_divergent"))
  expect_equal(result$max_Rhat,    1.003)
  expect_equal(result$n_divergent, 0)
})

test_that("check_convergence flags high Rhat", {
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")

  mock_fit  <- structure(list(), class = "brmsfit")
  mock_nuts <- data.frame(Parameter = rep("divergent__", 200), Value = rep(0, 200))

  local_mocked_bindings(
    posterior_summary = function(...) matrix(0, nrow = 2, ncol = 4),
    rhat              = function(...) c(b_Intercept = 1.10, b_PC1 = 1.08),
    neff_ratio        = function(...) c(b_Intercept = 0.1, b_PC1 = 0.08),
    as_draws_df       = function(...) posterior::as_draws_df(list(b_Intercept = rnorm(200))),
    nuts_params       = function(...) mock_nuts,
    .package = "brms"
  )

  result <- check_convergence(mock_fit)
  expect_gt(result$max_Rhat, 1.05)
})

test_that("check_convergence counts divergent transitions", {
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")

  mock_fit  <- structure(list(), class = "brmsfit")
  n_draws   <- 2000
  mock_nuts <- data.frame(
    Parameter = rep("divergent__", n_draws),
    Value     = c(rep(1, 50), rep(0, n_draws - 50))
  )

  local_mocked_bindings(
    posterior_summary = function(...) matrix(0, nrow = 1, ncol = 4),
    rhat              = function(...) c(b_Intercept = 1.001),
    neff_ratio        = function(...) c(b_Intercept = 0.5),
    as_draws_df       = function(...) posterior::as_draws_df(list(b_Intercept = rnorm(n_draws))),
    nuts_params       = function(...) mock_nuts,
    .package = "brms"
  )

  result <- check_convergence(mock_fit)
  expect_equal(result$n_divergent, 50)
})

test_that("check_convergence errors on non-brmsfit input", {
  expect_error(check_convergence(list()))
})


# =============================================================================
# 3. evaluate_bayes_model (mocked posterior_epred)
# =============================================================================

test_that("evaluate_bayes_model returns a confusionMatrix", {
  skip_if_not_installed("brms")
  skip_if_not_installed("caret")

  set.seed(10)
  n_obs   <- 40
  coords  <- data.frame(lon = runif(n_obs, -120, -115), lat = runif(n_obs, 35, 40))
  test_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  test_sf$occurrence <- rep(c(0L, 1L), n_obs / 2)
  test_sf$PC1        <- rnorm(n_obs)

  mock_fit <- structure(list(), class = "brmsfit")

  local_mocked_bindings(
    posterior_epred = function(object, newdata, ...) {
      # Return draws that perfectly match the occurrence column
      occ <- as.integer(newdata$occurrence)
      matrix(rep(ifelse(occ == 1L, 0.9, 0.1), each = 10L), nrow = 10L)
    },
    .package = "brms"
  )

  result <- evaluate_bayes_model(mock_fit, test_sf, pred_names = "PC1")

  expect_s3_class(result, "confusionMatrix")
  expect_equal(as.numeric(result$overall["Accuracy"]), 1.0)
})

# =============================================================================
# 4. create_bayes_spatial_predictions (mocked posterior_epred)
# =============================================================================

test_that("create_bayes_spatial_predictions returns mean, sd SpatRasters + pred_matrix", {
  skip_if_not_installed("brms")

  r   <- terra::rast(xmin = -120, xmax = -118, ymin = 35, ymax = 37,
                     nrows = 5, ncols = 5, crs = "EPSG:4326")
  PC1 <- r; terra::values(PC1) <- rnorm(terra::ncell(r))
  PC2 <- r; terra::values(PC2) <- rnorm(terra::ncell(r))
  predictors <- c(PC1, PC2); names(predictors) <- c("PC1", "PC2")

  mock_fit <- structure(list(), class = "brmsfit")

  local_mocked_bindings(
    posterior_epred = function(object, newdata, ndraws = 200, ...) {
      matrix(runif(ndraws * nrow(newdata), 0.3, 0.7), nrow = ndraws)
    },
    .package = "brms"
  )

  result <- create_bayes_spatial_predictions(
    fit               = mock_fit,
    predictors        = predictors,
    pred_names        = c("PC1", "PC2"),
    planar_projection = "EPSG:3857",
    iter              = 400
  )

  expect_named(result, c("mean", "sd", "pred_matrix"))
  expect_s4_class(result$mean, "SpatRaster")
  expect_s4_class(result$sd,   "SpatRaster")
  expect_equal(names(result$mean), "occurrence_prob_mean")
  expect_equal(names(result$sd),   "occurrence_prob_sd")
  expect_true(all(terra::values(result$mean)[!is.na(terra::values(result$mean))] >= 0))
  expect_true(all(terra::values(result$mean)[!is.na(terra::values(result$mean))] <= 1))
  expect_true(all(terra::values(result$sd)[!is.na(terra::values(result$sd))] >= 0))
  expect_true("gp_x" %in% names(result$pred_matrix))
  expect_true("gp_y" %in% names(result$pred_matrix))
})

test_that("create_bayes_spatial_predictions pred_matrix row count matches raster cells", {
  skip_if_not_installed("brms")

  r   <- terra::rast(xmin = -120, xmax = -118, ymin = 35, ymax = 37,
                     nrows = 4, ncols = 4, crs = "EPSG:4326")
  PC1 <- r; terra::values(PC1) <- runif(terra::ncell(r))
  predictors <- c(PC1); names(predictors) <- "PC1"

  mock_fit <- structure(list(), class = "brmsfit")

  local_mocked_bindings(
    posterior_epred = function(object, newdata, ndraws = 100, ...) {
      matrix(runif(ndraws * nrow(newdata), 0.2, 0.8), nrow = ndraws)
    },
    .package = "brms"
  )

  result <- create_bayes_spatial_predictions(
    fit               = mock_fit,
    predictors        = predictors,
    pred_names        = "PC1",
    planar_projection = "EPSG:3857",
    iter              = 200
  )

  expect_equal(nrow(result$pred_matrix), terra::ncell(r))
})

test_that("create_bayes_spatial_predictions uses ndraws = round(iter/2)", {
  skip_if_not_installed("brms")

  captured_ndraws <- NULL
  r   <- terra::rast(xmin = -120, xmax = -119, ymin = 35, ymax = 36,
                     nrows = 3, ncols = 3, crs = "EPSG:4326")
  PC1 <- r; terra::values(PC1) <- runif(9)
  predictors <- c(PC1); names(predictors) <- "PC1"
  mock_fit   <- structure(list(), class = "brmsfit")

  local_mocked_bindings(
    posterior_epred = function(object, newdata, ndraws = NULL, ...) {
      captured_ndraws <<- ndraws
      matrix(runif(ndraws * nrow(newdata)), nrow = ndraws)
    },
    .package = "brms"
  )

  create_bayes_spatial_predictions(
    fit = mock_fit, predictors = predictors,
    pred_names = "PC1", planar_projection = "EPSG:3857", iter = 600
  )

  expect_equal(captured_ndraws, 300L)
})


# =============================================================================
# 5. bayesianSDM branch tests (mocked brm + CAST)
# =============================================================================

test_that("PCAModel is NULL when pca_predictors = FALSE", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(2)

  with_standard_mocks({
    result <- bayesianSDM(
      x = pts, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
      chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 1
    )
    expect_null(result$PCAModel)
    expect_s4_class(result$Predictors, "SpatRaster")
  })
})

test_that("PCAModel is prcomp when pca_predictors = TRUE", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(4)

  local_mocked_bindings(
    brm = function(...) make_mock_brmsfit(), loo = function(...) make_mock_loo(),
    posterior_epred = stub_epred, posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.001), neff_ratio = function(...) c(a = 0.6),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.2, 0.8, -0.4), `Est.Error` = rep(0.1, 3))
      rownames(m) <- c("Intercept", "PC1", "PC2"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    .package = "CAST"
  )

  result <- bayesianSDM(
    x = pts, predictors = preds, planar_projection = 3857,
    pca_predictors = TRUE, pca_axes = 2,
    feature_selection = "none", vif = FALSE,
    chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 2
  )

  expect_s3_class(result$PCAModel, "prcomp")
  expect_true(all(grepl("^PC", names(result$Predictors))))
  expect_equal(terra::nlyr(result$Predictors), 2L)
})

test_that("min_ffs_var > nlyr is silently capped to 2 inside bayesianSDM", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(2)

  captured_minvar <- NULL

  local_mocked_bindings(
    brm = function(...) make_mock_brmsfit(), loo = function(...) make_mock_loo(),
    posterior_epred = stub_epred, posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.001), neff_ratio = function(...) c(a = 0.6),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.2, 0.8, 0.3), `Est.Error` = rep(0.1, 3))
      rownames(m) <- c("Intercept", "bio1", "bio2"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    ffs = function(predictors, minVar, ...) {
      captured_minvar <<- minVar
      list(selectedvars = colnames(predictors)[seq_len(min(2, ncol(predictors)))])
    },
    .package = "CAST"
  )

  expect_no_error(
    bayesianSDM(
      x = pts, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "ffs",
      min_ffs_var = 10, vif = FALSE,
      chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 3
    )
  )
  expect_equal(captured_minvar, 2L)
})

test_that("pseudo-absences are added when occurrence column is missing", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts_no_occ <- make_points(with_occurrence = FALSE)
  preds      <- make_predictors(2)

  with_standard_mocks({
    result <- bayesianSDM(
      x = pts_no_occ, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
      chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 4,
      fact = 1
    )
    all_occ <- c(result$TrainData$occurrence, result$TestData$occurrence)
    expect_true(0L %in% all_occ)
    expect_true(1L %in% all_occ)
  })
})

test_that("no pseudo-absences added when occurrence column is present", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts          <- make_points(n = 60, with_occurrence = TRUE)
  pts$occurrence <- 1L
  n_input      <- nrow(pts)
  preds        <- make_predictors(2)

  with_standard_mocks({
    result  <- bayesianSDM(
      x = pts, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
      chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 5
    )
    n_total <- nrow(result$TrainData) + nrow(result$TestData)
    expect_lte(n_total, n_input)
  })
})

test_that("feature_selection='ffs' subsets predictors to selected vars", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(4)

  local_mocked_bindings(
    brm = function(...) make_mock_brmsfit(), loo = function(...) make_mock_loo(),
    posterior_epred = stub_epred, posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.001), neff_ratio = function(...) c(a = 0.6),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.2, 0.7, -0.4), `Est.Error` = rep(0.1, 3))
      rownames(m) <- c("Intercept", "bio1", "bio2"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    ffs = function(...) list(selectedvars = c("bio1", "bio2")),
    .package = "CAST"
  )

  result <- bayesianSDM(
    x = pts, predictors = preds, planar_projection = 3857,
    pca_predictors = FALSE, feature_selection = "ffs",
    vif = FALSE, chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 6
  )

  expect_equal(terra::nlyr(result$Predictors), 2L)
  expect_setequal(names(result$Predictors), c("bio1", "bio2"))
})

test_that("warning issued and all predictors kept when ffs returns 0 vars", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(3)

  local_mocked_bindings(
    brm = function(...) make_mock_brmsfit(), loo = function(...) make_mock_loo(),
    posterior_epred = stub_epred, posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.001), neff_ratio = function(...) c(a = 0.6),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.2, 0.8, -0.3, 0.4), `Est.Error` = rep(0.1, 4))
      rownames(m) <- c("Intercept", "bio1", "bio2", "bio3"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    ffs = function(...) list(selectedvars = character(0)),
    .package = "CAST"
  )

  expect_warning(
    result <- bayesianSDM(
      x = pts, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "ffs",
      vif = FALSE, chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 7
    ),
    "Feature selection returned 0 variables"
  )
  expect_equal(terra::nlyr(result$Predictors), 3L)
})

test_that("vif=TRUE calls vifcor and subsets predictors", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")
  skip_if_not_installed("usdm")

  pts          <- make_points()
  preds        <- make_predictors(4)
  vifcor_called <- FALSE

  local_mocked_bindings(
    brm = function(...) make_mock_brmsfit(), loo = function(...) make_mock_loo(),
    posterior_epred = stub_epred, posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.001), neff_ratio = function(...) c(a = 0.6),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.2, 0.7, -0.3, 0.2), `Est.Error` = rep(0.1, 4))
      rownames(m) <- c("Intercept", "bio1", "bio2", "bio3"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    .package = "CAST"
  )
  local_mocked_bindings(
    vifcor = function(x, ...) {
      vifcor_called <<- TRUE
      if (!methods::isClass("MockVIF")) {
        methods::setClass("MockVIF", representation(results = "data.frame"))
      }
      methods::new("MockVIF", results = data.frame(Variables = c("bio1", "bio2", "bio3")))
    },
    .package = "usdm"
  )

  result <- bayesianSDM(
    x = pts, predictors = preds, planar_projection = 3857,
    pca_predictors = FALSE, feature_selection = "none", vif = TRUE,
    chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 8
  )

  expect_true(vifcor_called)
  expect_false("bio4" %in% names(result$Predictors))
  expect_equal(terra::nlyr(result$Predictors), 3L)
})

test_that("convergence warning issued when max Rhat > 1.05", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(2)

  local_mocked_bindings(
    brm = function(...) make_mock_brmsfit(), loo = function(...) make_mock_loo(),
    posterior_epred = stub_epred, posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.12, b = 1.08),
    neff_ratio = function(...) c(a = 0.1),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.1, 0.5, 0.3), `Est.Error` = rep(0.1, 3))
      rownames(m) <- c("Intercept", "bio1", "bio2"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    .package = "CAST"
  )

  expect_warning(
    bayesianSDM(
      x = pts, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
      chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 9
    ),
    "R-hat"
  )
})

test_that("prior_type='normal' is forwarded correctly to brm", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(2)
  brm_prior_captured <- NULL

  local_mocked_bindings(
    brm = function(prior, ...) { brm_prior_captured <<- prior; make_mock_brmsfit() },
    loo = function(...) make_mock_loo(), posterior_epred = stub_epred,
    posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.001), neff_ratio = function(...) c(a = 0.6),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.1, 0.5, 0.3), `Est.Error` = rep(0.1, 3))
      rownames(m) <- c("Intercept", "bio1", "bio2"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    .package = "CAST"
  )

  bayesianSDM(
    x = pts, predictors = preds, planar_projection = 3857,
    pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
    prior_type = "normal", prior_scale = 2,
    chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 10
  )

  expect_true(any(grepl("normal", brm_prior_captured$prior)))
})

test_that("prior_type='student' is forwarded correctly to brm", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(2)
  brm_prior_captured <- NULL

  local_mocked_bindings(
    brm = function(prior, ...) { brm_prior_captured <<- prior; make_mock_brmsfit() },
    loo = function(...) make_mock_loo(), posterior_epred = stub_epred,
    posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.001), neff_ratio = function(...) c(a = 0.6),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.1, 0.5, 0.3), `Est.Error` = rep(0.1, 3))
      rownames(m) <- c("Intercept", "bio1", "bio2"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    .package = "CAST"
  )

  bayesianSDM(
    x = pts, predictors = preds, planar_projection = 3857,
    pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
    prior_type = "student", prior_scale = 1.5,
    chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 11
  )

  expect_true(any(grepl("student_t", brm_prior_captured$prior)))
})

test_that("p0 dot-arg is forwarded and doesn't cause an error", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(3)

  local_mocked_bindings(
    brm = function(...) make_mock_brmsfit(), loo = function(...) make_mock_loo(),
    posterior_epred = stub_epred, posterior_summary = function(...) matrix(0, 2, 4),
    rhat = function(...) c(a = 1.001), neff_ratio = function(...) c(a = 0.6),
    as_draws_df = function(...) posterior::as_draws_df(list(x = rnorm(200))),
    nuts_params = function(...) data.frame(Parameter = rep("divergent__", 200),
                                           Value = rep(0, 200)),
    fixef = function(object, summary = TRUE, ...) {
      m <- cbind(Estimate = c(-0.2, 0.8, -0.3, 0.4), `Est.Error` = rep(0.1, 4))
      rownames(m) <- c("Intercept", "bio1", "bio2", "bio3"); m
    },
    .package = "brms"
  )
  local_mocked_bindings(
    aoa   = function(newdata, ...) make_mock_aoa(newdata[[1]]),
    knndm = function(tpoints, ...) {
      n <- nrow(tpoints)
      list(indx_train = list(seq_len(n)), indx_test = list(seq_len(n)))
    },
    .package = "CAST"
  )

  expect_no_error(
    bayesianSDM(
      x = pts, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none",
      prior_type = "horseshoe", vif = FALSE,
      chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 12,
      p0 = 1L
    )
  )
})

test_that("returned list has all required named slots with correct types", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(2)

  with_standard_mocks({
    result <- bayesianSDM(
      x = pts, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
      chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 13
    )

    expected_slots <- c(
      "RasterPredictions", "RasterPredictions_sd", "Predictors", "PCNM",
      "Model", "CVStructure", "LOO", "ConfusionMatrix", "TrainData",
      "TestData", "PredictMatrix", "AOA", "AOA_Diagnostics",
      "PCAModel", "Diagnostics"
    )
    expect_named(result, expected_slots, ignore.order = TRUE)
    expect_s4_class(result$RasterPredictions,    "SpatRaster")
    expect_s4_class(result$RasterPredictions_sd, "SpatRaster")
    expect_s4_class(result$Predictors,           "SpatRaster")
    expect_null(result$PCNM)
    expect_s3_class(result$Model,           "brmsfit")
    expect_type(result$CVStructure,         "list")
    expect_s3_class(result$LOO,             "loo")
    expect_s3_class(result$ConfusionMatrix, "confusionMatrix")
    expect_s3_class(result$TrainData,       "sf")
    expect_s3_class(result$TestData,        "sf")
    expect_s4_class(result$AOA,             "SpatRaster")
    expect_type(result$AOA_Diagnostics,     "list")
    expect_type(result$Diagnostics,         "list")
    expect_named(result$Diagnostics,
                 c("max_Rhat", "min_BulkESS", "min_TailESS", "n_divergent"))
  })
})

test_that("AOA slot is a 3-layer SpatRaster named AOA, DI, LPD", {
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")

  pts   <- make_points()
  preds <- make_predictors(2)

  with_standard_mocks({
    result <- bayesianSDM(
      x = pts, predictors = preds, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
      chains = 1, iter = 200, warmup = 100, cores = 1, k = 3, seed = 14
    )

    expect_equal(terra::nlyr(result$AOA), 3L)
    expect_named(result$AOA, c("AOA", "DI", "LPD"))
    expect_true(result$AOA_Diagnostics$AOA_coverage >= 0 &&
                  result$AOA_Diagnostics$AOA_coverage <= 1)
  })
})


# =============================================================================
# 6. Real integration tests (fit actual Stan model — slow, skip on CRAN)
# =============================================================================

test_that("bayesianSDM runs with minimal synthetic data and produces valid outputs", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  set.seed(42)
  n_pts  <- 100
  coords <- data.frame(lon = runif(n_pts, -120, -110), lat = runif(n_pts, 30, 40))
  x      <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  x$occurrence <- sample(0:1, n_pts, replace = TRUE, prob = c(0.3, 0.7))

  r <- terra::rast(xmin = -120, xmax = -110, ymin = 30, ymax = 40,
                   nrows = 20, ncols = 20, crs = "EPSG:4326")
  pred1 <- r; terra::values(pred1) <- runif(terra::ncell(r), 0, 100)
  pred2 <- r; terra::values(pred2) <- runif(terra::ncell(r), 10, 30)
  predictors <- c(pred1, pred2)
  names(predictors) <- c("bio1", "bio2")

  result <- suppressWarnings(
    bayesianSDM(
      x = x, predictors = predictors, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
      chains = 2, iter = 500, warmup = 250, cores = 1,
      k = 5, seed = 42, fact = 1.5, silent = 2
    )
  )

  expect_type(result, "list")
  expect_s4_class(result$RasterPredictions,    "SpatRaster")
  expect_s4_class(result$RasterPredictions_sd, "SpatRaster")
  expect_s4_class(result$Predictors,           "SpatRaster")
  expect_s3_class(result$Model,           "brmsfit")
  expect_type(result$CVStructure,         "list")
  expect_s3_class(result$LOO,             "loo")
  expect_s3_class(result$ConfusionMatrix, "confusionMatrix")
  expect_s3_class(result$TrainData,       "sf")
  expect_s3_class(result$TestData,        "sf")
  expect_s4_class(result$AOA,             "SpatRaster")
  expect_type(result$AOA_Diagnostics,     "list")
  expect_type(result$Diagnostics,         "list")
  expect_true(result$Diagnostics$max_Rhat    < 1.1)
  expect_true(result$Diagnostics$min_BulkESS > 100)

  # Predictions are valid probabilities
  pred_vals <- terra::values(result$RasterPredictions)
  pred_vals <- pred_vals[!is.na(pred_vals)]
  expect_true(all(pred_vals >= 0 & pred_vals <= 1))

  # SD is non-negative
  sd_vals <- terra::values(result$RasterPredictions_sd)
  sd_vals <- sd_vals[!is.na(sd_vals)]
  expect_true(all(sd_vals >= 0))
})

test_that("bayesianSDM works with PCA predictors", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  set.seed(123)
  n_pts  <- 100
  coords <- data.frame(lon = runif(n_pts, -120, -110), lat = runif(n_pts, 30, 40))
  x      <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

  r    <- terra::rast(xmin = -120, xmax = -110, ymin = 30, ymax = 40,
                      nrows = 15, ncols = 15, crs = "EPSG:4326")
  base <- matrix(runif(terra::ncell(r)), nrow = terra::nrow(r))
  pred_stack <- c(
    terra::rast(base),
    terra::rast(base * 0.8 + matrix(runif(terra::ncell(r), 0, 0.2), nrow = terra::nrow(r))),
    terra::rast(base * 0.6 + matrix(runif(terra::ncell(r), 0, 0.4), nrow = terra::nrow(r))),
    terra::rast(matrix(runif(terra::ncell(r)), nrow = terra::nrow(r)))
  )
  terra::ext(pred_stack) <- terra::ext(r)
  terra::crs(pred_stack) <- terra::crs(r)
  names(pred_stack) <- paste0("bio", 1:4)

  result <- suppressWarnings(
    bayesianSDM(
      x = x, predictors = pred_stack, planar_projection = 3857,
      pca_predictors = TRUE, pca_axes = 3,
      feature_selection = "none", vif = FALSE,
      chains = 2, iter = 400, warmup = 200, cores = 1,
      k = 5, seed = 99, silent = 2
    )
  )

  expect_s3_class(result$PCAModel, "prcomp")
  expect_equal(terra::nlyr(result$Predictors), 3)
  expect_match(names(result$Predictors)[1], "^PC")
  expect_true("occurrence" %in% names(result$TrainData))
})

test_that("bayesianSDM works with feature selection (ffs)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")
  skip_if_not_installed("cmdstanr")

  set.seed(456)
  n_pts  <- 100
  coords <- data.frame(lon = runif(n_pts, -120, -110), lat = runif(n_pts, 30, 40))
  x      <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  x$occurrence <- rbinom(n_pts, 1, 0.6)

  r <- terra::rast(xmin = -120, xmax = -110, ymin = 30, ymax = 40,
                   nrows = 12, ncols = 12, crs = "EPSG:4326")
  pred_list  <- lapply(1:5, function(i) { p <- r; terra::values(p) <- runif(terra::ncell(r)); p })
  predictors <- do.call(c, pred_list)
  names(predictors) <- paste0("var", 1:5)

  result <- suppressWarnings(
    bayesianSDM(
      x = x, predictors = predictors, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "ffs", min_ffs_var = 2,
      vif = FALSE, chains = 2, iter = 400, warmup = 200, cores = 1,
      k = 5, seed = 789, silent = 2
    )
  )

  expect_true(terra::nlyr(result$Predictors) <= 5)
  expect_true(terra::nlyr(result$Predictors) >= 2)
})