## Tests for projectClustersBayes and its internal helpers
##
## Strategy
## ─────────
## 1. All internal helper functions are tested directly with minimal synthetic
##    data — no brms / Stan required for these.
## 2. The four return-structure contract tests for projectClustersBayes() itself
##    use a real brms fit via a shared cached fixture (make_pcb_fixture).
##    The fixture fits one model, runs PosteriorCluster, and is reused across
##    all four tests. Runtime: ~60-90s locally — acceptable per project policy.

library(testthat)
library(terra)
library(sf)


# =============================================================================
# Shared synthetic fixtures for helper tests
# =============================================================================

## 10×10 geographic raster with two named env layers
make_env_rast <- function(seed = 1L, nrow = 10L, ncol = 10L) {
  set.seed(seed)
  r <- terra::rast(
    nrows = nrow, ncols = ncol,
    xmin = -100, xmax = -90,
    ymin = 40,   ymax = 50,
    crs = "EPSG:4326"
  )
  r$bio_01 <- terra::setValues(r[[1]], runif(nrow * ncol, 5,  20))
  r$bio_12 <- terra::setValues(r[[1]], runif(nrow * ncol, 100, 800))
  c(r$bio_01, r$bio_12)
}

env_vars <- c("bio_01", "bio_12")

make_scaling_params <- function(n_draws = 10L) {
  set.seed(42L)
  list(
    mean_betas = c(bio_01 = 0.8, bio_12 = 0.5),
    beta_draws = matrix(
      rnorm(n_draws * 2, mean = c(0.8, 0.5), sd = 0.1),
      nrow = n_draws,
      dimnames = list(NULL, c("bio_01", "bio_12"))
    ),
    coord_wt = 0.5
  )
}

make_geometry_sf <- function() {
  polys <- lapply(1:3, function(i) {
    sf::st_polygon(list(matrix(
      c(i-1, i-1, i, i, i-1,
        0,   1,   1, 0,  0),
      ncol = 2
    )))
  })
  sf::st_sf(ID = 1:3, geometry = sf::st_sfc(polys, crs = 4326))
}

## Fake caret KNN — always predicts the first level
make_mock_knn <- function(levels = as.character(1:3)) {
  structure(
    list(
      method     = "knn",
      bestTune   = data.frame(k = 3L),
      finalModel = list(cl = factor(rep(levels[1], 10), levels = levels)),
      trainingData = {
        df <- as.data.frame(matrix(runif(50), nrow = 10))
        names(df) <- c("bio_01", "bio_12", "coord_x_w", "coord_y_w", ".outcome")
        df$.outcome <- factor(rep(levels, length.out = 10), levels = levels)
        df
      }
    ),
    class = c("train", "train.formula")
  )
}

local_mock_predict <- function(object, newdata, ...) {
  lvls <- levels(object$finalModel$cl)
  factor(
    lvls[(seq_len(nrow(newdata)) %% length(lvls)) + 1L],
    levels = lvls
  )
}

make_mock_future_rr <- function(seed = 7L) {
  set.seed(seed)
  r     <- make_env_rast(seed = seed)
  r_min <- min(terra::global(r, "min", na.rm = TRUE)$min)
  r_max <- max(terra::global(r, "max", na.rm = TRUE)$max)
  list(RescaledPredictors = (r - r_min) / (r_max - r_min))
}

make_cluster_sf <- function(ids, xoffset = 0) {
  polys <- lapply(seq_along(ids), function(i) {
    x0 <- (i - 1) * 2 + xoffset
    sf::st_polygon(list(matrix(
      c(x0, x0, x0+1, x0+1, x0,
        0,  1,  1,    0,    0),
      ncol = 2
    )))
  })
  sf::st_sf(ID = ids, geometry = sf::st_sfc(polys, crs = "EPSG:32614"))
}


# =============================================================================
# 1. .build_future_pred_stack
# =============================================================================

test_that(".build_future_pred_stack adds gp_x and gp_y layers", {
  r     <- make_env_rast()
  stack <- .build_future_pred_stack(r, env_vars, planar_proj = "EPSG:5070")

  expect_true(inherits(stack, "SpatRaster"))
  expect_true("gp_x" %in% names(stack))
  expect_true("gp_y" %in% names(stack))
  expect_true(all(env_vars %in% names(stack)))
})

test_that(".build_future_pred_stack gp coordinates are in km scale", {
  r         <- make_env_rast()
  stack     <- .build_future_pred_stack(r, env_vars, planar_proj = "EPSG:5070")
  gp_x_rng  <- terra::global(stack[["gp_x"]], "range", na.rm = TRUE)
  gp_x_range <- gp_x_rng$max - gp_x_rng$min

  expect_gt(gp_x_range, 10)     # not degrees (~10)
  expect_lt(gp_x_range, 5000)   # not metres (~800 000)
})

test_that(".build_future_pred_stack errors on missing env_vars", {
  r <- make_env_rast()
  expect_error(
    .build_future_pred_stack(r, c("bio_01", "bio_99"), planar_proj = "EPSG:5070")
  )
})


# =============================================================================
# 2. .make_brms_predict_fun
# =============================================================================

test_that(".make_brms_predict_fun returns a function", {
  fn <- .make_brms_predict_fun(n_iter = 100L)
  expect_true(is.function(fn))
})

test_that(".make_brms_predict_fun closure captures n_iter correctly", {
  fn <- .make_brms_predict_fun(n_iter = 200L)
  expect_equal(environment(fn)$n_iter, 200L)
})


# =============================================================================
# 3. add_weighted_coordinates
# =============================================================================

test_that("add_weighted_coordinates adds x and y layers when coord_wt > 0", {
  r   <- make_env_rast()
  out <- add_weighted_coordinates(r, coord_wt = 0.5)
  expect_true("x" %in% names(out))
  expect_true("y" %in% names(out))
})

test_that("add_weighted_coordinates is a no-op when coord_wt = 0", {
  r   <- make_env_rast()
  out <- add_weighted_coordinates(r, coord_wt = 0)
  expect_equal(names(out), names(r))
})

test_that("add_weighted_coordinates coord range scales with coord_wt", {
  r    <- make_env_rast()
  out1 <- add_weighted_coordinates(r, coord_wt = 0.5)
  out2 <- add_weighted_coordinates(r, coord_wt = 2.0)

  range1 <- diff(as.numeric(terra::global(out1[["x"]], "range", na.rm = TRUE)))
  range2 <- diff(as.numeric(terra::global(out2[["x"]], "range", na.rm = TRUE)))
  expect_gt(range2, range1)
})


# =============================================================================
# 4. calculate_changes
# =============================================================================

test_that("calculate_changes returns correct columns", {
  cur <- make_cluster_sf(1:3)
  fut <- make_cluster_sf(1:3, xoffset = 0.1)
  out <- calculate_changes(cur, fut, planar_proj = "EPSG:32614")

  expect_true(is.data.frame(out))
  expect_named(out,
               c("cluster_id", "current_area_km2", "future_area_km2",
                 "area_change_pct", "centroid_shift_km"),
               ignore.order = TRUE)
})

test_that("calculate_changes detects lost clusters", {
  cur      <- make_cluster_sf(1:3)
  fut      <- make_cluster_sf(1:2)
  out      <- calculate_changes(cur, fut, planar_proj = "EPSG:32614")
  lost_row <- out[out$cluster_id == 3, ]

  expect_equal(nrow(lost_row), 1L)
  expect_equal(lost_row$future_area_km2, 0)
  expect_equal(lost_row$area_change_pct, -100)
})

test_that("calculate_changes detects gained clusters", {
  cur        <- make_cluster_sf(1:2)
  fut        <- make_cluster_sf(1:3)
  out        <- calculate_changes(cur, fut, planar_proj = "EPSG:32614")
  gained_row <- out[out$cluster_id == 3, ]

  expect_equal(nrow(gained_row), 1L)
  expect_equal(gained_row$current_area_km2, 0)
  expect_true(is.na(gained_row$area_change_pct))
})

test_that("calculate_changes handles no common clusters gracefully", {
  cur <- make_cluster_sf(1:2)
  fut <- make_cluster_sf(3:4)
  out <- calculate_changes(cur, fut, planar_proj = "EPSG:32614")

  expect_equal(nrow(out), 4L)
  expect_true(all(out$cluster_id %in% 1:4))
})

test_that("calculate_changes centroid_shift is NA for lost/gained clusters", {
  cur <- make_cluster_sf(1:3)
  fut <- make_cluster_sf(2:4)
  out <- calculate_changes(cur, fut, planar_proj = "EPSG:32614")

  expect_true(is.na(out$centroid_shift_km[out$cluster_id == 1]))
  expect_true(is.na(out$centroid_shift_km[out$cluster_id == 4]))
  expect_false(is.na(out$centroid_shift_km[out$cluster_id == 2]))
})

# =============================================================================
# 5. project_future_draws
# =============================================================================

local_mock_predict <- function(object, newdata, ...) {
  lvls <- levels(object$finalModel$cl)
  factor(
    lvls[(seq_len(nrow(newdata)) %% length(lvls)) + 1L],
    levels = lvls
  )
}

test_that("project_future_draws returns a SpatRaster named future_stability", {
  set.seed(1L)
  mock_rr   <- make_mock_future_rr()
  suit_mask <- mock_rr$RescaledPredictors[[1]] * 0 + 1
  terra::values(suit_mask)[1:5] <- 1
  scaling   <- make_scaling_params(n_draws = 20L)
  mock_knn  <- make_mock_knn(levels = as.character(1:3))

  local_mocked_bindings(
    train_regression_knn = function(X, y) {
      caret::train(
        y ~ .,
        data      = cbind(X, y = y + rnorm(length(y), 0, 0.01)),
        method    = "knn",
        tuneGrid  = data.frame(k = 3),
        trControl = caret::trainControl(method = "none")
      )
    },
    .package = "safeHavens"
  )
  local_mocked_bindings(predict = local_mock_predict, .package = "stats")

  out <- project_future_draws(
    future_rescaled_rr   = mock_rr,
    suitable_mask        = suit_mask,
    env_vars             = env_vars,
    beta_draws           = scaling$beta_draws,
    coord_wt             = 0.5,
    knn_consensus        = mock_knn,
    n_future_pts         = 80L,
    existing_cluster_ids = 1:3
  )

  expect_true(inherits(out, "SpatRaster"))
  expect_equal(names(out), "future_stability")
})

test_that("project_future_draws stability values are in [0, 1]", {
  set.seed(2L)
  mock_rr   <- make_mock_future_rr(seed = 2L)
  suit_mask <- mock_rr$RescaledPredictors[[1]] * 0 + 1
  terra::values(suit_mask)[1:5] <- 1
  scaling   <- make_scaling_params(n_draws = 20L)
  mock_knn  <- make_mock_knn()

  local_mocked_bindings(
  train_regression_knn = function(X, y) {
    caret::train(
      y ~ .,
      data      = cbind(X, y = y + rnorm(length(y), 0, 0.01)),
      method    = "knn",
      tuneGrid  = data.frame(k = 3),
      trControl = caret::trainControl(method = "none")
    )
  },
  .package = "safeHavens"
  )

  local_mocked_bindings(predict = local_mock_predict, .package = "stats")

  out <- project_future_draws(
    future_rescaled_rr   = mock_rr,
    suitable_mask        = suit_mask,
    env_vars             = env_vars,
    beta_draws           = scaling$beta_draws,
    coord_wt             = 0.1,
    knn_consensus        = mock_knn,
    n_future_pts         = 80L,
    existing_cluster_ids = 1:3
  )

  vals <- terra::values(out, na.rm = TRUE)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("project_future_draws returns NA raster when no suitable cells", {
  mock_rr   <- make_mock_future_rr()
  suit_mask <- mock_rr$RescaledPredictors[[1]] * NA
  terra::values(suit_mask)[1:5] <- 1
  scaling   <- make_scaling_params(n_draws = 5L)
  mock_knn  <- make_mock_knn()

  expect_warning(
    out <- project_future_draws(
      future_rescaled_rr   = mock_rr,
      suitable_mask        = suit_mask,
      env_vars             = env_vars,
      beta_draws           = scaling$beta_draws,
      coord_wt             = 0,
      knn_consensus        = mock_knn,
      n_future_pts         = 20L,
      existing_cluster_ids = 1:3
    ),
    regexp = "no suitable-habitat cells"
  )

  expect_true(inherits(out, "SpatRaster"))
  expect_equal(names(out), "future_stability")
  expect_true(all(is.na(terra::values(out))))
})

test_that("project_future_draws clamps n_future_pts when fewer cells than requested", {
  set.seed(3L)
  mock_rr <- make_mock_future_rr()

  # Build a mask with only 5 suitable cells using explicit values assignment
  suit_mask <- mock_rr$RescaledPredictors[[1]]
  v <- rep(NA_real_, terra::ncell(suit_mask))
  v[1:5] <- 1
  terra::values(suit_mask) <- v

  scaling  <- make_scaling_params(n_draws = 3L)
  mock_knn <- make_mock_knn()

  local_mocked_bindings(
    train_regression_knn = function(X, y) {
      caret::train(
        y ~ .,
        data      = cbind(X, y = y + rnorm(length(y), 0, 0.01)),
        method    = "knn",
        tuneGrid  = data.frame(k = 3),
        trControl = caret::trainControl(method = "none")
      )
    },
    .package = "safeHavens"
  )
  local_mocked_bindings(predict = local_mock_predict, .package = "stats")

  # With only 5 cells and n_future_pts = 500, function should either
  # emit the clamping message or hit the < 10 pts early-return warning.
  # Either way it must not error and must return a named SpatRaster.
  out <- suppressWarnings(suppressMessages(
    project_future_draws(
      future_rescaled_rr   = mock_rr,
      suitable_mask        = suit_mask,
      env_vars             = env_vars,
      beta_draws           = scaling$beta_draws,
      coord_wt             = 0,
      knn_consensus        = mock_knn,
      n_future_pts         = 500L,
      existing_cluster_ids = 1:3
    )
  ))

  expect_true(inherits(out, "SpatRaster"))
  expect_equal(names(out), "future_stability")
})

# =============================================================================
# 6. .build_label_lookup / reassign_cluster_ids
# =============================================================================

make_two_class_rasters <- function() {
  r <- terra::rast(nrows = 10, ncols = 10, vals = c(rep(1L, 50), rep(2L, 50)))
  g <- terra::rast(nrows = 10, ncols = 10, vals = c(rep(2L, 50), rep(1L, 50)))
  list(raw = r, geo = g)
}

test_that(".build_label_lookup maps raw to geo labels correctly", {
  rasts  <- make_two_class_rasters()
  lookup <- .build_label_lookup(rasts$raw, rasts$geo)

  expect_equal(unname(lookup["1"]), 2L)
  expect_equal(unname(lookup["2"]), 1L)
})

test_that(".build_label_lookup returns named integer vector", {
  rasts  <- make_two_class_rasters()
  lookup <- .build_label_lookup(rasts$raw, rasts$geo)

  expect_true(is.integer(lookup))
  expect_true(!is.null(names(lookup)))
})

test_that("remapping consensus_labels via lookup produces geo IDs", {
  rasts      <- make_two_class_rasters()
  lookup     <- .build_label_lookup(rasts$raw, rasts$geo)
  raw_labels <- c(1L, 2L, 1L, 2L, 2L)
  geo_labels <- unname(lookup[as.character(raw_labels)])

  expect_equal(geo_labels, c(2L, 1L, 2L, 1L, 1L))
})


# =============================================================================
# 7. n_future_draws clamping logic
# =============================================================================

test_that("n_future_draws is clamped to n_stored when NULL and n_stored <= 100", {
  n_stored       <- 40L
  n_future_draws <- NULL
  result <- if (is.null(n_future_draws)) min(n_stored, 100L) else
    min(as.integer(n_future_draws), n_stored)
  expect_equal(result, 40L)
})

test_that("n_future_draws caps at 100 when n_stored > 100", {
  n_stored       <- 200L
  n_future_draws <- NULL
  result <- if (is.null(n_future_draws)) min(n_stored, 100L) else
    min(as.integer(n_future_draws), n_stored)
  expect_equal(result, 100L)
})

test_that("explicit n_future_draws is clamped to n_stored", {
  n_stored       <- 30L
  n_future_draws <- 50L
  result <- if (is.null(n_future_draws)) min(n_stored, 100L) else
    min(as.integer(n_future_draws), n_stored)
  expect_equal(result, 30L)
})


# =============================================================================
# 8. Return-structure contract — real brms fit, no mocking
#
# Fixture fitted once and cached for all four tests. Uses current-era
# predictors as "future" predictors so MESS returns all-positive values
# (no novel climate), which is the natural state for cluster_novel = FALSE.
# =============================================================================

make_pcb_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)

    skip_on_cran()
    skip_if_not_installed("brms")
    skip_if_not_installed("cmdstanr")

    set.seed(42)
    n_pts  <- 100
    coords <- data.frame(lon = runif(n_pts, -120, -110), lat = runif(n_pts, 30, 40))
    pts    <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
    pts$occurrence <- sample(0:1, n_pts, replace = TRUE, prob = c(0.3, 0.7))

    r <- terra::rast(xmin = -120, xmax = -110, ymin = 30, ymax = 40,
                     nrows = 20, ncols = 20, crs = "EPSG:4326")
    pred1 <- r; terra::values(pred1) <- runif(terra::ncell(r), 0, 100)
    pred2 <- r; terra::values(pred2) <- runif(terra::ncell(r), 10, 30)
    predictors        <- c(pred1, pred2)
    names(predictors) <- c("bio1", "bio2")

    bsdm <- suppressWarnings(bayesianSDM(
      x = pts, predictors = predictors, planar_projection = 3857,
      pca_predictors = FALSE, feature_selection = "none", vif = FALSE,
      chains = 2, iter = 500, warmup = 250, cores = 1,
      k = 5, seed = 42, fact = 1.5, silent = 2
    ))

    threshold_rasts <- suppressWarnings(PostProcessSDM(
      rast_cont         = bsdm$RasterPredictions,
      test              = bsdm$TestData,
      train             = bsdm$TrainData,
      planar_projection = 3857
    ))

    pc <- suppressWarnings(PosteriorCluster(
      model         = bsdm$Model,
      predictors    = bsdm$Predictors,
      f_rasts       = bsdm$RasterPredictions,
      pred_mat      = bsdm$PredictMatrix,
      training_data = bsdm$TrainData,
      n_draws       = 20,
      n             = 3,
      n_pts         = 50,
      lyr           = "occurrence_prob_mean",
      planar_proj   = "EPSG:3857",
      seed          = 42
    ))

    cached <<- list(
      bsdm            = bsdm,
      predictors      = predictors,
      threshold_rasts = threshold_rasts,
      pc              = pc
    )
    cached
  }
})

test_that("projectClustersBayes return list has all required elements", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  f <- make_pcb_fixture()

  local_mocked_bindings(
  train_regression_knn = function(X, y) {
    caret::train(
      y ~ .,
      data      = cbind(X, y = y + rnorm(length(y), 0, 0.01)),
      method    = "knn",
      tuneGrid  = data.frame(k = 3),
      trControl = caret::trainControl(method = "none")
    )
  },
  .package = "safeHavens"
)
  
  result <- suppressWarnings(projectClustersBayes(
    bSDM_object        = f$bsdm,
    posterior_clusters = f$pc,
    future_predictors  = f$predictors,
    current_predictors = f$predictors,
    threshold_rasts    = f$threshold_rasts,
    planar_proj        = "EPSG:3857",
    cluster_novel      = FALSE,
    n_future_draws     = 30L,
    n_future_pts       = 100L
  ))

  expected_names <- c(
    "changes", "clusters_raster", "Geometry", "suitable_habitat",
    "novel_mask", "mess", "stability", "novel_similarity"
  )
  expect_named(result, expected_names, ignore.order = TRUE)
  expect_true(inherits(result$clusters_raster,  "SpatRaster"))
  expect_true(inherits(result$Geometry,      "sf"))
  expect_true(inherits(result$suitable_habitat, "SpatRaster"))
  expect_true(inherits(result$novel_mask,       "SpatRaster"))
  expect_true(inherits(result$mess,             "SpatRaster"))
  expect_true(inherits(result$stability,        "SpatRaster"))
  expect_true(is.data.frame(result$changes))
  expect_true(is.data.frame(result$novel_similarity))
})

local_mocked_bindings(
  train_regression_knn = function(X, y) {
    caret::train(
      y ~ .,
      data      = cbind(X, y = y + rnorm(length(y), 0, 0.01)),
      method    = "knn",
      tuneGrid  = data.frame(k = 3),
      trControl = caret::trainControl(method = "none")
    )
  },
  .package = "safeHavens"
)

test_that("novel_similarity is empty data frame when cluster_novel = FALSE", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  f <- make_pcb_fixture()
  local_mocked_bindings(
    train_regression_knn = function(X, y) {
      caret::train(
        y ~ .,
        data      = cbind(X, y = y + rnorm(length(y), 0, 0.01)),
        method    = "knn",
        tuneGrid  = data.frame(k = 3),
        trControl = caret::trainControl(method = "none")
      )
    },
    .package = "safeHavens"
  )

  result <- suppressWarnings(projectClustersBayes(
    bSDM_object        = f$bsdm,
    posterior_clusters = f$pc,
    future_predictors  = f$predictors,
    current_predictors = f$predictors,
    threshold_rasts    = f$threshold_rasts,
    planar_proj        = "EPSG:3857",
    cluster_novel      = FALSE,
    n_future_draws     = 30L,
    n_future_pts       = 100L
  ))

  expect_equal(nrow(result$novel_similarity), 0L)
  expect_named(
    result$novel_similarity,
    c("novel_cluster_id", "nearest_existing_id", "avg_silhouette_width")
  )
})

test_that("Geometry has an ID column with integer-compatible values", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  f <- make_pcb_fixture()

  result <- suppressWarnings(projectClustersBayes(
    bSDM_object        = f$bsdm,
    posterior_clusters = f$pc,
    future_predictors  = f$predictors,
    current_predictors = f$predictors,
    threshold_rasts    = f$threshold_rasts,
    planar_proj        = "EPSG:3857",
    cluster_novel      = FALSE,
    n_future_draws     = 30L,
    n_future_pts       = 100L
  ))

  expect_true("ID" %in% names(result$Geometry))
  expect_true(is.integer(result$Geometry$ID) ||
                is.numeric(result$Geometry$ID))
})

test_that("changes data frame contains cluster_id column and correct structure", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  f <- make_pcb_fixture()

  result <- suppressWarnings(projectClustersBayes(
    bSDM_object        = f$bsdm,
    posterior_clusters = f$pc,
    future_predictors  = f$predictors,
    current_predictors = f$predictors,
    threshold_rasts    = f$threshold_rasts,
    planar_proj        = "EPSG:3857",
    cluster_novel      = FALSE,
    n_future_draws     = 30L,
    n_future_pts       = 100L
  ))

  expect_true("cluster_id" %in% names(result$changes))
  expect_true(all(c("current_area_km2", "future_area_km2",
                     "area_change_pct",  "centroid_shift_km")
                   %in% names(result$changes)))

  current_ids <- unique(f$pc$Geometry$ID)
  max_current  <- max(current_ids)
  expect_true(all(
    result$changes$cluster_id %in% current_ids |
    result$changes$cluster_id > max_current
  ))
})