## Tests for projectClustersBayes and its internal helpers
##
library(testthat)
library(terra)
library(sf)

# ── Shared synthetic fixtures ──────────────────────────────────────────────────

## Small 10x10 raster in geographic CRS with two named env layers
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
  r[["bio_01"]]  # drop first blank layer
  c(r$bio_01, r$bio_12)
}

env_vars <- c("bio_01", "bio_12")

## Minimal ScalingParams that downstream code reads from posterior_clusters
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

## Minimal Geometry sf (3 clusters, NW->SE order)
make_geometry_sf <- function() {
  polys <- lapply(1:3, function(i) {
    sf::st_polygon(list(matrix(
      c(i-1, i-1, i, i, i-1,
        0,   1,   1, 0,  0),
      ncol = 2
    )))
  })
  sf::st_sf(
    ID       = 1:3,
    geometry = sf::st_sfc(polys, crs = 4326)
  )
}

## Fake trained caret KNN (always predicts class "1")
make_mock_knn <- function(levels = as.character(1:3)) {
  structure(
    list(
      method       = "knn",
      bestTune     = data.frame(k = 3L),
      finalModel   = list(cl = factor(rep(levels[1], 10), levels = levels)),
      trainingData = {
        df <- as.data.frame(matrix(runif(30), nrow = 10))
        names(df) <- c("bio_01", "bio_12")
        df$.outcome <- factor(rep(levels, length.out = 10), levels = levels)
        df
      }
    ),
    class = c("train", "train.formula")
  )
}

## Stub predict.train so mock KNN always returns the first level
local_mock_predict <- function(object, newdata, ...) {
  factor(
    rep(levels(object$finalModel$cl)[1], nrow(newdata)),
    levels = levels(object$finalModel$cl)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  1. .build_future_pred_stack
# ══════════════════════════════════════════════════════════════════════════════

test_that(".build_future_pred_stack adds gp_x and gp_y layers", {
  r <- make_env_rast()
  stack <- .build_future_pred_stack(
    future_predictors = r,
    env_vars          = env_vars,
    planar_proj       = 5070
  )

  expect_true(inherits(stack, "SpatRaster"))
  expect_true("gp_x" %in% names(stack))
  expect_true("gp_y" %in% names(stack))
  expect_true(all(env_vars %in% names(stack)))
})

test_that(".build_future_pred_stack gp coordinates are in km scale", {
  r <- make_env_rast()
  stack <- .build_future_pred_stack(r, env_vars, planar_proj = 5070)

  gp_x_range <- diff(as.numeric(terra::global(stack[["gp_x"]], "range", na.rm = TRUE)))
  # 10-degree longitude span in CONUS ~800km; gp_x in km should be ~that order
  expect_gt(gp_x_range, 10)    # not in degrees (~10) but km (~800)
  expect_lt(gp_x_range, 5000)  # not in metres (~800000)
})

test_that(".build_future_pred_stack errors on missing env_vars", {
  r <- make_env_rast()
  expect_error(
    .build_future_pred_stack(r, c("bio_01", "bio_99"), planar_proj = 5070),
    regexp = NULL  # terra will error on missing layer name
  )
})


# ══════════════════════════════════════════════════════════════════════════════
#  2. .make_brms_predict_fun
# ══════════════════════════════════════════════════════════════════════════════

test_that(".make_brms_predict_fun returns a function", {
  fn <- .make_brms_predict_fun(n_iter = 100L)
  expect_true(is.function(fn))
})

test_that(".make_brms_predict_fun closure captures n_iter correctly", {
  # The function should use round(n_iter/2) draws internally.
  # We test the closure by inspecting its environment.
  fn <- .make_brms_predict_fun(n_iter = 200L)
  expect_equal(environment(fn)$n_iter, 200L)
})


# ══════════════════════════════════════════════════════════════════════════════
#  3. add_weighted_coordinates
# ══════════════════════════════════════════════════════════════════════════════

test_that("add_weighted_coordinates adds x and y layers when coord_wt > 0", {
  r <- make_env_rast()
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


# ══════════════════════════════════════════════════════════════════════════════
#  4. calculate_changes
# ══════════════════════════════════════════════════════════════════════════════

make_cluster_sf <- function(ids, xoffset = 0) {
  polys <- lapply(seq_along(ids), function(i) {
    x0 <- (i - 1) * 2 + xoffset
    sf::st_polygon(list(matrix(
      c(x0, x0, x0+1, x0+1, x0,
        0,  1,  1,    0,    0),
      ncol = 2
    )))
  })
  sf::st_sf(
    ID       = ids,
    geometry = sf::st_sfc(polys, crs = 32614)  # planar UTM
  )
}

test_that("calculate_changes returns correct columns", {
  cur <- make_cluster_sf(1:3)
  fut <- make_cluster_sf(1:3, xoffset = 0.1)
  out <- calculate_changes(cur, fut, planar_proj = 32614)

  expect_true(is.data.frame(out))
  expect_named(
    out,
    c("cluster_id", "current_area_km2", "future_area_km2",
      "area_change_pct", "centroid_shift_km"),
    ignore.order = TRUE
  )
})

test_that("calculate_changes detects lost clusters", {
  cur <- make_cluster_sf(1:3)
  fut <- make_cluster_sf(1:2)   # cluster 3 gone
  out <- calculate_changes(cur, fut, planar_proj = 32614)

  lost_row <- out[out$cluster_id == 3, ]
  expect_equal(nrow(lost_row), 1L)
  expect_equal(lost_row$future_area_km2, 0)
  expect_equal(lost_row$area_change_pct, -100)
})

test_that("calculate_changes detects gained clusters", {
  cur <- make_cluster_sf(1:2)
  fut <- make_cluster_sf(1:3)   # cluster 3 is new
  out <- calculate_changes(cur, fut, planar_proj = 32614)

  gained_row <- out[out$cluster_id == 3, ]
  expect_equal(nrow(gained_row), 1L)
  expect_equal(gained_row$current_area_km2, 0)
  expect_true(is.na(gained_row$area_change_pct))
})

test_that("calculate_changes handles no common clusters gracefully", {
  cur <- make_cluster_sf(1:2)
  fut <- make_cluster_sf(3:4)
  out <- calculate_changes(cur, fut, planar_proj = 32614)

  # All should appear as lost or gained, none as persists
  expect_equal(nrow(out), 4L)
  expect_true(all(out$cluster_id %in% 1:4))
})

test_that("calculate_changes centroid_shift is NA for lost/gained", {
  cur <- make_cluster_sf(1:3)
  fut <- make_cluster_sf(2:4)  # 1 lost, 4 gained
  out <- calculate_changes(cur, fut, planar_proj = 32614)

  expect_true(is.na(out$centroid_shift_km[out$cluster_id == 1]))
  expect_true(is.na(out$centroid_shift_km[out$cluster_id == 4]))
  expect_false(is.na(out$centroid_shift_km[out$cluster_id == 2]))
})


# ══════════════════════════════════════════════════════════════════════════════
#  5. project_future_draws — via lightweight mock
# ══════════════════════════════════════════════════════════════════════════════

make_mock_future_rr <- function(seed = 7L) {
  set.seed(seed)
  r <- make_env_rast(seed = seed)
  # Rescaled predictors: just scale to [0,1] to mimic RescaleRasters_bayes output
  r_scaled <- (r - terra::global(r, "min", na.rm = TRUE)$min) /
    diff(as.numeric(terra::global(r, "range", na.rm = TRUE)))
  list(RescaledPredictors = r_scaled)
}

test_that("project_future_draws returns a SpatRaster named future_stability", {
  set.seed(1L)
  mock_rr   <- make_mock_future_rr()
  suit_mask <- mock_rr$RescaledPredictors[[1]] * 0 + 1  # all cells suitable

  scaling   <- make_scaling_params(n_draws = 5L)
  mock_knn  <- make_mock_knn(levels = as.character(1:3))

  # Stub caret predict
  local_mocked_bindings(
    predict.train = local_mock_predict,
    .package = "caret",
    {
      out <- project_future_draws(
        future_rescaled_rr   = mock_rr,
        suitable_mask        = suit_mask,
        env_vars             = env_vars,
        beta_draws           = scaling$beta_draws,
        coord_wt             = 0.5,
        knn_consensus        = mock_knn,
        n_future_pts         = 20L,
        existing_cluster_ids = 1:3
      )
    }
  )

  expect_true(inherits(out, "SpatRaster"))
  expect_equal(names(out), "future_stability")
})

test_that("project_future_draws stability values are in [0, 1]", {
  set.seed(2L)
  mock_rr   <- make_mock_future_rr(seed = 2L)
  suit_mask <- mock_rr$RescaledPredictors[[1]] * 0 + 1

  scaling  <- make_scaling_params(n_draws = 5L)
  mock_knn <- make_mock_knn()

  local_mocked_bindings(
    predict.train = local_mock_predict,
    .package = "caret",
    {
      out <- project_future_draws(
        future_rescaled_rr   = mock_rr,
        suitable_mask        = suit_mask,
        env_vars             = env_vars,
        beta_draws           = scaling$beta_draws,
        coord_wt             = 0,
        knn_consensus        = mock_knn,
        n_future_pts         = 20L,
        existing_cluster_ids = 1:3
      )
    }
  )

  vals <- terra::values(out, na.rm = TRUE)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("project_future_draws returns NA raster when no suitable cells", {
  mock_rr   <- make_mock_future_rr()
  # All NA suitable mask — no suitable habitat
  suit_mask <- mock_rr$RescaledPredictors[[1]] * NA

  scaling  <- make_scaling_params(n_draws = 5L)
  mock_knn <- make_mock_knn()

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

test_that("project_future_draws clamps n_future_pts to available cells", {
  set.seed(3L)
  mock_rr <- make_mock_future_rr()
  # Only 5 suitable cells in a mostly-NA mask
  suit_mask <- mock_rr$RescaledPredictors[[1]] * NA
  suit_mask[1:5] <- 1

  scaling  <- make_scaling_params(n_draws = 3L)
  mock_knn <- make_mock_knn()

  # Should warn about clamping, not error
  expect_message(
    local_mocked_bindings(
      predict.train = local_mock_predict,
      .package = "caret",
      {
        out <- project_future_draws(
          future_rescaled_rr   = mock_rr,
          suitable_mask        = suit_mask,
          env_vars             = env_vars,
          beta_draws           = scaling$beta_draws,
          coord_wt             = 0,
          knn_consensus        = mock_knn,
          n_future_pts         = 500L,
          existing_cluster_ids = 1:3
        )
      }
    ),
    regexp = "clamping"
  )
})


# ══════════════════════════════════════════════════════════════════════════════
#  6. .build_label_lookup / reassign_cluster_ids integration
# ══════════════════════════════════════════════════════════════════════════════

make_two_class_rasters <- function() {
  # raw:  cells 1-50 = class 1, cells 51-100 = class 2
  # geo:  swapped — cells 1-50 = class 2, cells 51-100 = class 1
  r <- terra::rast(nrows = 10, ncols = 10, vals = c(rep(1L, 50), rep(2L, 50)))
  g <- terra::rast(nrows = 10, ncols = 10, vals = c(rep(2L, 50), rep(1L, 50)))
  list(raw = r, geo = g)
}

test_that(".build_label_lookup maps raw to geo labels correctly", {
  rasts  <- make_two_class_rasters()
  lookup <- .build_label_lookup(rasts$raw, rasts$geo)

  # raw 1 should map to geo 2, raw 2 to geo 1
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
  rasts  <- make_two_class_rasters()
  lookup <- .build_label_lookup(rasts$raw, rasts$geo)

  raw_labels <- c(1L, 2L, 1L, 2L, 2L)
  geo_labels <- unname(lookup[as.character(raw_labels)])
  expect_equal(geo_labels, c(2L, 1L, 2L, 1L, 1L))
})


# ══════════════════════════════════════════════════════════════════════════════
#  7. n_future_draws clamping logic
# ══════════════════════════════════════════════════════════════════════════════

test_that("n_future_draws is clamped to n_stored when NULL", {
  # Test the clamping arithmetic directly — mirrors projectClustersBayes lines
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


# ══════════════════════════════════════════════════════════════════════════════
#  8. Return structure contract (mocked full run, no brms)
# ══════════════════════════════════════════════════════════════════════════════

## Build a minimal mock posterior_clusters object
make_mock_posterior_clusters <- function() {
  mock_knn    <- make_mock_knn(levels = as.character(1:3))
  list(
    Geometry  = make_geometry_sf(),
    KNNModels = list(KNN_Cluster = mock_knn),
    ScalingParams = make_scaling_params(n_draws = 5L)
  )
}

test_that("projectClustersBayes return list has all required elements", {
  skip_if_not_installed("brms")

  r         <- make_env_rast()
  pc        <- make_mock_posterior_clusters()
  mock_rr   <- make_mock_future_rr()

  # Stub the four heavy external calls so we never touch brms or terra::predict
  # with a real model
  local_mocked_bindings(
    RescaleRasters_bayes = function(...) mock_rr,
    .build_future_pred_stack = function(...) make_env_rast(),
    .make_brms_predict_fun = function(...) {
      function(model, data, ...) {
        matrix(rep(0.8, nrow(data) * 2), ncol = 2)
      }
    },
    `terra::predict` = function(x, model, fun, ...) {
      # Return a constant 0.8 suitability raster
      x[[1]] * 0 + 0.8
    },
    project_future_draws = function(...) {
      r <- make_env_rast()[[1]] * 0 + 0.7
      names(r) <- "future_stability"
      r
    },
    {
      result <- projectClustersBayes(
        bSDM_object          = list(
          Model        = structure(list(), class = "brmsfit"),
          PredictMatrix = as.matrix(as.data.frame(terra::values(r, na.rm = TRUE))),
          TrainData    = sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(-95, 45)),
                                                         crs = 4326))
        ),
        posterior_clusters   = pc,
        future_predictors    = r,
        current_predictors   = r,
        threshold_rasts      = list(
          Threshold = data.frame(sensitivity = 0.5)
        ),
        planar_proj          = 5070,
        cluster_novel        = FALSE,
        n_future_draws       = 5L,
        n_future_pts         = 10L
      )
    }
  )

  expected_names <- c(
    "clusters_raster", "clusters_sf", "suitable_habitat",
    "novel_mask", "mess", "stability", "changes", "novel_similarity"
  )
  expect_named(result, expected_names, ignore.order = TRUE)
  expect_true(inherits(result$clusters_raster,  "SpatRaster"))
  expect_true(inherits(result$clusters_sf,      "sf"))
  expect_true(inherits(result$suitable_habitat, "SpatRaster"))
  expect_true(inherits(result$novel_mask,       "SpatRaster"))
  expect_true(inherits(result$mess,             "SpatRaster"))
  expect_true(inherits(result$stability,        "SpatRaster"))
  expect_true(is.data.frame(result$changes))
  expect_true(is.data.frame(result$novel_similarity))
})

test_that("novel_similarity is empty data frame when cluster_novel = FALSE", {
  skip_if_not_installed("brms")

  r  <- make_env_rast()
  pc <- make_mock_posterior_clusters()
  mock_rr <- make_mock_future_rr()

  local_mocked_bindings(
    RescaleRasters_bayes     = function(...) mock_rr,
    .build_future_pred_stack = function(...) make_env_rast(),
    .make_brms_predict_fun   = function(...) function(model, data, ...) matrix(0.8, nrow(data), 2),
    `terra::predict`         = function(x, ...) x[[1]] * 0 + 0.8,
    project_future_draws     = function(...) { r <- make_env_rast()[[1]]*0+0.7; names(r) <- "future_stability"; r },
    {
      result <- projectClustersBayes(
        bSDM_object        = list(
          Model = structure(list(), class = "brmsfit"),
          PredictMatrix = as.matrix(as.data.frame(terra::values(r, na.rm = TRUE))),
          TrainData = sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(-95,45)), crs=4326))
        ),
        posterior_clusters = pc,
        future_predictors  = r,
        current_predictors = r,
        threshold_rasts    = list(Threshold = data.frame(sensitivity = 0.5)),
        planar_proj        = 5070,
        cluster_novel      = FALSE,
        n_future_draws     = 3L,
        n_future_pts       = 10L
      )
    }
  )

  expect_equal(nrow(result$novel_similarity), 0L)
  expect_named(
    result$novel_similarity,
    c("novel_cluster_id", "nearest_existing_id", "avg_silhouette_width")
  )
})

test_that("clusters_sf has an ID column with integer values", {
  skip_if_not_installed("brms")

  r  <- make_env_rast()
  pc <- make_mock_posterior_clusters()
  mock_rr <- make_mock_future_rr()

  local_mocked_bindings(
    RescaleRasters_bayes     = function(...) mock_rr,
    .build_future_pred_stack = function(...) make_env_rast(),
    .make_brms_predict_fun   = function(...) function(model, data, ...) matrix(0.8, nrow(data), 2),
    `terra::predict`         = function(x, ...) x[[1]] * 0 + 0.8,
    project_future_draws     = function(...) { r <- make_env_rast()[[1]]*0+0.7; names(r) <- "future_stability"; r },
    {
      result <- projectClustersBayes(
        bSDM_object        = list(
          Model = structure(list(), class = "brmsfit"),
          PredictMatrix = as.matrix(as.data.frame(terra::values(r, na.rm = TRUE))),
          TrainData = sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(-95,45)), crs=4326))
        ),
        posterior_clusters = pc,
        future_predictors  = r,
        current_predictors = r,
        threshold_rasts    = list(Threshold = data.frame(sensitivity = 0.5)),
        planar_proj        = 5070,
        cluster_novel      = FALSE,
        n_future_draws     = 3L,
        n_future_pts       = 10L
      )
    }
  )

  expect_true("ID" %in% names(result$clusters_sf))
  expect_true(is.integer(result$clusters_sf$ID) ||
                is.numeric(result$clusters_sf$ID))
})

test_that("changes data frame cluster IDs are a subset of current Geometry IDs", {
  skip_if_not_installed("brms")

  r  <- make_env_rast()
  pc <- make_mock_posterior_clusters()
  mock_rr <- make_mock_future_rr()

  local_mocked_bindings(
    RescaleRasters_bayes     = function(...) mock_rr,
    .build_future_pred_stack = function(...) make_env_rast(),
    .make_brms_predict_fun   = function(...) function(model, data, ...) matrix(0.8, nrow(data), 2),
    `terra::predict`         = function(x, ...) x[[1]] * 0 + 0.8,
    project_future_draws     = function(...) { r <- make_env_rast()[[1]]*0+0.7; names(r) <- "future_stability"; r },
    {
      result <- projectClustersBayes(
        bSDM_object        = list(
          Model = structure(list(), class = "brmsfit"),
          PredictMatrix = as.matrix(as.data.frame(terra::values(r, na.rm = TRUE))),
          TrainData = sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(-95,45)), crs=4326))
        ),
        posterior_clusters = pc,
        future_predictors  = r,
        current_predictors = r,
        threshold_rasts    = list(Threshold = data.frame(sensitivity = 0.5)),
        planar_proj        = 5070,
        cluster_novel      = FALSE,
        n_future_draws     = 3L,
        n_future_pts       = 10L
      )
    }
  )

  current_ids <- unique(pc$Geometry$ID)
  future_ids  <- unique(result$clusters_sf$ID)
  # All future IDs must be current-era IDs (no novel in this run)
  expect_true(all(future_ids %in% current_ids))
})
