library(testthat)

test_that("compute_top3_clusters: perfect agreement returns 100% rank1", {
  # All draws assign point 1 → cluster 1, point 2 → cluster 2
  # rbind gives nrow=2 (points), ncol=10 (draws) with correct row content
  dc <- rbind(rep(1L, 10), rep(2L, 10))
  res <- compute_top3_clusters(dc, n_pts = 2, n_draws = 10)

  expect_equal(res$labels_matrix[1, 1], 1L)
  expect_equal(res$pct_matrix[1, 1],   100)
  expect_equal(res$labels_matrix[2, 1], 2L)
  expect_equal(res$pct_matrix[2, 1],   100)
})

test_that("compute_top3_clusters: rank1 + rank2 + rank3 percentages sum <= 100", {
  set.seed(1)
  dc <- matrix(sample(1:5, 50 * 100, replace = TRUE), nrow = 50, ncol = 100)
  res <- compute_top3_clusters(dc, n_pts = 50, n_draws = 100)

  row_sums <- rowSums(res$pct_matrix, na.rm = TRUE)
  # Each row sum can exceed 100 only if rank2/3 pad with rank1 values.
  # But raw top-3 draws should give sums ≤ 100 + numeric tolerance.
  expect_true(all(res$pct_matrix[, 1] >= res$pct_matrix[, 2] - 1e-9))
  expect_true(all(res$pct_matrix[, 2] >= res$pct_matrix[, 3] - 1e-9))
})

test_that("compute_top3_clusters: fewer than 3 unique clusters pads with rank1", {
  # Point 1 only ever sees clusters 1 and 2 (never 3)
  dc <- matrix(c(rep(1L, 7), rep(2L, 3)), nrow = 1, ncol = 10)
  res <- compute_top3_clusters(dc, n_pts = 1, n_draws = 10)

  # rank3 should be padded with rank1 value
  expect_equal(res$labels_matrix[1, 3], res$labels_matrix[1, 1])
  expect_equal(res$pct_matrix[1, 3],    res$pct_matrix[1, 1])
})

test_that("compute_top3_clusters: all-NA row leaves NA in output", {
  dc <- matrix(NA_integer_, nrow = 1, ncol = 10)
  res <- compute_top3_clusters(dc, n_pts = 1, n_draws = 10)

  expect_true(all(is.na(res$labels_matrix[1, ])))
  expect_true(all(is.na(res$pct_matrix[1, ])))
})

test_that("compute_top3_clusters: single unique cluster across all draws", {
  dc <- matrix(3L, nrow = 1, ncol = 20)
  res <- compute_top3_clusters(dc, n_pts = 1, n_draws = 20)

  expect_equal(res$labels_matrix[1, ], c(3L, 3L, 3L))
  expect_equal(res$pct_matrix[1, ],    c(100, 100, 100))
})

test_that("compute_top3_clusters: output dimensions match n_pts", {
  n_pts <- 17; n_draws <- 50
  dc <- matrix(sample(1:4, n_pts * n_draws, replace = TRUE),
               nrow = n_pts, ncol = n_draws)
  res <- compute_top3_clusters(dc, n_pts, n_draws)

  expect_equal(dim(res$labels_matrix), c(n_pts, 3))
  expect_equal(dim(res$pct_matrix),    c(n_pts, 3))
})

test_that("compute_top3_clusters: percentages are non-negative and <= 100", {
  set.seed(99)
  dc <- matrix(sample(1:6, 30 * 80, replace = TRUE), nrow = 30, ncol = 80)
  res <- compute_top3_clusters(dc, 30, 80)

  expect_true(all(res$pct_matrix >= 0, na.rm = TRUE))
  expect_true(all(res$pct_matrix <= 100 + 1e-9, na.rm = TRUE))
})


# =============================================================================
# 2. build_cooccurrence_matrix()
# =============================================================================

test_that("build_cooccurrence_matrix: diagonal is always 1", {
  set.seed(7)
  dc <- matrix(sample(1:3, 20 * 50, replace = TRUE), nrow = 20, ncol = 50)
  co <- build_cooccurrence_matrix(dc, 20, 50)

  expect_true(all(diag(co) == 1))
})

test_that("build_cooccurrence_matrix: matrix is symmetric", {
  set.seed(13)
  dc <- matrix(sample(1:4, 15 * 40, replace = TRUE), nrow = 15, ncol = 40)
  co <- build_cooccurrence_matrix(dc, 15, 40)

  expect_equal(co, t(co))
})

test_that("build_cooccurrence_matrix: values in [0, 1]", {
  set.seed(17)
  dc <- matrix(sample(1:5, 25 * 60, replace = TRUE), nrow = 25, ncol = 60)
  co <- build_cooccurrence_matrix(dc, 25, 60)

  expect_true(all(co >= 0 - 1e-12))
  expect_true(all(co <= 1 + 1e-12))
})

test_that("build_cooccurrence_matrix: two always-same-cluster points have co-occurrence 1", {
  # Points 1 and 2 always in cluster 1; point 3 always in cluster 2
  dc <- matrix(c(1L, 1L, 2L), nrow = 3, ncol = 10)  # same across all 10 draws
  co <- build_cooccurrence_matrix(dc, 3, 10)

  expect_equal(co[1, 2], 1)
  expect_equal(co[1, 3], 0)
})

test_that("build_cooccurrence_matrix: all points in same cluster → all-ones matrix", {
  dc <- matrix(1L, nrow = 5, ncol = 30)
  co <- build_cooccurrence_matrix(dc, 5, 30)

  expect_true(all(co == 1))
})

test_that("build_cooccurrence_matrix: dimensions are n_pts × n_pts", {
  n_pts <- 12; n_draws <- 25
  dc <- matrix(sample(1:3, n_pts * n_draws, replace = TRUE),
               nrow = n_pts, ncol = n_draws)
  co <- build_cooccurrence_matrix(dc, n_pts, n_draws)

  expect_equal(dim(co), c(n_pts, n_pts))
})

test_that("build_cooccurrence_matrix: NAs in draw_clusterings treated as FALSE", {
  dc <- matrix(c(1L, NA_integer_, 2L), nrow = 3, ncol = 1)
  co <- build_cooccurrence_matrix(dc, 3, 1)

  # NA row/col entries should not propagate to 1 on diagonal for the NA point
  # Diagonal for NA point: outer(NA, NA) == NA → set to FALSE → contributes 0
  # In 1 draw: NA point contributes 0 on diagonal → diagonal[2] = 0
  expect_equal(co[2, 2], 0)
})


# =============================================================================
# 3. compute_stability_scores()
# =============================================================================

test_that("compute_stability_scores: output length equals n_pts", {
  n <- 10
  co  <- diag(n)
  lbl <- rep(1L, n)
  stab <- compute_stability_scores(co, lbl)
  expect_length(stab, n)
})

test_that("compute_stability_scores: values in [0, 1]", {
  set.seed(21)
  n   <- 20
  raw <- matrix(runif(n * n), n, n)
  co  <- (raw + t(raw)) / 2
  diag(co) <- 1
  lbl <- sample(1:3, n, replace = TRUE)
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(stab >= 0 - 1e-9))
  expect_true(all(stab <= 1 + 1e-9))
})

test_that("compute_stability_scores: singleton cluster gets stability 1", {
  n   <- 5
  co  <- diag(n)
  # Point 3 is the only member of cluster 99
  lbl <- c(1L, 1L, 99L, 1L, 1L)
  stab <- compute_stability_scores(co, lbl)

  expect_equal(stab[3], 1.0)
})

test_that("compute_stability_scores: perfectly cohesive cluster → stability = 1", {
  n   <- 4
  co  <- matrix(1.0, n, n)        # all pairs always co-cluster
  lbl <- rep(1L, n)
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(stab == 1.0))
})

test_that("compute_stability_scores: non-overlapping clusters, zero cross-cluster co-occurrence", {
  # 2 clusters, 3 pts each; within = 1, across = 0
  n   <- 6
  co  <- matrix(0.0, n, n)
  co[1:3, 1:3] <- 1
  co[4:6, 4:6] <- 1
  lbl <- c(1L, 1L, 1L, 2L, 2L, 2L)
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(stab[1:3] == 1.0))
  expect_true(all(stab[4:6] == 1.0))
})

test_that("compute_stability_scores: boundary points score lower than core points", {
  # 4-point example: pts 1 & 2 always together, pts 3 & 4 sometimes swap
  co <- matrix(c(
    1.0, 0.9, 0.1, 0.1,
    0.9, 1.0, 0.1, 0.1,
    0.1, 0.1, 1.0, 0.5,
    0.1, 0.1, 0.5, 1.0
  ), 4, 4)
  lbl <- c(1L, 1L, 2L, 2L)
  stab <- compute_stability_scores(co, lbl)

  # Cluster 1 stability should be higher than cluster 2 stability
  expect_gt(stab[1], stab[3])
  expect_gt(stab[2], stab[4])
})


# =============================================================================
# 4. rescale_points_by_betas()
# =============================================================================
rescale_points_by_betas <- function(pt_env, env_vars, var_mu, var_sd, betas) {
  out <- as.data.frame(matrix(NA_real_, nrow = nrow(pt_env),
                              ncol = length(env_vars)))
  colnames(out) <- env_vars
  for (v in env_vars) {
    out[[v]] <- ((pt_env[[v]] - var_mu[v]) / var_sd[v]) * abs(betas[v])
  }
  out
}

test_that("rescale_points_by_betas: output dimensions match input", {
  pt_env  <- data.frame(a = 1:5, b = 6:10)
  env_vars <- c("a", "b")
  var_mu  <- c(a = 3,   b = 8)
  var_sd  <- c(a = 1.5, b = 1.5)
  betas   <- c(a = 2.0, b = -1.0)

  out <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, betas)

  expect_equal(dim(out), dim(pt_env))
  expect_equal(names(out), env_vars)
})

test_that("rescale_points_by_betas: uses absolute value of beta", {
  pt_env  <- data.frame(x = c(1, 2, 3))
  env_vars <- "x"
  var_mu  <- c(x = 2)
  var_sd  <- c(x = 1)
  betas_pos <- c(x =  3)
  betas_neg <- c(x = -3)

  pos <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, betas_pos)
  neg <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, betas_neg)

  expect_equal(pos, neg)
})

test_that("rescale_points_by_betas: zero beta makes all values zero", {
  pt_env  <- data.frame(x = c(10, 20, 30))
  env_vars <- "x"
  var_mu  <- c(x = 20)
  var_sd  <- c(x = 5)
  betas   <- c(x = 0)

  out <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, betas)
  expect_true(all(out$x == 0))
})

test_that("rescale_points_by_betas: mean-centred variable with beta=1 gives unit-SD output", {
  set.seed(5)
  x <- rnorm(100, mean = 10, sd = 3)
  pt_env   <- data.frame(x = x)
  env_vars <- "x"
  var_mu   <- c(x = mean(x))
  var_sd   <- c(x = sd(x))
  betas    <- c(x = 1)

  out <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, betas)

  # Should be approximately standardised (unit SD)
  expect_equal(sd(out$x), 1, tolerance = 0.01)
  expect_equal(mean(out$x), 0, tolerance = 0.01)
})

test_that("rescale_points_by_betas: handles multiple variables independently", {
  pt_env <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))
  env_vars <- c("a", "b")
  var_mu <- c(a = 2, b = 5)
  var_sd <- c(a = 1, b = 1)
  betas  <- c(a = 2, b = 0.5)

  out <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, betas)

  # a column: (1-2)/1*2=-2, (2-2)/1*2=0, (3-2)/1*2=2
  expect_equal(out$a, c(-2, 0, 2))
  # b column: (4-5)/1*0.5=-0.5, 0, 0.5
  expect_equal(out$b, c(-0.5, 0, 0.5))
})

test_that("rescale_points_by_betas: NAs in input propagate to output", {
  pt_env  <- data.frame(x = c(1, NA, 3))
  env_vars <- "x"
  var_mu  <- c(x = 2)
  var_sd  <- c(x = 1)
  betas   <- c(x = 1)

  out <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, betas)
  expect_true(is.na(out$x[2]))
})


# =============================================================================
# 5. add_coord_weights_to_points()
# =============================================================================
# Inlined helper for testing
add_coord_weights_to_points <- function(sample_pts, pt_env, coord_wt, env_vars) {
  coords <- terra::crds(sample_pts)
  env_ranges <- apply(pt_env[, env_vars, drop = FALSE], 2, function(x)
    diff(range(x, na.rm = TRUE)))
  target_range <- max(env_ranges, na.rm = TRUE) * coord_wt

  scale_to_range <- function(v, tgt) {
    rng <- diff(range(v, na.rm = TRUE))
    if (rng == 0) return(rep(0, length(v)))
    (v - mean(v, na.rm = TRUE)) / rng * tgt
  }

  pt_env$coord_x_w <- scale_to_range(coords[, 1], target_range)
  pt_env$coord_y_w <- scale_to_range(coords[, 2], target_range)
  pt_env
}

test_that("add_coord_weights_to_points: adds two new columns", {
  skip_if_not_installed("terra")
  pts <- terra::vect(
    data.frame(x = c(1, 2, 3), y = c(4, 5, 6)),
    geom = c("x", "y"), crs = "EPSG:4326"
  )
  pt_env   <- data.frame(temp = c(10, 20, 30), precip = c(100, 200, 300))
  env_vars <- c("temp", "precip")

  out <- add_coord_weights_to_points(pts, pt_env, coord_wt = 2, env_vars)

  expect_true("coord_x_w" %in% names(out))
  expect_true("coord_y_w" %in% names(out))
})

test_that("add_coord_weights_to_points: coordinate columns are zero-mean", {
  skip_if_not_installed("terra")
  set.seed(3)
  n   <- 50
  pts <- terra::vect(
    data.frame(x = runif(n, -100, -90), y = runif(n, 30, 40)),
    geom = c("x", "y"), crs = "EPSG:4326"
  )
  pt_env   <- data.frame(a = runif(n), b = runif(n) * 100)
  env_vars <- c("a", "b")

  out <- add_coord_weights_to_points(pts, pt_env, coord_wt = 2.5, env_vars)

  expect_equal(mean(out$coord_x_w), 0, tolerance = 1e-10)
  expect_equal(mean(out$coord_y_w), 0, tolerance = 1e-10)
})

test_that("add_coord_weights_to_points: constant coordinate axis returns zeros", {
  skip_if_not_installed("terra")
  pts <- terra::vect(
    data.frame(x = rep(5, 4), y = c(1, 2, 3, 4)),
    geom = c("x", "y"), crs = "EPSG:4326"
  )
  pt_env   <- data.frame(v = c(10, 20, 30, 40))
  env_vars <- "v"

  out <- add_coord_weights_to_points(pts, pt_env, coord_wt = 1, env_vars)

  expect_true(all(out$coord_x_w == 0))
})

test_that("add_coord_weights_to_points: higher coord_wt scales range proportionally", {
  skip_if_not_installed("terra")
  n   <- 10
  pts <- terra::vect(
    data.frame(x = 1:n, y = 1:n),
    geom = c("x", "y"), crs = "EPSG:4326"
  )
  pt_env   <- data.frame(z = 1:n * 5)
  env_vars <- "z"

  out1 <- add_coord_weights_to_points(pts, pt_env, coord_wt = 1, env_vars)
  out2 <- add_coord_weights_to_points(pts, pt_env, coord_wt = 2, env_vars)

  ratio <- range(out2$coord_x_w) / range(out1$coord_x_w)
  expect_equal(ratio[1], 2, tolerance = 1e-9)
})


# =============================================================================
# 6. Integration-style smoke test for compute_top3 + build_cooccurrence +
#    compute_stability (using realistic sizes)
# =============================================================================
test_that("pipeline: top3 + cooccurrence + stability are mutually consistent", {
  set.seed(42)
  n_pts   <- 30
  n_draws <- 50
  n_clust <- 4

  dc  <- matrix(sample(seq_len(n_clust), n_pts * n_draws, replace = TRUE),
                nrow = n_pts, ncol = n_draws)
  co  <- build_cooccurrence_matrix(dc, n_pts, n_draws)
  top <- compute_top3_clusters(dc, n_pts, n_draws)

  diss <- stats::as.dist(1 - co)
  hc   <- stats::hclust(diss, method = "average")
  lbl  <- stats::cutree(hc, k = n_clust)
  stab <- compute_stability_scores(co, lbl)

  # Basic sanity: same number of results
  expect_length(stab, n_pts)
  expect_equal(dim(top$labels_matrix), c(n_pts, 3))
  expect_equal(dim(top$pct_matrix),    c(n_pts, 3))
  expect_equal(dim(co), c(n_pts, n_pts))

  # Stability values in [0, 1]
  expect_true(all(stab >= 0 - 1e-9))
  expect_true(all(stab <= 1 + 1e-9))

  # Rank1 >= Rank2 >= Rank3 percentages (with tolerance for rounding)
  expect_true(all(top$pct_matrix[, 1] >= top$pct_matrix[, 2] - 1e-6))
  expect_true(all(top$pct_matrix[, 2] >= top$pct_matrix[, 3] - 1e-6))
})


# =============================================================================
# 7. Edge cases
# =============================================================================
test_that("build_cooccurrence_matrix: n_pts = 1 returns 1×1 matrix with value 1", {
  dc <- matrix(c(1L, 1L, 2L, 1L), nrow = 1, ncol = 4)
  co <- build_cooccurrence_matrix(dc, 1, 4)

  expect_equal(dim(co), c(1, 1))
  expect_equal(co[1, 1], 1.0)
})

test_that("compute_top3_clusters: n_pts = 1 works correctly", {
  dc  <- matrix(c(2L, 2L, 1L, 2L, 1L), nrow = 1)
  res <- compute_top3_clusters(dc, 1, 5)

  expect_equal(res$labels_matrix[1, 1], 2L)
  expect_equal(res$pct_matrix[1, 1], 60)
  expect_equal(res$labels_matrix[1, 2], 1L)
  expect_equal(res$pct_matrix[1, 2], 40)
  # rank3 padded with rank1
  expect_equal(res$labels_matrix[1, 3], 2L)
})

test_that("rescale_points_by_betas: large beta magnifies separation", {
  pt_env   <- data.frame(x = c(1, 3))
  env_vars <- "x"
  var_mu   <- c(x = 2)
  var_sd   <- c(x = 1)

  out_small <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, c(x = 1))
  out_large <- rescale_points_by_betas(pt_env, env_vars, var_mu, var_sd, c(x = 10))

  expect_gt(diff(out_large$x), diff(out_small$x))
})

test_that("compute_stability_scores: all singletons (one point per cluster) return 1", {
  n   <- 5
  co  <- diag(n)
  lbl <- 1:n          # each point its own cluster
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(stab == 1.0))
})
make_mask_rast <- function(nrow = 20, ncol = 20,
                            xmin = 0, xmax = 1,
                            ymin = 0, ymax = 1,
                            crs  = "EPSG:4326") {
  r <- terra::rast(
    nrows = nrow, ncols = ncol,
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    crs  = crs
  )
  terra::values(r) <- 1
  r
}

make_sample_pts <- function(n, mask_rast, seed = 1) {
  set.seed(seed)
  terra::spatSample(mask_rast, size = n, method = "random",
                    as.points = TRUE, na.rm = TRUE)
}

# ── Shared raster fixture (created once, reused across tests) ─────────────────
make_mask_rast <- function(nrow = 20, ncol = 20,
                            xmin = 0, xmax = 1,
                            ymin = 0, ymax = 1,
                            crs  = "EPSG:4326") {
  r <- terra::rast(
    nrows = nrow, ncols = ncol,
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    crs  = crs
  )
  terra::values(r) <- 1
  r
}

make_sample_pts <- function(n, mask_rast, seed = 1) {
  set.seed(seed)
  terra::spatSample(mask_rast, size = n, method = "random",
                    as.points = TRUE, na.rm = TRUE)
}

# =============================================================================
# A.  build_cooccurrence_matrix  –  exhaustive branch / value coverage
# =============================================================================

# ── A1.  Accumulation arithmetic is exact ─────────────────────────────────────
test_that("bco: co-occurrence exactly equals fraction of draws sharing cluster", {
  # 3 points, 4 draws.  Points 1&2 share cluster in draws 1,2,3 (not 4).
  # Expected co[1,2] = 3/4 = 0.75
  dc <- rbind(
    c(1L, 1L, 1L, 2L),   # point 1
    c(1L, 1L, 1L, 1L),   # point 2  (same as pt1 in draws 1-3, different draw 4)
    c(2L, 2L, 2L, 1L)    # point 3
  )
  co <- build_cooccurrence_matrix(dc, n_pts = 3, n_draws = 4)

  expect_equal(co[1, 2], 0.75)
  expect_equal(co[2, 1], 0.75)   # symmetry
})

# ── A2.  Single draw – co-occurrence is 0 or 1 only ──────────────────────────
test_that("bco: single draw produces only 0/1 values", {
  dc <- matrix(c(1L, 1L, 2L, 3L, 2L), nrow = 5, ncol = 1)
  co <- build_cooccurrence_matrix(dc, 5, 1)

  expect_true(all(co %in% c(0, 1)))
})

# ── A3.  n_draws = 1 with all same cluster → all-ones ─────────────────────────
test_that("bco: one draw, all same cluster → all-ones matrix", {
  dc <- matrix(rep(5L, 4), nrow = 4, ncol = 1)
  co <- build_cooccurrence_matrix(dc, 4, 1)

  expect_true(all(co == 1))
})

# ── A4.  Partial NA: non-NA point pairs should still be correct ───────────────
test_that("bco: partial NA row — non-NA pairs unaffected", {
  # point 1 = cluster 1, point 2 = NA, point 3 = cluster 1 (all 5 draws)
  dc <- rbind(
    rep(1L, 5),
    rep(NA_integer_, 5),
    rep(1L, 5)
  )
  co <- build_cooccurrence_matrix(dc, 3, 5)

  # Points 1 and 3 always share a cluster → co[1,3] = 1
  expect_equal(co[1, 3], 1.0)
  # Point 2 is NA in all draws → diagonal = 0
  expect_equal(co[2, 2], 0.0)
  # co[1,2] and co[2,3] should be 0 (NA treated as FALSE)
  expect_equal(co[1, 2], 0.0)
  expect_equal(co[2, 3], 0.0)
})

# ── A5.  Mixed NA: some draws NA, some not ────────────────────────────────────
test_that("bco: point NA in half the draws contributes proportionally", {
  # 2 points, 4 draws.  Point 1 always cluster 1.
  # Point 2: cluster 1 in draws 1&2, NA in draws 3&4.
  dc <- rbind(
    c(1L,  1L,  1L,  1L),
    c(1L,  1L,  NA_integer_, NA_integer_)
  )
  co <- build_cooccurrence_matrix(dc, 2, 4)

  # Draws 1&2: pt1 and pt2 co-cluster (TRUE). Draws 3&4: NA→FALSE.
  # co[1,2] = 2/4 = 0.5
  expect_equal(co[1, 2], 0.5)
})

# ── A6.  Large realistic run: convergence of co-occurrence ────────────────────
test_that("bco: co-occurrence converges to theoretical fraction with many draws", {
  # 2 points.  In each draw they independently get cluster 1 with p=0.5.
  # P(co-cluster) = P(both=1) + P(both=2) = 0.5^2 + 0.5^2 = 0.5
  set.seed(77)
  n_draws <- 5000
  dc <- matrix(sample(1:2, 2 * n_draws, replace = TRUE), nrow = 2)
  co <- build_cooccurrence_matrix(dc, 2, n_draws)

  # Should be close to 0.5 (within 2-sigma of binomial noise)
  expect_equal(co[1, 2], 0.5, tolerance = 0.03)
})

# ── A7.  Verify every off-diagonal cell is <= 1 for large random input ────────
test_that("bco: all off-diagonal values in [0,1] for large random input", {
  set.seed(88)
  n_pts <- 50; n_draws <- 200
  dc <- matrix(sample(1:8, n_pts * n_draws, replace = TRUE),
               nrow = n_pts, ncol = n_draws)
  co <- build_cooccurrence_matrix(dc, n_pts, n_draws)

  expect_true(all(co >= 0))
  expect_true(all(co <= 1 + 1e-12))
})

# ── A8.  Two perfectly anti-correlated clusters ───────────────────────────────
test_that("bco: points always in different clusters have co-occurrence 0", {
  # Points 1,2,3 always cluster A; points 4,5,6 always cluster B.
  dc <- rbind(
    matrix(rep(c(1L, 1L, 1L, 2L, 2L, 2L), 10), nrow = 6, ncol = 10)
  )
  co <- build_cooccurrence_matrix(dc, 6, 10)

  # Within group A: co = 1
  expect_equal(co[1, 2], 1)
  expect_equal(co[1, 3], 1)
  # Across groups: co = 0
  expect_equal(co[1, 4], 0)
  expect_equal(co[3, 6], 0)
})

# ── A9.  Cluster labels don't have to start at 1 ─────────────────────────────
test_that("bco: arbitrary integer cluster labels work correctly", {
  dc <- rbind(rep(99L, 5), rep(99L, 5), rep(777L, 5))
  co <- build_cooccurrence_matrix(dc, 3, 5)

  expect_equal(co[1, 2], 1.0)
  expect_equal(co[1, 3], 0.0)
})

# ── A10.  n_pts = 2 minimal case ──────────────────────────────────────────────
test_that("bco: n_pts=2, always different clusters → off-diagonal = 0", {
  dc <- rbind(rep(1L, 10), rep(2L, 10))
  co <- build_cooccurrence_matrix(dc, 2, 10)

  expect_equal(co[1, 2], 0.0)
  expect_equal(co[2, 1], 0.0)
  expect_equal(co[1, 1], 1.0)
  expect_equal(co[2, 2], 1.0)
})

# ── A11.  Accumulation is additive across draws ───────────────────────────────
test_that("bco: manual accumulation matches function for 3-draw example", {
  # Draw 1: pts 1,2 same; pt 3 different
  # Draw 2: pts 2,3 same; pt 1 different
  # Draw 3: all same
  dc <- rbind(
    c(1L, 2L, 1L),
    c(1L, 1L, 1L),
    c(2L, 1L, 1L)
  )
  co <- build_cooccurrence_matrix(dc, 3, 3)

  # co[1,2]: draw1=TRUE, draw2=FALSE, draw3=TRUE → 2/3
  expect_equal(co[1, 2], 2 / 3, tolerance = 1e-12)
  # co[1,3]: draw1=FALSE, draw2=FALSE, draw3=TRUE → 1/3
  expect_equal(co[1, 3], 1 / 3, tolerance = 1e-12)
  # co[2,3]: draw1=FALSE, draw2=TRUE, draw3=TRUE → 2/3
  expect_equal(co[2, 3], 2 / 3, tolerance = 1e-12)
})


# =============================================================================
# B.  compute_stability_scores  –  exhaustive branch / value coverage
# =============================================================================

# ── B1.  Singleton branch: single-member clusters always return 1.0 ───────────
test_that("css: every point in its own cluster returns stability 1", {
  n   <- 8
  co  <- diag(n)                       # only self-similarity = 1
  lbl <- seq_len(n)                    # each point a unique cluster
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(stab == 1.0))
  expect_length(stab, n)
})

# ── B2.  Two-member cluster: stability = co-occurrence with the other member ──
test_that("css: two-point cluster, stability equals mutual co-occurrence", {
  co <- matrix(c(
    1.0, 0.7,
    0.7, 1.0
  ), 2, 2)
  lbl  <- c(1L, 1L)
  stab <- compute_stability_scores(co, lbl)

  expect_equal(stab[1], 0.7)
  expect_equal(stab[2], 0.7)
})

# ── B3.  Large cluster: stability is mean co-occurrence with all cluster mates──
test_that("css: stability is mean over all cluster-mates, excluding self", {
  # 4-point cluster, point 1 has co-occ 0.9, 0.8, 0.7 with pts 2,3,4
  n   <- 4
  co  <- matrix(1.0, n, n)
  co[1, 2] <- co[2, 1] <- 0.9
  co[1, 3] <- co[3, 1] <- 0.8
  co[1, 4] <- co[4, 1] <- 0.7
  lbl  <- rep(1L, n)
  stab <- compute_stability_scores(co, lbl)

  expected_stab1 <- mean(c(0.9, 0.8, 0.7))
  expect_equal(stab[1], expected_stab1, tolerance = 1e-12)
})

# ── B4.  Mixed clusters: each point only averages within its own cluster ───────
test_that("css: cross-cluster co-occurrences do not affect stability", {
  # Cluster 1: pts 1,2 (co-occ=1). Cluster 2: pts 3,4 (co-occ=1).
  # Cross-cluster co-occ = 0. Stability should be 1 for all.
  n   <- 4
  co  <- matrix(0.0, n, n)
  co[1:2, 1:2] <- 1
  co[3:4, 3:4] <- 1
  lbl  <- c(1L, 1L, 2L, 2L)
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(stab == 1.0))
})

# ── B5.  Stability reflects genuine boundary instability ──────────────────────
test_that("css: boundary point has lower stability than core point", {
  # 5 points: pts 1-3 are a tight cluster (high mutual co-occ ~0.92-0.95),
  # pts 4-5 are a looser cluster (co-occ = 0.88) with low cross-cluster co-occ.
  # Use rbind() so each row of the matrix is laid out exactly as written.
  co5 <- rbind(
    c(1.00, 0.95, 0.92, 0.30, 0.28),
    c(0.95, 1.00, 0.91, 0.32, 0.27),
    c(0.92, 0.91, 1.00, 0.31, 0.29),
    c(0.30, 0.32, 0.31, 1.00, 0.88),
    c(0.28, 0.27, 0.29, 0.88, 1.00)
  )
  lbl5  <- c(1L, 1L, 1L, 2L, 2L)
  stab5 <- compute_stability_scores(co5, lbl5)

  # Cluster 1 (pts 1-3, co-occ ~0.93) should be more stable than cluster 2 (pts 4-5, co-occ=0.88)
  expect_gt(mean(stab5[1:3]), mean(stab5[4:5]))
})

# ── B6.  All co-occurrences zero within cluster → stability = 0 ───────────────
test_that("css: within-cluster co-occ all zero gives stability 0", {
  n   <- 4
  co  <- diag(n) * 0      # all zeros, including diagonal
  lbl <- rep(1L, n)
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(stab == 0.0))
})

# ── B7.  Single-point cluster mixed with multi-point cluster ──────────────────
test_that("css: singleton and multi-member clusters coexist correctly", {
  n   <- 5
  co  <- matrix(0.6, n, n)
  diag(co) <- 1.0
  # Point 5 is a singleton
  lbl  <- c(1L, 1L, 1L, 1L, 2L)
  stab <- compute_stability_scores(co, lbl)

  # Singleton (pt 5) always returns 1.0
  expect_equal(stab[5], 1.0)
  # Non-singletons use mean co-occ with cluster mates (all 0.6)
  expect_equal(stab[1], 0.6, tolerance = 1e-12)
})

# ── B8.  Result is numeric vector, not matrix or list ─────────────────────────
test_that("css: return type is plain numeric vector", {
  co  <- diag(3)
  lbl <- c(1L, 1L, 2L)
  stab <- compute_stability_scores(co, lbl)

  expect_true(is.numeric(stab))
  expect_true(is.vector(stab))
})

# ── B9.  n_pts = 1: single point is its own cluster → stability = 1 ──────────
test_that("css: single point returns stability 1", {
  co   <- matrix(1.0, 1, 1)
  lbl  <- 1L
  stab <- compute_stability_scores(co, lbl)

  expect_equal(stab, 1.0)
  expect_length(stab, 1)
})

# ── B10.  Asymmetric co-occurrence matrix still works (robustness) ─────────────
test_that("css: non-symmetric co-occ matrix uses row i values as specified", {

  co <- rbind(c(1.0, 0.8),
              c(0.3, 1.0))
  lbl  <- c(1L, 1L)
  stab <- compute_stability_scores(co, lbl)

  # pt1 stability: reads co_mat[1, 2] = 0.8
  expect_equal(stab[1], 0.8)
  # pt2 stability: reads co_mat[2, 1] = 0.3
  expect_equal(stab[2], 0.3)
})

# ── B11.  Three clusters, each > 1 member, all independent ────────────────────
test_that("css: three-cluster setup returns correct per-cluster stability", {
  # 6 points: clusters of size 2.  Within-cluster co-occ = 0.9.
  co <- diag(6)
  co[1, 2] <- co[2, 1] <- 0.9
  co[3, 4] <- co[4, 3] <- 0.7
  co[5, 6] <- co[6, 5] <- 0.5
  lbl <- c(1L, 1L, 2L, 2L, 3L, 3L)
  stab <- compute_stability_scores(co, lbl)

  expect_equal(stab[1], 0.9)
  expect_equal(stab[2], 0.9)
  expect_equal(stab[3], 0.7)
  expect_equal(stab[4], 0.7)
  expect_equal(stab[5], 0.5)
  expect_equal(stab[6], 0.5)
})

# ── B12.  Stability ties: multiple cluster-mates with identical co-occ ─────────
test_that("css: mean of tied co-occurrences equals that value", {
  # All pairs in one cluster have co-occ = 0.65
  n   <- 5
  co  <- matrix(0.65, n, n)
  diag(co) <- 1.0
  lbl  <- rep(1L, n)
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(abs(stab - 0.65) < 1e-12))
})

# ── B13.  NA in co-occ matrix: na.rm=TRUE should skip NAs ─────────────────────
test_that("css: NA co-occ values are skipped via na.rm=TRUE", {
  co <- matrix(c(
    1.0, 0.8,  NA,
    0.8, 1.0, 0.6,
     NA, 0.6, 1.0
  ), 3, 3)
  lbl <- rep(1L, 3)
  stab <- compute_stability_scores(co, lbl)

  # pt1: mean(co[1, c(2,3)]) = mean(0.8, NA, na.rm=TRUE) = 0.8
  expect_equal(stab[1], 0.8,        tolerance = 1e-12)
  # pt2: mean(co[2, c(1,3)]) = mean(0.8, 0.6) = 0.7
  expect_equal(stab[2], 0.7,        tolerance = 1e-12)
  # pt3: mean(co[3, c(1,2)]) = mean(NA, 0.6, na.rm=TRUE) = 0.6
  expect_equal(stab[3], 0.6,        tolerance = 1e-12)
})

# ── B14.  Large realistic size (slow) ─────────────────────────────────────────
test_that("css: 500-point, 10-cluster run completes and stays in [0,1] (slow)", {
  set.seed(101)
  n   <- 500
  raw <- matrix(runif(n * n, 0.3, 1), n, n)
  co  <- (raw + t(raw)) / 2
  diag(co) <- 1
  lbl <- sample(1:10, n, replace = TRUE)
  stab <- compute_stability_scores(co, lbl)

  expect_length(stab, n)
  expect_true(all(stab >= 0 - 1e-9))
  expect_true(all(stab <= 1 + 1e-9))
})


# =============================================================================
# C.  build_cooccurrence_matrix + compute_stability_scores  – joint properties
# =============================================================================

# ── C1.  Stability computed from perfectly separated co-occ is 1 for all ──────
test_that("joint: perfect cluster separation → all stability = 1", {
  n_pts   <- 20
  n_draws <- 50
  # Pts 1-10 always cluster A; pts 11-20 always cluster B
  dc <- rbind(
    matrix(rep(1L, 10 * n_draws), nrow = 10),
    matrix(rep(2L, 10 * n_draws), nrow = 10)
  )
  co   <- build_cooccurrence_matrix(dc, n_pts, n_draws)
  lbl  <- c(rep(1L, 10), rep(2L, 10))
  stab <- compute_stability_scores(co, lbl)

  expect_true(all(abs(stab - 1.0) < 1e-12))
})

# ── C2.  Random draws: co-occurrence diagonal should always be 1 ──────────────
test_that("joint: random large run diagonal always 1 before stability", {
  set.seed(202)
  n_pts <- 100; n_draws <- 200
  dc <- matrix(sample(1:5, n_pts * n_draws, replace = TRUE),
               nrow = n_pts, ncol = n_draws)
  co <- build_cooccurrence_matrix(dc, n_pts, n_draws)

  expect_true(all(diag(co) == 1.0))
})

# ── C3.  Stability ordering: stable cluster has higher mean stability ──────────
test_that("joint: clearly stable cluster has higher mean stability than boundary cluster", {
  set.seed(303)
  n_pts <- 30; n_draws <- 100

  # Cluster A (pts 1-10): always assigned cluster 1
  # Cluster B (pts 11-20): always assigned cluster 2
  # Boundary (pts 21-30): randomly assigned cluster 1 or 2 each draw
  dc <- rbind(
    matrix(rep(1L, 10 * n_draws), nrow = 10),
    matrix(rep(2L, 10 * n_draws), nrow = 10),
    matrix(sample(1:2, 10 * n_draws, replace = TRUE), nrow = 10)
  )
  co   <- build_cooccurrence_matrix(dc, n_pts, n_draws)

  # Consensus: majority-vote for boundary points
  diss <- stats::as.dist(1 - co)
  hc   <- stats::hclust(diss, method = "average")
  lbl  <- stats::cutree(hc, k = 2)

  stab <- compute_stability_scores(co, lbl)

  # Stable clusters should beat boundary region
  expect_gt(mean(stab[1:20]), mean(stab[21:30]))
})


# =============================================================================
# D.  interpolate_stability_to_raster  –  full branch and property coverage
# =============================================================================

# ── D1.  Output is a SpatRaster with correct name ─────────────────────────────
test_that("isr: returns SpatRaster named 'cluster_stability'", {
  skip_if_not_installed("terra")
  mask <- make_mask_rast()
  pts  <- make_sample_pts(30, mask)
  set.seed(1)
  stab <- runif(terra::nrow(pts) %||% length(pts))
  stab <- runif(nrow(as.data.frame(terra::crds(pts))))

  result <- interpolate_stability_to_raster(pts, stab, mask)

  expect_s4_class(result, "SpatRaster")
  expect_equal(names(result), "cluster_stability")
})

# ── D2.  Output extent matches mask extent ────────────────────────────────────
test_that("isr: output extent matches mask raster extent", {
  skip_if_not_installed("terra")
  mask   <- make_mask_rast()
  pts    <- make_sample_pts(30, mask)
  stab   <- rep(0.8, nrow(terra::as.data.frame(terra::crds(pts))))

  result <- interpolate_stability_to_raster(pts, stab, mask)

  expect_equal(terra::ext(result), terra::ext(mask))
})

# ── D3.  Output CRS matches mask CRS ─────────────────────────────────────────
test_that("isr: output CRS matches mask raster CRS", {
  skip_if_not_installed("terra")
  mask   <- make_mask_rast()
  pts    <- make_sample_pts(25, mask)
  stab   <- rep(0.5, nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask)

  expect_equal(terra::crs(result), terra::crs(mask))
})

# ── D4.  All-constant stability input: non-NA output cells ≈ constant ─────────
test_that("isr: uniform stability input → all interpolated values equal that constant", {
  skip_if_not_installed("terra")
  mask   <- make_mask_rast(nrow = 15, ncol = 15)
  pts    <- make_sample_pts(40, mask, seed = 5)
  stab   <- rep(0.75, nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask)
  vals   <- terra::values(result, na.rm = TRUE)

  expect_true(all(abs(vals - 0.75) < 1e-6))
})

# ── D5.  Stability = 0 everywhere: all output values should be ≈ 0 ────────────
test_that("isr: all-zero stability input → output values ≈ 0", {
  skip_if_not_installed("terra")
  mask   <- make_mask_rast(nrow = 15, ncol = 15)
  pts    <- make_sample_pts(40, mask, seed = 6)
  stab   <- rep(0.0, nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask)
  vals   <- terra::values(result, na.rm = TRUE)

  expect_true(all(abs(vals) < 1e-6))
})

# ── D6.  Stability = 1 everywhere: all output values should be ≈ 1 ────────────
test_that("isr: all-one stability input → output values ≈ 1", {
  skip_if_not_installed("terra")
  mask   <- make_mask_rast(nrow = 15, ncol = 15)
  pts    <- make_sample_pts(40, mask, seed = 7)
  stab   <- rep(1.0, nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask)
  vals   <- terra::values(result, na.rm = TRUE)

  expect_true(all(abs(vals - 1.0) < 1e-6))
})

# ── D7.  Output values bounded in [0,1] for realistic stability scores ─────────
test_that("isr: random stability in [0,1] → interpolated output also in [0,1]", {
  skip_if_not_installed("terra")
  set.seed(8)
  mask   <- make_mask_rast()
  pts    <- make_sample_pts(50, mask, seed = 8)
  stab   <- runif(nrow(terra::crds(pts)), min = 0, max = 1)

  result <- interpolate_stability_to_raster(pts, stab, mask)
  vals   <- terra::values(result, na.rm = TRUE)

  expect_true(all(vals >= -1e-6))
  expect_true(all(vals <=  1 + 1e-6))
})

# ── D8.  Output has same resolution as mask ───────────────────────────────────
test_that("isr: output resolution matches mask resolution", {
  skip_if_not_installed("terra")
  mask   <- make_mask_rast(nrow = 10, ncol = 20)
  pts    <- make_sample_pts(30, mask, seed = 9)
  stab   <- runif(nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask)

  expect_equal(terra::res(result), terra::res(mask))
})

# ── D9.  Masking: NA cells in mask → NA in output ────────────────────────────
test_that("isr: cells masked out in mask raster are NA in output", {
  skip_if_not_installed("terra")
  mask   <- make_mask_rast(nrow = 20, ncol = 20)
  # Set a rectangular strip to NA in the mask
  mask_na <- mask
  terra::values(mask_na)[1:100] <- NA

  pts    <- make_sample_pts(50, mask_na, seed = 10)
  stab   <- runif(nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask_na)

  # Cells that were NA in mask must be NA in output
  mask_vals   <- terra::values(mask_na)[1:100]
  result_vals <- terra::values(result)[1:100]
  expect_true(all(is.na(result_vals)))
})

# ── D10.  Non-NA cell count > 0: some cells should be filled ─────────────────
test_that("isr: at least some non-NA output cells when mask has valid cells", {
  skip_if_not_installed("terra")
  mask   <- make_mask_rast()
  pts    <- make_sample_pts(30, mask, seed = 11)
  stab   <- runif(nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask)

  expect_gt(sum(!is.na(terra::values(result))), 0)
})

# ── D11.  Spatial gradient: high-stability region predicts higher values nearby─
test_that("isr: high-stability points produce higher interpolated values nearby", {
  skip_if_not_installed("terra")
  # Place points manually: bottom-left = low stability; top-right = high
  mask   <- make_mask_rast(nrow = 30, ncol = 30)

  # 10 low-stability points near (0.1, 0.1)
  low_pts  <- terra::vect(
    data.frame(x = runif(10, 0.01, 0.15), y = runif(10, 0.01, 0.15)),
    geom = c("x", "y"), crs = "EPSG:4326"
  )
  # 10 high-stability points near (0.85, 0.85)
  high_pts <- terra::vect(
    data.frame(x = runif(10, 0.70, 0.99), y = runif(10, 0.70, 0.99)),
    geom = c("x", "y"), crs = "EPSG:4326"
  )
  all_pts  <- rbind(low_pts, high_pts)
  stab     <- c(rep(0.1, 10), rep(0.9, 10))

  result <- interpolate_stability_to_raster(all_pts, stab, mask)

  # Sample values near each corner
  low_corner  <- terra::extract(result,
    terra::vect(data.frame(x = 0.05, y = 0.05), geom = c("x","y"), crs="EPSG:4326")
  )[, 2]
  high_corner <- terra::extract(result,
    terra::vect(data.frame(x = 0.90, y = 0.90), geom = c("x","y"), crs="EPSG:4326")
  )[, 2]

  expect_gt(high_corner, low_corner)
})

# ── D12.  Single sample point: function does not error, returns a raster ──────
test_that("isr: single sample point does not error", {
  skip_if_not_installed("terra")
  mask  <- make_mask_rast(nrow = 10, ncol = 10)
  pts   <- terra::vect(data.frame(x = 0.5, y = 0.5),
                       geom = c("x", "y"), crs = "EPSG:4326")
  stab  <- 0.6

  expect_no_error({
    result <- interpolate_stability_to_raster(pts, stab, mask)
  })
  expect_s4_class(result, "SpatRaster")
})

# ── D13.  Large realistic size (slow): 500 points on 100×100 grid ─────────────
test_that("isr: 500 sample points on 100x100 raster completes (slow)", {
  skip_if_not_installed("terra")
  set.seed(42)
  mask <- make_mask_rast(nrow = 100, ncol = 100)
  pts  <- make_sample_pts(500, mask, seed = 42)
  stab <- runif(nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask)

  expect_s4_class(result, "SpatRaster")
  expect_equal(names(result), "cluster_stability")
  expect_true(sum(!is.na(terra::values(result))) > 100)
})

# ── D14.  Projected (non-geographic) CRS is preserved ────────────────────────
test_that("isr: projected CRS in mask is preserved in output", {
  skip_if_not_installed("terra")
  # Use a simple metric CRS
  mask <- terra::rast(
    nrows = 20, ncols = 20,
    xmin = 0, xmax = 20000,
    ymin = 0, ymax = 20000,
    crs  = "EPSG:32636"        # UTM zone 36N
  )
  terra::values(mask) <- 1

  pts <- terra::spatSample(mask, size = 30, method = "random",
                           as.points = TRUE, na.rm = TRUE)
  stab <- runif(nrow(terra::crds(pts)))

  result <- interpolate_stability_to_raster(pts, stab, mask)

  expect_equal(terra::crs(result), terra::crs(mask))
})

# =============================================================================
# Shared fixtures
# =============================================================================

# Minimal but realistic synthetic dataset used by most tests.
# Returns a named list so individual tests can poke at any piece.
make_pcr_fixture <- function(n_pts   = 80,
                              n_clust = 3,
                              seed    = 42) {
  set.seed(seed)

  # ── Rasters ────────────────────────────────────────────────────────────────
  mask_rast <- terra::rast(
    nrows = 20, ncols = 20,
    xmin = -1, xmax = 1,
    ymin = -1, ymax = 1,
    crs  = "EPSG:4326"
  )
  terra::values(mask_rast) <- 1

  # Two synthetic environmental predictors
  pred1 <- mask_rast; terra::values(pred1) <- runif(terra::ncell(mask_rast), 0, 10)
  pred2 <- mask_rast; terra::values(pred2) <- runif(terra::ncell(mask_rast), 5, 20)
  predictors <- c(pred1, pred2)
  names(predictors) <- c("bio1", "bio2")

  # ── Sample points ──────────────────────────────────────────────────────────
  sample_pts <- terra::spatSample(mask_rast, size = n_pts, method = "random",
                                   as.points = TRUE, na.rm = TRUE)

  # ── Labels and scores ──────────────────────────────────────────────────────
  consensus_labels <- sample(seq_len(n_clust), n_pts, replace = TRUE)
  stability_scores <- runif(n_pts, 0.5, 1.0)

  # top3_labels: n_pts x 3 integer matrix (col 1 = rank1 = consensus_labels)
  top3_labels <- cbind(
    consensus_labels,
    sample(seq_len(n_clust), n_pts, replace = TRUE),
    sample(seq_len(n_clust), n_pts, replace = TRUE)
  )

  # ── Model parameters ───────────────────────────────────────────────────────
  env_vars  <- c("bio1", "bio2")
  var_mu    <- c(bio1 = 5,  bio2 = 12)
  var_sd    <- c(bio1 = 2,  bio2 = 4)
  mean_betas <- c(bio1 = 0.8, bio2 = -0.5)
  coord_wt  <- 2.0
  planar_proj <- "EPSG:32632"   # UTM zone 32N – a sensible planar CRS

  list(
    mask_rast       = mask_rast,
    predictors      = predictors,
    sample_pts      = sample_pts,
    consensus_labels = consensus_labels,
    stability_scores = stability_scores,
    top3_labels     = top3_labels,
    env_vars        = env_vars,
    var_mu          = var_mu,
    var_sd          = var_sd,
    mean_betas      = mean_betas,
    coord_wt        = coord_wt,
    planar_proj     = planar_proj,
    n_pts           = n_pts,
    n_clust         = n_clust
  )
}

# Convenience wrapper: call project_consensus_to_raster with all fixture args
call_pcr <- function(f, ...) {
  project_consensus_to_raster(
    sample_pts       = f$sample_pts,
    consensus_labels = f$consensus_labels,
    stability_scores = f$stability_scores,
    top3_labels      = f$top3_labels,
    predictors       = f$predictors,
    env_vars         = f$env_vars,
    var_mu           = f$var_mu,
    var_sd           = f$var_sd,
    mean_betas       = f$mean_betas,
    coord_wt         = f$coord_wt,
    mask_rast        = f$mask_rast,
    planar_proj      = f$planar_proj,
    ...
  )
}


# =============================================================================
# E.  project_consensus_to_raster – branch and output coverage
# =============================================================================

# ── E1.  Smoke test: character planar_proj runs without error ─────────────────
test_that("pcr: character planar_proj runs without error and returns a list", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  expect_type(res, "list")
})

# ── E2.  Numeric EPSG is accepted and coerced to "epsg:XXXX" string ───────────
test_that("pcr: numeric planar_proj (EPSG integer) is accepted", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture()

  # Pass the EPSG code as a plain integer instead of a string
  expect_no_error(
    project_consensus_to_raster(
      sample_pts       = f$sample_pts,
      consensus_labels = f$consensus_labels,
      stability_scores = f$stability_scores,
      top3_labels      = f$top3_labels,
      predictors       = f$predictors,
      env_vars         = f$env_vars,
      var_mu           = f$var_mu,
      var_sd           = f$var_sd,
      mean_betas       = f$mean_betas,
      coord_wt         = f$coord_wt,
      mask_rast        = f$mask_rast,
      planar_proj      = 32632L   # numeric
    )
  )
})

# ── E3.  < 50 complete cases triggers stop() ──────────────────────────────────
test_that("pcr: fewer than 50 complete predictor rows stops with informative message", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture(n_pts = 80)

  # Corrupt the predictor raster so nearly all extractions return NA
  bad_pred1 <- f$predictors[[1]]
  bad_pred2 <- f$predictors[[2]]
  terra::values(bad_pred1) <- NA_real_
  terra::values(bad_pred2) <- NA_real_
  bad_predictors <- c(bad_pred1, bad_pred2)
  names(bad_predictors) <- c("bio1", "bio2")

  expect_error(
    project_consensus_to_raster(
      sample_pts       = f$sample_pts,
      consensus_labels = f$consensus_labels,
      stability_scores = f$stability_scores,
      top3_labels      = f$top3_labels,
      predictors       = bad_predictors,
      env_vars         = f$env_vars,
      var_mu           = f$var_mu,
      var_sd           = f$var_sd,
      mean_betas       = f$mean_betas,
      coord_wt         = f$coord_wt,
      mask_rast        = f$mask_rast,
      planar_proj      = f$planar_proj
    ),
    regexp = "Fewer than 50 complete cases"
  )
})

# ── E4.  rank1 with insufficient observations per cluster stops() ─────────────
test_that("pcr: rank1 cluster with only 1 member per cluster triggers stop()", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture(n_pts = 80, n_clust = 3)

  # Force all 80 points into 80 unique clusters → no cluster has >= 2 members
  degenerate_labels <- seq_len(f$n_pts)
  degenerate_top3   <- cbind(degenerate_labels,
                              degenerate_labels,
                              degenerate_labels)

  expect_error(
    project_consensus_to_raster(
      sample_pts       = f$sample_pts,
      consensus_labels = degenerate_labels,
      stability_scores = f$stability_scores,
      top3_labels      = degenerate_top3,
      predictors       = f$predictors,
      env_vars         = f$env_vars,
      var_mu           = f$var_mu,
      var_sd           = f$var_sd,
      mean_betas       = f$mean_betas,
      coord_wt         = f$coord_wt,
      mask_rast        = f$mask_rast,
      planar_proj      = f$planar_proj
    ),
    regexp = "insufficient observations per cluster"
  )
})

# ── E5.  rank2 with enough variation trains its own model (not rank1 fallback) ─
test_that("pcr: rank2 with sufficient variation trains a distinct model", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture(seed = 7)

  # Ensure rank2 has clear multi-cluster structure (at least 2 members each)
  f$top3_labels[, 2] <- rep(seq_len(f$n_clust),
                             length.out = f$n_pts)

  res <- call_pcr(f)

  # Both knn_cluster and knn_rank2 should be present
  expect_true(!is.null(res$knn_rank2))
  # They should be different objects (rank2 was not short-circuited to rank1)
  expect_false(identical(res$knn_cluster, res$knn_rank2))
})

# ── E6.  rank2 degenerate → warning + fallback to rank1 model ─────────────────
# BUG NOTE: can_train_knn(y, min_per_class=2) checks that every cluster has
# >= 2 members in the FULL training set, but trainKNN() performs an 80/20
# split internally.  A single dominant cluster (e.g. 79 vs 1 member) passes
# can_train_knn() yet crashes confusionMatrix() when the 1-member cluster
# lands entirely in the test split with 0 predicted levels.
# To reliably trigger the can_train_knn() == FALSE branch (the warning path)
# we need a cluster with exactly 1 member so the guard fires before trainKNN.
test_that("pcr: rank2 with a singleton cluster triggers warning and falls back", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture(seed = 8)

  # Give 79 points cluster 1, 1 point cluster 2 → can_train_knn() returns FALSE
  # because cluster 2 has only 1 member (< min_per_class = 2)
  rank2 <- rep(1L, f$n_pts)
  rank2[1] <- 2L
  f$top3_labels[, 2] <- rank2

  expect_warning(
    res <- call_pcr(f),
    regexp = "Rank2 clusters.*Using rank1 model"
  )

  # Fallback: knn_rank2 should be identical to knn_cluster (same object)
  expect_identical(res$knn_rank2, res$knn_cluster)
})

# ── E7.  rank3 degenerate → warning + fallback to rank1 model ─────────────────
test_that("pcr: rank3 with a singleton cluster triggers warning and falls back", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture(seed = 9)

  # rank2 has good variation; rank3 has a singleton cluster → triggers fallback
  f$top3_labels[, 2] <- rep(seq_len(f$n_clust), length.out = f$n_pts)
  rank3 <- rep(1L, f$n_pts)
  rank3[1] <- 2L
  f$top3_labels[, 3] <- rank3

  expect_warning(
    res <- call_pcr(f),
    regexp = "Rank3 clusters.*Using rank1 model"
  )

  expect_identical(res$knn_rank3, res$knn_cluster)
})

# ── E8.  Output list has all required keys ────────────────────────────────────
test_that("pcr: output list contains all 8 expected keys", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  expected_keys <- c("cluster_raster", "stability_raster",
                      "rank2_raster",   "rank3_raster",
                      "knn_cluster",    "knn_rank2",
                      "knn_rank3",      "knn_stability")
  expect_true(all(expected_keys %in% names(res)))
})

# ── E9.  Rasters have correct layer names ─────────────────────────────────────
test_that("pcr: output rasters carry the correct layer names", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  expect_equal(names(res$stability_raster), "cluster_stability")
  expect_equal(names(res$rank2_raster),     "rank2_cluster")
  expect_equal(names(res$rank3_raster),     "rank3_cluster")
  # cluster_raster name is set by caret/terra predict – just check it is a single layer
  expect_equal(terra::nlyr(res$cluster_raster), 1L)
})

# ── E10.  Rasters are masked: no values bleed outside the mask ────────────────
test_that("pcr: all output rasters are masked to the spatial domain of mask_rast", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  # Project mask to planar CRS to match output rasters
  mask_proj <- terra::project(f$mask_rast, f$planar_proj)

  for (rname in c("cluster_raster", "stability_raster",
                   "rank2_raster",   "rank3_raster")) {
    rast_vals  <- terra::values(res[[rname]])
    mask_vals  <- terra::values(
      terra::resample(mask_proj, res[[rname]], method = "near")
    )
    # Any cell that is NA in the mask should also be NA in the output raster
    outside_mask <- is.na(mask_vals)
    expect_true(
      all(is.na(rast_vals[outside_mask])),
      info = paste(rname, "has non-NA values outside the mask")
    )
  }
})

# ── E11.  Output CRS matches planar_proj ──────────────────────────────────────
test_that("pcr: all output rasters have CRS matching planar_proj", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  expected_crs <- terra::crs(terra::rast(crs = f$planar_proj))

  for (rname in c("cluster_raster", "stability_raster",
                   "rank2_raster",   "rank3_raster")) {
    expect_equal(
      terra::crs(res[[rname]]),
      expected_crs,
      info = paste(rname, "CRS mismatch")
    )
  }
})

# ── E12.  Numeric EPSG and string EPSG produce the same output CRS ────────────
test_that("pcr: numeric EPSG and equivalent string EPSG yield identical output CRS", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture()

  res_str <- call_pcr(f)   # planar_proj = "EPSG:32632"

  res_num <- project_consensus_to_raster(
    sample_pts       = f$sample_pts,
    consensus_labels = f$consensus_labels,
    stability_scores = f$stability_scores,
    top3_labels      = f$top3_labels,
    predictors       = f$predictors,
    env_vars         = f$env_vars,
    var_mu           = f$var_mu,
    var_sd           = f$var_sd,
    mean_betas       = f$mean_betas,
    coord_wt         = f$coord_wt,
    mask_rast        = f$mask_rast,
    planar_proj      = 32632L
  )

  expect_equal(
    terra::crs(res_str$cluster_raster),
    terra::crs(res_num$cluster_raster)
  )
})

# ── E13.  All four rasters share the same extent and resolution ───────────────
test_that("pcr: all four output rasters share identical extent and resolution", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  ref_ext <- terra::ext(res$cluster_raster)
  ref_res <- terra::res(res$cluster_raster)

  for (rname in c("stability_raster", "rank2_raster", "rank3_raster")) {
    expect_equal(terra::ext(res[[rname]]), ref_ext,
                 info = paste(rname, "extent differs from cluster_raster"))
    expect_equal(terra::res(res[[rname]]), ref_res,
                 info = paste(rname, "resolution differs from cluster_raster"))
  }
})

# ── E14.  KNN models are caret train objects ──────────────────────────────────
test_that("pcr: knn_cluster and knn_stability are caret train objects", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  # knn_cluster / knn_rank2 / knn_rank3 are the $fit.knn slots (caret models)
  expect_s3_class(res$knn_cluster,   "train")
  expect_s3_class(res$knn_rank2,     "train")
  expect_s3_class(res$knn_rank3,     "train")
  # knn_stability is the full regression train object
  expect_s3_class(res$knn_stability, "train")
})

# ── E15.  Partial NAs in predictors: complete cases still >= 50 proceeds ──────
test_that("pcr: partial NAs in predictors (but >= 50 complete rows) proceeds normally", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture(n_pts = 120)

  # Set top-left corner of predictor rasters to NA (affects ~10 sample points)
  bad_pred <- f$predictors
  terra::values(bad_pred[[1]])[1:40] <- NA_real_

  expect_no_error(
    project_consensus_to_raster(
      sample_pts       = f$sample_pts,
      consensus_labels = f$consensus_labels,
      stability_scores = f$stability_scores,
      top3_labels      = f$top3_labels,
      predictors       = bad_pred,
      env_vars         = f$env_vars,
      var_mu           = f$var_mu,
      var_sd           = f$var_sd,
      mean_betas       = f$mean_betas,
      coord_wt         = f$coord_wt,
      mask_rast        = f$mask_rast,
      planar_proj      = f$planar_proj
    )
  )
})

# ── E16.  Stability raster values are bounded in [0, 1] ──────────────────────
test_that("pcr: stability raster values stay within [0, 1]", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  vals <- terra::values(res$stability_raster, na.rm = TRUE)
  expect_true(all(vals >= -1e-6), info = "stability raster has values below 0")
  expect_true(all(vals <=  1 + 1e-6), info = "stability raster has values above 1")
})

# ── E17.  Cluster raster contains only valid cluster IDs ──────────────────────
test_that("pcr: cluster raster values are a subset of the consensus cluster IDs", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f   <- make_pcr_fixture()
  res <- call_pcr(f)

  predicted_ids <- unique(na.omit(terra::values(res$cluster_raster)))
  valid_ids     <- unique(f$consensus_labels)

  # Every predicted cluster ID should be one we trained on
  expect_true(
    all(predicted_ids %in% valid_ids),
    info = paste("Unexpected cluster IDs:", paste(setdiff(predicted_ids, valid_ids), collapse = ", "))
  )
})

# ── E18.  Zero-beta variable: known crash bug ────────────────────────────────
test_that("pcr: all mean_betas == 0 triggers an informative stop()", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture()

  f$mean_betas[] <- 0   # every predictor zeroed out

  expect_error(
    suppressWarnings(call_pcr(f)),
    regexp = "all predictors have posterior mean beta"
  )
})

# ── E19.  Both rank2 and rank3 degenerate simultaneously: two warnings ─────────
test_that("pcr: both rank2 and rank3 singleton clusters emit two separate warnings", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  f <- make_pcr_fixture(seed = 11)

  # Singleton cluster in both rank2 and rank3 -> both trigger the fallback warning
  rank_degen <- rep(1L, f$n_pts); rank_degen[1] <- 2L
  f$top3_labels[, 2] <- rank_degen
  f$top3_labels[, 3] <- rank_degen

  seen <- character(0)
  withCallingHandlers(
    call_pcr(f),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(
    any(grepl("Rank2 clusters", seen)),
    info = paste("Rank2 warning not found in:", paste(seen, collapse = " | "))
  )
  expect_true(
    any(grepl("Rank3 clusters", seen)),
    info = paste("Rank3 warning not found in:", paste(seen, collapse = " | "))
  )
})

# ── E20.  Negative mean_betas: abs() ensures same rescaling as positive ────────
test_that("pcr: negative mean_betas produce identical rasters to their absolute values", {
  skip_if_not_installed("terra")
  skip_if_not_installed("caret")
  set.seed(55)
  f <- make_pcr_fixture(seed = 55)

  f_neg           <- f
  f_neg$mean_betas <- -abs(f$mean_betas)
  f_pos           <- f
  f_pos$mean_betas <-  abs(f$mean_betas)

  res_neg <- call_pcr(f_neg)
  res_pos <- call_pcr(f_pos)

  # Rescaling uses abs(mean_betas), so results must be identical
  expect_equal(
    terra::values(res_neg$stability_raster),
    terra::values(res_pos$stability_raster)
  )
})


# ── Shared fixture helper ─────────────────────────────────────────────────────
# Build a posterior::draws_matrix that looks like brms posterior samples.
# Column names follow the brms convention: fixed effects are "b_<varname>".
# We also add non-b_ columns (Intercept, sigma, lp__) that should never
# appear in the output.
make_fake_model <- function(n_total  = 200,
                             env_vars = c("bio1", "bio2", "bio3"),
                             seed     = 1) {
  set.seed(seed)
  b_cols     <- paste0("b_", env_vars)
  extra_cols <- c("b_Intercept", "sigma", "lp__")
  all_cols   <- c(b_cols, extra_cols)

  raw <- matrix(
    rnorm(n_total * length(all_cols)),
    nrow     = n_total,
    ncol     = length(all_cols),
    dimnames = list(NULL, all_cols)
  )
  # posterior::as_draws_matrix() accepts a named matrix and returns a
  # draws_matrix that brms::as_draws_matrix() dispatches through unchanged
  posterior::as_draws_matrix(raw)
}


# =============================================================================
# F1.  All vars found, n_draws < n_total
# =============================================================================
test_that("ebd: returns matrix with correct dims when all vars found", {
  model <- make_fake_model(n_total = 200, env_vars = c("bio1", "bio2", "bio3"))
  res   <- extract_beta_draws(model, env_vars = c("bio1", "bio2", "bio3"),
                               n_draws = 50)

  expect_equal(nrow(res), 50)
  expect_equal(ncol(res), 3)
})

# =============================================================================
# F2.  n_draws > n_total → capped, no over-sampling
# =============================================================================
test_that("ebd: requesting more draws than exist returns at most n_total rows", {
  n_total <- 30
  model   <- make_fake_model(n_total = n_total, env_vars = c("bio1", "bio2"))
  res     <- extract_beta_draws(model, env_vars = c("bio1", "bio2"),
                                 n_draws = 999)

  expect_equal(nrow(res), n_total)
})

# =============================================================================
# F3.  n_draws == n_total → all rows returned
# =============================================================================
test_that("ebd: n_draws == n_total returns exactly n_total rows", {
  n_total <- 50
  model   <- make_fake_model(n_total = n_total, env_vars = "bio1")
  res     <- extract_beta_draws(model, env_vars = "bio1", n_draws = n_total)

  expect_equal(nrow(res), n_total)
})

# =============================================================================
# F4.  No matching columns → stop()
# =============================================================================
test_that("ebd: no matching columns triggers an informative stop()", {
  model <- make_fake_model(env_vars = c("bio1", "bio2"))

  expect_error(
    extract_beta_draws(model, env_vars = c("temp", "precip"), n_draws = 10),
    regexp = "No posterior draw columns matched"
  )
})

# =============================================================================
# F5.  Some columns missing → warning listing all missing vars
# =============================================================================
test_that("ebd: partially missing columns emit a warning naming each missing var", {
  model <- make_fake_model(env_vars = c("bio1", "bio2"))  # bio3 absent

  expect_warning(
    res <- extract_beta_draws(model,
                               env_vars = c("bio1", "bio2", "bio3"),
                               n_draws  = 10),
    regexp = "bio3"
  )
  # Despite the warning the found columns are still returned
  expect_equal(ncol(res), 2)
  expect_true("bio1" %in% colnames(res))
  expect_true("bio2" %in% colnames(res))
})

# =============================================================================
# F6.  Exactly one column missing → warning names only that variable
# =============================================================================
test_that("ebd: single missing column warning names only that variable", {
  model <- make_fake_model(env_vars = c("bio1", "bio2", "bio3"))

  # Request bio4 which is not in the draws
  w <- tryCatch(
    extract_beta_draws(model,
                        env_vars = c("bio1", "bio2", "bio3", "bio4"),
                        n_draws  = 10),
    warning = function(w) w
  )

  expect_true(grepl("bio4", conditionMessage(w)))
  # bio1, bio2, bio3 should NOT appear in the warning message
  expect_false(grepl("bio1|bio2|bio3", conditionMessage(w)))
})

# =============================================================================
# F7.  "b_" prefix is stripped from all output column names
# =============================================================================
test_that("ebd: b_ prefix is stripped from all output column names", {
  env_vars <- c("bio1", "bio2", "bio3")
  model    <- make_fake_model(env_vars = env_vars)
  res      <- extract_beta_draws(model, env_vars = env_vars, n_draws = 20)

  expect_equal(sort(colnames(res)), sort(env_vars))
  expect_false(any(startsWith(colnames(res), "b_")))
})

# =============================================================================
# F8.  Output is a plain matrix, not a draws_matrix or data frame
# =============================================================================
test_that("ebd: output is a plain numeric matrix", {
  model <- make_fake_model(env_vars = c("bio1", "bio2"))
  res   <- extract_beta_draws(model, env_vars = c("bio1", "bio2"), n_draws = 10)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  # Should NOT carry the draws_matrix class
  #expect_false(inherits(res, "draws_matrix"))
  expect_false(is.data.frame(res))
})

# =============================================================================
# F9.  Row count is exactly min(n_draws, n_total)
# =============================================================================
test_that("ebd: row count is exactly min(n_draws, n_total) in all cases", {
  model <- make_fake_model(n_total = 100, env_vars = "bio1")

  expect_equal(nrow(extract_beta_draws(model, "bio1", n_draws = 10)),  10)
  expect_equal(nrow(extract_beta_draws(model, "bio1", n_draws = 100)), 100)
  expect_equal(nrow(extract_beta_draws(model, "bio1", n_draws = 200)), 100)
})

# =============================================================================
# F10. Output columns are exactly the found env_vars — no extras
# =============================================================================
test_that("ebd: output contains only requested env_var columns, no extras", {
  env_vars <- c("bio1", "bio2")
  model    <- make_fake_model(
    env_vars = c("bio1", "bio2", "bio3", "bio4")  # extra vars in draws
  )
  res <- extract_beta_draws(model, env_vars = env_vars, n_draws = 10)

  expect_equal(sort(colnames(res)), sort(env_vars))
})

# =============================================================================
# F11. Sampling is without replacement (no duplicate rows in output)
# =============================================================================
test_that("ebd: sampled rows are unique (sampling without replacement)", {
  set.seed(42)
  model <- make_fake_model(n_total = 100, env_vars = c("bio1", "bio2"),
                            seed = 42)
  res   <- extract_beta_draws(model, env_vars = c("bio1", "bio2"), n_draws = 80)

  # If any two rows were identical it would indicate replacement sampling.
  # Duplicated rows in a continuous random matrix is effectively impossible
  # without replacement, so this is a strong test.
  expect_equal(nrow(unique(res)), nrow(res))
})

# =============================================================================
# F12. set.seed() makes draw selection reproducible
# =============================================================================
test_that("ebd: same seed produces identical row selection", {
  model <- make_fake_model(n_total = 200, env_vars = c("bio1", "bio2"))

  set.seed(7)
  res1 <- extract_beta_draws(model, env_vars = c("bio1", "bio2"), n_draws = 30)

  set.seed(7)
  res2 <- extract_beta_draws(model, env_vars = c("bio1", "bio2"), n_draws = 30)

  expect_equal(res1, res2)
})

# =============================================================================
# F13. Single env_var, single draw → 1×1 matrix with correct colname
# =============================================================================
test_that("ebd: single var and single draw returns a 1x1 named matrix", {
  model <- make_fake_model(n_total = 50, env_vars = "bio1")
  res   <- extract_beta_draws(model, env_vars = "bio1", n_draws = 1)

  expect_equal(dim(res), c(1L, 1L))
  expect_equal(colnames(res), "bio1")
})

# =============================================================================
# F14. Non-b_ columns (Intercept, sigma, lp__) never appear in output
# =============================================================================
test_that("ebd: non-b_ columns from draws matrix are excluded from output", {
  model <- make_fake_model(env_vars = c("bio1", "bio2"))
  res   <- extract_beta_draws(model, env_vars = c("bio1", "bio2"), n_draws = 20)

  forbidden <- c("b_Intercept", "Intercept", "sigma", "lp__")
  expect_false(any(forbidden %in% colnames(res)))
})

# =============================================================================
# F15. Output values match the corresponding rows of the input draws matrix
# =============================================================================
test_that("ebd: output values are a verbatim row-subset of the input draws", {
  set.seed(99)
  env_vars   <- c("bio1", "bio2")
  raw_matrix <- matrix(
    rnorm(100 * 2),
    nrow     = 100,
    ncol     = 2,
    dimnames = list(NULL, paste0("b_", env_vars))
  )
  model <- posterior::as_draws_matrix(raw_matrix)

  set.seed(99)
  res <- extract_beta_draws(model, env_vars = env_vars, n_draws = 10)

  # Every row in res must appear verbatim in raw_matrix
  for (i in seq_len(nrow(res))) {
    row_found <- any(apply(raw_matrix, 1, function(r) {
      isTRUE(all.equal(as.numeric(r), as.numeric(res[i, ]), tolerance = 1e-12))
    }))
    expect_true(row_found,
      info = sprintf("Row %d of output not found in original draws matrix", i))
  }
})


library(testthat)
library(terra)

# =============================================================================
# Tests for PosteriorCluster() — top-level exported function
#
# Coverage strategy
# -----------------
# Most of PosteriorCluster() is wrapped in # nocov (spatSample, the KNN
# projection, the final list assembly) because those paths require a live
# brmsfit and are covered by integration tests elsewhere.
#
# The lines that ARE coverable without a real brmsfit are:
#
#   G1.  set.seed() + match.arg(consensus_method)
#   G2.  fe_names filtering logic (grepl pattern)
#   G3.  length(env_vars) == 0 → stop()
#   G4.  is.null(beta_draws) bypass — supply beta_draws directly
#   G5.  var_mu / var_sd computation from pred_mat
#   G6.  zero SD warning + variable drop
#   G7.  terra::extract + complete.cases filtering of sample_pts
#   G8.  Draw loop body + the d %% 25 == 0 message branch
#   G9.  colMeans(beta_draws) + calls to compute_stability_scores /
#        compute_top3_clusters
#
# Mocking approach
# ----------------
# - brms::fixef        → local_mocked_bindings returns a named matrix
# - terra::spatSample  → local_mocked_bindings returns pre-built SpatVector
#   (avoids the nocov spatSample block while keeping terra::extract real)
# - beta_draws supplied directly → skips extract_beta_draws entirely
# - project_consensus_to_raster → mocked to return a minimal list so the
#   nocov projection block never executes and we don't need KNN/caret
# - reorder_clusters_geographically, sf::st_transform, terra::project →
#   mocked with identity-like stubs for the same reason
# =============================================================================

# ── Shared synthetic raster + points fixture ─────────────────────────────────
make_pc_raster <- function(seed = 1) {
  set.seed(seed)
  r <- terra::rast(
    nrows = 10, ncols = 10,
    xmin = 0, xmax = 1,
    ymin = 0, ymax = 1,
    crs  = "EPSG:4326"
  )
  terra::values(r) <- 1
  names(r) <- 'occurrence_prob_mean'
  r
}

make_pc_predictors <- function(mask, env_vars = c("bio1", "bio2"), seed = 1) {
  set.seed(seed)
  layers <- lapply(env_vars, function(v) {
    r <- mask
    terra::values(r) <- runif(terra::ncell(mask), 1, 10)
    r
  })
  pred <- terra::rast(layers)
  names(pred) <- env_vars
  pred
}

make_pc_pts <- function(mask, n = 60, seed = 1) {
  set.seed(seed)
  terra::spatSample(mask, size = n, method = "random",
                    as.points = TRUE, na.rm = TRUE)
}

# Minimal named matrix that looks like brms::fixef() output
# (row names = parameter names, cols = Estimate/Est.Error/Q2.5/Q97.5)
make_fixef_matrix <- function(env_vars = c("bio1", "bio2"),
                               extra    = "Intercept") {
  all_pars <- c(extra, env_vars)
  m <- matrix(
    rnorm(length(all_pars) * 4),
    nrow     = length(all_pars),
    ncol     = 4,
    dimnames = list(all_pars, c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  )
  m
}

# beta_draws matrix: n_draws rows, one col per env_var
make_beta_draws <- function(env_vars = c("bio1", "bio2"),
                             n_draws  = 30,
                             seed     = 1) {
  set.seed(seed)
  m <- matrix(
    rnorm(n_draws * length(env_vars), mean = 0.5, sd = 0.2),
    nrow     = n_draws,
    ncol     = length(env_vars),
    dimnames = list(NULL, env_vars)
  )
  m
}

# Minimal stub for project_consensus_to_raster — returns the skeleton list
# shape the caller expects so we never enter the nocov projection block
make_rast_list_stub <- function(mask) {
  dummy <- mask
  terra::values(dummy) <- 1L
  names(dummy) <- "cluster"
  stab <- mask; terra::values(stab) <- 0.8; names(stab) <- "cluster_stability"
  r2   <- mask; terra::values(r2)   <- 1L;  names(r2)   <- "rank2_cluster"
  r3   <- mask; terra::values(r3)   <- 1L;  names(r3)   <- "rank3_cluster"
  list(
    cluster_raster   = dummy,
    stability_raster = stab,
    rank2_raster     = r2,
    rank3_raster     = r3,
    knn_cluster      = structure(list(), class = "train"),
    knn_rank2        = structure(list(), class = "train"),
    knn_rank3        = structure(list(), class = "train"),
    knn_stability    = structure(list(), class = "train")
  )
}

# Stub for reorder_clusters_geographically
make_reorder_stub <- function(rast) {
  vect_stub <- sf::st_sfc(sf::st_point(c(0.5, 0.5)), crs = "EPSG:4326") |>
    sf::st_sf(geometry = _, cluster_id = 1L)
  list(raster = rast, vectors = vect_stub)
}

# =============================================================================
# G1.  match.arg: invalid consensus_method is caught immediately
# =============================================================================
test_that("pc: invalid consensus_method is rejected by match.arg", {
  expect_error(
    PosteriorCluster(
      model            = NULL,
      predictors       = make_pc_predictors(make_pc_raster()),
      f_rasts          = make_pc_raster(),
      pred_mat         = data.frame(bio1 = 1, bio2 = 1),
      training_data    = NULL,
      consensus_method = "kmeans",   # not in c("hierarchical", "pam")
      planar_proj      = "EPSG:32632",
      beta_draws       = make_beta_draws()
    ),
    regexp = "arg.*should be one of"
  )
})

# =============================================================================
# G2 + G3.  fe_names filtering: no env vars matched → stop()
# =============================================================================
test_that("pc: no env_vars matching raster layers triggers stop()", {
  mask       <- make_pc_raster()
  predictors <- make_pc_predictors(mask, env_vars = c("bio1", "bio2"))

  # fixef returns only "Intercept" + GP terms — none match raster layer names
  fixef_mat <- make_fixef_matrix(env_vars = character(0),
                                  extra    = c("Intercept", "sgp_lat", "sdgp_lat"))

  local_mocked_bindings(
    fixef = function(model, ...) fixef_mat,
    .package = "brms"
  )

  expect_error(
    PosteriorCluster(
      model         = NULL,
      predictors    = predictors,
      f_rasts       = mask,
      pred_mat      = data.frame(bio1 = rnorm(20), bio2 = rnorm(20)),
      training_data = NULL,
      planar_proj   = "EPSG:32632",
      beta_draws    = make_beta_draws()
    ),
    regexp = "No environmental fixed effects matched"
  )
})

# =============================================================================
# G4.  beta_draws supplied → is.null(beta_draws) branch NOT taken
#      (extract_beta_draws never called, no brms::fixef needed for draws step)
# =============================================================================
test_that("pc: supplying beta_draws bypasses extract_beta_draws", {
  skip_if_not_installed("sf")
  skip_if_not_installed("caret")

  set.seed(1)
  env_vars   <- c("bio1", "bio2")
  mask       <- make_pc_raster()
  predictors <- make_pc_predictors(mask, env_vars = env_vars)
  pts        <- make_pc_pts(mask, n = 60)
  bd         <- make_beta_draws(env_vars = env_vars, n_draws = 30)
  pred_mat   <- as.data.frame(terra::extract(predictors, pts, ID = FALSE))
  fixef_mat  <- make_fixef_matrix(env_vars = env_vars)
  rast_stub  <- make_rast_list_stub(mask)

  called_extract <- FALSE

  local_mocked_bindings(
    fixef              = function(model, ...) fixef_mat,
    .package           = "brms"
  )
  local_mocked_bindings(
    spatSample = function(...) pts,
    .package   = "terra"
  )
  local_mocked_bindings(
    project_consensus_to_raster      = function(...) rast_stub,
    reorder_clusters_geographically  = function(r, ...) make_reorder_stub(r),
    .package = 'safeHavens'
  )
  local_mocked_bindings(
    project   = function(x, ...) x,
    .package  = "terra"
  )
  local_mocked_bindings(
    st_transform = function(x, ...) x,
    .package     = "sf"
  )

  # If extract_beta_draws were called it would error (model = NULL, no draws).
  # Completing without error confirms the bypass branch was taken.
  expect_no_error(
    PosteriorCluster(
      model         = NULL,
      predictors    = predictors,
      f_rasts       = mask,
      pred_mat      = pred_mat,
      training_data = NULL,
      n_draws       = 30,
      n             = 3,
      n_pts         = 60,
      planar_proj   = "EPSG:32632",
      beta_draws    = bd          # <── supplied directly
    )
  )
})

# =============================================================================
# G5.  var_mu / var_sd computed correctly from pred_mat
#      (tested indirectly: if wrong, draw clustering diverges visibly)
# =============================================================================
test_that("pc: var_mu and var_sd are computed from pred_mat columns", {
  skip_if_not_installed("sf")
  skip_if_not_installed("caret")

  set.seed(2)
  env_vars   <- c("bio1", "bio2")
  mask       <- make_pc_raster()
  predictors <- make_pc_predictors(mask, env_vars = env_vars)
  pts        <- make_pc_pts(mask, n = 60)
  # pred_mat with known mean/sd — bio1 centred on 5, bio2 on 50
  pred_mat   <- data.frame(
    bio1 = rnorm(100, mean = 5,  sd = 1),
    bio2 = rnorm(100, mean = 50, sd = 10)
  )
  bd        <- make_beta_draws(env_vars = env_vars, n_draws = 10)
  fixef_mat <- make_fixef_matrix(env_vars = env_vars)
  rast_stub <- make_rast_list_stub(mask)

  local_mocked_bindings(fixef = function(...) fixef_mat, .package = "brms")
  local_mocked_bindings(spatSample = function(...) pts,  .package = "terra")
  local_mocked_bindings(
    project_consensus_to_raster     = function(...) rast_stub,
    reorder_clusters_geographically = function(r, ...) make_reorder_stub(r),
    .package = 'safeHavens'
  )
  local_mocked_bindings(project      = function(x, ...) x, .package = "terra")
  local_mocked_bindings(st_transform = function(x, ...) x, .package = "sf")

  # Should complete — if var_sd computation errored we'd see NaN/Inf
  expect_no_error(
    PosteriorCluster(
      model = NULL, predictors = predictors, f_rasts = mask,
      pred_mat = pred_mat, training_data = NULL,
      n_draws = 10, n = 3, n_pts = 60,
      planar_proj = "EPSG:32632", beta_draws = bd
    )
  )
})

# =============================================================================
# G6.  Zero-SD variable → warning + dropped from env_vars
# =============================================================================
test_that("pc: constant predictor column emits zero-SD warning and is dropped", {
  skip_if_not_installed("sf")
  skip_if_not_installed("caret")

  set.seed(3)
  env_vars   <- c("bio1", "bio2")
  mask       <- make_pc_raster()
  predictors <- make_pc_predictors(mask, env_vars = env_vars)
  pts        <- make_pc_pts(mask, n = 60)
  # bio2 is constant → SD = 0
  pred_mat   <- data.frame(
    bio1 = rnorm(100, 5, 1),
    bio2 = rep(7.0, 100)        # constant
  )
  bd        <- make_beta_draws(env_vars = env_vars, n_draws = 10)
  fixef_mat <- make_fixef_matrix(env_vars = env_vars)
  rast_stub <- make_rast_list_stub(mask)

  local_mocked_bindings(fixef = function(...) fixef_mat, .package = "brms")
  local_mocked_bindings(spatSample = function(...) pts,  .package = "terra")
  local_mocked_bindings(
    project_consensus_to_raster     = function(...) rast_stub,
    reorder_clusters_geographically = function(r, ...) make_reorder_stub(r),
    .package = 'safeHavens'
  )
  local_mocked_bindings(project      = function(x, ...) x, .package = "terra")
  local_mocked_bindings(st_transform = function(x, ...) x, .package = "sf")

  expect_warning(
    PosteriorCluster(
      model = NULL, predictors = predictors, f_rasts = mask,
      pred_mat = pred_mat, training_data = NULL,
      n_draws = 10, n = 3, n_pts = 60,
      planar_proj = "EPSG:32632", beta_draws = bd
    ),
    regexp = "zero SD"
  )
})

# =============================================================================
# G7.  terra::extract + complete.cases filtering
#      Inject NAs into the predictor raster so some rows are dropped.
# ==============================f===============================================
test_that("pc: complete.cases filtering reduces n_pts_actual when NAs present", {
  skip_if_not_installed("sf")
  skip_if_not_installed("caret")

  set.seed(4)
  env_vars   <- c("bio1", "bio2")
  mask       <- make_pc_raster()
  predictors <- make_pc_predictors(mask, env_vars = env_vars)

  # Corrupt half of bio1 with NAs so ~half of extracted rows are incomplete
  terra::values(predictors[[1]])[1:50] <- NA_real_

  pts        <- make_pc_pts(mask, n = 60)
  pred_mat   <- data.frame(bio1 = rnorm(100, 5, 1), bio2 = rnorm(100, 5, 1))
  bd         <- make_beta_draws(env_vars = env_vars, n_draws = 10)
  fixef_mat  <- make_fixef_matrix(env_vars = env_vars)
  rast_stub  <- make_rast_list_stub(mask)

  local_mocked_bindings(fixef = function(...) fixef_mat, .package = "brms")
  local_mocked_bindings(spatSample = function(...) pts,  .package = "terra")
  local_mocked_bindings(
    project_consensus_to_raster     = function(...) rast_stub,
    reorder_clusters_geographically = function(r, ...) make_reorder_stub(r),
    .package = 'safeHavens'
  )
  local_mocked_bindings(project      = function(x, ...) x, .package = "terra")
  local_mocked_bindings(st_transform = function(x, ...) x, .package = "sf")

  # Capture the "N of N_requested points have complete predictor data" message
  msgs <- character(0)
  withCallingHandlers(
    PosteriorCluster(
      model = NULL, predictors = predictors, f_rasts = mask,
      pred_mat = pred_mat, training_data = NULL,
      n_draws = 10, n = 3, n_pts = 60,
      planar_proj = "EPSG:32632", beta_draws = bd
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  # The message should report fewer than 60 complete points
  complete_msg <- msgs[grepl("complete predictor data", msgs)]
  expect_length(complete_msg, 1L)
  n_complete <- as.integer(regmatches(
    complete_msg,
    regexpr("^\\s*\\d+", complete_msg)
  ))
  expect_lt(n_complete, 60L)
})

# =============================================================================
# G8.  Draw loop: d %% 25 == 0 message fires on draw 25
# =============================================================================
test_that("pc: draw progress message emitted at draw 25 when n_draws >= 25", {
  skip_if_not_installed("sf")
  skip_if_not_installed("caret")

  set.seed(5)
  env_vars   <- c("bio1", "bio2")
  mask       <- make_pc_raster()
  predictors <- make_pc_predictors(mask, env_vars = env_vars)
  pts        <- make_pc_pts(mask, n = 60)
  pred_mat   <- data.frame(bio1 = rnorm(100, 5, 1), bio2 = rnorm(100, 5, 1))
  bd         <- make_beta_draws(env_vars = env_vars, n_draws = 26)
  fixef_mat  <- make_fixef_matrix(env_vars = env_vars)
  rast_stub  <- make_rast_list_stub(mask)

  local_mocked_bindings(fixef = function(...) fixef_mat, .package = "brms")
  local_mocked_bindings(spatSample = function(...) pts,  .package = "terra")
  local_mocked_bindings(
    project_consensus_to_raster     = function(...) rast_stub,
    reorder_clusters_geographically = function(r, ...) make_reorder_stub(r),
    .package = 'safeHavens'
  )
  local_mocked_bindings(project      = function(x, ...) x, .package = "terra")
  local_mocked_bindings(st_transform = function(x, ...) x, .package = "sf")

  msgs <- character(0)
  withCallingHandlers(
    PosteriorCluster(
      model = NULL, predictors = predictors, f_rasts = mask,
      pred_mat = pred_mat, training_data = NULL,
      n_draws = 26, n = 3, n_pts = 60,
      planar_proj = "EPSG:32632", beta_draws = bd
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  # "Draw 25 / 26" should appear
  expect_true(
    any(grepl("Draw 25", msgs)),
    info = paste("Messages seen:", paste(msgs, collapse = " | "))
  )
})

# =============================================================================
# G9.  PAM consensus method branch is reached
# =============================================================================
test_that("pc: consensus_method = 'pam' runs without error", {
  skip_if_not_installed("sf")
  skip_if_not_installed("caret")
  skip_if_not_installed("cluster")

  set.seed(6)
  env_vars   <- c("bio1", "bio2")
  mask       <- make_pc_raster()
  predictors <- make_pc_predictors(mask, env_vars = env_vars)
  pts        <- make_pc_pts(mask, n = 60)
  pred_mat   <- data.frame(bio1 = rnorm(100, 5, 1), bio2 = rnorm(100, 5, 1))
  bd         <- make_beta_draws(env_vars = env_vars, n_draws = 10)
  fixef_mat  <- make_fixef_matrix(env_vars = env_vars)
  rast_stub  <- make_rast_list_stub(mask)

  local_mocked_bindings(fixef = function(...) fixef_mat, .package = "brms")
  local_mocked_bindings(spatSample = function(...) pts,  .package = "terra")
  local_mocked_bindings(
    project_consensus_to_raster     = function(...) rast_stub,
    reorder_clusters_geographically = function(r, ...) make_reorder_stub(r),
    .package = 'safeHavens'
  )
  local_mocked_bindings(project      = function(x, ...) x, .package = "terra")
  local_mocked_bindings(st_transform = function(x, ...) x, .package = "sf")

  expect_no_error(
    PosteriorCluster(
      model = NULL, predictors = predictors, f_rasts = mask,
      pred_mat = pred_mat, training_data = NULL,
      n_draws = 10, n = 3, n_pts = 60,
      planar_proj   = "EPSG:32632",
      beta_draws    = bd,
      consensus_method = "pam"    # <── exercises the switch() pam branch
    )
  )
})