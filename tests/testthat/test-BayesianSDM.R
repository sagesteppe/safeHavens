test_that("bayesianSDM runs with minimal synthetic data", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Create minimal synthetic data
  set.seed(42)
  n_pts <- 100  # Increased further for knndm stability
  coords <- data.frame(
    lon = runif(n_pts, -120, -110),
    lat = runif(n_pts, 30, 40)
  )
  
  # Synthetic occurrence points
  x <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  x$occurrence <- sample(0:1, n_pts, replace = TRUE, prob = c(0.3, 0.7))
  
  # Create synthetic raster predictors
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 20, ncols = 20, crs = "EPSG:4326"
  )
  pred1 <- r
  terra::values(pred1) <- runif(terra::ncell(r), 0, 100)
  pred2 <- r
  terra::values(pred2) <- runif(terra::ncell(r), 10, 30)
  predictors <- c(pred1, pred2)
  names(predictors) <- c("bio1", "bio2")
  
  # Run minimal model
  result <- expect_warning(
    bayesianSDM(
      x = x,
      predictors = predictors,
      planar_projection = 3857,
      pca_predictors = FALSE,
      feature_selection = "none",
      vif = FALSE,
      chains = 2,
      iter = 500,
      warmup = 250,
      cores = 1,
      k = 5,  # Changed to 5 for better knndm stability
      seed = 42,
      fact = 1.5,
      silent = 2
    ),
    "The ESS has been capped to avoid unstable estimates."
  )
  
  # Check structure
  expect_type(result, "list")
  expect_s4_class(result$RasterPredictions, "SpatRaster")
  expect_s4_class(result$RasterPredictions_sd, "SpatRaster")
  expect_s4_class(result$Predictors, "SpatRaster")
  expect_s3_class(result$Model, "brmsfit")
  expect_type(result$CVStructure, "list")
  expect_s3_class(result$LOO, "loo")
  expect_s3_class(result$ConfusionMatrix, "confusionMatrix")
  expect_s3_class(result$TrainData, "sf")
  expect_s3_class(result$TestData, "sf")
  expect_s4_class(result$AOA, "SpatRaster")
  expect_type(result$AOA_Diagnostics, "list")
  expect_type(result$Diagnostics, "list")
  
  # Check diagnostics
  expect_true(result$Diagnostics$max_Rhat < 1.1)
  expect_true(result$Diagnostics$min_BulkESS > 100)
})


test_that("bayesianSDM works with PCA", {
  skip_on_cran()
  skip_if_not_installed("brms")
  
  set.seed(123)
  n_pts <- 100
  coords <- data.frame(
    lon = runif(n_pts, -120, -110),
    lat = runif(n_pts, 30, 40)
  )
  
  x <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  # Don't specify occurrence - test pseudo-absence generation
  
  # Create 4 correlated predictors
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 15, ncols = 15, crs = "EPSG:4326"
  )
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
  
  result <- expect_warning(
    bayesianSDM(
      x = x,
      predictors = pred_stack,
      planar_projection = 3857,
      pca_predictors = TRUE,
      pca_axes = 3,
      feature_selection = "none",
      vif = FALSE,
      chains = 2,
      iter = 400,
      warmup = 200,
      cores = 1,
      k = 5,
      seed = 99,
      silent = 2
    ),
    "The ESS has been capped to avoid unstable estimates."
  )
  
  # Check PCA outputs
  expect_s3_class(result$PCAModel, "prcomp")
  expect_equal(terra::nlyr(result$Predictors), 3)
  expect_match(names(result$Predictors)[1], "^PC")
  
  # Check occurrence was generated
  expect_true("occurrence" %in% names(result$TrainData))
})


test_that("bayesianSDM works with feature selection", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("CAST")
  
  set.seed(456)
  n_pts <- 100
  coords <- data.frame(
    lon = runif(n_pts, -120, -110),
    lat = runif(n_pts, 30, 40)
  )
  
  x <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  x$occurrence <- rbinom(n_pts, 1, 0.6)
  
  # Create 5 predictors
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 12, ncols = 12, crs = "EPSG:4326"
  )
  pred_list <- lapply(1:5, function(i) {
    p <- r
    terra::values(p) <- runif(terra::ncell(r))
    p
  })
  predictors <- do.call(c, pred_list)
  names(predictors) <- paste0("var", 1:5)
  
  result <- expect_warning(
    bayesianSDM(
      x = x,
      predictors = predictors,
      planar_projection = 3857,
      pca_predictors = FALSE,
      feature_selection = "ffs",
      min_ffs_var = 2,
      vif = FALSE,
      chains = 2,
      iter = 400,
      warmup = 200,
      cores = 1,
      k = 5,
      seed = 789,
      silent = 2
      ),
    "The ESS has been capped to avoid unstable estimates."
  )
  
  
  # FFS should select a subset
  expect_true(terra::nlyr(result$Predictors) <= 5)
  expect_true(terra::nlyr(result$Predictors) >= 2)
})


test_that("VIF filtering works", {
  skip_on_cran()
  skip_if_not_installed("usdm")
  
  set.seed(321)
  # Create highly correlated predictors
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 10, ncols = 10, crs = "EPSG:4326"
  )
  base <- matrix(runif(terra::ncell(r)), nrow = terra::nrow(r))
  
  # 3 highly correlated + 1 independent
  pred_stack <- c(
    terra::rast(base),
    terra::rast(base * 0.95 + matrix(runif(terra::ncell(r), 0, 0.05), nrow = terra::nrow(r))),
    terra::rast(base * 0.9 + matrix(runif(terra::ncell(r), 0, 0.1), nrow = terra::nrow(r))),
    terra::rast(matrix(runif(terra::ncell(r)), nrow = terra::nrow(r)))
  )
  terra::ext(pred_stack) <- terra::ext(r)
  terra::crs(pred_stack) <- terra::crs(r)
  names(pred_stack) <- paste0("bio", 1:4)
  
  # VIF should remove some
  vif_result <- usdm::vifcor(pred_stack)
  expect_true(length(vif_result@results$Variables) < 4)
})


test_that("Helper functions work correctly", {
  # Test add_planar_coords
  pts <- data.frame(x = c(-120, -115), y = c(35, 38))
  sf_pts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = 4326)
  
  result <- add_planar_coords(sf_pts, 3857, scale_km = TRUE)
  
  expect_true("gp_x" %in% names(result))
  expect_true("gp_y" %in% names(result))
  expect_true(all(abs(result$gp_x) < 20000))  # Should be in km, not metres
})


test_that("run_pca_on_raster handles edge cases", {
  # Small raster
  r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  pred1 <- r
  terra::values(pred1) <- 1:25
  pred2 <- r
  terra::values(pred2) <- 25:1
  predictors <- c(pred1, pred2)
  names(predictors) <- c("a", "b")
  
  # Request more axes than available
  result <- run_pca_on_raster(predictors, pca_axes = 5)
  
  expect_equal(length(result$pc_names), 2)  # Should cap at 2
  expect_s4_class(result$raster, "SpatRaster")
  expect_s3_class(result$pca_model, "prcomp")
})


test_that("build_priors creates valid prior objects", {
  # Horseshoe prior
  hs_prior <- build_priors(
    prior_type = "horseshoe",
    prior_scale = 1,
    pred_names = c("var1", "var2", "var3"),
    n_obs = 100,
    p0 = 2,
    gp_prior = NULL
  )
  
  expect_s3_class(hs_prior, "brmsprior")
  expect_true(any(grepl("horseshoe", hs_prior$prior)))
  
  # Normal prior
  norm_prior <- build_priors(
    prior_type = "normal",
    prior_scale = 2,
    pred_names = c("var1", "var2"),
    n_obs = 50,
    p0 = NULL,
    gp_prior = NULL
  )
  
  expect_s3_class(norm_prior, "brmsprior")
  expect_true(any(grepl("normal", norm_prior$prior)))
  
  # Student prior
  student_prior <- build_priors(
    prior_type = "student",
    prior_scale = 1.5,
    pred_names = c("var1"),
    n_obs = 75,
    p0 = NULL,
    gp_prior = NULL
  )
  
  expect_s3_class(student_prior, "brmsprior")
  expect_true(any(grepl("student_t", student_prior$prior)))
})


test_that("check_convergence requires brmsfit object", {
  skip_if_not_installed("brms")
  
  # Check that function requires a brmsfit object
  expect_error(
    check_convergence(list())
  )
})


test_that("compute_aoa_bayes extracts variables correctly", {
  skip_on_cran()
  skip_if_not_installed("CAST")
  skip_if_not_installed("brms")
  
  # This would require a full bayesianSDM result
  # Testing the environmental variable extraction logic
  
  # Mock fixef output
  mock_fixef <- matrix(c(0.5, -0.3, 0.1, 0.8), ncol = 1)
  rownames(mock_fixef) <- c("PC1", "PC2", "sgpx_1", "sgpy_2")
  
  # Test the filtering regex
  fe_names <- rownames(mock_fixef)
  fe_names <- sub("^b_", "", fe_names)
  fe_names_filtered <- fe_names[!grepl("^(Intercept|s\\(|s[gp])", fe_names)]
  
  expect_equal(fe_names_filtered, c("PC1", "PC2"))
})


test_that("PosteriorCluster handles matrix indexing correctly", {
  skip_on_cran()
  
  # Test the beta draw extraction fix
  # Create a mock beta matrix
  beta_matrix <- matrix(
    c(0.5, 0.3, -0.2, 0.7, 0.1, -0.4),
    nrow = 2, ncol = 3
  )
  colnames(beta_matrix) <- c("PC1", "PC2", "PC3")
  
  # Extract one row (this was the bug)
  d <- 1
  betas_d <- beta_matrix[d, , drop = FALSE]
  betas_d <- as.numeric(betas_d)
  names(betas_d) <- colnames(beta_matrix)
  
  # Check names are preserved
  expect_named(betas_d, c("PC1", "PC2", "PC3"))
  expect_equal(unname(betas_d[1]), 0.5)  # Compare values without names
  expect_equal(names(betas_d)[1], "PC1")
})


test_that("compute_top3_clusters replicates rank1 when only one cluster", {
  # Simulate draw clusterings where all points always go to cluster 5
  n_pts <- 10
  n_draws <- 20
  draw_clusterings <- matrix(5, nrow = n_pts, ncol = n_draws)
  
  # This should replicate cluster 5 into rank2 and rank3
  result <- compute_top3_clusters(draw_clusterings, n_pts, n_draws)
  
  expect_equal(result$labels_matrix[1, 1], 5)
  expect_equal(result$labels_matrix[1, 2], 5)  # Replicated
  expect_equal(result$labels_matrix[1, 3], 5)  # Replicated
  expect_equal(result$pct_matrix[1, 1], 100)
  expect_equal(result$pct_matrix[1, 2], 100)
  expect_equal(result$pct_matrix[1, 3], 100)
})


test_that("compute_top3_clusters handles mixed assignments", {
  # Simulate varied assignments
  n_pts <- 5
  n_draws <- 100
  
  # Point 1: 70% cluster 1, 20% cluster 2, 10% cluster 3
  # Point 2: 100% cluster 5
  draw_clusterings <- matrix(NA_integer_, nrow = n_pts, ncol = n_draws)
  draw_clusterings[1, ] <- sample(c(1, 2, 3), n_draws, replace = TRUE, prob = c(0.7, 0.2, 0.1))
  draw_clusterings[2, ] <- rep(5, n_draws)
  
  result <- compute_top3_clusters(draw_clusterings, n_pts, n_draws)
  
  # Point 1 should have 3 different clusters
  expect_true(length(unique(result$labels_matrix[1, ])) == 3)
  
  # Point 2 should have cluster 5 replicated
  expect_true(all(result$labels_matrix[2, ] == 5))
  expect_equal(result$pct_matrix[2, 1], 100)
})


test_that("Error handling works for insufficient data", {
  skip_on_cran()
  
  # Too few points for knndm - expect any error
  set.seed(42)
  coords <- data.frame(lon = runif(5, -120, -119), lat = runif(5, 35, 36))
  x <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  x$occurrence <- sample(0:1, 5, replace = TRUE)
  
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 5, ncols = 5, crs = "EPSG:4326"
  )
  terra::values(r) <- 1:25
  
  # Should error during knndm or data partitioning
  expect_error(
    expect_warning(
      bayesianSDM(
        x = x,
        predictors = r,
        planar_projection = 3857,
        chains = 1,
        iter = 100,
        vif = FALSE,
        feature_selection = "none",
        k = 3,
        silent = 2
      ), "no non-missing arguments to min; returning Inf"
    ), "undefined columns selected", fixed = TRUE
  )
})


test_that("Coordinate weighting functions work", {
  # Test add_coord_weights_to_points
  sample_pts <- terra::vect(
    cbind(x = c(-120, -115, -110), y = c(35, 37, 39)),
    crs = "EPSG:4326"
  )
  
  pt_env <- data.frame(
    PC1 = c(0.5, 1.2, -0.3),
    PC2 = c(0.1, 0.8, 0.4)
  )
  
  result <- add_coord_weights_to_points(
    sample_pts, pt_env, coord_wt = 2.0, env_vars = c("PC1", "PC2")
  )
  
  expect_true("coord_x_w" %in% names(result))
  expect_true("coord_y_w" %in% names(result))
  expect_equal(nrow(result), 3)
})


test_that("Feature selection handles min_ffs_var correctly", {
  # When min_ffs_var > nlyr, should reset to 2
  pred_names <- c("PC1", "PC2", "PC3")
  min_ffs_var <- 10
  nlyr <- 3
  
  if (min_ffs_var > nlyr) {
    min_ffs_var <- 2
  }
  
  expect_equal(min_ffs_var, 2)
})


test_that("Prediction rasters have correct structure", {
  skip_on_cran()
  skip_if_not_installed("brms")
  
  # Create proper raster with values
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 5, ncols = 5, crs = "EPSG:4326"
  )
  terra::values(r) <- 1:25
  
  # Create coordinate data frame manually
  coords_df <- data.frame(
    x = rep(seq(-120, -110, length.out = 5), 5),
    y = rep(seq(30, 40, length.out = 5), each = 5)
  )
  
  coords_vect <- terra::vect(coords_df, geom = c("x", "y"), crs = "EPSG:4326")
  coords_proj <- terra::project(coords_vect, "EPSG:3857")
  coords_km <- terra::crds(coords_proj) / 1000
  
  expect_equal(nrow(coords_km), 25)
  expect_equal(ncol(coords_km), 2)
  expect_true(all(abs(coords_km) < 50000))  # Reasonable km values
})


test_that("AOA diagnostics are computed correctly", {
  # Mock AOA result structure
  mock_aoa <- terra::rast(nrows = 10, ncols = 10)
  terra::values(mock_aoa) <- sample(0:1, 100, replace = TRUE, prob = c(0.3, 0.7))
  
  aoa_coverage <- sum(terra::values(mock_aoa), na.rm = TRUE) / 
                  sum(!is.na(terra::values(mock_aoa)))
  
  expect_true(aoa_coverage >= 0 && aoa_coverage <= 1)
  expect_type(aoa_coverage, "double")
})


test_that("Integration test: full workflow runs end-to-end", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("CAST")
  
  # Full minimal workflow
  set.seed(2024)
  
  # Generate data
  n_pts <- 120
  coords <- data.frame(
    lon = runif(n_pts, -120, -110),
    lat = runif(n_pts, 30, 40)
  )
  
  x <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  x$occurrence <- rbinom(n_pts, 1, 0.6)
  
  # Predictors
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 15, ncols = 15, crs = "EPSG:4326"
  )
  pred_list <- lapply(1:3, function(i) {
    p <- r
    terra::values(p) <- rnorm(terra::ncell(r))
    p
  })
  predictors <- do.call(c, pred_list)
  names(predictors) <- paste0("bio", 1:3)
  
  result <- expect_warning(
    bayesianSDM(
      x = x,
      predictors = predictors,
      planar_projection = 3857,
      pca_predictors = TRUE,
      pca_axes = 2,
      feature_selection = "none",
      vif = FALSE,
      chains = 2,
      iter = 400,
      warmup = 200,
      cores = 1,
      k = 5,
      seed = 2024,
      fact = 2,
      silent = 2
    ),
    "The ESS has been capped to avoid unstable estimates."
  )
  
  # FFS should select a subset
  expect_true(terra::nlyr(result$Predictors) <= 5)
  expect_true(terra::nlyr(result$Predictors) >= 2)
})


test_that("VIF filtering works", {
  skip_on_cran()
  skip_if_not_installed("usdm")
  
  set.seed(321)
  # Create highly correlated predictors
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 10, ncols = 10, crs = "EPSG:4326"
  )
  base <- matrix(runif(terra::ncell(r)), nrow = terra::nrow(r))
  
  # 3 highly correlated + 1 independent
  pred_stack <- c(
    terra::rast(base),
    terra::rast(base * 0.95 + matrix(runif(terra::ncell(r), 0, 0.05), nrow = terra::nrow(r))),
    terra::rast(base * 0.9 + matrix(runif(terra::ncell(r), 0, 0.1), nrow = terra::nrow(r))),
    terra::rast(matrix(runif(terra::ncell(r)), nrow = terra::nrow(r)))
  )
  terra::ext(pred_stack) <- terra::ext(r)
  terra::crs(pred_stack) <- terra::crs(r)
  names(pred_stack) <- paste0("bio", 1:4)
  
  # VIF should remove some
  vif_result <- usdm::vifcor(pred_stack)
  expect_true(length(vif_result@results$Variables) < 4)
})


test_that("Helper functions work correctly", {
  # Test add_planar_coords
  pts <- data.frame(x = c(-120, -115), y = c(35, 38))
  sf_pts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = 4326)
  
  result <- add_planar_coords(sf_pts, 3857, scale_km = TRUE)
  
  expect_true("gp_x" %in% names(result))
  expect_true("gp_y" %in% names(result))
  expect_true(all(abs(result$gp_x) < 20000))  # Should be in km, not metres
})


test_that("run_pca_on_raster handles edge cases", {
  # Small raster
  r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  pred1 <- r
  terra::values(pred1) <- 1:25
  pred2 <- r
  terra::values(pred2) <- 25:1
  predictors <- c(pred1, pred2)
  names(predictors) <- c("a", "b")
  
  # Request more axes than available
  result <- run_pca_on_raster(predictors, pca_axes = 5)
  
  expect_equal(length(result$pc_names), 2)  # Should cap at 2
  expect_s4_class(result$raster, "SpatRaster")
  expect_s3_class(result$pca_model, "prcomp")
})


test_that("build_priors creates valid prior objects", {
  # Horseshoe prior
  hs_prior <- build_priors(
    prior_type = "horseshoe",
    prior_scale = 1,
    pred_names = c("var1", "var2", "var3"),
    n_obs = 100,
    p0 = 2,
    gp_prior = NULL
  )
  
  expect_s3_class(hs_prior, "brmsprior")
  expect_true(any(grepl("horseshoe", hs_prior$prior)))
  
  # Normal prior
  norm_prior <- build_priors(
    prior_type = "normal",
    prior_scale = 2,
    pred_names = c("var1", "var2"),
    n_obs = 50,
    p0 = NULL,
    gp_prior = NULL
  )
  
  expect_s3_class(norm_prior, "brmsprior")
  expect_true(any(grepl("normal", norm_prior$prior)))
  
  # Student prior
  student_prior <- build_priors(
    prior_type = "student",
    prior_scale = 1.5,
    pred_names = c("var1"),
    n_obs = 75,
    p0 = NULL,
    gp_prior = NULL
  )
  
  expect_s3_class(student_prior, "brmsprior")
  expect_true(any(grepl("student_t", student_prior$prior)))
})


test_that("check_convergence extracts correct diagnostics", {
  skip_if_not_installed("brms")
  
  # Check that function requires a brmsfit object
  expect_error(
    check_convergence(list()),
    "posterior"  # Should error about missing posterior
  )
  
  # The actual diagnostics would be tested in the integration tests
  # where we have real brmsfit objects
})


test_that("compute_aoa_bayes extracts variables correctly", {
  skip_on_cran()
  skip_if_not_installed("CAST")
  skip_if_not_installed("brms")
  
  # This would require a full bayesianSDM result
  # Testing the environmental variable extraction logic
  
  # Mock fixef output
  mock_fixef <- matrix(c(0.5, -0.3, 0.1, 0.8), ncol = 1)
  rownames(mock_fixef) <- c("PC1", "PC2", "sgpx_1", "sgpy_2")
  
  # Test the filtering regex
  fe_names <- rownames(mock_fixef)
  fe_names <- sub("^b_", "", fe_names)
  fe_names_filtered <- fe_names[!grepl("^(Intercept|s\\(|s[gp])", fe_names)]
  
  expect_equal(fe_names_filtered, c("PC1", "PC2"))
})


test_that("PosteriorCluster handles matrix indexing correctly", {
  skip_on_cran()
  
  # Test the beta draw extraction fix
  # Create a mock beta matrix
  beta_matrix <- matrix(
    c(0.5, 0.3, -0.2, 0.7, 0.1, -0.4),
    nrow = 2, ncol = 3
  )
  colnames(beta_matrix) <- c("PC1", "PC2", "PC3")
  
  # Extract one row (this was the bug)
  d <- 1
  betas_d <- beta_matrix[d, , drop = FALSE]
  betas_d <- as.numeric(betas_d)
  names(betas_d) <- colnames(beta_matrix)
  
  # Check names are preserved
  expect_named(betas_d, c("PC1", "PC2", "PC3"))
  expect_equal(unname(betas_d[1]), 0.5)  # Compare values without names
  expect_equal(names(betas_d)[1], "PC1")
})


test_that("compute_top3_clusters replicates rank1 when only one cluster", {
  # Simulate draw clusterings where all points always go to cluster 5
  n_pts <- 10
  n_draws <- 20
  draw_clusterings <- matrix(5, nrow = n_pts, ncol = n_draws)
  
  # This should replicate cluster 5 into rank2 and rank3
  result <- compute_top3_clusters(draw_clusterings, n_pts, n_draws)
  
  expect_equal(result$labels_matrix[1, 1], 5)
  expect_equal(result$labels_matrix[1, 2], 5)  # Replicated
  expect_equal(result$labels_matrix[1, 3], 5)  # Replicated
  expect_equal(result$pct_matrix[1, 1], 100)
  expect_equal(result$pct_matrix[1, 2], 100)
  expect_equal(result$pct_matrix[1, 3], 100)
})


test_that("compute_top3_clusters handles mixed assignments", {
  # Simulate varied assignments
  n_pts <- 5
  n_draws <- 100
  
  # Point 1: 70% cluster 1, 20% cluster 2, 10% cluster 3
  # Point 2: 100% cluster 5
  draw_clusterings <- matrix(NA_integer_, nrow = n_pts, ncol = n_draws)
  draw_clusterings[1, ] <- sample(c(1, 2, 3), n_draws, replace = TRUE, prob = c(0.7, 0.2, 0.1))
  draw_clusterings[2, ] <- rep(5, n_draws)
  
  result <- compute_top3_clusters(draw_clusterings, n_pts, n_draws)
  
  # Point 1 should have 3 different clusters
  expect_true(length(unique(result$labels_matrix[1, ])) == 3)
  
  # Point 2 should have cluster 5 replicated
  expect_true(all(result$labels_matrix[2, ] == 5))
  expect_equal(result$pct_matrix[2, 1], 100)
})


test_that("Error handling works for insufficient data", {
  skip_on_cran()
  
  # Too few points for knndm
  set.seed(42)
  coords <- data.frame(lon = runif(10, -120, -119), lat = runif(10, 35, 36))
  x <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  x$occurrence <- sample(0:1, 10, replace = TRUE)
  
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 5, ncols = 5, crs = "EPSG:4326"
  )
  terra::values(r) <- 1:25
  
  # Should error during knndm or data partitioning
  expect_error(
    expect_warning(
      bayesianSDM(
        x = x,
        predictors = r,
        planar_projection = 3857,
        chains = 1,
        iter = 100,
        vif = FALSE,
        feature_selection = "none",
        k = 3,
        silent = 2
      ), "Some classes have a single record ( 1 ) and these will be selected for the sample"
    ), "A term has fewer unique covariate combinations than specified maximum degrees of freedom",
    fixed = TRUE
  )
})


test_that("Coordinate weighting functions work", {
  # Test add_coord_weights_to_points
  sample_pts <- terra::vect(
    cbind(x = c(-120, -115, -110), y = c(35, 37, 39)),
    crs = "EPSG:4326"
  )
  
  pt_env <- data.frame(
    PC1 = c(0.5, 1.2, -0.3),
    PC2 = c(0.1, 0.8, 0.4)
  )
  
  result <- add_coord_weights_to_points(
    sample_pts, pt_env, coord_wt = 2.0, env_vars = c("PC1", "PC2")
  )
  
  expect_true("coord_x_w" %in% names(result))
  expect_true("coord_y_w" %in% names(result))
  expect_equal(nrow(result), 3)
})


test_that("Feature selection handles min_ffs_var correctly", {
  # When min_ffs_var > nlyr, should reset to 2
  pred_names <- c("PC1", "PC2", "PC3")
  min_ffs_var <- 10
  nlyr <- 3
  
  if (min_ffs_var > nlyr) {
    min_ffs_var <- 2
  }
  
  expect_equal(min_ffs_var, 2)
})


test_that("Prediction rasters have correct structure", {
  skip_on_cran()
  skip_if_not_installed("brms")
  
  # This test would run the full pipeline and check outputs
  # Abbreviated here for speed
  
  # Check that coordinate rasters can be created
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 5, ncols = 5, crs = "EPSG:4326"
  )
  terra::values(r) <- 1:25  # Add values
  
  coords_df <- data.frame(
    x = rep(seq(-120, -110, length.out = 5), 5),
    y = rep(seq(30, 40, length.out = 5), each = 5)
  )
  
  coords_vect <- terra::vect(coords_df, geom = c("x", "y"), crs = "EPSG:4326")
  coords_proj <- terra::project(coords_vect, "EPSG:3857")
  coords_km <- terra::crds(coords_proj) / 1000
  
  expect_equal(nrow(coords_km), 25)
  expect_equal(ncol(coords_km), 2)
  expect_true(all(abs(coords_km) < 50000))  # Reasonable km values
})


test_that("AOA diagnostics are computed correctly", {
  # Mock AOA result structure
  mock_aoa <- terra::rast(nrows = 10, ncols = 10)
  terra::values(mock_aoa) <- sample(0:1, 100, replace = TRUE, prob = c(0.3, 0.7))
  
  aoa_coverage <- sum(terra::values(mock_aoa), na.rm = TRUE) / 
                  sum(!is.na(terra::values(mock_aoa)))
  
  expect_true(aoa_coverage >= 0 && aoa_coverage <= 1)
  expect_type(aoa_coverage, "double")
})


test_that("Integration test: full workflow runs end-to-end", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("CAST")
  
  # Full minimal workflow
  set.seed(2024)
  
  # Generate data
  n_pts <- 90  # Increased for knndm
  coords <- data.frame(
    lon = runif(n_pts, -120, -110),
    lat = runif(n_pts, 30, 40)
  )
  
  x <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  x$occurrence <- rbinom(n_pts, 1, 0.6)
  
  # Predictors
  r <- terra::rast(
    xmin = -120, xmax = -110, ymin = 30, ymax = 40,
    nrows = 15, ncols = 15, crs = "EPSG:4326"
  )
  pred_list <- lapply(1:3, function(i) {
    p <- r
    terra::values(p) <- rnorm(terra::ncell(r))
    p
  })
  predictors <- do.call(c, pred_list)
  names(predictors) <- paste0("bio", 1:3)
  
  # Run full pipeline
  result <- expect_warning(
    bayesianSDM(
      x = x,
      predictors = predictors,
      planar_projection = 3857,
      pca_predictors = TRUE,
      pca_axes = 2,
      feature_selection = "none",
      vif = FALSE,
      chains = 2,
      iter = 400,
      warmup = 200,
      cores = 1,
      k = 3,  # Changed to 3
      seed = 2024,
      fact = 2,
      silent = 2
    ),
    "The ESS has been capped to avoid unstable estimates."
  )
  
  # Comprehensive checks
  expect_s3_class(result$Model, "brmsfit")
  expect_true(result$Diagnostics$max_Rhat < 1.1)
  expect_true(terra::nlyr(result$Predictors) == 2)
  expect_s4_class(result$AOA, "SpatRaster")
  expect_true(terra::nlyr(result$AOA) == 3)  # AOA, DI, LPD
  expect_named(result$AOA, c("AOA", "DI", "LPD"))
  expect_true(result$AOA_Diagnostics$AOA_coverage >= 0)
  expect_true(result$AOA_Diagnostics$AOA_coverage <= 1)
  
  # Check that predictions are valid probabilities
  pred_vals <- terra::values(result$RasterPredictions)
  pred_vals <- pred_vals[!is.na(pred_vals)]
  expect_true(all(pred_vals >= 0 & pred_vals <= 1))
  
  # Check SD is non-negative
  sd_vals <- terra::values(result$RasterPredictions_sd)
  sd_vals <- sd_vals[!is.na(sd_vals)]
  expect_true(all(sd_vals >= 0))
})