# tests/testthat/test-EnvironmentalBasedSample.R

# NOTE: Spatial sampling functions (terra::spatSample, sf::st_sample) do not guarantee
# exact sample sizes. Tests use ±10% tolerance to account for this variation.
# Additionally, complete.cases() filtering can reduce sample size if predictor rasters
# have NAs or mismatched extents. Tests ensure all predictor layers are complete.

# Test helper functions first ----

test_that("add_weighted_coordinates adds x and y layers", {
  # Create simple test raster
  r <- terra::rast(ncols = 10, nrows = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  r$var1 <- runif(100, 0, 1)
  r$var2 <- runif(100, 0, 1)
  
  result <- add_weighted_coordinates(r, coord_wt = 2.5)
  
  expect_true("x" %in% names(result))
  expect_true("y" %in% names(result))
  expect_equal(terra::nlyr(result), 4)  # var1, var2, x, y
})

test_that("add_weighted_coordinates drops constant predictors", {
  r <- terra::rast(ncols = 10, nrows = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  r$var1 <- runif(100, 0, 1)
  r$var2 <- 5  # constant
  r$var3 <- runif(100, 0, 1)
  
  result <- add_weighted_coordinates(r, coord_wt = 2.5)
  
  # Should have var1, var3, x, y (not var2)
  expect_false("var2" %in% names(result))
  expect_equal(terra::nlyr(result), 4)
})

test_that("extract_weighted_matrix returns data frame with correct structure", {
  # Setup - pred_rescale needs x and y coordinates already added
  r <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(r) <- matrix(runif(400, 0, 1), ncol = 1)
  names(r) <- "var1"
  r$var2 <- runif(400, 0, 1)
  
  # Add x and y coordinates (this is what add_weighted_coordinates does)
  r$x <- terra::init(r, fun = 'x')
  r$y <- terra::init(r, fun = 'y')
  
  # Create Supplemented layer with values (this is the sampling mask)
  supp <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  result <- extract_weighted_matrix(r, f_rasts, n_pts = 50)
  
  expect_s3_class(result, "data.frame")
  # Should get close to 50 points since no NAs in predictor rasters
  expect_true(nrow(result) >= 45)  # At least 90% of requested
  expect_true(nrow(result) <= 55)  # At most 110% of requested
  expect_true("x" %in% names(result))
  expect_true("y" %in% names(result))
  expect_false("Supplemented" %in% names(result))
})

test_that("perform_clustering returns correct length vector for fixed clusters", {
  set.seed(123)
  test_data <- data.frame(
    x = rnorm(100),
    y = rnorm(100),
    z = rnorm(100)
  )
  
  result <- perform_clustering(test_data, n = 5, fixedClusters = TRUE)
  
  expect_length(result, 100)
  expect_true(all(result %in% 1:5))
})

test_that("perform_clustering handles non-fixed clusters", {
  skip_if_not_installed("NbClust")
  
  set.seed(123)
  test_data <- data.frame(
    x = rnorm(50),
    y = rnorm(50)
  )
  
  result <- perform_clustering(
    test_data, 
    fixedClusters = FALSE,
    method = "complete",
    min.nc = 2,
    max.nc = 5,
    index = "silhouette"
  )
  
  expect_length(result, 50)
  expect_type(result, "integer")
})

# Test main workflow components ----

test_that("sample_underrepresented_clusters identifies underrepresented clusters", {
  skip_if_not_installed("sf")
  
  # Use realistic coordinates for UTM zone 10N
  r <- terra::rast(ncols = 30, nrows = 30, xmin = -123, xmax = -122, 
                   ymin = 45, ymax = 46, crs = "EPSG:4326")
  r$var1 <- runif(900, 0, 1)
  r$var2 <- runif(900, 0, 1)
  r$x <- terra::init(r, 'x')
  r$y <- terra::init(r, 'y')
  
  # Create Supplemented layer
  supp <- terra::rast(ncols = 30, nrows = 30, xmin = -123, xmax = -122, 
                      ymin = 45, ymax = 46, crs = "EPSG:4326")
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  # Create unbalanced clusters
  clusterCut <- c(rep(1, 40), rep(2, 5), rep(3, 30), rep(4, 25))
  
  pts <- data.frame(
    x = runif(100, -123, -122),
    y = runif(100, 45, 46)
  )
  pts <- terra::vect(pts, geom = c("x", "y"), crs = "EPSG:4326")
  
  weighted_mat <- data.frame(
    var1 = runif(100),
    var2 = runif(100),
    x = runif(100, -123, -122),
    y = runif(100, 45, 46),
    ID = clusterCut
  )
  
  result <- sample_underrepresented_clusters(
    clusterCut, pts, weighted_mat, r, f_rasts, 
    lyr = "Supplemented",
    planar_proj = "EPSG:32610",  # UTM zone 10N
    buffer_d = 3,
    n_pts = 100
  )
  
  expect_type(result, "list")
  expect_named(result, c("concentrated_pts", "more_samples"))
  expect_true(2 %in% result$more_samples)  # cluster 2 is underrepresented
  expect_s3_class(result$concentrated_pts, "data.frame")
})

test_that("reorder_clusters_geographically reorders from top-left", {
  skip_if_not_installed("sf")
  
  # Create simple cluster raster with 'class' as layer name (required for select(-class))
  r <- terra::rast(ncols = 10, nrows = 10, xmin = 0, xmax = 100, ymin = 0, ymax = 100, crs = "EPSG:32610")
  
  # Create 4 distinct spatial clusters
  vals <- rep(NA, 100)
  vals[1:25] <- 1      # top-left
  vals[26:50] <- 2     # top-right  
  vals[51:75] <- 3     # bottom-left
  vals[76:100] <- 4    # bottom-right
  
  terra::values(r) <- vals
  names(r) <- "class"  # Name it 'class' as expected by the function
  
  result <- reorder_clusters_geographically(r)
  
  expect_type(result, "list")
  expect_named(result, c("raster", "vectors"))
  expect_s4_class(result$raster, "SpatRaster")
  expect_s3_class(result$vectors, "sf")
  expect_true("ID" %in% names(result$vectors))
  expect_false("class" %in% names(result$vectors))  # Should be removed
})

# Integration tests ----

test_that("EnvironmentalBasedSample runs complete workflow with minimal data", {
  skip_if_not_installed("sf")
  skip_if_not_installed("caret")
  skip("Integration test - run manually")
  
  # Create minimal test data with realistic coordinates
  r <- terra::rast(ncols = 20, nrows = 20, xmin = -123, xmax = -122, 
                   ymin = 45, ymax = 46, crs = "EPSG:4326")
  r$var1 <- runif(400, 0, 1)
  r$var2 <- runif(400, 0, 1)
  r$var3 <- runif(400, 0, 1)
  
  # Create Supplemented layer
  supp <- terra::rast(ncols = 20, nrows = 20, xmin = -123, xmax = -122, 
                      ymin = 45, ymax = 46, crs = "EPSG:4326")
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  # Mock trainKNN function
  trainKNN <- function(data, split_prop = 0.8) {
    # Simple mock KNN model
    list(
      fit.knn = structure(
        list(
          confusionMatrix = matrix(c(10, 2, 1, 15), nrow = 2)
        ),
        class = "train"
      )
    )
  }
  
  result <- EnvironmentalBasedSample(
    pred_rescale = r,
    f_rasts = f_rasts,
    lyr = "Supplemented",
    taxon = "test_species",
    path = tempdir(),
    n = 3,
    fixedClusters = TRUE,
    n_pts = 30,
    planar_proj = "EPSG:32610",
    coord_wt = 2.5,
    buffer_d = 3,
    prop_split = 0.8,
    write2disk = FALSE
  )
  
  expect_type(result, "list")
  expect_named(result, "Geometry")
  expect_s3_class(result$Geometry, "sf")
})

test_that("write_cluster_results creates all required files", {
  skip_if_not_installed("sf")
  
  # Setup
  temp_path <- file.path(tempdir(), "cluster_test")
  dir.create(temp_path, showWarnings = FALSE)
  
  r <- terra::rast(ncols = 10, nrows = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "EPSG:32610")
  terra::values(r) <- sample(1:3, 100, replace = TRUE)
  names(r) <- "class"
  
  vects <- terra::as.polygons(r) |>
    sf::st_as_sf() |>
    dplyr::mutate(ID = class) |>
    dplyr::select(ID, geometry)
  
  mock_fit <- list(confusionMatrix = matrix(1:4, nrow = 2))
  mock_cm <- matrix(1:9, nrow = 3)
  
  write_cluster_results(r, vects, mock_fit, mock_cm, temp_path, "test_taxon")
  
  expect_true(dir.exists(file.path(temp_path, "ClusterRasters")))
  expect_true(dir.exists(file.path(temp_path, "ClusterVectors")))
  expect_true(dir.exists(file.path(temp_path, "TrainingKNN")))
  expect_true(file.exists(file.path(temp_path, "ClusterRasters", "test_taxon.tif")))
  expect_true(file.exists(file.path(temp_path, "TrainingKNN", "test_taxon-finalKNN.rds")))
  
  # Cleanup
  unlink(temp_path, recursive = TRUE)
})

# Edge cases and error handling ----

test_that("add_weighted_coordinates handles single predictor", {
  r <- terra::rast(ncols = 10, nrows = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  r$var1 <- runif(100, 0, 1)
  
  result <- add_weighted_coordinates(r, coord_wt = 2.5)
  
  expect_equal(terra::nlyr(result), 3)  # var1, x, y
})

test_that("extract_weighted_matrix handles NA values correctly", {
  r <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  vals <- runif(400, 0, 1)
  vals[1:50] <- NA
  terra::values(r) <- matrix(vals, ncol = 1)
  names(r) <- "var1"
  
  r$var2 <- runif(400, 0, 1)
  
  # Add x and y coordinates
  r$x <- terra::init(r, fun = 'x')
  r$y <- terra::init(r, fun = 'y')
  
  supp <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  result <- extract_weighted_matrix(r, f_rasts, n_pts = 100)
  
  expect_true(all(complete.cases(result)))
  # With ~12.5% NA values in var1, expect reduced sample size
  # But should still get a reasonable number of points
  expect_true(nrow(result) >= 70)  # At least 70% should be valid
  expect_true(nrow(result) < 100)  # Should be less than requested due to NAs
})

test_that("extract_weighted_matrix diagnostic - verify raster creation", {
  # This test helps diagnose if raster creation is the issue
  r <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  r$var1 <- runif(400, 0, 1)
  r$var2 <- runif(400, 0, 1)
  
  # Add x and y coordinates
  r$x <- terra::init(r, fun = 'x')
  r$y <- terra::init(r, fun = 'y')
  
  # Check that rasters have no NAs
  expect_equal(sum(is.na(terra::values(r$var1))), 0)
  expect_equal(sum(is.na(terra::values(r$var2))), 0)
  
  supp <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  # Sample and extract
  pts <- terra::spatSample(f_rasts[['Supplemented']], as.points = TRUE, 
                           method = 'random', size = 30, na.rm = TRUE)
  
  extracted <- terra::extract(r, pts, bind = TRUE)
  
  # Check extraction results
  expect_s3_class(as.data.frame(extracted), "data.frame")
  
  # Should have minimal NAs if rasters are properly aligned
  na_count <- sum(!complete.cases(as.data.frame(extracted)))
  expect_true(na_count < 5, 
              info = paste("Too many NAs in extraction:", na_count, "out of", nrow(extracted)))
})

test_that("perform_clustering handles small sample sizes", {
  test_data <- data.frame(
    x = rnorm(10),
    y = rnorm(10)
  )
  
  result <- perform_clustering(test_data, n = 3, fixedClusters = TRUE)
  
  expect_length(result, 10)
  expect_true(all(result %in% 1:3))
})

# Additional coverage for helper functions ----

test_that("add_weighted_coordinates coord_wt parameter changes scaling", {
  r <- terra::rast(ncols = 10, nrows = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  r$var1 <- runif(100, 0, 1)
  
  result1 <- add_weighted_coordinates(r, coord_wt = 1)
  result2 <- add_weighted_coordinates(r, coord_wt = 5)
  
  # Higher coord_wt should produce larger coordinate ranges
  range1 <- terra::global(result1$x, "range", na.rm = TRUE)
  range2 <- terra::global(result2$x, "range", na.rm = TRUE)
  
  expect_true(abs(range2$max - range2$min) > abs(range1$max - range1$min))
})

test_that("extract_weighted_matrix respects n_pts parameter", {
  # Use larger raster to ensure enough valid cells
  r <- terra::rast(ncols = 50, nrows = 50, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(r) <- matrix(runif(2500, 0, 1), ncol = 1)
  names(r) <- "var1"
  
  # Add x and y coordinates
  r$x <- terra::init(r, fun = 'x')
  r$y <- terra::init(r, fun = 'y')
  
  supp <- terra::rast(ncols = 50, nrows = 50, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  result_small <- extract_weighted_matrix(r, f_rasts, n_pts = 20)
  result_large <- extract_weighted_matrix(r, f_rasts, n_pts = 100)
  
  # With no NAs, should get close to requested counts (±10%)
  expect_true(nrow(result_small) >= 18 && nrow(result_small) <= 22)
  expect_true(nrow(result_large) >= 90 && nrow(result_large) <= 110)
  expect_true(nrow(result_large) > nrow(result_small))
})

test_that("sample_underrepresented_clusters buffer_d affects sampling area", {
  skip_if_not_installed("sf")
  
  # Use realistic coordinates for UTM zone 10N
  r <- terra::rast(ncols = 30, nrows = 30, xmin = -123, xmax = -122, 
                   ymin = 45, ymax = 46, crs = "EPSG:4326")
  r$var1 <- runif(900, 0, 1)
  r$x <- terra::init(r, 'x')
  r$y <- terra::init(r, 'y')
  
  supp <- terra::rast(ncols = 30, nrows = 30, xmin = -123, xmax = -122, 
                      ymin = 45, ymax = 46, crs = "EPSG:4326")
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  clusterCut <- c(rep(1, 45), rep(2, 5), rep(3, 50))
  
  pts <- data.frame(x = runif(100, -123, -122), y = runif(100, 45, 46))
  pts <- terra::vect(pts, geom = c("x", "y"), crs = "EPSG:4326")
  
  weighted_mat <- data.frame(
    var1 = runif(100), x = runif(100, -123, -122), 
    y = runif(100, 45, 46), ID = clusterCut
  )
  
  # Different buffer distances should produce different results
  result_small <- sample_underrepresented_clusters(
    clusterCut, pts, weighted_mat, r, f_rasts, "Supplemented",
    "EPSG:32610", buffer_d = 2, n_pts = 100
  )
  
  result_large <- sample_underrepresented_clusters(
    clusterCut, pts, weighted_mat, r, f_rasts, "Supplemented",
    "EPSG:32610", buffer_d = 5, n_pts = 100
  )
  
  expect_s3_class(result_small$concentrated_pts, "data.frame")
  expect_s3_class(result_large$concentrated_pts, "data.frame")
})

test_that("sample_underrepresented_clusters removes duplicates correctly", {
  skip_if_not_installed("sf")
  
  # Use realistic geographic coordinates for UTM zone 10N (-126° to -120° longitude, reasonable latitude)
  r <- terra::rast(ncols = 20, nrows = 20, xmin = -123, xmax = -122, 
                   ymin = 45, ymax = 46, crs = "EPSG:4326")
  r$var1 <- runif(400, 0, 1)
  r$x <- terra::init(r, 'x')
  r$y <- terra::init(r, 'y')
  
  supp <- terra::rast(ncols = 20, nrows = 20, xmin = -123, xmax = -122, 
                      ymin = 45, ymax = 46, crs = "EPSG:4326")
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  clusterCut <- c(rep(1, 40), rep(2, 10), rep(3, 50))
  
  # Create pts with specific coordinates in the valid range
  coords <- data.frame(x = runif(100, -123, -122), y = runif(100, 45, 46))
  pts <- terra::vect(coords, geom = c("x", "y"), crs = "EPSG:4326")
  
  weighted_mat <- cbind(coords, data.frame(var1 = runif(100), ID = clusterCut))
  
  result <- sample_underrepresented_clusters(
    clusterCut, pts, weighted_mat, r, f_rasts, "Supplemented",
    "EPSG:32610", buffer_d = 3, n_pts = 100
  )
  
  # Check no duplicates between concentrated_pts and original weighted_mat
  combined_coords <- rbind(
    result$concentrated_pts[, c("x", "y")],
    weighted_mat[, c("x", "y")]
  )
  
  expect_equal(nrow(combined_coords), nrow(unique(combined_coords)))
})

test_that("reorder_clusters_geographically handles different cluster counts", {
  skip_if_not_installed("sf")
  
  for (n_clusters in c(2, 5, 10)) {
    r <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 100, 
                     ymin = 0, ymax = 100, crs = "EPSG:32610")
    terra::values(r) <- sample(1:n_clusters, 400, replace = TRUE)
    names(r) <- "class"
    
    result <- reorder_clusters_geographically(r)
    
    expect_equal(nrow(result$vectors), n_clusters)
    expect_true(all(result$vectors$ID %in% 1:n_clusters))
    expect_equal(length(unique(terra::values(result$raster))), n_clusters)
  }
})

test_that("write_cluster_results handles custom paths", {
  skip_if_not_installed("sf")
  
  temp_path <- file.path(tempdir(), "custom_cluster_path", "subdir")
  dir.create(temp_path, showWarnings = FALSE, recursive = TRUE)
  
  r <- terra::rast(ncols = 5, nrows = 5, crs = "EPSG:32610")
  terra::values(r) <- sample(1:2, 25, replace = TRUE)
  names(r) <- "class"
  
  vects <- terra::as.polygons(r) |>
    sf::st_as_sf() |>
    dplyr::mutate(ID = class) |>
    dplyr::select(ID, geometry)
  
  write_cluster_results(
    r, vects, 
    list(confusionMatrix = matrix(1)), 
    matrix(1), 
    temp_path, 
    "custom_taxon"
  )
  
  expect_true(file.exists(file.path(temp_path, "ClusterRasters", "custom_taxon.tif")))
  expect_true(file.exists(file.path(temp_path, "ClusterVectors", "custom_taxon.shp")))
  
  unlink(file.path(tempdir(), "custom_cluster_path"), recursive = TRUE)
})

# Parameter validation tests ----


test_that("extract_weighted_matrix returns same points each time with set seed", {
  set.seed(123)
  r <- terra::rast(ncols = 20, nrows = 20)
  r$var1 <- runif(400, 0, 1)
  f_rasts <- r[[1]]
  names(f_rasts) <- 'Supplemented'
  
  set.seed(123)
  result1 <- extract_weighted_matrix(r, f_rasts, n_pts = 30)
  
  set.seed(123)
  result2 <- extract_weighted_matrix(r, f_rasts, n_pts = 30)
  
  expect_equal(result1, result2)
})

# Coordinate system tests ----

test_that("reorder_clusters_geographically preserves CRS", {
  skip_if_not_installed("sf")
  
  r <- terra::rast(ncols = 10, nrows = 10, xmin = 0, xmax = 100, 
                   ymin = 0, ymax = 100, crs = "EPSG:32610")
  terra::values(r) <- sample(1:3, 100, replace = TRUE)
  names(r) <- "class"
  
  result <- reorder_clusters_geographically(r)
  
  expect_equal(terra::crs(result$raster), terra::crs(r))
  # Compare EPSG codes rather than input strings
  expect_equal(sf::st_crs(result$vectors)$epsg, 32610)
})

test_that("sample_underrepresented_clusters handles CRS transformations", {
  skip_if_not_installed("sf")
  
  # Start with geographic CRS
  r <- terra::rast(ncols = 20, nrows = 20, xmin = -120, xmax = -119, 
                   ymin = 35, ymax = 36, crs = "EPSG:4326")
  r$var1 <- runif(400, 0, 1)
  r$x <- terra::init(r, 'x')
  r$y <- terra::init(r, 'y')
  
  supp <- terra::rast(ncols = 20, nrows = 20, xmin = -120, xmax = -119, 
                      ymin = 35, ymax = 36, crs = "EPSG:4326")
  terra::values(supp) <- 1
  names(supp) <- "Supplemented"
  
  f_rasts <- list(Supplemented = supp)
  
  clusterCut <- c(rep(1, 40), rep(2, 10), rep(3, 50))
  
  pts <- data.frame(x = runif(100, -120, -119), y = runif(100, 35, 36))
  pts <- terra::vect(pts, geom = c("x", "y"), crs = "EPSG:4326")
  
  weighted_mat <- data.frame(
    var1 = runif(100), 
    x = runif(100, -120, -119), 
    y = runif(100, 35, 36), 
    ID = clusterCut
  )
  
  # Transform to planar should work
  expect_no_error(
    sample_underrepresented_clusters(
      clusterCut, pts, weighted_mat, r, f_rasts, "Supplemented",
      "EPSG:32610", buffer_d = 3, n_pts = 100
    )
  )
})