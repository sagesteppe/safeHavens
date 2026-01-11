# Test suite for EnvironmentalBasedSample helper functions
# Uses dismo package example data (Bradypus occurrence data)

library(testthat)
library(terra)
library(sf)
library(dplyr)

# Setup: Load example data used across tests
setup_test_data <- function() {
  # Load Bradypus occurrence data
  x <- read.csv(file.path(system.file(package = "dismo"), 'ex', 'bradypus.csv'))
  x <- x[, c('lon', 'lat')]
  x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)
  
  # Load predictor rasters
  files <- list.files(
    path = file.path(system.file(package = "dismo"), 'ex'), 
    pattern = 'grd', 
    full.names = TRUE
  )
  predictors <- terra::rast(files)
  
  # Create a simple rescaled version (simulate elastic net output)
  # In real use, these would be rescaled by beta coefficients
  pred_rescale <- predictors[[1:4]]  # Use first 4 layers
  
  list(
    occurrences = x,
    predictors = predictors,
    pred_rescale = pred_rescale,
    planar_proj = 3857  # Web Mercator
  )
}

# ==============================================================================
# Tests for add_weighted_coordinates()
# ==============================================================================

test_that("add_weighted_coordinates adds x and y layers", {
  data <- setup_test_data()
  
  result <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  # Check that x and y layers were added
  expect_true("x" %in% names(result))
  expect_true("y" %in% names(result))
  
  # Check that we have original layers plus x and y
  expect_equal(terra::nlyr(result), terra::nlyr(data$pred_rescale) + 2)
})

test_that("add_weighted_coordinates scales coordinates appropriately", {
  data <- setup_test_data()
  
  result <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  # Extract coordinate ranges
  x_range <- terra::global(result$x, fun = 'range', na.rm = TRUE)
  y_range <- terra::global(result$y, fun = 'range', na.rm = TRUE)
  
  # Coordinates should have non-zero range
  expect_true(abs(x_range$max - x_range$min) > 0)
  expect_true(abs(y_range$max - y_range$min) > 0)
  
  # Check that coordinate range is scaled relative to environmental variables
  env_ranges <- terra::global(result[[1:4]], fun = 'range', na.rm = TRUE)
  max_env_range <- max(abs(env_ranges$max - env_ranges$min))
  coord_range <- abs(x_range$max - x_range$min)
  
  # Coordinate range should be approximately coord_wt times the max env range
  # Allow some tolerance for rounding
  expect_true(coord_range > max_env_range * 2)
  expect_true(coord_range < max_env_range * 3)
})

test_that("add_weighted_coordinates removes zero-variance predictors", {
  data <- setup_test_data()
  
  # Add a constant layer (zero variance)
  data$pred_rescale$constant <- 5
  
  result <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  # Constant layer should be removed (except x and y which are added)
  expect_false("constant" %in% names(result))
})

test_that("add_weighted_coordinates handles different coord_wt values", {
  data <- setup_test_data()
  
  result_low <- add_weighted_coordinates(data$pred_rescale, coord_wt = 1.0)
  result_high <- add_weighted_coordinates(data$pred_rescale, coord_wt = 5.0)
  
  # Extract coordinate ranges
  x_range_low <- terra::global(result_low$x, fun = 'range', na.rm = TRUE)
  x_range_high <- terra::global(result_high$x, fun = 'range', na.rm = TRUE)
  
  coord_range_low <- abs(x_range_low$max - x_range_low$min)
  coord_range_high <- abs(x_range_high$max - x_range_high$min)
  
  # Higher weight should produce larger coordinate range
  expect_true(coord_range_high > coord_range_low)
})

# ==============================================================================
# Tests for extract_weighted_matrix()
# ==============================================================================

test_that("extract_weighted_matrix returns correct structure", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  # Create a simple f_rasts structure
  f_rasts <- list(Supplemented = pred_weighted[[1]])
  
  result <- extract_weighted_matrix(pred_weighted, f_rasts, n_pts = 100)
  
  # Should return a data frame
  expect_s3_class(result, "data.frame")
  
  # Should have the expected number of columns (all predictor layers + x + y)
  expect_equal(ncol(result), terra::nlyr(pred_weighted))
  
  # Should have no more than n_pts rows (may be fewer due to NA removal)
  expect_true(nrow(result) <= 100)
  
  # Should have no NA values
  expect_true(all(complete.cases(result)))
})

test_that("extract_weighted_matrix samples correct number of points", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  f_rasts <- list(Supplemented = pred_weighted[[1]])
  
  # Test different sample sizes
  for (n in c(50, 100, 200)) {
    result <- extract_weighted_matrix(pred_weighted, f_rasts, n_pts = n)
    # May be slightly fewer due to NA removal, but should be close
    expect_true(nrow(result) >= n * 0.8)
    expect_true(nrow(result) <= n)
  }
})

test_that("extract_weighted_matrix includes coordinate columns", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  f_rasts <- list(Supplemented = pred_weighted[[1]])
  
  result <- extract_weighted_matrix(pred_weighted, f_rasts, n_pts = 50)
  
  # Should have x and y columns
  expect_true("x" %in% names(result))
  expect_true("y" %in% names(result))
})

# ==============================================================================
# Tests for perform_clustering()
# ==============================================================================

test_that("perform_clustering returns correct number of clusters with fixedClusters=TRUE", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  f_rasts <- list(Supplemented = pred_weighted[[1]])
  weighted_mat <- extract_weighted_matrix(pred_weighted, f_rasts, n_pts = 100)
  
  for (n_clusters in c(3, 5, 10)) {
    result <- perform_clustering(weighted_mat, n = n_clusters, fixedClusters = TRUE)
    
    # Should return a vector of cluster assignments
    expect_type(result, "integer")
    expect_length(result, nrow(weighted_mat))
    
    # Should have exactly n_clusters unique values
    expect_equal(length(unique(result)), n_clusters)
    
    # All values should be between 1 and n_clusters
    expect_true(all(result >= 1 & result <= n_clusters))
  }
})

test_that("perform_clustering with fixedClusters=FALSE uses NbClust", {
  skip_if_not_installed("NbClust")
  
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  f_rasts <- list(Supplemented = pred_weighted[[1]])
  weighted_mat <- extract_weighted_matrix(pred_weighted, f_rasts, n_pts = 100)
  
  # This may take a while, so use small dataset
  result <- perform_clustering(
    weighted_mat[1:50, ], 
    fixedClusters = FALSE,
    min.nc = 2,
    max.nc = 5,
    method = 'complete'
  )
  
  # Should return cluster assignments
  expect_type(result, "integer")
  expect_length(result, 50)
  
  # Should have between min.nc and max.nc clusters
  expect_true(length(unique(result)) >= 2)
  expect_true(length(unique(result)) <= 5)
})

test_that("perform_clustering produces consistent results with same data", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  f_rasts <- list(Supplemented = pred_weighted[[1]])
  weighted_mat <- extract_weighted_matrix(pred_weighted, f_rasts, n_pts = 100)
  
  # Run clustering twice with same parameters
  result1 <- perform_clustering(weighted_mat, n = 5, fixedClusters = TRUE)
  result2 <- perform_clustering(weighted_mat, n = 5, fixedClusters = TRUE)
  
  # Results should be identical (deterministic)
  expect_identical(result1, result2)
})

# ==============================================================================
# Tests for sample_underrepresented_clusters()
# ==============================================================================

test_that("sample_underrepresented_clusters identifies underrepresented clusters", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  f_rasts <- list(Supplemented = pred_weighted[[1]])
  weighted_mat <- extract_weighted_matrix(pred_weighted, f_rasts, n_pts = 100)
  
  # Get sample points
  pts <- terra::spatSample(
    f_rasts[['Supplemented']], 
    as.points = TRUE,
    method = 'random', 
    size = 100, 
    na.rm = TRUE
  )
  
  # Perform clustering
  clusterCut <- perform_clustering(weighted_mat, n = 5, fixedClusters = TRUE)
  weighted_mat$ID <- factor(clusterCut)
  
  result <- sample_underrepresented_clusters(
    clusterCut, pts, weighted_mat, pred_weighted, 
    list(Supplemented = f_rasts$Supplemented),
    lyr = 'Supplemented',
    planar_proj = data$planar_proj,
    buffer_d = 3,
    n_pts = 100
  )
  
  # Should return a list with two elements
  expect_type(result, "list")
  expect_named(result, c("concentrated_pts", "more_samples"))
  
  # concentrated_pts should be a data frame
  expect_s3_class(result$concentrated_pts, "data.frame")
  
  # more_samples should be numeric vector
  expect_type(result$more_samples, "double")
  
  # more_samples should identify clusters below median count
  cluster_counts <- table(clusterCut)
  median_count <- median(cluster_counts)
  expect_true(all(cluster_counts[result$more_samples] < median_count))
})

test_that("sample_underrepresented_clusters removes duplicates", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  f_rasts <- list(Supplemented = pred_weighted[[1]])
  weighted_mat <- extract_weighted_matrix(pred_weighted, f_rasts, n_pts = 100)
  
  pts <- terra::spatSample(
    f_rasts[['Supplemented']], 
    as.points = TRUE,
    method = 'random', 
    size = 100, 
    na.rm = TRUE
  )
  
  clusterCut <- perform_clustering(weighted_mat, n = 5, fixedClusters = TRUE)
  weighted_mat$ID <- factor(clusterCut)
  
  result <- sample_underrepresented_clusters(
    clusterCut, pts, weighted_mat, pred_weighted, 
    list(Supplemented = f_rasts$Supplemented),
    lyr = 'Supplemented',
    planar_proj = data$planar_proj,
    buffer_d = 3,
    n_pts = 100
  )
  
  # Check that no x,y coordinates from concentrated_pts match weighted_mat
  original_coords <- weighted_mat[, c('x', 'y')]
  new_coords <- result$concentrated_pts[, c('x', 'y')]
  
  # Use anti_join to find coordinates that don't overlap
  non_duplicates <- dplyr::anti_join(new_coords, original_coords, by = c('x', 'y'))
  
  # All new coordinates should be non-duplicates
  expect_equal(nrow(non_duplicates), nrow(new_coords))
})

# ==============================================================================
# Tests for reorder_clusters_geographically()
# ==============================================================================

test_that("reorder_clusters_geographically returns correct structure", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  # Create a simple cluster raster (5 clusters)
  cluster_rast <- pred_weighted[[1]]
  values(cluster_rast) <- sample(1:5, ncell(cluster_rast), replace = TRUE)
  cluster_rast <- terra::mask(cluster_rast, pred_weighted[[1]])
  
  result <- reorder_clusters_geographically(cluster_rast)
  
  # Should return a list with raster and vectors
  expect_type(result, "list")
  expect_named(result, c("raster", "vectors"))
  
  # Raster should be SpatRaster
  expect_s4_class(result$raster, "SpatRaster")
  
  # Vectors should be sf object
  expect_s3_class(result$vectors, "sf")
})

test_that("reorder_clusters_geographically preserves number of clusters", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  n_clusters <- 7
  cluster_rast <- pred_weighted[[1]]
  values(cluster_rast) <- sample(1:n_clusters, ncell(cluster_rast), replace = TRUE)
  cluster_rast <- terra::mask(cluster_rast, pred_weighted[[1]])
  
  result <- reorder_clusters_geographically(cluster_rast)
  
  # Check number of unique values in raster
  unique_vals <- unique(terra::values(result$raster, na.rm = TRUE))
  expect_equal(length(unique_vals), n_clusters)
  
  # Check number of polygons
  expect_equal(nrow(result$vectors), n_clusters)
})

test_that("reorder_clusters_geographically orders from north to south", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  cluster_rast <- pred_weighted[[1]]
  values(cluster_rast) <- sample(1:5, ncell(cluster_rast), replace = TRUE)
  cluster_rast <- terra::mask(cluster_rast, pred_weighted[[1]])
  
  result <- reorder_clusters_geographically(cluster_rast)
  
  # Extract Y coordinates of cluster centroids
  centroids <- sf::st_point_on_surface(result$vectors)
  y_coords <- sf::st_coordinates(centroids)[, 2]
  
  # Y coordinates should generally decrease as ID increases (north to south)
  # Allow some variation but check overall trend
  expect_true(cor(result$vectors$ID, y_coords, method = "spearman") < 0)
})

# ==============================================================================
# Tests for write_cluster_results()
# ==============================================================================

test_that("write_cluster_results creates all expected directories and files", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  # Create test cluster data
  cluster_rast <- pred_weighted[[1]]
  values(cluster_rast) <- sample(1:5, ncell(cluster_rast), replace = TRUE)
  cluster_rast <- terra::mask(cluster_rast, pred_weighted[[1]])
  
  reordered <- reorder_clusters_geographically(cluster_rast)
  
  # Create dummy KNN model and confusion matrix
  mock_fit <- list(confusionMatrix = matrix(1:25, 5, 5))
  mock_cm <- matrix(1:25, 5, 5)
  
  # Use temporary directory
  temp_dir <- tempdir()
  taxon <- "test_species"
  
  write_cluster_results(
    reordered$raster, 
    reordered$vectors, 
    mock_fit, 
    mock_cm,
    temp_dir, 
    taxon
  )
  
  # Check that directories were created
  expect_true(dir.exists(file.path(temp_dir, 'ClusterRasters')))
  expect_true(dir.exists(file.path(temp_dir, 'ClusterVectors')))
  expect_true(dir.exists(file.path(temp_dir, 'TrainingKNN')))
  
  # Check that files were created
  expect_true(file.exists(file.path(temp_dir, 'ClusterRasters', paste0(taxon, '.tif'))))
  expect_true(file.exists(file.path(temp_dir, 'TrainingKNN', paste0(taxon, '-finalKNN.rds'))))
  expect_true(file.exists(file.path(temp_dir, 'TrainingKNN', paste0(taxon, '-finalCM.rds'))))
  expect_true(file.exists(file.path(temp_dir, 'ClusterVectors', paste0(taxon, '.shp'))))
  
  # Clean up
  unlink(file.path(temp_dir, 'ClusterRasters'), recursive = TRUE)
  unlink(file.path(temp_dir, 'ClusterVectors'), recursive = TRUE)
  unlink(file.path(temp_dir, 'TrainingKNN'), recursive = TRUE)
})

test_that("write_cluster_results can be read back correctly", {
  data <- setup_test_data()
  pred_weighted <- add_weighted_coordinates(data$pred_rescale, coord_wt = 2.5)
  
  cluster_rast <- pred_weighted[[1]]
  values(cluster_rast) <- sample(1:5, ncell(cluster_rast), replace = TRUE)
  cluster_rast <- terra::mask(cluster_rast, pred_weighted[[1]])
  
  reordered <- reorder_clusters_geographically(cluster_rast)
  
  mock_fit <- list(confusionMatrix = matrix(1:25, 5, 5))
  mock_cm <- matrix(1:25, 5, 5)
  
  temp_dir <- tempdir()
  taxon <- "test_species"
  
  write_cluster_results(
    reordered$raster, 
    reordered$vectors, 
    mock_fit, 
    mock_cm,
    temp_dir, 
    taxon
  )
  
  # Read back the raster
  rast_read <- terra::rast(file.path(temp_dir, 'ClusterRasters', paste0(taxon, '.tif')))
  expect_s4_class(rast_read, "SpatRaster")
  
  # Read back the vectors
  vect_read <- sf::st_read(
    file.path(temp_dir, 'ClusterVectors', paste0(taxon, '.shp')),
    quiet = TRUE
  )
  expect_s3_class(vect_read, "sf")
  expect_equal(nrow(vect_read), nrow(reordered$vectors))
  
  # Read back the RDS files
  fit_read <- readRDS(file.path(temp_dir, 'TrainingKNN', paste0(taxon, '-finalKNN.rds')))
  cm_read <- readRDS(file.path(temp_dir, 'TrainingKNN', paste0(taxon, '-finalCM.rds')))
  expect_identical(fit_read, mock_fit)
  expect_identical(cm_read, mock_cm)
  
  # Clean up
  unlink(file.path(temp_dir, 'ClusterRasters'), recursive = TRUE)
  unlink(file.path(temp_dir, 'ClusterVectors'), recursive = TRUE)
  unlink(file.path(temp_dir, 'TrainingKNN'), recursive = TRUE)
})

# ==============================================================================
# Integration test for the full workflow
# ==============================================================================

test_that("full workflow integration test", {
  skip_if_not_installed("class")  # Required for KNN
  skip("Integration test - requires trainKNN function")
  
  data <- setup_test_data()
  
  # This would test the full EnvironmentalBasedSample function
  # Skipped because it requires the trainKNN function which we don't have
  # But this is where you'd test the complete workflow
})
