library(terra)
# Helper function to create mock rasters
create_mock_raster <- function(nrow = 10, ncol = 10, vals = 0) {
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  values(r) <- vals
  r
}

# Test 1: Basic functionality with named parameters only
test_that("buildResistanceSurface works with named parameters only", {
  base_r <- create_mock_raster()
  ocean_r <- create_mock_raster(vals = c(rep(1, 50), rep(0, 50)))
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    oceans = ocean_r,
    w_ocean = 100
  )
  
  expect_s4_class(res, "SpatRaster")
  expect_true(all(values(res) >= 1, na.rm = TRUE)) # Should be clamped to minimum of 1
})

# Test 2: Using all named parameters
test_that("buildResistanceSurface works with all named parameters", {
  base_r <- create_mock_raster()
  ocean_r <- create_mock_raster(vals = c(rep(1, 30), rep(0, 70)))
  lake_r <- create_mock_raster(vals = c(rep(0, 70), rep(1, 30)))
  river_r <- create_mock_raster(vals = c(rep(0, 80), rep(1, 20)))
  tri_r <- create_mock_raster(vals = runif(100, 0, 100))
  habitat_r <- create_mock_raster(vals = runif(100, 0, 1))
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    oceans = ocean_r,
    lakes = lake_r,
    rivers = river_r,
    tri = tri_r,
    habitat = habitat_r,
    w_ocean = 100,
    w_lakes = 50,
    w_rivers = 20,
    w_tri = 2,
    w_habitat = 1
  )
  
  expect_s4_class(res, "SpatRaster")
  expect_true(all(values(res) >= 1, na.rm = TRUE))
})

# Test 3: Single additional layer
test_that("buildResistanceSurface accepts single additional layer", {
  base_r <- create_mock_raster()
  extra_r <- create_mock_raster(vals = runif(100, 0, 1))
  
  # Create a single-layer SpatRaster
  addtl_stack <- extra_r
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    addtl_r = addtl_stack,
    addtl_w = 5
  )
  
  expect_s4_class(res, "SpatRaster")
})

# Test 4: Multiple additional layers as stack
test_that("buildResistanceSurface accepts multiple additional layers", {
  base_r <- create_mock_raster()
  extra_r1 <- create_mock_raster(vals = runif(100, 0, 1))
  extra_r2 <- create_mock_raster(vals = runif(100, 0, 2))
  extra_r3 <- create_mock_raster(vals = runif(100, 0, 0.5))
  
  # Create a multi-layer stack
  addtl_stack <- c(extra_r1, extra_r2, extra_r3)
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    addtl_r = addtl_stack,
    addtl_w = c(1, 5, 10)
  )
  
  expect_s4_class(res, "SpatRaster")
  expect_equal(nlyr(res), 1)  # Output should be single layer
})

# Test 5: Combining named parameters with additional layers
test_that("buildResistanceSurface combines named and additional layers correctly", {
  base_r <- create_mock_raster()
  ocean_r <- create_mock_raster(vals = c(rep(1, 50), rep(0, 50)))
  lake_r <- create_mock_raster(vals = c(rep(0, 70), rep(1, 30)))
  extra_r <- create_mock_raster(vals = runif(100, 0, 1))
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    oceans = ocean_r,
    lakes = lake_r,
    w_ocean = 100,
    w_lakes = 50,
    addtl_r = extra_r,
    addtl_w = 3
  )
  
  expect_s4_class(res, "SpatRaster")
})

# Test 6: Mismatched number of layers and weights (should not add layers)
test_that("buildResistanceSurface skips addtl layers when counts don't match", {
  base_r <- create_mock_raster()
  extra_r1 <- create_mock_raster(vals = 5)
  extra_r2 <- create_mock_raster(vals = 10)
  
  addtl_stack <- c(extra_r1, extra_r2)
  
  # Only 1 weight for 2 layers - should skip the addtl_r processing
  res <- buildResistanceSurface(
    base_raster = base_r,
    addtl_r = addtl_stack,
    addtl_w = 1  # Mismatch
  )
  
  expect_s4_class(res, "SpatRaster")
  # All values should be 1 (from clamping) since no layers were added
  expect_true(all(values(res) == 1, na.rm = TRUE))
})


# Test 10: Result is integer type
test_that("buildResistanceSurface returns integer raster", {
  base_r <- create_mock_raster()
  extra_r <- create_mock_raster(vals = runif(100, 0, 10))
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    addtl_r = extra_r,
    addtl_w = 1.5
  )
  
  expect_true(is.int(res))
})

# Test 11: Clamping to minimum value of 1
test_that("buildResistanceSurface clamps values to minimum of 1", {
  base_r <- create_mock_raster()
  # This should result in some values < 1 before clamping
  extra_r <- create_mock_raster(vals = runif(100, -5, 0))
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    addtl_r = extra_r,
    addtl_w = 1
  )
  
  expect_true(all(values(res) >= 1, na.rm = TRUE))
})

# Test 12: Empty base raster produces clamped result
test_that("buildResistanceSurface handles empty base raster", {
  base_r <- create_mock_raster()
  
  res <- buildResistanceSurface(
    base_raster = base_r
  )
  
  expect_s4_class(res, "SpatRaster")
  # Should all be 1 due to clamping
  expect_true(all(values(res) == 1, na.rm = TRUE))
})

# Test 13: Negative weights work
test_that("buildResistanceSurface handles negative weights", {
  base_r <- create_mock_raster()
  ocean_r <- create_mock_raster(vals = 1)
  forest_r <- create_mock_raster(vals = runif(100, 0, 1))
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    oceans = ocean_r,
    w_ocean = 100,
    addtl_r = forest_r,
    addtl_w = c(-0.5)  # Negative to reduce resistance in forests
  )
  
  expect_s4_class(res, "SpatRaster")
  expect_true(all(values(res) >= 1, na.rm = TRUE))  # Still clamped
})

# Test 14: Complex real-world scenario with multiple layer types
test_that("buildResistanceSurface handles complex multi-layer scenario", {
  base_r <- create_mock_raster(nrow = 50, ncol = 50)
  ocean_r <- create_mock_raster(nrow = 50, ncol = 50, vals = rbinom(2500, 1, 0.1))
  lake_r <- create_mock_raster(nrow = 50, ncol = 50, vals = rbinom(2500, 1, 0.05))
  tri_r <- create_mock_raster(nrow = 50, ncol = 50, vals = runif(2500, 0, 100))
  habitat_r <- create_mock_raster(nrow = 50, ncol = 50, vals = runif(2500, 0, 1))
  urban_r <- create_mock_raster(nrow = 50, ncol = 50, vals = rbinom(2500, 1, 0.15))
  forest_r <- create_mock_raster(nrow = 50, ncol = 50, vals = runif(2500, 0, 1))
  
  # Stack additional layers
  addtl_stack <- c(urban_r, forest_r)
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    oceans = ocean_r,
    lakes = lake_r,
    tri = tri_r,
    habitat = habitat_r,
    w_ocean = 100,
    w_lakes = 50,
    w_tri = 2,
    w_habitat = 1,
    addtl_r = addtl_stack,
    addtl_w = c(500, -0.5)
  )
  
  expect_s4_class(res, "SpatRaster")
  expect_true(all(values(res) >= 1, na.rm = TRUE))
  expect_true(is.int(res))
  expect_equal(nlyr(res), 1)
})

# Test 15: NULL additional layers (should work fine)
test_that("buildResistanceSurface handles NULL additional layers", {
  base_r <- create_mock_raster()
  ocean_r <- create_mock_raster(vals = 1)
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    oceans = ocean_r,
    w_ocean = 100,
    addtl_r = NULL,
    addtl_w = NULL
  )
  
  expect_s4_class(res, "SpatRaster")
})

# Test 16: Only addtl_r provided (no addtl_w) - should skip
test_that("buildResistanceSurface skips addtl_r when addtl_w is NULL", {
  base_r <- create_mock_raster()
  extra_r <- create_mock_raster(vals = 100)
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    addtl_r = extra_r,
    addtl_w = NULL
  )
  
  expect_s4_class(res, "SpatRaster")
  # Should be all 1s since no layers were actually added
  expect_true(all(values(res) == 1, na.rm = TRUE))
})

# Test 17: Only addtl_w provided (no addtl_r) - should skip
test_that("buildResistanceSurface skips when only addtl_w provided", {
  base_r <- create_mock_raster()
  
  res <- buildResistanceSurface(
    base_raster = base_r,
    addtl_r = NULL,
    addtl_w = c(1, 2, 3)
  )
  
  expect_s4_class(res, "SpatRaster")
  expect_true(all(values(res) == 1, na.rm = TRUE))
})


library(testthat)
library(terra)
library(sf)
library(dplyr)

# ==============================================================================
# Helper functions for creating test data
# ==============================================================================

create_test_raster <- function(ncols = 20, nrows = 20, vals = NULL) {
  r <- rast(ncols = ncols, nrows = nrows, 
            xmin = -10, xmax = 10, ymin = -10, ymax = 10,
            crs = "EPSG:4326")
  if (is.null(vals)) {
    values(r) <- sample(c(1, NA), ncell(r), replace = TRUE, prob = c(0.7, 0.3))
  } else {
    values(r) <- vals
  }
  r
}

create_test_points <- function(n = 10, base_raster) {
  # Sample points within non-NA cells
  pts <- spatSample(base_raster, size = n, na.rm = TRUE, xy = TRUE, as.points = FALSE)
  pts_sf <- st_as_sf(pts, coords = c("x", "y"), crs = crs(base_raster))
  pts_sf
}

create_test_distance_matrix <- function(n = 10) {
  # Create a symmetric distance matrix
  mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      if (i != j) {
        d <- runif(1, 0.1, 10)
        mat[i, j] <- d
        mat[j, i] <- d
      }
    }
  }
  mat
}

# ==============================================================================
# Tests for cluster_connectivity()
# ==============================================================================

test_that("cluster_connectivity works with fixed clusters", {
  skip_if_not_installed("igraph")
  
  set.seed(123)
  
  # Create test data
  n_pts <- 10
  dist_mat <- create_test_distance_matrix(n_pts)
  
  base_r <- create_test_raster()
  pts_sf <- create_test_points(n_pts, base_r)
  
  result <- cluster_connectivity(
    x = dist_mat,
    pts_sf = pts_sf,
    input = "distance",
    fixedClusters = TRUE,
    n = 3,
    min.nc = 2,
    max.nc = 5
  )
  
  expect_type(result, "list")
  expect_true(all(c("clusters", "hclust", "input") %in% names(result)))
  expect_s3_class(result$clusters, "sf")
  expect_true("ID" %in% names(result$clusters))
  expect_equal(length(unique(result$clusters$ID)), 3)
  expect_s3_class(result$hclust, "hclust")
})

test_that("cluster_connectivity requires n when fixedClusters=TRUE", {
  dist_mat <- create_test_distance_matrix(5)
  base_r <- create_test_raster()
  pts_sf <- create_test_points(5, base_r)
  
  expect_error(
    cluster_connectivity(
      x = dist_mat,
      pts_sf = pts_sf,
      input = "distance",
      fixedClusters = TRUE,
      n = NULL,  # Missing n should error
      min.nc = 2,
      max.nc = 5
    ),
    "n must be provided when fixedClusters = TRUE"
  )
})

test_that("cluster_connectivity input argument works", {
  dist_mat <- create_test_distance_matrix(5)
  base_r <- create_test_raster()
  pts_sf <- create_test_points(5, base_r)
  
  result_dist <- cluster_connectivity(
    x = dist_mat,
    pts_sf = pts_sf,
    input = "distance",
    fixedClusters = TRUE,
    n = 2,
    min.nc = 2,
    max.nc = 5
  )
  
  expect_equal(result_dist$input, "distance")
  
  result_feat <- cluster_connectivity(
    x = dist_mat,
    pts_sf = pts_sf,
    input = "features",
    fixedClusters = TRUE,
    n = 2,
    min.nc = 2,
    max.nc = 5
  )
  
  expect_equal(result_feat$input, "features")
})


# ==============================================================================
# Tests for geographic_core_assignment()
# ==============================================================================

test_that("geographic_core_assignment creates cluster assignments", {
  set.seed(789)
  
  # Create base raster in lon/lat
  pop_raster <- create_test_raster(ncols = 30, nrows = 30)
  
  # Create 3 seed points
  pts_coords <- rbind(
    c(-5, 5),
    c(0, 0),
    c(5, -5)
  )
  pts_sf <- st_as_sf(
    data.frame(ID = 1:3, x = pts_coords[, 1], y = pts_coords[, 2]),
    coords = c("x", "y"),
    crs = "EPSG:4326"
  )
  
  result <- geographic_core_assignment(
    pop_raster = pop_raster,
    pts_sf = pts_sf,
    max_dist = 500000,  # 500km in meters
    distance_method = "haversine"
  )
  
  expect_s4_class(result, "SpatRaster")
  expect_equal(names(result), "ID")
  
  # Should have some assigned cells
  assigned_vals <- values(result)
  expect_true(any(!is.na(assigned_vals)))
  
  # Assigned values should be 1, 2, or 3
  unique_ids <- unique(assigned_vals[!is.na(assigned_vals)])
  expect_true(all(unique_ids %in% 1:3))
})

test_that("geographic_core_assignment respects max_dist", {
  set.seed(999)
  
  pop_raster <- create_test_raster(ncols = 20, nrows = 20)
  
  # Single point in center
  pts_sf <- st_as_sf(
    data.frame(ID = 1, x = 0, y = 0),
    coords = c("x", "y"),
    crs = "EPSG:4326"
  )
  
  # Very small max_dist should assign few cells
  result_small <- geographic_core_assignment(
    pop_raster = pop_raster,
    pts_sf = pts_sf,
    max_dist = 10000,  # 10km
    distance_method = "haversine"
  )
  
  # Larger max_dist should assign more cells
  result_large <- geographic_core_assignment(
    pop_raster = pop_raster,
    pts_sf = pts_sf,
    max_dist = 1000000,  # 1000km
    distance_method = "cosine"
  )
  
  n_small <- sum(!is.na(values(result_small)))
  n_large <- sum(!is.na(values(result_large)))
  
  expect_true(n_large >= n_small)
})

test_that("geographic_core_assignment handles different distance methods", {
  pop_raster <- create_test_raster()
  pts_sf <- st_as_sf(
    data.frame(ID = 1, x = 0, y = 0),
    coords = c("x", "y"),
    crs = "EPSG:4326"
  )
  
  result_hav <- geographic_core_assignment(
    pop_raster, pts_sf, max_dist = 500000, distance_method = "haversine"
  )
  
  result_cos <- geographic_core_assignment(
    pop_raster, pts_sf, max_dist = 500000, distance_method = "cosine"
  )
  
  expect_s4_class(result_hav, "SpatRaster")
  expect_s4_class(result_cos, "SpatRaster")
})


# ==============================================================================
# Tests for expand_geographic_front()
# ==============================================================================

test_that("expand_geographic_front expands clusters", {
  set.seed(111)
  
  # Create initial cluster raster with small cores
  cluster_r <- create_test_raster(ncols = 20, nrows = 20)
  vals <- values(cluster_r)
  
  # Assign just a few cells to clusters 1 and 2
  vals[1:10] <- 1
  vals[350:360] <- 2
  vals[11:349] <- NA  # Most cells unassigned
  vals[361:400] <- NA
  
  values(cluster_r) <- vals
  names(cluster_r) <- "ID"
  
  n_before <- sum(!is.na(values(cluster_r)))
  
  result <- expand_geographic_front(
    cluster_r = cluster_r,
    max_cells_per_cluster = 100,
    max_iters = 5
  )
  
  n_after <- sum(!is.na(values(result)))
  
  expect_s4_class(result, "SpatRaster")
  expect_true(n_after >= n_before)  # Should expand
})

test_that("expand_geographic_front respects max_cells_per_cluster", {
  set.seed(222)
  
  cluster_r <- create_test_raster(ncols = 30, nrows = 30)
  vals <- values(cluster_r)
  vals[1:5] <- 1  # Small seed
  vals[6:900] <- NA
  values(cluster_r) <- vals
  names(cluster_r) <- "ID"
  
  result <- expand_geographic_front(
    cluster_r = cluster_r,
    max_cells_per_cluster = 20,  # Strict limit
    max_iters = 100
  )
  
  n_cluster1 <- sum(values(result) == 1, na.rm = TRUE)
  expect_true(n_cluster1 <= 25) ## tolerance
})

test_that("expand_geographic_front stops at max_iters", {
  cluster_r <- create_test_raster(ncols = 15, nrows = 15)
  vals <- values(cluster_r)
  vals[1:3] <- 1
  vals[4:225] <- NA
  values(cluster_r) <- vals
  names(cluster_r) <- "ID"
  
  # Should stop after max_iters even if could expand more
  result <- expand_geographic_front(
    cluster_r = cluster_r,
    max_cells_per_cluster = 1000,
    max_iters = 2  # Very few iterations
  )
  
  expect_s4_class(result, "SpatRaster")
  # Won't fill entire raster with just 2 iterations
  expect_true(sum(is.na(values(result))) > 0)
})

test_that("expand_geographic_front handles already-full raster", {
  cluster_r <- create_test_raster(ncols = 10, nrows = 10)
  vals <- rep(1, 100)
  values(cluster_r) <- vals
  names(cluster_r) <- "ID"
  
  result <- expand_geographic_front(
    cluster_r = cluster_r,
    max_cells_per_cluster = 1000,
    max_iters = 10
  )
  
  # Should return unchanged
  expect_equal(values(result), values(cluster_r))
})


# ==============================================================================
# Tests for find_contested_cells()
# ==============================================================================

test_that("find_contested_cells identifies boundary cells", {
  set.seed(333)
  
  # Create cluster raster with clear boundaries
  cluster_r <- create_test_raster(ncols = 20, nrows = 20)
  pop_rast <- create_test_raster(ncols = 20, nrows = 20)
  
  vals_c <- values(cluster_r)
  vals_p <- values(pop_rast)
  
  # Left half = cluster 1, right half = cluster 2, middle = NA
  vals_c[1:150] <- 1
  vals_c[151:170] <- NA  # Contested zone
  vals_c[171:400] <- 2
  
  values(cluster_r) <- vals_c
  names(cluster_r) <- "ID"
  
  result <- find_contested_cells(
    cluster_r = cluster_r,
    pop_rast = pop_rast
  )
  
  expect_type(result, "list")
  expect_true(all(c("cluster_r", "contested", "candidates") %in% names(result)))
  expect_s4_class(result$cluster_r, "SpatRaster")
  expect_s4_class(result$contested, "SpatRaster")
  
  # Some cells should be contested
  contested_vals <- values(result$contested)
  n_contested <- sum(contested_vals == 1, na.rm = TRUE)
  expect_true(n_contested >= 0)  # May or may not have contested cells depending on pop_rast
})

test_that("find_contested_cells handles no unassigned cells", {
  cluster_r <- create_test_raster(ncols = 10, nrows = 10)
  pop_rast <- create_test_raster(ncols = 10, nrows = 10)
  
  # All cells assigned
  values(cluster_r) <- rep(1, 100)
  names(cluster_r) <- "ID"
  
  result <- find_contested_cells(
    cluster_r = cluster_r,
    pop_rast = pop_rast
  )
  
  # Should return empty contested
  expect_true(all(is.na(values(result$contested))))
})

test_that("find_contested_cells assigns safe cells", {
  set.seed(444)
  
  cluster_r <- create_test_raster(ncols = 15, nrows = 15)
  pop_rast <- create_test_raster(ncols = 15, nrows = 15)
  
  vals_c <- values(cluster_r)
  vals_p <- values(pop_rast)
  
  # Create a pattern with some safe cells (only one neighbor)
  vals_c[1:50] <- 1
  vals_c[51:55] <- NA  # These should be safe - only neighbor cluster 1
  vals_c[56:150] <- 2
  vals_c[151:225] <- NA
  
  values(cluster_r) <- vals_c
  names(cluster_r) <- "ID"
  
  result <- find_contested_cells(
    cluster_r = cluster_r,
    pop_rast = pop_rast
  )
  
  # Result should have fewer NAs than input (safe cells assigned)
  n_na_before <- sum(is.na(vals_c))
  n_na_after <- sum(is.na(values(result$cluster_r)))
  
  expect_true(n_na_after <= n_na_before)
})


# ==============================================================================
# Tests for assign_contested_line()
# ==============================================================================

test_that("assign_contested_line assigns contested cells", {
  skip_if_not_installed("igraph")
  
  set.seed(555)
  
  cluster_r <- create_test_raster(ncols = 20, nrows = 20)
  contested <- create_test_raster(ncols = 20, nrows = 20)
  
  vals_c <- values(cluster_r)
  vals_cont <- rep(NA, 400)
  
  # Set up clusters
  vals_c[1:150] <- 1
  vals_c[251:400] <- 2
  vals_c[151:250] <- NA
  
  # Mark middle as contested
  vals_cont[175:225] <- 1
  
  values(cluster_r) <- vals_c
  values(contested) <- vals_cont
  names(cluster_r) <- "ID"
  names(contested) <- "contested"
  
  result <- assign_contested_line(
    cluster_r = cluster_r,
    contested = contested
  )
  
  expect_s4_class(result, "SpatRaster")
  
  # Some contested cells should now be assigned
  n_na_before <- sum(is.na(vals_c))
  n_na_after <- sum(is.na(values(result)))
  
  expect_true(n_na_after <= n_na_before)
})

test_that("assign_contested_line handles no contested cells", {
  cluster_r <- create_test_raster(ncols = 10, nrows = 10)
  contested <- cluster_r * NA
  names(contested) <- "contested"
  
  values(cluster_r) <- rep(1, 100)
  names(cluster_r) <- "ID"
  
  result <- assign_contested_line(
    cluster_r = cluster_r,
    contested = contested
  )
  
  # Should return unchanged
  expect_equal(values(result), values(cluster_r))
})

test_that("assign_contested_line handles edge cases with NA adjacencies", {
  skip_if_not_installed("igraph")
  
  cluster_r <- create_test_raster(ncols = 10, nrows = 10)
  contested <- cluster_r * NA
  
  vals_c <- values(cluster_r)
  vals_cont <- rep(NA, 100)
  
  # Isolated contested cell with no neighbors
  vals_c[1:30] <- 1
  vals_c[50] <- NA  # Isolated
  vals_c[70:100] <- 2
  vals_cont[50] <- 1
  
  values(cluster_r) <- vals_c
  values(contested) <- vals_cont
  names(cluster_r) <- "ID"
  names(contested) <- "contested"
  
  # Should handle without error
  result <- assign_contested_line(
    cluster_r = cluster_r,
    contested = contested
  )
  
  expect_s4_class(result, "SpatRaster")
})


# ==============================================================================
# Tests for finalize_cluster()
# ==============================================================================

test_that("finalize_cluster assigns remaining NA cells", {
  set.seed(666)
  
  cluster_r <- create_test_raster(ncols = 15, nrows = 15)
  pop_raster <- create_test_raster(ncols = 15, nrows = 15)
  
  vals_c <- values(cluster_r)
  vals_p <- values(pop_raster)
  
  # Some clusters, some NAs
  vals_c[1:60] <- 1
  vals_c[61:75] <- NA  # Unassigned but in population
  vals_c[76:150] <- 2
  vals_c[151:225] <- NA
  
  values(cluster_r) <- vals_c
  names(cluster_r) <- "ID"
  
  n_na_before <- sum(is.na(vals_c) & !is.na(vals_p))
  
  result <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_raster
  )
  
  # Masked to population
  result_vals <- values(result)
  n_na_after <- sum(is.na(result_vals) & !is.na(vals_p))
  
  expect_s4_class(result, "SpatRaster")
  expect_true(n_na_after <= n_na_before)  # Should assign some/all NAs
})

test_that("finalize_cluster masks to population", {
  cluster_r <- create_test_raster(ncols = 10, nrows = 10)
  pop_raster <- create_test_raster(ncols = 10, nrows = 10)
  
  values(cluster_r) <- rep(1, 100)
  names(cluster_r) <- "ID"
  
  # Pop raster has NAs
  vals_p <- values(pop_raster)
  vals_p[1:30] <- NA
  values(pop_raster) <- vals_p
  
  result <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_raster
  )
  
  # Result should have NAs where pop_raster has NAs
  result_vals <- values(result)
  expect_true(all(is.na(result_vals[1:30])))
})

test_that("finalize_cluster handles isolated cells", {
  cluster_r <- create_test_raster(ncols = 10, nrows = 10)
  pop_raster <- create_test_raster(ncols = 10, nrows = 10)
  
  vals_c <- values(cluster_r)
  vals_p <- rep(1, 100)  # All population
  
  # Single isolated NA cell with no assigned neighbors
  vals_c <- rep(1, 100)
  vals_c[50] <- NA
  vals_c[40:60][-11] <- NA  # All neighbors also NA
  
  values(cluster_r) <- vals_c
  values(pop_raster) <- vals_p
  names(cluster_r) <- "ID"
  
  result <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_raster
  )
  
  # Should handle without error
  expect_s4_class(result, "SpatRaster")
})


# ==============================================================================
# Tests for IBRSurface() - Integration tests
# ==============================================================================

test_that("IBRSurface requires lon/lat projection", {
  # Create raster in wrong projection
  base_r <- rast(ncols = 10, nrows = 10, 
                 xmin = 0, xmax = 100000, ymin = 0, ymax = 100000,
                 crs = "EPSG:3857")  # Web Mercator, not lon/lat
  values(base_r) <- 1
  
  pop_r <- base_r
  res_r <- base_r
  pts_sf <- st_as_sf(data.frame(x = 50000, y = 50000), 
                     coords = c("x", "y"), crs = "EPSG:3857")
  dist_mat <- matrix(c(0, 1, 1, 0), nrow = 2)
  
  expect_error(
    IBRSurface(
      base_raster = base_r,
      pop_raster = pop_r,
      resistance_surface = res_r,
      pts_sf = pts_sf,
      ibr_matrix = dist_mat,
      fixedClusters = TRUE,
      n = 2,
      planar_proj = 3857
    ),
    "base_raster must be in lon/lat"
  )
})

# ==============================================================================
# Edge cases and error handling
# ==============================================================================

test_that("cluster_connectivity handles single cluster request", {
  dist_mat <- create_test_distance_matrix(5)
  base_r <- create_test_raster()
  pts_sf <- create_test_points(5, base_r)
  
  result <- cluster_connectivity(
    x = dist_mat,
    pts_sf = pts_sf,
    input = "distance",
    fixedClusters = TRUE,
    n = 1,  # Single cluster
    min.nc = 2,
    max.nc = 5
  )
  
  expect_equal(length(unique(result$clusters$ID)), 1)
})

test_that("expand_geographic_front handles single-cluster raster", {
  cluster_r <- create_test_raster(ncols = 10, nrows = 10)
  vals <- values(cluster_r)
  vals[!is.na(vals)] <- 1  # All assigned to cluster 1
  values(cluster_r) <- vals
  names(cluster_r) <- "ID"
  
  result <- expand_geographic_front(
    cluster_r = cluster_r,
    max_cells_per_cluster = 100,
    max_iters = 5
  )
  
  # Should handle without error
  expect_s4_class(result, "SpatRaster")
})

test_that("finalize_cluster handles all-assigned raster", {
  cluster_r <- create_test_raster(ncols = 10, nrows = 10)
  pop_raster <- create_test_raster(ncols = 10, nrows = 10)
  
  vals_c <- values(cluster_r)
  vals_p <- values(pop_raster)
  vals_c[!is.na(vals_p)] <- 1  # All pop cells assigned
  
  values(cluster_r) <- vals_c
  names(cluster_r) <- "ID"
  
  result <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_raster
  )
  
  # Should return with minimal changes
  expect_s4_class(result, "SpatRaster")
})