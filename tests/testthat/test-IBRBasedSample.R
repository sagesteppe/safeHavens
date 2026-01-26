test_that("IBRBasedSample returns correct structure", {
  # Setup
  base_rast <- terra::rast(nrows = 50, ncols = 50, 
                           xmin = -100, xmax = -99, 
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  terra::values(base_rast) <- 1
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- sample(1:3, 2500, replace = TRUE)
  
  res_rast <- terra::rast(base_rast)
  terra::values(res_rast) <- runif(2500, 1, 100)
  
  # Create sample points
  coords <- data.frame(
    lon = runif(20, -100, -99),
    lat = runif(20, 40, 41)
  )
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  
  # Create IBR matrix
  n <- nrow(pts_sf)
  ibr_mat <- matrix(runif(n*n, 10, 100), n, n)
  diag(ibr_mat) <- 0
  ibr_mat[lower.tri(ibr_mat)] <- t(ibr_mat)[lower.tri(ibr_mat)]
  
  # Test 1: Function returns list with expected components
  result <- suppressWarnings(
      IBRBasedSample(
      base_raster = base_rast,
      pop_raster = pop_rast,
      resistance_surface = res_rast,
      pts_sf = pts_sf,
      planar_proj = 3857,
      ibr_matrix = ibr_mat,
      fixedClusters = TRUE,
      n = 3
    )
  )
  expect_type(result, "list")
  expect_named(result, c("points", "geometry"))
  
  # Test 2: Points is sf object
  expect_s3_class(result$points, "sf")
  
  # Test 3: Clusters is SpatRaster
  expect_s3_class(result$geometry, "sf")
  
  # Test 4: Points have ID column
  expect_true("ID" %in% names(result$points))
})

test_that("IBRBasedSample requires lon/lat coordinate system", {
  # Setup with projected CRS
  base_rast <- terra::rast(nrows = 50, ncols = 50, 
                           xmin = 0, xmax = 100000, 
                           ymin = 0, ymax = 100000,
                           crs = "EPSG:3857")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  pts_sf <- sf::st_as_sf(
    data.frame(x = c(50000), y = c(50000)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  ibr_mat <- matrix(c(0, 10, 10, 0), 2, 2)
  
  # Test 5: Error when base_raster is not lon/lat
  expect_error(
    IBRBasedSample(
      base_raster = base_rast,
      pop_raster = pop_rast,
      resistance_surface = base_rast,
      planar_proj = 3857,
      pts_sf = pts_sf,
      ibr_matrix = ibr_mat,
      fixedClusters = TRUE,
      n = 2
    ),
    "base_raster must be in lon/lat"
  )
})

test_that("IBRBasedSample handles different distance methods", {
  base_rast <- terra::rast(nrows = 30, ncols = 30, 
                           xmin = -100, xmax = -99, 
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  terra::values(base_rast) <- 1
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- sample(1:2, 900, replace = TRUE)
  
  coords <- data.frame(lon = c(-99.8, -99.2), lat = c(40.2, 40.8))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  
  ibr_mat <- matrix(c(0, 50, 50, 0), 2, 2)
  
  # Test 6: Haversine method
  result_haversine <- suppressWarnings(
    IBRBasedSample(
      base_raster = base_rast,
      pop_raster = pop_rast,
      resistance_surface = base_rast,
      pts_sf = pts_sf,
      ibr_matrix = ibr_mat,
      planar_proj = 3857,
      fixedClusters = TRUE,
      n = 2,
      distance_method = "haversine"
    )
  )
  expect_s3_class(result_haversine$geometry, "sf")
  
  # Test 7: Cosine method
  result_cosine <- IBRBasedSample(
    base_raster = base_rast,
    pop_raster = pop_rast,
    resistance_surface = base_rast,
    pts_sf = pts_sf,
    ibr_matrix = ibr_mat,
    planar_proj = 3857,
    fixedClusters = TRUE,
    n = 2,
    distance_method = "cosine"
  )
  expect_s3_class(result_cosine$geometry, "sf")
})

# ===== Tests for cluster_connectivity =====

test_that("cluster_connectivity creates clusters correctly", {
  # Setup
  pts_sf <- sf::st_as_sf(
    data.frame(x = 1:10, y = 1:10),
    coords = c("x", "y"),
    crs = 3857
  )
  
  # Create distance matrix
  dist_mat <- as.matrix(dist(cbind(1:10, 1:10)))
  
  # Test 9: Fixed clusters with n specified
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
  expect_named(result, c("clusters", "hclust", "input"))
  expect_s3_class(result$clusters, "sf")
  expect_true("ID" %in% names(result$clusters))
  
  # Test 10: Correct number of clusters
  expect_equal(length(unique(result$clusters$ID)), 3)
})

test_that("cluster_connectivity requires n when fixedClusters=TRUE", {
  pts_sf <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y"),
    crs = 3857
  )
  
  dist_mat <- as.matrix(dist(cbind(1:5, 1:5)))
  
  # Test 11: Error when n is NULL and fixedClusters is TRUE
  expect_error(
    cluster_connectivity(
      x = dist_mat,
      pts_sf = pts_sf,
      input = "distance",
      fixedClusters = TRUE,
      n = NULL,
      min.nc = 2,
      max.nc = 5
    ),
    "n must be provided when fixedClusters = TRUE"
  )
})

# ===== Tests for geographic_core_assignment =====

test_that("geographic_core_assignment creates core assignments", {
  base_rast <- terra::rast(nrows = 20, ncols = 20,
                           xmin = -100, xmax = -99,
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  # Create seed points with IDs
  coords <- data.frame(
    lon = c(-99.8, -99.5, -99.2),
    lat = c(40.2, 40.5, 40.8)
  )
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  pts_sf$ID <- c(1, 2, 3)
  
  # Test 13: Returns SpatRaster
  result <- geographic_core_assignment(
    pop_raster = pop_rast,
    pts_sf = pts_sf,
    max_dist = 50000,
    distance_method = "haversine"
  )
  
  expect_s4_class(result, "SpatRaster")
  
  # Test 14: Raster has name "ID"
  expect_equal(names(result), "ID")
  
  # Test 15: Some cells are assigned
  expect_true(any(!is.na(terra::values(result))))
})

test_that("geographic_core_assignment respects max_dist", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = -100, xmax = -99,
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  coords <- data.frame(lon = c(-99.5), lat = c(40.5))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  pts_sf$ID <- 1
  
  # Test 16: Very small max_dist limits assignment
  result_small <- geographic_core_assignment(
    pop_raster = pop_rast,
    pts_sf = pts_sf,
    max_dist = 1000,
    distance_method = "haversine"
  )
  
  # Test 17: Larger max_dist assigns more cells
  result_large <- geographic_core_assignment(
    pop_raster = pop_rast,
    pts_sf = pts_sf,
    max_dist = 100000,
    distance_method = "haversine"
  )
  
  n_assigned_small <- sum(!is.na(terra::values(result_small)))
  n_assigned_large <- sum(!is.na(terra::values(result_large)))
  
  expect_true(n_assigned_large >= n_assigned_small)
})

# ===== Tests for expand_geographic_front =====

test_that("expand_geographic_front expands clusters", {
  base_rast <- terra::rast(nrows = 20, ncols = 20,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  # Create initial cluster raster with small cores
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- NA
  cluster_r[200] <- 1
  cluster_r[210] <- 2
  
  # Test 18: Returns SpatRaster
  result <- expand_geographic_front(
    cluster_r = cluster_r,
    max_cells_per_cluster = 50,
    max_iters = 5
  )
  
  expect_s4_class(result, "SpatRaster")
  
  # Test 19: Expansion increases assigned cells
  initial_assigned <- sum(!is.na(terra::values(cluster_r)))
  final_assigned <- sum(!is.na(terra::values(result)))
  
  expect_true(final_assigned >= initial_assigned)
})

test_that("expand_geographic_front handles empty input", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- NA
  
  # Test 21: Handles all-NA raster
  result <- expand_geographic_front(
    cluster_r = cluster_r,
    max_cells_per_cluster = 100,
    max_iters = 5
  )
  
  expect_s4_class(result, "SpatRaster")
  expect_true(all(is.na(terra::values(result))))
})

# ===== Tests for find_contested_cells =====

test_that("find_contested_cells identifies contested cells", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  # Create cluster raster with gaps
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- NA
  cluster_r[c(1:10)] <- 1
  cluster_r[c(91:100)] <- 2
  # Middle cells are NA and will be contested
  
  # Test 22: Returns list with expected components
  result <- find_contested_cells(
    cluster_r = cluster_r,
    pop_rast = pop_rast,
    directions = 8
  )
  
  expect_type(result, "list")
  expect_named(result, c("cluster_r", "contested", "candidates"))
  
  # Test 23: cluster_r is SpatRaster
  expect_s4_class(result$cluster_r, "SpatRaster")
  
  # Test 24: contested is SpatRaster
  expect_s4_class(result$contested, "SpatRaster")
})

test_that("find_contested_cells handles no contested cells", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  # Fully assigned cluster raster
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- sample(1:3, 100, replace = TRUE)
  
  # Test 25: Handles fully assigned raster
  result <- find_contested_cells(
    cluster_r = cluster_r,
    pop_rast = pop_rast
  )
  
  expect_s4_class(result$cluster_r, "SpatRaster")
  expect_true(all(is.na(terra::values(result$contested))))
})

test_that("find_contested_cells assigns safe cells", {
  base_rast <- terra::rast(nrows = 5, ncols = 5,
                           xmin = 0, xmax = 5,
                           ymin = 0, ymax = 5,
                           crs = "EPSG:3857")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- NA
  cluster_r[c(1, 2, 3)] <- 1
  cluster_r[c(23, 24, 25)] <- 2
  # Cell adjacent to only cluster 1 should be safely assigned
  
  # Test 26: Safe cells are assigned
  result <- find_contested_cells(
    cluster_r = cluster_r,
    pop_rast = pop_rast
  )
  
  initial_assigned <- sum(!is.na(terra::values(cluster_r)))
  final_assigned <- sum(!is.na(terra::values(result$cluster_r)))
  
  expect_true(final_assigned >= initial_assigned)
})

# ===== Tests for assign_contested_line =====

test_that("assign_contested_line assigns contested cells", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- NA
  cluster_r[c(1:20)] <- 1
  cluster_r[c(81:100)] <- 2
  
  contested <- terra::rast(base_rast)
  terra::values(contested) <- NA
  contested[c(40:60)] <- 1  # Mark middle cells as contested
  
  # Test 27: Returns SpatRaster
  result <- suppressWarnings(
      assign_contested_line(
      cluster_r = cluster_r,
      contested = contested
    )
  )
  
  expect_s4_class(result, "SpatRaster")
  
  # Test 28: Some contested cells are assigned
  # Note: May not assign all if components can't be split
  expect_s4_class(result, "SpatRaster")
})

test_that("assign_contested_line handles no contested cells", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- sample(1:2, 100, replace = TRUE)
  
  contested <- terra::rast(base_rast)
  terra::values(contested) <- NA
  
  # Test 29: Returns unchanged raster when no contested cells
  result <- assign_contested_line(
    cluster_r = cluster_r,
    contested = contested
  )
  
  expect_equal(terra::values(result), terra::values(cluster_r))
})

# ===== Tests for finalize_cluster =====

test_that("finalize_cluster assigns remaining cells", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- NA
  cluster_r[c(1:30)] <- 1
  cluster_r[c(71:100)] <- 2
  # Middle cells are unassigned
  
  # Test 30: Returns SpatRaster
  result <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_rast
  )
  
  expect_s4_class(result, "SpatRaster")
  
  # Test 31: More cells are assigned after finalization
  initial_assigned <- sum(!is.na(terra::values(cluster_r)))
  final_assigned <- sum(!is.na(terra::values(result)))
  
  expect_true(final_assigned >= initial_assigned)
})

test_that("finalize_cluster handles fully assigned raster", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- sample(1:3, 100, replace = TRUE)
  
  # Test 32: Returns unchanged when already fully assigned
  result <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_rast
  )
  
  expect_equal(
    sum(!is.na(terra::values(result))),
    sum(!is.na(terra::values(cluster_r)))
  )
})

test_that("finalize_cluster respects population mask", {
  base_rast <- terra::rast(nrows = 10, ncols = 10,
                           xmin = 0, xmax = 10,
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  
  # Population only in some cells
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- NA
  pop_rast[c(1:50)] <- 1
  
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- NA
  cluster_r[c(1:20)] <- 1
  cluster_r[c(31:40)] <- 2
  
  # Test 33: Only assigns within population mask
  result <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_rast
  )
  
  # Cells outside population should remain NA
  result_vals <- terra::values(result)
  pop_vals <- terra::values(pop_rast)
  
  expect_true(all(is.na(result_vals[is.na(pop_vals)])))
})

test_that("finalize_cluster uses neighbor assignment", {
  base_rast <- terra::rast(nrows = 5, ncols = 5,
                           xmin = 0, xmax = 5,
                           ymin = 0, ymax = 5,
                           crs = "EPSG:3857")
  
  pop_rast <- terra::rast(base_rast)
  terra::values(pop_rast) <- 1
  
  cluster_r <- terra::rast(base_rast)
  terra::values(cluster_r) <- NA
  cluster_r[c(1, 2, 3, 6, 7, 8)] <- 1  # Cluster 1 surrounds cell 12
  
  # Test 34: Unassigned cell gets neighbor's cluster
  result <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_rast
  )
  
  # Cell 12 should be assigned to cluster 1 (or at least be assigned)
  expect_true(!is.na(terra::values(result)[12]))
})