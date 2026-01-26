test_that("populationResistance returns correct structure", {
  # Setup
  base_rast <- terra::rast(nrows = 50, ncols = 50, 
                           xmin = -100, xmax = -99, 
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  terra::values(base_rast) <- 1
  
  # Create sample populations
  pops_coords <- data.frame(
    lon = c(-99.8, -99.5, -99.2),
    lat = c(40.2, 40.5, 40.8)
  )
  pops_sf <- sf::st_as_sf(pops_coords, coords = c("lon", "lat"), crs = 4326)
  
  # Test 1: Function returns a list with expected components
  result <- populationResistance(
    populations_sf = pops_sf,
    base_raster = base_rast,
    buffer_dist = 5000,
    planar_proj = 32614,
    n_points = 10,
    graph_method = "delaunay"
  )
  
  expect_type(result, "list")
  expect_named(result, c("pop_raster", "sampled_points", "spatial_graph", 
                         "edge_list", "ibr_matrix"))
  
  # Test 2: pop_raster is a SpatRaster
  expect_s4_class(result$pop_raster, "SpatRaster")
  
  # Test 3: sampled_points is sf object
  expect_s3_class(result$sampled_points, "sf")
  
  # Test 4: spatial_graph is igraph object
  expect_s3_class(result$spatial_graph, "igraph")
  
  # Test 5: edge_list is data.frame
  expect_s3_class(result$edge_list, "data.frame")
  
  # Test 6: ibr_matrix is numeric matrix
  expect_true(is.matrix(result$ibr_matrix))
  expect_type(result$ibr_matrix, "double")
})

test_that("populationResistance requires planar_proj", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, 
                           xmin = -100, xmax = -99, 
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  
  pops_coords <- data.frame(lon = c(-99.5), lat = c(40.5))
  pops_sf <- sf::st_as_sf(pops_coords, coords = c("lon", "lat"), crs = 4326)
  
  # Test 7: Error when planar_proj is NULL
  expect_error(
    populationResistance(
      populations_sf = pops_sf,
      base_raster = base_rast,
      planar_proj = NULL
    ),
    "arg to `planar proj` is required"
  )
})

test_that("populationResistance handles different graph methods", {
  base_rast <- terra::rast(nrows = 30, ncols = 30, 
                           xmin = -100, xmax = -99, 
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  terra::values(base_rast) <- 1
  
  pops_coords <- data.frame(
    lon = c(-99.8, -99.5, -99.2),
    lat = c(40.2, 40.5, 40.8)
  )
  pops_sf <- sf::st_as_sf(pops_coords, coords = c("lon", "lat"), crs = 4326)
  
  # Test 8: Delaunay graph method
  result_delaunay <- populationResistance(
    populations_sf = pops_sf,
    base_raster = base_rast,
    planar_proj = 32614,
    n_points = 10,
    graph_method = "delaunay"
  )
  expect_s3_class(result_delaunay$spatial_graph, "igraph")
  
  # Test 9: Complete graph method
  result_complete <- populationResistance(
    populations_sf = pops_sf,
    base_raster = base_rast,
    planar_proj = 32614,
    n_points = 10,
    graph_method = "complete"
  )
  expect_s3_class(result_complete$spatial_graph, "igraph")
  
  # Test 10: Complete graph has more edges than Delaunay
  n_edges_delaunay <- igraph::ecount(result_delaunay$spatial_graph)
  n_edges_complete <- igraph::ecount(result_complete$spatial_graph)
  expect_true(n_edges_complete >= n_edges_delaunay)
})

test_that("populationResistance uses pre-computed resistance surface", {
  base_rast <- terra::rast(nrows = 30, ncols = 30, 
                           xmin = -100, xmax = -99, 
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  terra::values(base_rast) <- 1
  
  # Create pre-computed resistance
  res_surf <- terra::rast(base_rast)
  terra::values(res_surf) <- 50
  
  pops_coords <- data.frame(lon = c(-99.5, -99.3), lat = c(40.5, 40.7))
  pops_sf <- sf::st_as_sf(pops_coords, coords = c("lon", "lat"), crs = 4326)
  
  # Test 11: Pre-computed surface is used
  result <- populationResistance(
    populations_sf = pops_sf,
    base_raster = base_rast,
    planar_proj = 32614,
    n_points = 10,
    resistance_surface = res_surf,
    graph_method = "delaunay"
  )
  
  expect_s4_class(result$pop_raster, "SpatRaster")
  expect_true(nrow(result$ibr_matrix) > 0)
})

test_that("populationResistance creates symmetric IBR matrix", {
  base_rast <- terra::rast(nrows = 30, ncols = 30, 
                           xmin = -100, xmax = -99, 
                           ymin = 40, ymax = 41,
                           crs = "EPSG:4326")
  terra::values(base_rast) <- 1
  
  pops_coords <- data.frame(
    lon = c(-99.8, -99.5, -99.2),
    lat = c(40.2, 40.5, 40.8)
  )
  pops_sf <- sf::st_as_sf(pops_coords, coords = c("lon", "lat"), crs = 4326)
  
  result <- populationResistance(
    populations_sf = pops_sf,
    base_raster = base_rast,
    planar_proj = 32614,
    n_points = 15,
    graph_method = "complete"
  )
  
  # Test 12: Matrix is symmetric
  mat <- result$ibr_matrix
  expect_true(isSymmetric(mat))
  
  # Test 13: Diagonal is zero
  expect_true(all(diag(mat) == 0))
})

# ===== Tests for sample_population_cells =====

test_that("sample_population_cells returns correct structure", {
  base_rast <- terra::rast(nrows = 20, ncols = 20, 
                           xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10)
  terra::values(base_rast) <- 1:400
  
  # Test 14: Returns sf object
  result <- sample_population_cells(
    pop_raster = base_rast,
    n_total = 10,
    method = "spread"
  )
  
  expect_s3_class(result, "sf")
  expect_true(all(sf::st_geometry_type(result) == "POINT"))
})

test_that("sample_population_cells respects n_total", {
  base_rast <- terra::rast(nrows = 20, ncols = 20, 
                           xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10)
  terra::values(base_rast) <- 1:400
  
  # Test 15: Correct number of points
  result <- sample_population_cells(
    pop_raster = base_rast,
    n_total = 25,
    method = "spread"
  )
  
  expect_equal(nrow(result), 25)
})

# ===== Tests for build_spatial_graph =====

test_that("build_spatial_graph validates input", {
  # Test 16: Error on non-sf input
  expect_error(
    build_spatial_graph(data.frame(x = 1, y = 1), method = "delaunay"),
    "pts_sf must be an sf object"
  )
})

test_that("build_spatial_graph creates correct graph types", {
  coords <- data.frame(
    x = c(0, 1, 0, 1),
    y = c(0, 0, 1, 1)
  )
  pts_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 3857)
  
  # Test 18: Delaunay graph
  result_del <- build_spatial_graph(pts_sf, method = "delaunay", great_circle = FALSE)
  expect_s3_class(result_del$graph, "igraph")
  expect_true(igraph::ecount(result_del$graph) > 0)
  
  # Test 19: Complete graph
  result_comp <- build_spatial_graph(pts_sf, method = "complete", great_circle = FALSE)
  n <- nrow(pts_sf)
  expected_edges <- choose(n, 2)
  expect_equal(igraph::ecount(result_comp$graph), expected_edges)
})

test_that("build_spatial_graph assigns edge weights", {
  coords <- data.frame(x = c(0, 1, 2), y = c(0, 0, 0))
  pts_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 3857)
  
  result <- build_spatial_graph(pts_sf, method = "complete", great_circle = FALSE)
  
  # Test 20: Edge weights exist
  expect_true(!is.null(igraph::E(result$graph)$weight))
  expect_true(all(igraph::E(result$graph)$weight > 0))
  
  # Test 21: Edges data frame has weight column
  expect_true("weight" %in% names(result$edges))
})

test_that("build_spatial_graph handles great circle vs planar distance", {
  coords <- data.frame(lon = c(-100, -99), lat = c(40, 41))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  
  # Test 22: Great circle distances
  result_gc <- build_spatial_graph(pts_sf, method = "complete", great_circle = TRUE)
  expect_s3_class(result_gc$graph, "igraph")
  
  # Test 23: Planar distances (after projection)
  pts_proj <- sf::st_transform(pts_sf, 3857)
  result_planar <- build_spatial_graph(pts_proj, method = "complete", great_circle = FALSE)
  expect_s3_class(result_planar$graph, "igraph")
})

test_that("build_spatial_graph creates undirected graphs", {
  coords <- data.frame(x = c(0, 1, 2), y = c(0, 1, 0))
  pts_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 3857)
  
  result <- build_spatial_graph(pts_sf, method = "delaunay", great_circle = FALSE)
  
  # Test 24: Graph is undirected
  expect_false(igraph::is_directed(result$graph))
})

# ===== Tests for compute_ibr_edges =====

test_that("compute_ibr_edges returns symmetric matrix", {
  base_rast <- terra::rast(nrows = 20, ncols = 20, 
                           xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  terra::values(base_rast) <- runif(400, 1, 100)
  
  coords <- data.frame(x = c(2, 5, 8), y = c(2, 5, 8))
  pts_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 3857)
  
  edges <- data.frame(from = c(1, 1, 2), to = c(2, 3, 3))
  
  result <- compute_ibr_edges(
    resistance_raster = base_rast,
    pts_sf = pts_sf,
    edges = edges,
    method = "leastcost"
  )
  
  # Test 25: Returns matrix
  expect_true(is.matrix(result))
  
  # Test 26: Matrix is symmetric
  expect_true(isSymmetric(result, check.attributes = FALSE))
})

test_that("compute_ibr_edges handles edge cases", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, 
                           xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10,
                           crs = "EPSG:3857")
  terra::values(base_rast) <- 10
  
  coords <- data.frame(x = c(2, 8), y = c(2, 8))
  pts_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 3857)
  
  edges <- data.frame(from = 1, to = 2)
  
  # Test 27: Single edge
  result <- compute_ibr_edges(
    resistance_raster = base_rast,
    pts_sf = pts_sf,
    edges = edges,
    method = "leastcost"
  )
  
  expect_equal(dim(result), c(2, 2))
  expect_true(!is.na(result[1, 2]))
  expect_equal(result[1, 2], result[2, 1])
})

# ===== Tests for inflate_to_dense =====

test_that("inflate_to_dense creates correct dense matrix", {
  edges <- data.frame(
    from = c(1, 1, 2),
    to = c(2, 3, 3),
    weight = c(10, 20, 15)
  )
  
  result <- inflate_to_dense(edges)
  
  # Test 28: Correct dimensions
  expect_equal(dim(result), c(3, 3))
  
  # Test 29: Diagonal is zero
  expect_true(all(diag(result) == 0))
  
  # Test 30: Symmetric
  expect_true(isSymmetric(result))
  
  # Test 31: Edge weights are preserved
  expect_equal(result[1, 2], 10)
  expect_equal(result[1, 3], 20)
  expect_equal(result[2, 3], 15)
})

test_that("inflate_to_dense handles missing edges", {
  edges <- data.frame(
    from = c(1, 2),
    to = c(2, 4),
    weight = c(5, 10)
  )
  
  result <- inflate_to_dense(edges)
  
  # Test 32: Size determined by max index
  expect_equal(dim(result), c(4, 4))
  
  # Test 33: Unconnected pairs have large finite value
  finite_max <- max(edges$weight) * 2
  expect_equal(result[1, 3], finite_max)
  expect_equal(result[1, 4], finite_max)
  expect_equal(result[3, 4], finite_max)
})

test_that("inflate_to_dense handles single edge", {
  edges <- data.frame(from = 1, to = 2, weight = 42)
  
  result <- inflate_to_dense(edges)
  
  # Test 34: Correct size
  expect_equal(dim(result), c(2, 2))
  
  # Test 35: Value preserved
  expect_equal(result[1, 2], 42)
  expect_equal(result[2, 1], 42)
})

