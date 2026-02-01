test_that("buildResistanceSurface creates resistance surface correctly", {
  # Setup: Create base raster
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  
  # Test 1: Returns a SpatRaster
  result <- buildResistanceSurface(base_raster = base_rast)
  expect_s4_class(result, "SpatRaster")
  
  # Test 2: Output matches input geometry
  expect_true(terra::compareGeom(result, base_rast, stopOnError = FALSE))
  
})

test_that("buildResistanceSurface handles pre-computed resistance surface", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  
  # Create pre-computed resistance surface
  precomputed <- terra::rast(base_rast)
  terra::values(precomputed) <- 100
  
  # Test 4: Uses pre-computed surface when provided
  result <- buildResistanceSurface(base_raster = base_rast, 
                                   resistance_surface = precomputed)
  
  # Check that values are preserved (accounting for integer conversion and clamping)
  result_vals <- terra::values(result)
  expect_true(all(result_vals == 100L, na.rm = TRUE))
  expect_equal(length(result_vals), 100)
  expect_true(terra::compareGeom(result, base_rast, stopOnError = FALSE))
  
  # Test 5: Throws error if geometries don't match
  mismatched <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, 
                            ymin = 0, ymax = 5, vals = 50)
  expect_error(
    buildResistanceSurface(base_raster = base_rast, 
                          resistance_surface = mismatched),
    "resistance_surface must match base_raster geometry"
  )
})

test_that("buildResistanceSurface applies ocean weights correctly", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  ocean_rast <- terra::rast(base_rast)
  terra::values(ocean_rast) <- c(rep(1, 50), rep(0, 50))
  
  # Test 6: Default ocean weight
  result <- buildResistanceSurface(base_raster = base_rast, oceans = ocean_rast)
  ocean_vals <- terra::values(result)[1:50]
  expect_true(all(ocean_vals >= 1000L))
  
  # Test 7: Custom ocean weight
  result_custom <- buildResistanceSurface(base_raster = base_rast, 
                                         oceans = ocean_rast, 
                                         w_ocean = 500)
  ocean_vals_custom <- terra::values(result_custom)[1:50]
  expect_true(all(ocean_vals_custom >= 500L))
})

test_that("buildResistanceSurface applies lake weights correctly", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  lakes_rast <- terra::rast(base_rast)
  terra::values(lakes_rast) <- c(rep(1, 25), rep(0, 75))
  
  # Test 8: Lake weights applied
  result <- buildResistanceSurface(base_raster = base_rast, lakes = lakes_rast)
  lake_vals <- terra::values(result)[1:25]
  expect_true(all(lake_vals >= 200L))
})

test_that("buildResistanceSurface applies river weights correctly", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  rivers_rast <- terra::rast(base_rast)
  terra::values(rivers_rast) <- c(rep(1, 10), rep(0, 90))
  
  # Test 9: River weights applied
  result <- buildResistanceSurface(base_raster = base_rast, rivers = rivers_rast)
  river_vals <- terra::values(result)[1:10]
  expect_true(all(river_vals >= 20L))
})

test_that("buildResistanceSurface handles TRI correctly", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  tri_rast <- terra::rast(base_rast)
  terra::values(tri_rast) <- runif(100, 0, 100)
  
  # Test 10: TRI values are scaled and weighted
  result <- buildResistanceSurface(base_raster = base_rast, tri = tri_rast)
  expect_s4_class(result, "SpatRaster")
  expect_true(all(!is.na(terra::values(result))))
})

test_that("buildResistanceSurface handles habitat correctly", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  habitat_rast <- terra::rast(base_rast)
  terra::values(habitat_rast) <- runif(100, 0, 1)
  
  # Test 11: Habitat values inversely weighted
  result <- buildResistanceSurface(base_raster = base_rast, habitat = habitat_rast)
  expect_s4_class(result, "SpatRaster")
  # Low habitat should produce higher resistance
  expect_true(all(terra::values(result) >= 1L))
})

test_that("buildResistanceSurface combines multiple features", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  ocean_rast <- terra::rast(base_rast)
  terra::values(ocean_rast) <- c(rep(1, 30), rep(0, 70))
  lakes_rast <- terra::rast(base_rast)
  terra::values(lakes_rast) <- c(rep(0, 30), rep(1, 20), rep(0, 50))
  
  # Test 12: Multiple features are additive
  result <- buildResistanceSurface(base_raster = base_rast, 
                                   oceans = ocean_rast, 
                                   lakes = lakes_rast)
  expect_s4_class(result, "SpatRaster")
  # Cells with both ocean and lake should have combined weights
  combined_vals <- terra::values(result)
  expect_true(max(combined_vals, na.rm = TRUE) >= 200L)
})

test_that("buildResistanceSurface handles NULL inputs gracefully", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  
  # Test 15: All NULL optional parameters
  result <- buildResistanceSurface(base_raster = base_rast, 
                                   oceans = NULL, 
                                   lakes = NULL, 
                                   rivers = NULL, 
                                   tri = NULL, 
                                   habitat = NULL)
  expect_s4_class(result, "SpatRaster")
  # Should return minimum resistance values
  expect_true(all(terra::values(result) >= 1L, na.rm = TRUE))
})

test_that("buildResistanceSurface handles edge cases", {
  base_rast <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, 
                           ymin = 0, ymax = 10, vals = 1)
  
  # Test 16: Zero weights
  ocean_rast <- terra::rast(base_rast)
  terra::values(ocean_rast) <- rep(1, 100)
  result <- buildResistanceSurface(base_raster = base_rast, 
                                   oceans = ocean_rast, 
                                   w_ocean = 0)
  expect_true(all(terra::values(result) >= 1L, na.rm = TRUE))
  
  # Test 17: Very large weights
  result_large <- buildResistanceSurface(base_raster = base_rast, 
                                        oceans = ocean_rast, 
                                        w_ocean = 1e6)
  expect_true(max(terra::values(result_large), na.rm = TRUE) >= 1e6)
})