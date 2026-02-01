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