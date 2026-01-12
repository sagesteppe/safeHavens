library(testthat)
library(sf)
library(dplyr)

test_that("reduceFinalGrids reduces grids to top 20 when >20 provided", {
  # Create 25 grids of varying sizes
  grids <- lapply(1:25, function(i) {
    size <- 100 - i * 2  # Decreasing sizes
    sf::st_polygon(list(rbind(
      c(i * 150, 0), c(i * 150 + size, 0),
      c(i * 150 + size, size), c(i * 150, size),
      c(i * 150, 0)
    )))
  })
  
  final_grids <- sf::st_sf(
    Assigned = 1:25,
    geometry = sf::st_sfc(grids, crs = 3857)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  # Should reduce to 20 grids
  expect_equal(nrow(result), 20)
  expect_s3_class(result, "sf")
})

test_that("reduceFinalGrids handles exactly 20 grids without merging", {
  # Create exactly 20 grids
  grids <- lapply(1:20, function(i) {
    sf::st_polygon(list(rbind(
      c(i*100, 0), c(i*100+80, 0),
      c(i*100+80, 80), c(i*100, 80),
      c(i*100, 0)
    )))
  })
  
  final_grids <- sf::st_sf(
    Assigned = 1:20,
    geometry = sf::st_sfc(grids, crs = 3857)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  # Should still have 20 grids (no merging)
  expect_equal(nrow(result), 20)
  expect_equal(result$Assigned, 1:20)
})

test_that("reduceFinalGrids handles fewer than 20 grids", {
  # Create only 15 grids
  grids <- lapply(1:15, function(i) {
    sf::st_polygon(list(rbind(
      c(i*100, 0), c(i*100+80, 0),
      c(i*100+80, 80), c(i*100, 80),
      c(i*100, 0)
    )))
  })
  
  final_grids <- sf::st_sf(
    Assigned = 1:15,
    geometry = sf::st_sfc(grids, crs = 3857)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  # Should keep all 15 grids
  expect_equal(nrow(result), 15)
  expect_equal(result$Assigned, 1:15)
})

test_that("reduceFinalGrids assigns sequential IDs", {
  # Create 25 simple grids
  grids <- lapply(1:25, function(i) {
    sf::st_polygon(list(rbind(
      c(i*100, 0), c(i*100+80, 0),
      c(i*100+80, 80), c(i*100, 80),
      c(i*100, 0)
    )))
  })
  
  final_grids <- sf::st_sf(
    Assigned = 1:25,
    geometry = sf::st_sfc(grids, crs = 3857)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  # IDs should be sequential from 1 to 20
  expect_equal(result$Assigned, 1:20)
  expect_true(all(diff(result$Assigned) == 1))
})

test_that("reduceFinalGrids orders by Y then X coordinates", {
  # Create grids in a known pattern
  grids <- list()
  idx <- 1
  for (y in c(200, 100, 0)) {  # Descending Y
    for (x in seq(0, 700, by = 100)) {  # Ascending X
      grids[[idx]] <- sf::st_polygon(list(rbind(
        c(x, y), c(x+80, y),
        c(x+80, y+80), c(x, y+80),
        c(x, y)
      )))
      idx <- idx + 1
      if (idx > 25) break
    }
    if (idx > 25) break
  }
  
  final_grids <- sf::st_sf(
    Assigned = seq_along(length(grids)),
    geometry = sf::st_sfc(grids, crs = 3857)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  # Check ordering: should be arranged by -Y (descending), then X (ascending)
  sf::st_agr(result) <- 'constant'
  coords <- sf::st_coordinates(sf::st_point_on_surface(result))
  
  # First grid should have highest Y
  expect_true(coords[1, "Y"] >= coords[nrow(coords), "Y"])
})

test_that("reduceFinalGrids merges smallest grids into neighbors", {
  # Create 20 large grids
  large_grids <- lapply(1:20, function(i) {
    sf::st_polygon(list(rbind(
      c(i*200, 0), c(i*200+150, 0),
      c(i*200+150, 150), c(i*200, 150),
      c(i*200, 0)
    )))
  })
  
  # Add 5 tiny slivers adjacent to large grids
  small_grids <- lapply(1:5, function(i) {
    sf::st_polygon(list(rbind(
      c(i*200+155, 10), c(i*200+160, 10),
      c(i*200+160, 15), c(i*200+155, 15),
      c(i*200+155, 10)
    )))
  })
  
  all_grids <- c(large_grids, small_grids)
  
  final_grids <- sf::st_sf(
    Assigned = 1:25,
    geometry = sf::st_sfc(all_grids, crs = 3857)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  # Should merge down to 20 grids
  expect_equal(nrow(result), 20)
  
  # All geometries should be valid
  expect_true(all(sf::st_is_valid(result)))
})

test_that("reduceFinalGrids removes X and Y temporary columns", {
  grids <- lapply(1:25, function(i) {
    sf::st_polygon(list(rbind(
      c(i*100, 0), c(i*100+80, 0),
      c(i*100+80, 80), c(i*100, 80),
      c(i*100, 0)
    )))
  })
  
  final_grids <- sf::st_sf(
    Assigned = 1:25,
    geometry = sf::st_sfc(grids, crs = 3857)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  # X and Y should not be in final output
  expect_false("X" %in% colnames(result))
  expect_false("Y" %in% colnames(result))
  
  # Should only have Assigned and geometry
  expect_equal(colnames(result), c("Assigned", "geometry"))
})

test_that("reduceFinalGrids preserves CRS", {
  grids <- lapply(1:25, function(i) {
    sf::st_polygon(list(rbind(
      c(i*100, 0), c(i*100+80, 0),
      c(i*100+80, 80), c(i*100, 80),
      c(i*100, 0)
    )))
  })
  
  final_grids <- sf::st_sf(
    Assigned = 1:25,
    geometry = sf::st_sfc(grids, crs = 32617)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  expect_equal(sf::st_crs(result), sf::st_crs(32617))
})

test_that("reduceFinalGrids with 30 grids reduces to 20", {
  # Create 30 grids
  grids <- lapply(1:30, function(i) {
    sf::st_polygon(list(rbind(
      c(i*100, 0), c(i*100+80, 0),
      c(i*100+80, 80), c(i*100, 80),
      c(i*100, 0)
    )))
  })
  
  final_grids <- sf::st_sf(
    Assigned = 1:30,
    geometry = sf::st_sfc(grids, crs = 3857)
  )
  
  result <- reduceFinalGrids(final_grids)
  
  # Should reduce to 20
  expect_equal(nrow(result), 20)
  expect_equal(result$Assigned, 1:20)
})