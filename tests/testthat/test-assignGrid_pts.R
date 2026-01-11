test_that("assignGrid_pts handles disconnected focal grids", {
  library(sf)
  library(dplyr)
  library(spData)

  # Rhode Island example
  ri <- spData::us_states |>
    filter(NAME == "Rhode Island") |>
    sf::st_transform(32615)

  # Split RI into disconnected parts for focal grid
  polys <- sf::st_cast(sf::st_geometry(ri), "POLYGON")
  focal_grid <- sf::st_sf(polys, ID = seq_along(polys))

  # Neighbor grid: slightly coarser (combine some polygons)
  neighb_grid <- focal_grid[1:3, ] # just first 3 polygons as neighbors

  # Example props (for assignGrid_pts)
  props <- setNames(rep(100 / nrow(focal_grid), nrow(focal_grid)), focal_grid$ID)
  nf_pct <- 100

  pts <- assignGrid_pts(neighb_grid, focal_grid, props, nf_pct)

  # Basic output checks
  expect_s3_class(pts, "sf")
  expect_true(all(c("Assigned", "geometry") %in% colnames(pts)))
  expect_type(pts$Assigned, "double")
  expect_equal(nrow(pts), 100)
  expect_false(any(sf::st_is_empty(pts)))

  # All points assigned to a neighbor
  expect_true(all(pts$Assigned %in% neighb_grid$ID))

  # At least one point assigned to each neighbor (basic coverage test)
  assigned_counts <- table(pts$Assigned)
  expect_true(all(assigned_counts > 0))
})

test_that("assignGrid_pts handles single neighbor grid", {
  # Create a simple focal grid
  focal_grid <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(c(0,0), c(100,0), c(100,100), c(0,100), c(0,0))))
    ),
    crs = 3857
  )
  
  # Single neighbor
  neighb_grid <- sf::st_sf(
    ID = 1,
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(-10,-10), c(110,-10), c(110,110), c(-10,110), c(-10,-10))))
    ),
    crs = 3857
  )
  
  props <- c("1" = 100)
  nf_pct <- 100
  
  pts <- assignGrid_pts(neighb_grid, focal_grid, props, nf_pct)
  
  # With single neighbor, all points assigned to that neighbor
  expect_equal(nrow(pts), 100)
  expect_true(all(pts$Assigned == 1))
})

test_that("assignGrid_pts handles when initial sample is under 100", {
  # Create a very small/narrow polygon that may not yield 100 points initially
  focal_grid <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
    ),
    crs = 3857
  )
  
  neighb_grid <- sf::st_sf(
    ID = c(1, 2),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(-1,-1), c(2,-1), c(2,2), c(-1,2), c(-1,-1)))),
      sf::st_polygon(list(rbind(c(-1,-1), c(2,-1), c(2,2), c(-1,2), c(-1,-1))))
    ),
    crs = 3857
  )
  
  props <- c("1" = 50, "2" = 50)
  nf_pct <- 50
  
  pts <- assignGrid_pts(neighb_grid, focal_grid, props, nf_pct)
  
  # Should eventually get 100 points via the while loop
  expect_equal(nrow(pts), 100)
})

test_that("assignGrid_pts assigns points based on nearest feature and proportions", {
  # Create a rectangular focal grid
  focal_grid <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(c(0,0), c(200,0), c(200,100), c(0,100), c(0,0))))
    ),
    crs = 3857
  )
  
  # Two neighbors on left and right sides
  neighb_grid <- sf::st_sf(
    ID = c(1, 2),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(-50,0), c(50,0), c(50,100), c(-50,100), c(-50,0)))),
      sf::st_polygon(list(rbind(c(150,0), c(250,0), c(250,100), c(150,100), c(150,0))))
    ),
    crs = 3857
  )
  
  props <- c("1" = 60, "2" = 40)
  nf_pct <- 70
  
  pts <- assignGrid_pts(neighb_grid, focal_grid, props, nf_pct)
  
  expect_equal(nrow(pts), 100)
  expect_true(all(pts$Assigned %in% c(1, 2)))
  
  # Check that assignments roughly follow proportions (with some tolerance)
  counts <- table(pts$Assigned)
  expect_true(counts["1"] > 40 && counts["1"] < 80)
  expect_true(counts["2"] > 20 && counts["2"] < 60)
})

test_that("assignGrid_pts handles unassigned points in center", {
  # Create focal grid with multiple neighbors around it
  focal_grid <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(c(0,0), c(100,0), c(100,100), c(0,100), c(0,0))))
    ),
    crs = 3857
  )
  
  # Four neighbors around the edges
  neighb_grid <- sf::st_sf(
    ID = 1:4,
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(-20,0), c(30,0), c(30,100), c(-20,100), c(-20,0)))),    # left
      sf::st_polygon(list(rbind(c(70,0), c(120,0), c(120,100), c(70,100), c(70,0)))),     # right
      sf::st_polygon(list(rbind(c(0,-20), c(100,-20), c(100,30), c(0,30), c(0,-20)))),    # bottom
      sf::st_polygon(list(rbind(c(0,70), c(100,70), c(100,120), c(0,120), c(0,70))))      # top
    ),
    crs = 3857
  )
  
  props <- c("1" = 25, "2" = 25, "3" = 25, "4" = 25)
  nf_pct <- 40
  
  pts <- assignGrid_pts(neighb_grid, focal_grid, props, nf_pct)
  
  # All points should be assigned (none NA)
  expect_equal(nrow(pts), 100)
  expect_false(any(is.na(pts$Assigned)))
  expect_true(all(pts$Assigned %in% 1:4))
})

test_that("assignGrid_pts handles disconnected points correctly", {
  # Elongated focal grid
  focal_grid <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(c(0,0), c(500,0), c(500,50), c(0,50), c(0,0))))
    ),
    crs = 3857
  )
  
  # Neighbors at opposite ends
  neighb_grid <- sf::st_sf(
    ID = c(1, 2),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(-50,-50), c(100,-50), c(100,100), c(-50,100), c(-50,-50)))),
      sf::st_polygon(list(rbind(c(400,-50), c(550,-50), c(550,100), c(400,100), c(400,-50))))
    ),
    crs = 3857
  )
  
  props <- c("1" = 50, "2" = 50)
  nf_pct <- 60
  
  pts <- assignGrid_pts(neighb_grid, focal_grid, props, nf_pct)
  
  expect_equal(nrow(pts), 100)
  expect_true(all(pts$Assigned %in% c(1, 2)))
  
  # Should handle disconnected assignments
  counts <- table(pts$Assigned)
  expect_true(all(counts > 0))
})

test_that("assignGrid_pts output has correct structure", {
  focal_grid <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(c(0,0), c(100,0), c(100,100), c(0,100), c(0,0))))
    ),
    crs = 3857
  )
  
  neighb_grid <- sf::st_sf(
    ID = c(1, 2),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(-10,-10), c(60,-10), c(60,110), c(-10,110), c(-10,-10)))),
      sf::st_polygon(list(rbind(c(40,-10), c(110,-10), c(110,110), c(40,110), c(40,-10))))
    ),
    crs = 3857
  )
  
  props <- c("1" = 50, "2" = 50)
  nf_pct <- 60
  
  pts <- assignGrid_pts(neighb_grid, focal_grid, props, nf_pct)
  
  # Check structure
  expect_s3_class(pts, "sf")
  expect_equal(colnames(pts), c("Assigned", "geometry"))
  expect_type(pts$Assigned, "double")
  expect_equal(nrow(pts), 100)
  expect_false(any(sf::st_is_empty(pts)))
})