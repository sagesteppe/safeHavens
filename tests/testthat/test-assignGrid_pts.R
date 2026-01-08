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
  focal_grid <- sf::st_sf(polys, ID = seq_len(length(polys)))

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

