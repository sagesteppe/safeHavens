test_that("snapGrids merges focal grid into neighbors correctly", {
  # Neighbor grids
  neighb_grid <- sf::st_sf(
    ID = c(1, 2),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0,0), c(2,0), c(2,2), c(0,2), c(0,0)
      ))),
      sf::st_polygon(list(rbind(
        c(4,0), c(6,0), c(6,2), c(4,2), c(4,0)
      )))
    ),
    crs = 3857
  )

  # Focal grid between them
  focal_grid <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(2,0), c(4,0), c(4,2), c(2,2), c(2,0)
      )))
    ),
    crs = 3857
  )

  # Fake assignGrid_pts output
  x <- sf::st_sf(
    Assigned = c(1, 2),
    geometry = sf::st_sfc(
      sf::st_point(c(2.2, 1)),
      sf::st_point(c(3.8, 1))
    ),
    crs = 3857
  )

  res <- snapGrids(x, neighb_grid, focal_grid)

  expect_s3_class(res, "sf")
  expect_true(all(sf::st_is_valid(res)))
  expect_true(all(sf::st_is(res, c("POLYGON", "MULTIPOLYGON"))))
})

