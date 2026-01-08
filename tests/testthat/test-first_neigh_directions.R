test_that("first_neigh_directions returns logical vector of correct length", {
  from <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 3857)

  destinations <- sf::st_sf(
    ID = 1:2,
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(5, -1), c(6, -1), c(6, 1), c(5, 1), c(5, -1)
      ))),
      sf::st_polygon(list(rbind(
        c(8, -1), c(9, -1), c(9, 1), c(8, 1), c(8, -1)
      )))
    ),
    crs = 3857
  )

  res <- first_neigh_directions(from, destinations)

  expect_type(res, "logical")
  expect_length(res, nrow(destinations))
})

