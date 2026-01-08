test_that("healPolygons preserves rows, IDs, and validity", {
  poly1 <- sf::st_polygon(list(rbind(
    c(0,0), c(2,0), c(2,2), c(0,2), c(0,0)
  )))
  poly2 <- sf::st_polygon(list(rbind(
    c(3,0), c(5,0), c(5,2), c(3,2), c(3,0)
  )))

  x <- sf::st_sf(
    Assigned = c(1, 2),
    geometry = sf::st_sfc(poly1, poly2),
    crs = 3857
  )

  res <- healPolygons(x)

  expect_s3_class(res, "sf")
  expect_equal(nrow(res), nrow(x))
  expect_equal(sort(res$Assigned), sort(x$Assigned))
  expect_true(all(sf::st_is_valid(res)))
  expect_true(all(sf::st_is(res, c("POLYGON", "MULTIPOLYGON"))))
})

