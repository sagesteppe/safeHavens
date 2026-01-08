test_that("snapR unions geometry and preserves Assigned", {
  poly1 <- sf::st_polygon(list(rbind(
    c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)
  )))
  poly2 <- sf::st_polygon(list(rbind(
    c(1,0), c(2,0), c(2,1), c(1,1), c(1,0)
  )))

  x <- sf::st_sf(
    Assigned = 5,
    geometry = sf::st_sfc(poly1, poly2),
    crs = 3857
  )

  res <- snapR(x)

  expect_s3_class(res, "sf")
  expect_equal(nrow(res), 1)
  expect_equal(res$Assigned, 5)
  expect_true(sf::st_is_valid(res))
  expect_true(sf::st_is(res, c("POLYGON", "MULTIPOLYGON")))
})

