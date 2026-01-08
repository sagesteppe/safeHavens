square_poly <- sf::st_as_sf(
  sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(0,0), c(1000,0), c(1000,1000), c(0,1000), c(0,0)
    )))
  ),
  crs = 3857
)

test_that("PointBasedSample returns expected structure", {

  set.seed(1)

  out <- PointBasedSample(
    polygon = square_poly,
    n = 5,
    reps = 10,
    BS.reps = 50
  )

  expect_type(out, "list")
  expect_named(out, c("SummaryData", "Geometry"))

  expect_s3_class(out$Geometry, "sf")
  expect_equal(nrow(out$Geometry), 5)
})


test_that("PointBasedSample preserves CRS", {

  set.seed(1)

  out <- PointBasedSample(
    polygon = square_poly,
    n = 5,
    reps = 10,
    BS.reps = 50
  )

  expect_equal(
    sf::st_crs(out$Geometry),
    sf::st_crs(square_poly)
  )
})

test_that("PointBasedSample assigns sequential IDs", {

  set.seed(1)

  out <- PointBasedSample(
    polygon = square_poly,
    n = 5,
    reps = 10,
    BS.reps = 50
  )

  expect_true("ID" %in% names(out$Geometry))
  expect_equal(out$Geometry$ID, seq_len(5))
})

test_that("PointBasedSample returns valid SummaryData", {

  set.seed(1)

  out <- PointBasedSample(
    polygon = square_poly,
    n = 5,
    reps = 10,
    BS.reps = 50
  )

  sd <- out$SummaryData

  expect_s3_class(sd, "data.frame")
  expect_named(sd, c("Metric", "Value"))

  expect_true(
    all(c(
      "variance.observed",
      "quantile.0.001",
      "lwr.95.CI",
      "upr.95.CI",
      "Voronoi.reps.asked",
      "Voronoi.reps.received",
      "BS.reps"
    ) %in% sd$Metric)
  )
})

