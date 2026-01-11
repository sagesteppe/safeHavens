library(testthat)
library(sf)
library(dplyr)

# Create a simple square polygon for testing
square_poly <- sf::st_as_sf(
  sf::st_sfc(
    sf::st_polygon(list(rbind(
      c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)
    )))
  ),
  crs = 3857
)

test_that("OpportunisticSample returns expected structure", {
  set.seed(123)
  
  # Create some existing collections
  existing <- sf::st_as_sf(
    data.frame(
      x = c(250, 750),
      y = c(250, 750)
    ),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 5,
    collections = existing,
    reps = 10,
    BS.reps = 50
  )
  
  expect_type(out, "list")
  expect_named(out, c("SummaryData", "Geometry"))
  expect_s3_class(out$Geometry, "sf")
  expect_s3_class(out$SummaryData, "data.frame")
})

test_that("OpportunisticSample creates correct number of polygons", {
  set.seed(123)
  
  existing <- sf::st_as_sf(
    data.frame(x = c(500), y = c(500)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 8,
    collections = existing,
    reps = 10,
    BS.reps = 50
  )
  
  expect_equal(nrow(out$Geometry), 8)
})

test_that("OpportunisticSample assigns sequential IDs", {
  set.seed(42)
  
  existing <- sf::st_as_sf(
    data.frame(x = c(200, 800), y = c(200, 800)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 6,
    collections = existing,
    reps = 10,
    BS.reps = 50
  )
  
  expect_true("ID" %in% colnames(out$Geometry))
  expect_equal(out$Geometry$ID, seq_len(6))
  expect_true(all(diff(out$Geometry$ID) == 1))
})

test_that("OpportunisticSample SummaryData has correct metrics", {
  set.seed(123)
  
  existing <- sf::st_as_sf(
    data.frame(x = c(300, 700), y = c(300, 700)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 5,
    collections = existing,
    reps = 15,
    BS.reps = 100
  )
  
  sd <- out$SummaryData
  expect_named(sd, c("Metric", "Value"))
  
  expected_metrics <- c(
    "variance.observed",
    "quantile.0.001",
    "lwr.95.CI",
    "upr.95.CI",
    "Voronoi.reps.asked",
    "Voronoi.reps.received",
    "BS.reps"
  )
  
  expect_true(all(expected_metrics %in% sd$Metric))
  expect_equal(nrow(sd), 7)
})

test_that("OpportunisticSample records correct rep counts", {
  set.seed(123)
  
  existing <- sf::st_as_sf(
    data.frame(x = c(500), y = c(500)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  reps_asked <- 20
  bs_reps_asked <- 200
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 5,
    collections = existing,
    reps = reps_asked,
    BS.reps = bs_reps_asked
  )
  
  sd <- out$SummaryData
  
  expect_equal(
    sd$Value[sd$Metric == "Voronoi.reps.asked"],
    reps_asked
  )
  expect_equal(
    sd$Value[sd$Metric == "BS.reps"],
    bs_reps_asked
  )
  # Received may be <= asked due to filtering
  expect_true(
    sd$Value[sd$Metric == "Voronoi.reps.received"] <= reps_asked
  )
})

test_that("OpportunisticSample filters voronoi reps correctly", {
  set.seed(456)
  
  existing <- sf::st_as_sf(
    data.frame(x = c(250, 750), y = c(250, 750)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 5,
    collections = existing,
    reps = 25,
    BS.reps = 50
  )
  
  sd <- out$SummaryData
  received <- sd$Value[sd$Metric == "Voronoi.reps.received"]
  
  # Should have some successful reps
  expect_gt(received, 0)
  # May be less than asked due to filtering
  expect_true(received <= 25)
})

test_that("OpportunisticSample preserves CRS", {
  set.seed(123)
  
  existing <- sf::st_as_sf(
    data.frame(x = c(400, 600), y = c(400, 600)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 5,
    collections = existing,
    reps = 10,
    BS.reps = 50
  )
  
  expect_equal(
    sf::st_crs(out$Geometry),
    sf::st_crs(square_poly)
  )
})

test_that("OpportunisticSample works with default n=20", {
  set.seed(789)
  
  # Use larger polygon for n=20
  large_poly <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(5000, 0), c(5000, 5000), c(0, 5000), c(0, 0)
      )))
    ),
    crs = 3857
  )
  
  existing <- sf::st_as_sf(
    data.frame(x = c(1000, 4000), y = c(1000, 4000)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = large_poly,
    collections = existing,
    # n defaults to 20
    reps = 10,
    BS.reps = 50
  )
  
  expect_equal(nrow(out$Geometry), 20)
})

test_that("OpportunisticSample variance values are sensible", {
  set.seed(123)
  
  existing <- sf::st_as_sf(
    data.frame(x = c(300, 700), y = c(300, 700)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 5,
    collections = existing,
    reps = 15,
    BS.reps = 100
  )
  
  sd <- out$SummaryData
  
  var_obs <- sd$Value[sd$Metric == "variance.observed"]
  quantile_001 <- sd$Value[sd$Metric == "quantile.0.001"]
  lwr_ci <- sd$Value[sd$Metric == "lwr.95.CI"]
  upr_ci <- sd$Value[sd$Metric == "upr.95.CI"]
  
  # All should be numeric and non-negative
  expect_true(is.numeric(var_obs) && var_obs >= 0)
  expect_true(is.numeric(quantile_001) && quantile_001 >= 0)
  expect_true(is.numeric(lwr_ci) && lwr_ci >= 0)
  expect_true(is.numeric(upr_ci) && upr_ci >= 0)
  
  # CI bounds should make sense
  expect_true(lwr_ci <= upr_ci)
})

test_that("OpportunisticSample geometry column is properly named", {
  set.seed(123)
  
  existing <- sf::st_as_sf(
    data.frame(x = c(500), y = c(500)),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 5,
    collections = existing,
    reps = 10,
    BS.reps = 50
  )
  
  # Geometry column should be named "geometry"
  expect_true("geometry" %in% colnames(out$Geometry))
  expect_s3_class(out$Geometry$geometry, "sfc")
})

test_that("OpportunisticSample with multiple existing collections", {
  set.seed(999)
  
  # Create several existing collection points
  existing <- sf::st_as_sf(
    data.frame(
      x = c(200, 400, 600, 800),
      y = c(200, 400, 600, 800)
    ),
    coords = c("x", "y"),
    crs = 3857
  )
  
  out <- OpportunisticSample(
    polygon = square_poly,
    n = 10,
    collections = existing,
    reps = 10,
    BS.reps = 50
  )
  
  expect_equal(nrow(out$Geometry), 10)
  expect_true(all(!is.na(out$Geometry$ID)))
})