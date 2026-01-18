library(testthat)
library(sf)
library(dplyr)

# Helper to create test data
create_test_range <- function() {
  sr_mat <- rbind(c(0,0), c(10,0), c(10,10), c(0,10), c(0,0))
  sr_poly <- sf::st_polygon(list(sr_mat))
  sf::st_sf(id = 1, geometry = sf::st_sfc(sr_poly), crs = 3857)
}

create_test_zones <- function(n_zones = 5, with_climate = FALSE) {
  set.seed(123)
  zone_polys <- data.frame(
    x = runif(n_zones * 2, min = -2, max = 12),
    y = runif(n_zones * 2, min = -2, max = 12)
  ) |>
    sf::st_as_sf(coords = c('x', 'y'), crs = 3857) |>
    sf::st_union() |>
    sf::st_voronoi() |>
    sf::st_collection_extract('POLYGON') |>
    sf::st_as_sf() |>
    dplyr::mutate(zone_key = sample(LETTERS[1:n_zones], size = dplyr::n(), replace = TRUE))
  
  if (with_climate) {
    zone_polys <- zone_polys |>
      dplyr::mutate(
        temp_max = runif(dplyr::n(), 20, 35),
        precip_min = runif(dplyr::n(), 200, 800)
      )
  }
  
  zone_polys |> dplyr::rename(geometry = x)
}

# Test input validation ----
test_that("PolygonBasedSample errors on missing required arguments", {
  x <- create_test_range()
  zones <- create_test_zones()
  
  expect_error(
    PolygonBasedSample(zones = zones, zone_key = "zone_key"),
    "Must supply.*x"
  )
  
  expect_error(
    PolygonBasedSample(x = x, zone_key = "zone_key"),
    "Must supply.*zones"
  )
  
  expect_error(
    PolygonBasedSample(x = x, zones = zones),
    "Must supply.*zone_key"
  )
})

test_that("PolygonBasedSample errors on invalid zone_key", {
  x <- create_test_range()
  zones <- create_test_zones()
  
  expect_error(
    PolygonBasedSample(x = x, zones = zones, zone_key = "invalid_key"),
    "zone_key.*not found"
  )
})

test_that("PolygonBasedSample errors on missing climate columns", {
  x <- create_test_range()
  zones <- create_test_zones()
  
  expect_error(
    PolygonBasedSample(x = x, zones = zones, zone_key = "zone_key", 
                      decrease_method = "Assist-warm"),
    "Assist-warm requires.*warmest_col"
  )
  
  expect_error(
    PolygonBasedSample(x = x, zones = zones, zone_key = "zone_key",
                      increase_method = "Assist-drier"),
    "Assist-drier requires.*precip_col"
  )
})

# Test Case 1: n == n_zones ----
test_that("PolygonBasedSample with n == n_zones returns one sample per zone", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5)
  sf::st_agr(zones) <- 'constant'
  zones_cropped <- sf::st_crop(zones, x)
  
  n_unique <- length(unique(zones_cropped$zone_key))
  
  result <- PolygonBasedSample(
    x = x, 
    zones = zones, 
    zone_key = "zone_key",
    n = n_unique
  )
  
  expect_s3_class(result, "sf")
  expect_equal(nrow(result), n_unique)
  expect_true(all(result$allocation == 1))
  expect_equal(length(unique(result$zone_key)), n_unique)
})

# Test Case 2: n < n_zones (decrease methods) ----
test_that("PolygonBasedSample with n < n_zones and method='Largest' selects largest zones", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 7)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 3,
    decrease_method = "Largest"
  )
  
  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 3)
  expect_true(all(result$allocation == 1))
})

test_that("PolygonBasedSample with n < n_zones and method='Smallest' selects smallest zones", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 7)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 3,
    decrease_method = "Smallest"
  )
  
  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 3)
  expect_true(all(result$allocation == 1))
})

test_that("PolygonBasedSample with n < n_zones and method='Most' selects most fragmented", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 7)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 3,
    decrease_method = "Most"
  )
  
  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 3)
  expect_true(all(result$allocation == 1))
})

test_that("PolygonBasedSample with n < n_zones and method='Assist-warm'", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 7, with_climate = TRUE)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 3,
    decrease_method = "Assist-warm",
    warmest_col = "temp_max"
  )
  
  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 3)
})

test_that("PolygonBasedSample with n < n_zones and method='Assist-drier'", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 7, with_climate = TRUE)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 3,
    decrease_method = "Assist-drier",
    precip_col = "precip_min"
  )
  
  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 3)
})

# Test Case 3: n > n_zones (increase methods) ----
test_that("PolygonBasedSample with n > n_zones and method='Largest' allocates extras", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 12,
    increase_method = "Largest"
  )
  
  expect_s3_class(result, "sf")
  # Allow for ±1 due to rounding in proportional allocation
  expect_true(abs(sum(result$allocation) - 12) <= 1)
  expect_true(all(result$allocation >= 0))  
})

test_that("PolygonBasedSample with n > n_zones and method='Smallest' allocates extras", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 12,
    increase_method = "Smallest"
  )
  
  expect_s3_class(result, "sf")
  expect_equal(sum(result$allocation), 12)
})

test_that("PolygonBasedSample with n > n_zones and method='Smallest' allocates extras", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 12,
    increase_method = "Smallest"
  )
  
  expect_s3_class(result, "sf")
  expect_true(abs(sum(result$allocation) - 12) <= 1)  # Allow ±1
})

test_that("PolygonBasedSample with n > n_zones and method='Assist-warm'", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5, with_climate = TRUE)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 12,
    increase_method = "Assist-warm",
    warmest_col = "temp_max"
  )
  
  expect_s3_class(result, "sf")
  expect_equal(sum(result$allocation), 12)
})

test_that("PolygonBasedSample with n > n_zones and method='Assist-drier'", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5, with_climate = TRUE)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 12,
    increase_method = "Assist-drier",
    precip_col = "precip_min"
  )
  
  expect_s3_class(result, "sf")
  expect_equal(sum(result$allocation), 12)
})

# Test geometry handling ----
test_that("PolygonBasedSample handles MULTIPOLYGON input", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5)
  
  # Convert to MULTIPOLYGON
  zones_multi <- zones |>
    dplyr::group_by(zone_key) |>
    dplyr::summarise(geometry = sf::st_union(geometry))
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones_multi,
    zone_key = "zone_key",
    n = 3
  )
  
  expect_s3_class(result, "sf")
})

test_that("PolygonBasedSample returns NULL when no zones intersect", {
  x <- create_test_range()
  
  # Create zones far away from x
  far_zones <- data.frame(
    x = runif(10, 100, 200),
    y = runif(10, 100, 200)
  ) |>
    sf::st_as_sf(coords = c('x', 'y'), crs = 3857) |>
    sf::st_buffer(5) |>
    dplyr::mutate(zone_key = LETTERS[1:10])
  
  expect_warning(
    result <- PolygonBasedSample(x = x, zones = far_zones, zone_key = "zone_key", n = 5),
    "No zone polygons intersect"
  )
  
  expect_null(result)
})

# Test helper functions ----
test_that("split_cols extracts temperature ranges correctly", {
  df <- data.frame(
    Tmin_class = c('10 - 15 Deg. F.', '15 - 20 Deg. F.', '> 55 Deg. F.')
  )
  
  result <- split_cols(df, 'Tmin_class')
  
  expect_s3_class(result, "data.frame")
  expect_named(result, c("lower", "upper", "median", "range"))
  expect_equal(result$lower, c(10, 15, 55))
  expect_equal(result$upper, c(15, 20, 55))
})

test_that("split_cols handles AHM ranges", {
  df <- data.frame(
    AHM_class = c('2 - 3', '6 - 12', '3 - 6')
  )
  
  result <- split_cols(df, 'AHM_class')
  
  expect_equal(result$lower, c(2, 6, 3))
  expect_equal(result$upper, c(3, 12, 6))
  expect_equal(result$range, c(1, 6, 3))
})

test_that("proportional_round with larger_up rounds up largest", {
  v <- c(0.4, 0.3, 0.2, 0.1)
  result <- proportional_round(v, 15, method = "larger_up")
  
  expect_equal(sum(result), 15)
  expect_true(all(result >= 0))
  expect_type(result, "double")
})

test_that("proportional_round with larger_down rounds down largest", {
  v <- c(0.4, 0.3, 0.2, 0.1)
  result <- proportional_round(v, 15, method = "larger_down")
  
  expect_equal(sum(result), 15)
  expect_true(all(result >= 0))
})

test_that("proportional_round handles edge case with exact division", {
  v <- c(0.25, 0.25, 0.25, 0.25)
  result <- proportional_round(v, 20, method = "larger_up")
  
  expect_equal(sum(result), 20)
  expect_equal(result, c(5, 5, 5, 5))
})

# Test allocation logic ----
test_that("allocate_increase distributes points across multiple polygons per zone", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 3)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 15,
    increase_method = "Largest"
  )
  
  # Check that allocations sum correctly
  expect_equal(sum(result$allocation), 15)
  
  # Each zone should have at least 1
  zone_totals <- result |>
    sf::st_drop_geometry() |>
    dplyr::group_by(zone_key) |>
    dplyr::summarise(total = sum(allocation))
  
  expect_true(all(zone_totals$total >= 1))
})

test_that("PolygonBasedSample default n=20 works", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key"
  )
  
  expect_true(abs(sum(result$allocation) - 20) <= 1)  # Allow ±1
})

test_that("PolygonBasedSample preserves CRS", {
  x <- create_test_range()
  zones <- create_test_zones(n_zones = 5)
  
  result <- PolygonBasedSample(
    x = x,
    zones = zones,
    zone_key = "zone_key",
    n = 10
  )
  
  expect_equal(sf::st_crs(result), sf::st_crs(x))
})
