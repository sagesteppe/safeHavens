library(testthat)
library(sf)

test_that("PrioritizeSample returns sf object with correct columns", {

  # Create a simple test polygon
  poly_g <- sf::st_as_sf(
    data.frame(ID = 1),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0,0), c(0,10), c(10,10), c(10,0), c(0,0)
      )))
    ),
    crs = 32617
  )

  res <- PrioritizeSample(poly_g, n_breaks = 3, verbose = FALSE, metric = "var")

  expect_s3_class(res$Geometry, "sf")
  expect_true(all(c("ID", "SampleOrder", "Level", "geometry") %in% colnames(res$Geometry)))
  expect_true(all(res$Geometry$Level %in% 1:3))
  expect_equal(sf::st_crs(res$Geometry), sf::st_crs(poly_g))
})

test_that("PrioritizeSample handles multiple polygons", {

  polys <- sf::st_as_sf(
    data.frame(ID = 1:2),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(0,0),c(0,10),c(10,10),c(10,0),c(0,0)))),
      sf::st_polygon(list(rbind(c(20,0),c(20,10),c(30,10),c(30,0),c(20,0))))
    ),
    crs = 32617
  )

  res <- PrioritizeSample(polys, n_breaks = 2, metric = "sd", verbose = FALSE)
  #expect_equal(length(unique(res$Geometry$ID)), 2)
  expect_true(all(res$Geometry$Level %in% 1:2))
})

test_that("PrioritizeSample stops for longlat coordinates", {
  library(sf)
  poly <- sf::st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE) |>
        dplyr::select(NAME) |>
        sf::st_transform(5070)  #
  poly_ll <- sf::st_transform(poly, 4326)
  expect_error(PrioritizeSample(poly_ll), "requires planar")
})

