library(testthat)
library(sf)

test_that("EqualAreaSample basic functionality", {
  
  nc <- sf::st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE) |>
    dplyr::select(NAME)
  
  set.seed(123)
  
  res <- EqualAreaSample(
    x = nc, 
    n = 5, 
    pts = 100, 
    planar_proj = 5070, 
    returnProjected = FALSE,
    reps = 5,
    BS.reps = 50
  )
  
  # --- Structure checks ---
  expect_true(is.list(res))
  expect_true(all(c("SummaryData", "Geometry") %in% names(res)))
  
  expect_s3_class(res$Geometry, "sf")
  expect_s3_class(res$SummaryData, "data.frame")
  
  expect_equal(colnames(res$SummaryData), c("Metric", "Value"))
  expect_true(all(sapply(res$SummaryData$Value, is.numeric)))
  
  # --- CRS checks ---
  expect_equal(sf::st_crs(res$Geometry), sf::st_crs(nc))
  
  # --- Number of polygons ---
  expect_equal(nrow(res$Geometry), 5) # n polygons requested
  
  # --- Valid geometries ---
  expect_true(all(sf::st_is_valid(res$Geometry)))
  
  # --- Bootstrap / variance reflected correctly ---
  expect_equal(res$SummaryData$Value[res$SummaryData$Metric == "Voronoi.reps.asked"], 5)
  expect_equal(res$SummaryData$Value[res$SummaryData$Metric == "BS.reps"], 50)
  
  # --- Message trigger ---
  expect_error(EqualAreaSample(nc, n = 3), "Argument to `planar_proj` is required")
  
})

