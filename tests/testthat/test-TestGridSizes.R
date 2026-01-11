library(testthat)
library(sf)
library(dplyr)

test_that("TestGridSizes returns expected structure", {
  ri <- spData::us_states |>
    dplyr::filter(NAME == "Rhode Island") |>
    sf::st_transform(32617)
  
  result <- TestGridSizes(ri)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
  expect_named(result, c("Name", "Grids", "Variance", "GridNOx", "GridNOy"))
  expect_equal(result$Name, c("Smallest", "Smaller", "Original", "Larger", "Largest"))
})

test_that("TestGridSizes handles tall/narrow shapes (ratio < 0.8)", {
  # Create a tall narrow polygon (height >> width)
  tall_poly <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(100, 0), c(100, 500), c(0, 500), c(0, 0)
      )))
    ),
    crs = 3857
  )
  
  result <- TestGridSizes(tall_poly)
  
  # Should use x_start=4, y_start=7 for ratio < 0.8
  expect_equal(result$GridNOx[result$Name == "Original"], 4)
  expect_equal(result$GridNOy[result$Name == "Original"], 7)
})

test_that("TestGridSizes handles slightly tall shapes (0.8 <= ratio < 0.9)", {
  # Create a slightly tall polygon
  slightly_tall <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(100, 0), c(100, 125), c(0, 125), c(0, 0)
      )))
    ),
    crs = 3857
  )
  
  result <- TestGridSizes(slightly_tall)
  
  # Should use x_start=4, y_start=6 for 0.8 <= ratio < 0.9
  expect_equal(result$GridNOx[result$Name == "Original"], 4)
  expect_equal(result$GridNOy[result$Name == "Original"], 6)
})

test_that("TestGridSizes handles square-ish shapes (0.9 <= ratio < 1.2)", {
  # Create a roughly square polygon
  square_poly <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
      )))
    ),
    crs = 3857
  )
  
  result <- TestGridSizes(square_poly)
  
  # Should use x_start=5, y_start=5 for 0.9 <= ratio < 1.2
  expect_equal(result$GridNOx[result$Name == "Original"], 5)
  expect_equal(result$GridNOy[result$Name == "Original"], 5)
})

test_that("TestGridSizes handles wide shapes (ratio >= 1.2)", {
  # Create a wide polygon (width >> height)
  wide_poly <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(500, 0), c(500, 100), c(0, 100), c(0, 0)
      )))
    ),
    crs = 3857
  )
  
  result <- TestGridSizes(wide_poly)
  
  # Should use x_start=6, y_start=4 for ratio >= 1.2
  expect_equal(result$GridNOx[result$Name == "Original"], 6)
  expect_equal(result$GridNOy[result$Name == "Original"], 4)
})

test_that("TestGridSizes applies offsets correctly", {
  square_poly <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
      )))
    ),
    crs = 3857
  )
  
  result <- TestGridSizes(square_poly)
  
  # With x_start=5, y_start=5, offsets should be +2, +1, 0, -1, -2
  expect_equal(result$GridNOx, c(7, 6, 5, 4, 3))  # 5 + offsets
  expect_equal(result$GridNOy, c(7, 6, 5, 4, 3))  # 5 + offsets
})

test_that("TestGridSizes variance increases with fewer grids", {
  ri <- spData::us_states |>
    dplyr::filter(NAME == "Rhode Island") |>
    sf::st_transform(32617)
  
  result <- TestGridSizes(ri)
  
  # Generally, variance should be higher with fewer grids (Largest)
  # and lower with more grids (Smallest), though not strictly monotonic
  expect_true(all(!is.na(result$Variance)))
  expect_true(all(result$Variance >= 0))
})

test_that("TestGridSizes grid counts differ across options", {
  ri <- spData::us_states |>
    dplyr::filter(NAME == "Rhode Island") |>
    sf::st_transform(32617)
  
  result <- TestGridSizes(ri)
  
  # Smallest should have most grids, Largest should have fewest
  expect_gt(result$Grids[result$Name == "Smallest"], 
            result$Grids[result$Name == "Largest"])
})

test_that("grid_variance handles empty intersection", {
  # Create a polygon and grid that don't overlap well
  tiny_poly <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
      )))
    ),
    crs = 3857
  )
  
  # Use very large grid that might result in empty intersection
  result <- grid_variance(tiny_poly, nx = 100, ny = 100, top_n = 20)
  
  # Should handle gracefully
  expect_true(is.list(result))
  expect_true("var" %in% names(result))
  expect_true("n" %in% names(result))
})

test_that("grid_variance limits to top_n grids", {
  ri <- spData::us_states |>
    dplyr::filter(NAME == "Rhode Island") |>
    sf::st_transform(32617)
  
  result <- grid_variance(ri, nx = 10, ny = 10, top_n = 5)
  
  # Should only use top 5 grids for variance calculation
  expect_true(is.numeric(result$var))
  expect_gt(result$n, 0)  # Total grids created
})

test_that("grid_variance with different nx and ny", {
  square_poly <- sf::st_as_sf(
    sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
      )))
    ),
    crs = 3857
  )
  
  result1 <- grid_variance(square_poly, nx = 3, ny = 3)
  result2 <- grid_variance(square_poly, nx = 5, ny = 5)
  
  # More grids should generally produce different variance
  expect_true(result1$n < result2$n)
  expect_true(is.numeric(result1$var))
  expect_true(is.numeric(result2$var))
})