skip_if_not_installed("terra")
skip_if_not_installed("sf")
skip_if_not_installed("caret")
skip_if_not_installed("NbClust")

test_that("IBDBasedSample runs with minimal inputs", {

  set.seed(1)

  r <- terra::rast(
    ncols = 10, nrows = 10,
    xmin = 0, xmax = 10,
    ymin = 0, ymax = 10,
    crs = "EPSG:3857"
  )
  terra::values(r) <- 1

  template <- r

  expect_no_error(
    res <- IBDBasedSample(
      x = r,
      n = 3,
      fixedClusters = TRUE,
      n_pts = 50,
      template = template
    )
  )

  expect_true(is.list(res))
  expect_true("Geometry" %in% names(res))
})

test_that("IBDBasedSample returns sf geometry with IDs", {

  set.seed(1)

  r <- terra::rast(ncols = 10, nrows = 10, crs = "EPSG:3857")
  terra::values(r) <- 1

  res <- IBDBasedSample(
    x = r,
    n = 4,
    fixedClusters = TRUE,
    n_pts = 50,
    template = r
  )

  geom <- res$Geometry

  expect_s3_class(geom, "sf")
  expect_true("ID" %in% names(geom))
  expect_true(all(!is.na(geom$ID)))
})

test_that("IBDBasedSample creates exactly n clusters when fixedClusters = TRUE", {

  set.seed(1)

  r <- terra::rast(ncols = 10, nrows = 10, crs = "EPSG:3857")
  terra::values(r) <- 1

  n_clusters <- 5

  res <- IBDBasedSample(
    x = r,
    n = n_clusters,
    fixedClusters = TRUE,
    n_pts = 60,
    template = r
  )

  expect_equal(
    length(unique(res$Geometry$ID)),
    n_clusters
  )
})

test_that("IBDBasedSample respects planar_proj and returns CRS of x", {

  set.seed(1)

  r_ll <- terra::rast(
    ncols = 10, nrows = 10,
    crs = "EPSG:4326"
  )
  terra::values(r_ll) <- 1

  # Project raster to planar CRS for safe sampling
  r <- terra::project(r_ll, "EPSG:3857")

  planar_proj <- "EPSG:3857"

  res <- IBDBasedSample(
    x = r,
    n = 3,
    fixedClusters = TRUE,
    n_pts = 50,
    template = r,
    planar_proj = planar_proj
  )

  expect_equal(
    sf::st_crs(res$Geometry)$epsg,
    sf::st_crs(r)$epsg
  )
})


test_that("IBDBasedSample runs when fixedClusters = FALSE", {

  skip_if_not_installed("NbClust")

  set.seed(1)

  r <- terra::rast(ncols = 10, nrows = 10, crs = "EPSG:3857")
  terra::values(r) <- 1

  expect_no_error(
    IBDBasedSample(
      x = r,
      n = 3,
      fixedClusters = FALSE,
      n_pts = 50,
      template = r,
      min.nc = 2,
      max.nc = 4
    )
  )
})

test_that("IBDBasedSample errors on empty raster", {

  r <- terra::rast(ncols = 10, nrows = 10)

  expect_error(
    IBDBasedSample(
      x = r,
      n = 3,
      fixedClusters = TRUE,
      n_pts = 50,
      template = r
    )
  )
})

