test_that("OpportunisticSample runs on a small example", {

  skip_if_not_installed("sf")
  skip_if_not_installed("spData")
  skip_if_not_installed("nptest")
  skip_on_cran() 

  set.seed(42)
  samples = 15
  # small polygon: Rhode Island
  ri <- spData::us_states |>
    dplyr::filter(NAME == "Rhode Island") |>
    sf::st_transform(32617)

  existing_collections <- sf::st_sample(ri, size = 5) |>
   sf::st_as_sf() |>
   dplyr::rename(geometry = x)

  # call OpportunisticSample with small reps / BS.reps for speed
  res <- OpportunisticSample(
    polygon = ri,
    n = samples,
    collections = existing_collections,
    reps = 10,      # small number for test speed
    BS.reps = 250    # small bootstrap reps
  )

  # top-level structure
  expect_type(res, "list")
  expect_named(res, c("SummaryData", "Geometry"))

  # SummaryData
  expect_s3_class(res$SummaryData, "data.frame")
  expect_true(all(c("Metric", "Value") %in% colnames(res$SummaryData)))
  expect_gt(nrow(res$SummaryData), 0)

  # Geometry
  expect_s3_class(res$Geometry, "sf")
  expect_true(nrow(res$Geometry) == samples)  
  expect_true("geometry" %in% colnames(res$Geometry))

  # IDs are sequential
  expect_equal(res$Geometry$ID, seq_len(nrow(res$Geometry)))
})

