test_that("grid variance helper returns sensible output", {
  ri <- spData::us_states |>
    dplyr::filter(NAME == "Rhode Island") |>
    sf::st_transform(32617)

  out <- grid_variance(ri, 5, 5)

  expect_true(is.numeric(out$var))
  expect_gt(out$n, 0)
})

