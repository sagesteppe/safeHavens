#' Get an estimate for how many grids to draw over a species range
#' 
#' @description This function uses the dimensions of a species grid to estimate how many grids
#' would need to be added in the x and y directions to cover it with 20 grid cells
#' of roughly equal areas. 
#' @param target a species range as a simple feature (`sf`) object. 
#' @examples
#' ri <- spData::us_states |>
#' dplyr::select(NAME) |>
#'    dplyr::filter(NAME == 'Rhode Island') |>
#'    sf::st_transform(32617)
#'
#' sizeOptions <- TestGridSizes(ri)
#' head(sizeOptions)
#' @return A dataframe with testing results for each grid combination. 
#' A user needs to select the optimal grid size based on a tradeoff with minimizing variance, without creating too many grids which will need to be erased.
#' In the Rhode Island example I would use the 'Original' option which asks for 4 x grids and 7 y grids. 
#' @export
TestGridSizes <- function(target){
  
  bound <- sf::st_bbox(target)
  ratio <- (bound["xmax"] - bound["xmin"]) /
           (bound["ymax"] - bound["ymin"])

  if (ratio < 0.8) {
    x_start <- 4; y_start <- 7
  } else if (ratio < 0.9) {
    x_start <- 4; y_start <- 6
  } else if (ratio < 1.2) {
    x_start <- 5; y_start <- 5
  } else {
    x_start <- 6; y_start <- 4
  }

  offsets <- c(+2, +1, 0, -1, -2)
  names   <- c("Smallest", "Smaller", "Original", "Larger", "Largest")

  res <- lapply(offsets, function(o) {
    grid_variance(
      target,
      nx = x_start + o,
      ny = y_start + o
    )
  })

  data.frame(
    Name     = names,
    Grids    = vapply(res, `[[`, numeric(1), "n"),
    Variance = vapply(res, `[[`, numeric(1), "var"),
    GridNOx  = x_start + offsets,
    GridNOy  = y_start + offsets
  )
}

#' pass grid combinations through 
#' @keywords internal
#' @noRd
grid_variance <- function(target, nx, ny, top_n = 20) {
  gr <- sf::st_make_grid(target, n = c(nx, ny), square = FALSE)
  gr <- sf::st_intersection(gr, target)
  
  areas <- as.numeric(sf::st_area(gr))
  if (length(areas) == 0) {
    return(list(var = NA_real_, n = 0))
  }

  areas <- sort(areas / max(areas) * 100, decreasing = TRUE)
  areas <- head(areas, top_n)

  list(
    var = stats::var(areas, na.rm = TRUE),
    n   = length(gr)
  )
}