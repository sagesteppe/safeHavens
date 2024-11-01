#' Get an estimate for how many grids to draw over a species range
#' 
#' @description This function uses the dimensions of a species grid to estimate how many grids
#' would need to be added in the x and y directions to cover it with 20 grid cells
#' of roughly equal areas. 
#' @param target a species range as a simple feature (`sf`) object. 
#' @examples
#' RI <- spData::us_states |>
#' dplyr::select(NAME) |>
#'    dplyr::filter(NAME == 'Rhode Island') |>
#'    sf::st_transform(32617)
#'
#' sizeOptions <- TestGridSizes(RI)
#' head(sizeOptions)
#' @return A dataframe with testing results for each grid combination. 
#' A user needs to select the optimal grid size based on a tradeoff with minimizing variance, without creating too many grids which will need to be erased.
#' In the Rhode Island example I would use the 'Original' option which asks for 4 x grids and 7 y grids. 
#' @export
TestGridSizes <- function(target){
  
  bound <- sf::st_bbox(target)
  x_dist <- bound['xmax'] - bound['xmin']
  y_dist <- bound['ymax'] - bound['ymin']
  
  ratio <- x_dist/y_dist
  rm(bound, x_dist, y_dist)
  # values < 0.8 indicate y is considerable greater (longer) than x
  # values near 1 indicate a natural symmetry between x and y, both values can start at 5. 
  # values < 0.9 > 0.8 indicate x is greater (wider) than y
  # values > 1.4 indicate x is much greater (wider) than y
  
  if(ratio < 0.8){ # these areas are very long
    x_start = 4; y_start = 7} else if(ratio > 0.8 & ratio < 0.9) {
      x_start = 4; y_start = 6} else if(ratio > 0.9 & ratio < 1.2) { # equilateral
        x_start = 5; y_start = 5} else {
          x_start = 6; y_start = 4} 
  
  # intitial casting of recommended grid sizes. 
  gr <- sf::st_make_grid(target, n = c(x_start, y_start), square = FALSE)
  gr <- sf::st_intersection(gr, target) 
  areas <- as.numeric(sf::st_area(gr))
  areas <- sort(areas / max(areas) * 100, decreasing = TRUE)[1:20]
  var_Original <- var(areas,  na.rm = TRUE)
  
  # try with bigger grids
  gr_larger <- sf::st_make_grid(target, n = c(x_start-1, y_start-1), square = FALSE)
  gr_larger <- sf::st_intersection(gr_larger, target) 
  gr_larger_area <- as.numeric(sf::st_area(gr_larger))
  gr_larger_area <- sort(gr_larger_area / max(gr_larger_area) * 100, decreasing = TRUE)[1:20]
  var_larger <- var(gr_larger_area, na.rm = TRUE)
  
  # try with biggest grids
  gr_largest <- sf::st_make_grid(target, n = c(x_start-2, y_start-2), square = FALSE)
  gr_largest <- sf::st_intersection(gr_largest, target) 
  gr_largest_area <- as.numeric(sf::st_area(gr_largest))
  gr_largest_area <- sort(gr_largest_area / max(gr_largest_area) * 100, decreasing = TRUE)[1:20]
  var_largest <- var(gr_largest_area, na.rm = TRUE)
  
  # try with smaller grids
  gr_smaller <- sf::st_make_grid(target, n = c(x_start+1, y_start+1), square = FALSE)
  gr_smaller <- sf::st_intersection(gr_smaller, target) 
  gr_smaller_area <- as.numeric(sf::st_area(gr_smaller))
  gr_smaller_area <- sort(gr_smaller_area / max(gr_smaller_area) * 100, decreasing = TRUE)[1:20]
  var_smaller <- var(gr_smaller_area, na.rm = TRUE)
  
  # try with smallest grids
  gr_smallest <- sf::st_make_grid(target, n = c(x_start+2, y_start+2), square = FALSE)
  gr_smallest <- sf::st_intersection(gr_smallest, target) 
  gr_smallest_area <- as.numeric(sf::st_area(gr_smallest))
  gr_smallest_area <- sort(gr_smallest_area / max(gr_smallest_area) * 100, decreasing = TRUE)[1:20]
  var_smallest <- var(gr_smallest_area, na.rm = TRUE)
  
  # results for the user. 
  results <- data.frame(
    Name = c('Smallest', 'Smaller', 'Original',   'Larger', 'Largest'),
    Grids = c(
      length(gr_smallest), length(gr_smaller), 
      length(gr), length(gr_larger), length(gr_largest)), 
    Variance = c(var_smallest, var_smaller, var_Original, var_larger, var_largest),
    GridNOx = c(x_start+2, x_start+1, x_start, x_start-1, x_start-2), 
    GridNOy = c(y_start+2, y_start+1, y_start, y_start-1, y_start-2)
  )
  return(results)
}