#' Clean up unioned geometries - part 1
#' 
#' this function uses sf::st_snap to remove small lines and other artifacts associated
#' with the unioning of polygons. This is ran within `snapGrids`
#' @param x most of the output of `snapgrids`
#' @keywords internal
healPolygons <- function(x){
  
  healR <- function(x){
    Assigned <- x$Assigned
    sf::st_agr(x) = "constant" 
    
    x <- x |> sf::st_buffer(0.0001) |> 
      sf::st_union() |> 
      sf::st_combine() |>
      sf::st_as_sf() |> 
      dplyr::mutate(Assigned = Assigned) |> 
      dplyr::rename(geometry = x) 
    x
  }
  
  rows <- split(x, f = seq_len(nrow(x)))
  rows <- lapply(rows, healR)
  rows <- dplyr::bind_rows(rows)
  
} 

#' Clean up unioned geometries - part 2
#' 
#' this function uses sf::st_snap to remove small lines and other artifacts associated
#' with the unioning of polygons
#' @param x the output of healPolygons
#' @keywords internal
#' @noRd
snapR <- function(x){
  sf::st_agr(x) = 'constant'
  Assigned <- sf::st_drop_geometry(x$Assigned)[1]
  x <- sf::st_snap(x = x, y = x, tolerance = 0.0001)|>
    sf::st_union() |>
    sf::st_as_sf() |> 
    dplyr::mutate(Assigned = Assigned) |> 
    dplyr::rename(geometry = x) 
  
  x
}
