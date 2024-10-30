#' We use spdep::knearneigh to ensure we obtain 4 of the nearest neighbors to a polygon
#' 
#' However, this system has limitations, if an island is ringed by land, it would work adequately.
#' But when the island is off a coast, it may return neighbors which are 'behind' another neighbor. 
#' We draw lines from the focal portion of the grid we will be combining to each neighbor, 
#' if the lines cross any other neighbor they are discarded. 
#' 
#' This is preferable to simply using the nearest feature (e.g. sf:;st_nearest_feature) when the chain of islands is elongated and may more appropriately be split across multiple grids downstream from here. 
#' @param from a point on surface of the grid which will be merged. POS will ensure it lands on a feature.
#' @param destinations the set of kearneigh which are accepting polygon merges 
first_neigh_directions <- function(from, destinations){
  
  dest_POS <- sf::st_point_on_surface(destinations)
  positions <- sf::st_union(from, dest_POS) |> 
    sf::st_cast('LINESTRING') |>
    sf::st_intersects(destinations) |>
    lengths() == 1
  
  return(positions)
}
