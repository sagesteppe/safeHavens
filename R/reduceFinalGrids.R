#' More sliver fixing
#' 
#' in some instances tiny little grids are tagged along from processing. 
#' we will just give these to an arbitrary nearest feature. They would be 
#' annoying if by chance (... they are about a millionth of the area...)
#' an arbitrary random point was in them and missed, 
#' but the areas themselves tend to be remarkable inconsequential, and worth randomly
#' reassigning to a neighbor
#' @param final_grids truthfully nearly final at this point. 
reduceFinalGrids <- function(final_grids){
  
  grid_sz_order <- order(
    sf::st_area(final_grids)/sum(sf::st_area(final_grids)) * 100, decreasing = TRUE)
  grids2receive <- spdep::knearneigh(
    sf::st_point_on_surface(final_grids), k=1
  )[['nn']][grid_sz_order[21:length(grid_sz_order)]]
  
  final_grids[grid_sz_order[21:length(grid_sz_order)],'Assigned'] <- grids2receive 
  
  final_grids <- final_grids |>
    dplyr::group_by(Assigned) |>
    dplyr::summarise(geometry = sf::st_union(geometry))
  
  final_grids <- final_grids |>
    dplyr::ungroup() |>
    mutate(
      X = sf::st_coordinates(sf::st_point_on_surface(final_grids))[,1],
      Y = sf::st_coordinates(sf::st_point_on_surface(final_grids))[,2]
    ) |>
    dplyr::arrange(-Y, X) |>
    dplyr::mutate(Assigned = 1:dplyr::n()) |>
    dplyr::select(-X, -Y)
  
  return(final_grids)
}
