#' More sliver fixing
#' 
#' in some instances tiny little grids are tagged along from processing. 
#' we will just give these to an arbitrary nearest feature. They would be 
#' annoying if by chance (... they are about a millionth of the area...)
#' an arbitrary random point was in them and missed, 
#' but the areas themselves tend to be remarkable inconsequential, and worth randomly
#' reassigning to a neighbor
#' @param final_grids truthfully nearly final at this point. 
#' @keywords internal
#' @noRd
reduceFinalGrids <- function(final_grids){
  
  n_grids <- nrow(final_grids)
  
  # Only reduce if we have more than 20 grids
  if (n_grids <= 20) {
    sf::st_agr(final_grids) <- 'constant'
    # Just renumber if 20 or fewer
    final_grids <- final_grids |>
      dplyr::ungroup() |>
      dplyr::mutate(
        X = sf::st_coordinates(sf::st_point_on_surface(final_grids))[,1],
        Y = sf::st_coordinates(sf::st_point_on_surface(final_grids))[,2]
      ) |>
      dplyr::arrange(-Y, X) |>
      dplyr::mutate(Assigned = seq_len(dplyr::n())) |>
      dplyr::select(-X, -Y)
    return(final_grids)
  }
  
  # Calculate grid sizes and order
  grid_sz_order <- order(
    sf::st_area(final_grids) / sum(sf::st_area(final_grids)) * 100, 
    decreasing = TRUE
  )
  
  # Keep the 20 largest grids
  keep_indices <- grid_sz_order[1:20]
  small_grid_indices <- grid_sz_order[21:n_grids]
  
  # Find nearest neighbor among the KEPT grids only
  keep_grids <- final_grids[keep_indices, ]
  small_grids <- final_grids[small_grid_indices, ]
  
  sf::st_agr(keep_grids) <- 'constant'
  sf::st_agr(small_grids) <- 'constant'

  # Find which of the kept grids each small grid is nearest to
  nearest_kept <- spdep::knearneigh(
    sf::st_point_on_surface(small_grids),
    k = 1 )[['nn']]
  
  # Map these back to the original Assigned values of kept grids
  final_grids$NewAssigned <- final_grids$Assigned
  final_grids$NewAssigned[small_grid_indices] <- keep_grids$Assigned[nearest_kept]
  
  # Merge grids by NewAssigned value
  final_grids <- final_grids |>
    dplyr::group_by(NewAssigned) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = 'drop') |>
    dplyr::rename(Assigned = NewAssigned)
  
  sf::st_agr(final_grids) <- 'constant'
  # Renumber based on spatial position
  final_grids <- final_grids |>
    dplyr::mutate(
      X = sf::st_coordinates(sf::st_point_on_surface(final_grids))[,1],
      Y = sf::st_coordinates(sf::st_point_on_surface(final_grids))[,2]
    ) |>
    dplyr::arrange(-Y, X) |>
    dplyr::mutate(Assigned = seq_len(dplyr::n())) |>
    dplyr::select(-X, -Y)
  
  final_grids
}
