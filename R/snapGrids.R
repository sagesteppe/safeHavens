#' turn the point grid suggestions made by assignGrid_pts into polygons
#' 
#' This function is part of the grid based sampling process to turn small grid cells, 
#' which are to be broken up, into a larger existing grid cells. 
#' @param x output of `assignGrid_pts`
#' @param neighb_grid the neighboring grid options.
#' @param focal_grid the grid to reassign the area of. 
#' @keywords internal
snapGrids <- function(x, neighb_grid, focal_grid){
  
  # combine the output of `assignGrid_pts` together and create a small
  # polygon around there extent. 
  
  focal_grid <- sf::st_make_valid(focal_grid)
  neighb_grid <- sf::st_make_valid(neighb_grid)
  
  sf::st_agr(focal_grid) = "constant"
  sf::st_agr(x) = "constant"
  sf::st_agr(neighb_grid) = "constant"
  
  x <- x |>
    dplyr::group_by(Assigned) |> 
    dplyr::summarise(geometry = sf::st_combine(geometry)) |> 
    geos::geos_concave_hull(ratio = 0.01) |> 
    sf::st_difference() |> 
    sf::st_make_valid() |>
    sf::st_intersection(focal_grid)
  
  # now place points all across the grid to be divided. 
  pts <- sf::st_sample(focal_grid, 10000) |>
    sf::st_as_sf()
  sf::st_agr(pts) = "constant"
  #assign each of the newly generated points to the nearest polygon 
  # created above. The closest polygon will become the points
  # identity. 
  pts$Assigned <- sf::st_drop_geometry(
    x$Assigned)[ sf::st_nearest_feature(pts, x)] 
  pts <- pts |>
    dplyr::group_by(Assigned) |>
    dplyr::summarise(geometry = sf::st_union(x)) |>
    sf::st_convex_hull()
  
  # these points were cast to polygons, and now we remove overlapping areas
  all_pts_surface <- sf::st_difference(pts)
  sf::st_agr(all_pts_surface) = "constant"
  
  # we will determine snap and interval distances by drawing points which should land into the remaining slivers. This will remove any major
  # holes left by the above process. 
  slivers <- rmapshaper::ms_erase(focal_grid, all_pts_surface)
  sliver_pts <- sf::st_sample(slivers, size = 250)
  snap_dist <- ceiling(
    as.numeric(
      max(
        sf::st_distance(
          sliver_pts, 
          all_pts_surface[
            sf::st_nearest_feature(sliver_pts, all_pts_surface),], 
          by_element = TRUE)
      )
    )
  ) + 1
  
  # now fill in the gaps across the spatial data set. 
  all_pts_surface <- rmapshaper::ms_simplify(
      all_pts_surface,  keep_shapes = TRUE, snap = TRUE,
      keep = 1, weighting = 0, snap_interval = snap_dist) |>
    sf::st_make_valid() |>
    sf::st_difference() |>
    sf::st_buffer(snap_dist) |>
    sf::st_difference()
  
  sf::st_agr(all_pts_surface) <- 'constant'
  
  all_pts_surface <- sf::st_difference(
    all_pts_surface,  sf::st_make_valid(sf::st_union(sf::st_combine(neighb_grid)))
    ) |>
    sf::st_make_valid()
  
  # join the target grids onto this output. 
  final_grids <- neighb_grid |> 
    dplyr::select(Assigned = ID) |> 
    dplyr::bind_rows(all_pts_surface) |> 
    dplyr::group_by(Assigned) |> 
    dplyr::summarize(geometry = sf::st_union(geometry)) |> 
    sf::st_make_valid()
  
  # remove small internal holes which may arise from when the
  # created geometries were joined back to the original grid cells. 
  final_grids <- nngeo::st_remove_holes(final_grids, max_area = 1000)

  # now we determine which MULTIPOLYGONS are touching - if two polygons are contiguous, than we want to convert them to a single 
  # polygon. These exist because of very minor/short distance between
  # the existing polygons. 
  
  final_grids <- healPolygons(final_grids)
}
