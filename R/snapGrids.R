#' turn the point grid suggestions made by assignGrid_pts into polygons
#' 
#' This function is part of the grid based sampling process to turn small grid cells, 
#' which are to be broken up, into a larger existing grid cells. 
#' @param x output of `assignGrid_pts`
#' @param neighb_grid the neighboring grid options.
#' @param focal_grid the grid to reassign the area of. 
#' @keywords internal
snapGrids <- function(x, neighb_grid, focal_grid){

  focal_grid  <- sf::st_make_valid(focal_grid)
  neighb_grid <- sf::st_make_valid(neighb_grid)

  sf::st_agr(focal_grid)  <- "constant"
  sf::st_agr(neighb_grid) <- "constant"
  sf::st_agr(x)           <- "constant"

  # ---- build assignment polygons ----
  x <- x |>
    dplyr::group_by(Assigned) |> 
    dplyr::summarise(geometry = sf::st_combine(geometry), .groups = "drop") |> 
    sf::st_concave_hull(ratio = 0.01) |> 
    sf::st_make_valid() |>
    sf::st_intersection(focal_grid)

  # ---- sample points across focal grid ----
  pts <- sf::st_sample(focal_grid, 10000) |> sf::st_as_sf()
  sf::st_agr(pts) <- "constant"

  pts$Assigned <- x$Assigned[
    sf::st_nearest_feature(pts, x)
  ]

  pts <- pts <- pts |>
  dplyr::group_by(Assigned) |>
  dplyr::summarise(.groups = "drop") |>
  sf::st_convex_hull() |>
    sf::st_make_valid()

  all_pts_surface <- sf::st_as_sf(
    sf::st_cast(
      sf::st_make_valid(
        sf::st_union(
          sf::st_combine(pts)
        )
      ),
      "POLYGON",
      warn = FALSE
    )
  )
  sf::st_agr(all_pts_surface) <- "constant"

  slivers <- rmapshaper::ms_erase(focal_grid, all_pts_surface)
  sliver_pts <- sf::st_sample(slivers, size = 250)

  snap_dist <- as.numeric(
    ceiling(
      max(
        sf::st_distance(
          sliver_pts,
          all_pts_surface[ sf::st_nearest_feature(sliver_pts, all_pts_surface), ],
          by_element = TRUE
        )
      )
    )
  ) + 1

  all_pts_surface <- rmapshaper::ms_simplify(
    all_pts_surface,
    keep_shapes = TRUE,
    snap = TRUE,
    keep = 1,
    weighting = 0,
    snap_interval = snap_dist
  ) |>
    sf::st_make_valid() |>
    sf::st_buffer(snap_dist) |>
    sf::st_make_valid()

  all_pts_surface <- sf::st_difference(
    all_pts_surface,
    sf::st_make_valid(sf::st_union(sf::st_combine(neighb_grid)))
  ) |>
    sf::st_make_valid()

  # ---- rejoin with neighbor grids ----
  final_grids <- final_grids <- neighb_grid |> 
    dplyr::select(Assigned = ID) |> 
    dplyr::bind_rows(all_pts_surface) |> 
    dplyr::group_by(Assigned) |> 
    dplyr::summarise(.groups = "drop") |> 
    sf::st_make_valid()

  final_grids <- sf::st_collection_extract(final_grids, "POLYGON")
  final_grids <- final_grids[sf::st_is(final_grids, c("POLYGON", "MULTIPOLYGON")), ]
  final_grids <- nngeo::st_remove_holes(final_grids, max_area = 1000)
  final_grids <- healPolygons(final_grids)

  final_grids
}
