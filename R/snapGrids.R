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
    sf::st_set_agr('constant') |>
    sf::st_intersection(focal_grid) |>
    sf::st_set_agr('constant') 

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
  sf::st_agr(pts) <- "constant"

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
  sf::st_agr(slivers) <- "constant"
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

  sf::st_agr(all_pts_surface)  <- "constant"
  all_pts_surface <- sf::st_difference(
    all_pts_surface,
    sf::st_make_valid(sf::st_union(sf::st_combine(neighb_grid)))
  ) |>
    sf::st_make_valid()
  sf::st_agr(all_pts_surface) <- "constant"
  
  # ---- rejoin with neighbor grids ----
  # Extract geometry explicitly
  geom1 <- sf::st_geometry(neighb_grid)
  geom2 <- sf::st_geometry(all_pts_surface)

  # Bind attributes only
  attrs <- dplyr::bind_rows(
    sf::st_drop_geometry(neighb_grid) |> dplyr::select(Assigned = ID),
  sf::st_drop_geometry(all_pts_surface)
  )

  # Combine geometries safely
  geom <- c(geom1, geom2)

  # Rebuild sf object explicitly
  final_grids <- sf::st_sf(
    attrs,
    geometry = geom,
    crs = sf::st_crs(neighb_grid)
  ) |>
    dplyr::group_by(Assigned) |>
    dplyr::summarise(.groups = "drop") |>
    sf::st_make_valid()

  geom_types <- sf::st_geometry_type(final_grids)

  if (any(geom_types %in% c("GEOMETRYCOLLECTION", "MULTIPOLYGON"))) {
    final_grids <- sf::st_collection_extract(
      final_grids,
      "POLYGON",
      warn = FALSE
    )
  }

  final_grids <- final_grids[sf::st_is(final_grids, c("POLYGON", "MULTIPOLYGON")), ]
  final_grids <- nngeo::st_remove_holes(final_grids, max_area = 1000)
  final_grids <- healPolygons(final_grids)

  final_grids
}

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
