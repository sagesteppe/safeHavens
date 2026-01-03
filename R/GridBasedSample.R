#' Create hexagonal grid based polygons over a geographic range
#' 
#' @description This function creates 20 grid cells over a geographic area (`x`), typically a species range.  
#' @param x An SF object or terra spatraster. the range over which to generate the clusters.
#' @param planar_projection Numeric, or character vector. An EPSG code, or a proj4 string, for a planar coordinate projection, in meters, for use with the function. For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise a continental scale projection like 5070 See https://projectionwizard.org/ for more information on CRS. The value is simply passed to sf::st_transform if you need to experiment. 
#' @param gridDimensions A single row form the ouput of `TestGridSizes` with the optimal number of grids to generate. 
#' @examples \dontrun{ # not ran to bypass CRAN check time limits. ~6 seconds to treat Rhode Island. 
#' ri <- spData::us_states |> 
#' dplyr::filter(NAME == 'Rhode Island') |>
#'   sf::st_transform(32615)
#'
#' sizeOptions <- TestGridSizes(ri)
#' head(sizeOptions) # in this case let's shoot for 33 and see what happens
#' sizeOptions <- sizeOptions[sizeOptions$Name == 'Original',]
#' 
#' output <- GridBasedSample(ri, 5070, gridDimensions = sizeOptions)
#' plot(output$geometry)
#' }
#' @return An simple features (sf) object containing the final grids for saving to computer. See the vignette for questions about saving the two main types of spatial data models (vector - used here, and raster). 
#' @export
GridBasedSample <- function(x, planar_projection, gridDimensions){
  
  if(missing(planar_projection)){planar_projection = 5070}
  if(sf::st_crs(x) != planar_projection){x <-  sf::st_transform(x, planar_projection)}
  
  # Determine the size of each grid. 
  gr <- sf::st_make_grid(x, n = c(gridDimensions$GridNOx, gridDimensions$GridNOy), square = FALSE) 
  
  gr <- sf::st_intersection(gr, x) 
  gr <- sf::st_collection_extract(gr, 'POLYGON')
  
  grid_areas <- sf::st_as_sf(gr) |> 
    sf::st_make_valid() |>
    dplyr::mutate(
      ID   = seq_len(dplyr::n()),
      Area = as.numeric(sf::st_area(gr))
    )
  
  # order polygons by size, all polygons > 20 will be merged with a neighboring polygon
  indices <- grid_areas$Area >= sort(grid_areas$Area, decreasing = TRUE)[20]
  to_merge_sf <- gr[!indices,]
  merge_into_sf <- gr[indices,] |> sf::st_as_sf()
  # before calculating neighbors, we will union, adjacent polygons which we will 
  # end up merging to the polygons we will keep. This SHOULD allow for better 
  # distribution of there areas into the remaining polygons, making the kept polygons more equal in size. 
  
  gr <- to_merge_sf |> 
    sf::st_union() |> 
    sf::st_cast('POLYGON') |> 
    sf::st_as_sf()  %>% # gotta use pipe to set position of tibbles 
    dplyr::bind_rows(
      gr[indices,] |> sf::st_as_sf(), .
    )
  
  # Determine neighboring polygons
  sf::st_agr(gr) <- 'constant'
  neighbors <- spdep::poly2nb(gr, queen = FALSE)[21:nrow(gr)]
  
  full_sized_neighbors <- which( # consider these to be full sized grids
    grid_areas$Area[1:20] / max(grid_areas$Area) >= 0.975) 
  
  # identify neighboring polygons
  to_merge_sf <- gr[21:nrow(gr),]
  merge_into_sf <- gr[1:20,] 
  
  # if their are no neighbors, this implies that the focal grid area is isolated, e.g. an island
  # we will use distance based neighbors for that focal area.  
  for (i in seq_along(neighbors)){
    if(sum(neighbors[[i]])==0){
      need_distance_neighbs <- which(sum(neighbors[[i]])==0)
    #  sf::st_agr(need_distance_neighbs) <- 'constant'
      gr_pos <- sf::st_point_on_surface(gr)
      ob <- spdep::knearneigh(gr_pos, k=4)[['nn']]
      neighbors[[i]] <- ob[20+need_distance_neighbs,]
      neighbors[[i]] <- neighbors[[i]][neighbors[[i]]<=20] # remove neighbors which will be dropped. 
      
      from <- sf::st_point_on_surface(to_merge_sf[need_distance_neighbs,])
      destinations <- merge_into_sf[neighbors[[i]],]
      neighbors[[i]] <- neighbors[[i]][first_neigh_directions(from, destinations = destinations)]
      
    }
  }
  
  area2be_reassigned <- vector(length = length(neighbors))
  areas <- vector(mode = 'list', length = length(neighbors))
  prop_areas <- vector(mode = 'list', length = length(neighbors))
  prop_donor <- numeric(length = length(neighbors))
  for (i in seq_along(neighbors)){
    
    area2be_reassigned[i] <- sf::st_area(to_merge_sf[i,])
    areas[[i]] <- as.numeric(sf::st_area(gr[neighbors[[i]],]))
    
    prop_donor[i] <- area2be_reassigned[i] / (sum(areas[[i]]) + area2be_reassigned[i]) 
    prop_areas[[i]] <- areas[[i]] / (sum(areas[[i]]) + area2be_reassigned[i]) 
  }
  
  area2be_reassigned <- vector(length = length(neighbors))
  areas <- vector(mode = 'list', length = length(neighbors))
  prop_areas <- vector(mode = 'list', length = length(neighbors))
  prop_donor <- numeric(length = length(neighbors))
  
  for (i in seq_along(neighbors)){
    
    area2be_reassigned[i] <- sf::st_area(to_merge_sf[i,])
    areas[[i]] <- as.numeric(sf::st_area(gr[neighbors[[i]],]))
    
    prop_donor[i] <- area2be_reassigned[i] / (sum(areas[[i]]) + area2be_reassigned[i]) 
    prop_areas[[i]] <- areas[[i]] / (sum(areas[[i]]) + area2be_reassigned[i]) 
  }
  
  rm(area2be_reassigned)
  
  prop_target <- vector(mode = 'list', length = length(prop_areas))
  area_sort <- vector(mode = 'list', length = length(prop_areas))
  area_des <- vector(mode = 'list', length = length(prop_areas))
  nf_pct <- vector(mode = 'list', length = length(prop_areas))
  props <- vector(mode = 'list', length = length(prop_areas))
  # Using the polygons which will be merged, try to make the following polygons
  # as equally sized as possible - without ever removing area from an existing grid. 
  for (i in seq_along(area_sort)){
    
    area_des <- (sum(prop_areas[[i]]) + prop_donor[i]) / length(prop_areas[[i]])
    
    if(all(prop_areas[[i]] < area_des)==TRUE){
      
      # these polygons will all be the same size!... roughly... 
      prop_target[[i]] <- rep(area_des, times = length(prop_areas[[i]]))
      
    } else if(any(prop_areas[[i]] > area_des)==TRUE){
      
      prop_target[[i]] <- numeric(length(prop_areas[[i]]))
      kp <- prop_areas[[i]] < area_des # determine which grids are to big
      area_des <- (sum(prop_areas[[i]][kp]) + prop_donor[i]) / 
        length(prop_areas[[i]][kp])
      
      # make grids smaller than the goal threshold size the threshold, 
      # return grids larger than the threshold size as they are. 
      prop_target[[i]][kp] <- area_des
      prop_target[[i]][!kp] <-prop_areas[[i]][!kp]
    }
    
    nf_pct[[i]] <- stats::setNames( # the existing cover for each grid. 
      prop_areas[[i]] * 100, 
      neighbors[[i]]
    )
    
    props[[i]] <- stats::setNames( # the desired cover for each grid 
      prop_target[[i]] * 100, 
      neighbors[[i]]
    )
  }

  rm(area_des, area_sort, i, prop_donor, prop_target, areas)
  
  gr <- dplyr::mutate(gr, ID = seq_len(dplyr::n()),  .before = x) |>
    dplyr::rename(geometry = x)
  neighb_grid <- vector(mode = 'list', length = length(prop_areas))
  for (i in seq_along(neighbors)){
    neighb_grid[[i]] <- gr[neighbors[[i]], ]
  }

  # place points throughout the grids which need to be merged to determine
  # how they will be reallocated into larger grids. 
  to_merge_sf <- dplyr::rename(to_merge_sf, geometry = x)
  to_merge_sf <- split(to_merge_sf, f = seq_len(nrow(to_merge_sf)))
  
  out <- vector(mode = 'list', length = length(prop_areas))
  for (i in seq_along(out)){
    out[[i]] <- assignGrid_pts(
      neighb_grid =  neighb_grid[[i]], 
      focal_grid = to_merge_sf[[i]], 
      props = props[[i]], 
      nf_pct = nf_pct[[i]]
    )
  }
  
  # finally create polygons from the point samples
  final <- vector(mode = 'list', length = length(out))
  for (i in seq_along(final)){
    final[[i]] <- snapGrids(
      x = out[[i]],
      neighb_grid =  neighb_grid[[i]], 
      focal_grid = to_merge_sf[[i]]
    )
  } 
  final_grids <- dplyr::bind_rows(final)
  
  groups <- split(final_grids, f = final_grids$Assigned)
  final_grids <- lapply(groups, snapR) |> dplyr::bind_rows()
  final_grids <- sf::st_make_valid(final_grids)
  sf::st_agr(final_grids) <- 'constant'
  
  # reconstitute all original input grids, i.e. those without neighbors, 
  # with all reassigned grids. 
  gr2 <- sf::st_difference(
    gr,  sf::st_make_valid(sf::st_union(sf::st_combine(final_grids)))
  ) |>
    sf::st_make_valid()
  gr2 <- gr2[as.numeric(sf::st_area(gr2))/1000 > 10,]
  final_grids <- gr2 |>
    dplyr::rename(Assigned = ID) |>
    dplyr::bind_rows(final_grids)  |>
    sf::st_as_sf()
  
  # now number the grids in a uniform fashion
  cents <- sf::st_point_on_surface(final_grids)
  sf::st_agr(cents) = "constant"
  cents <- cents |>
    dplyr::mutate(
      X = sf::st_coordinates(cents)[,1],
      Y = sf::st_coordinates(cents)[,2]
    ) |>
    dplyr::arrange(-Y, X) |>
    dplyr::mutate(NEWID = seq_len(dplyr::n())) |>
    sf::st_drop_geometry() |>
    dplyr::select(Assigned, NEWID)
  
  final_grids <- dplyr::left_join(final_grids, cents, by = 'Assigned') |>
    dplyr::select(-Assigned, Assigned = NEWID) 
  
  if(max(final_grids$Assigned)>20){
    final_grids <- reduceFinalGrids(final_grids)
    if(max(final_grids$Assigned)>20){
      final_grids <- reduceFinalGrids(final_grids)}
  }
  final_grids <- sf::st_as_sf(final_grids) |>
    dplyr::rename(dplyr::any_of( c(ID = 'Assigned'))) |>
    dplyr::arrange(ID) |>
    dplyr::relocate(ID)
  
  list(Geometry = final_grids)
}
