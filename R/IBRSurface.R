#' Sample a species based on Isolation by Resistance Distance (IBR)
#' 
#' @description Create `n` seed collection areas based on the habitat, and geographic, resistance between points. 
#' @param base_raster a Raster surface to sample the points within, .... 
#' @param pop_raster Spatiall buffered population raster from `populationResistance`. 
#' @param resistance_surface a resistance surface raster from  `build_resistance_surface` or `populationResistance`.
#' @param ibr_matrix Distance matrix, from `populationResistance` slot $Rmat
#' @param pts_sf sf/tibble/dataframe of point locations from `populationResistance` $pts_sf 
#' @param n Numeric. the number of clusters desired. 
#' @param planar_proj For numbering clusters in return object. 
#' @param fixedClusters Boolean. Defaults to FALSE, which will create n clusters. If False then use NbClust::NbClust to determine the optimal number of clusters. TRUE not yet suported. 
#' @param min.nc Numeric. Minimum number of clusters to test if fixedClusters=FALSE, defaults to 5. 
#' @param max.nc Numeric. Maximum number of clusters to test if fixedClusters=FALSE, defaults to 20. 
#' @param geo_iters Numeric. Maximum number of iterations to grow hulls by. 
#' @param boundary_width Numeric. Maximum boundary width before switching to final assignment method. 
#' @param max_geo_cells Numeric. Maximum number of cells to grow cheap geographic front growth of cluster centers by. 
#' @param distance_method Great circle distance calculation method passed onto ot terra. One of 'haversine' (default), or 'cosine'.
#' @examples  # see package vignette
#' @export
IBRSurface <- function(
  base_raster,
  pop_raster,
  resistance_surface,
  pts_sf, 
  ibr_matrix, 
  fixedClusters = FALSE,
  n = NULL,
  min.nc = 5, 
  max.nc = 20,
  planar_proj, 
  max_geo_cells = NULL,
  geo_iters = 10,
  boundary_width = 2,
  distance_method = c("haversine", "cosine")
) {

  distance_method <- match.arg(distance_method)

  ## ---- Setup ---------------------------------------------------------------
  if (!terra::is.lonlat(base_raster)) {
    stop("base_raster must be in lon/lat (EPSG:4326) for geodesic distances")
  }

  pts_sf <- sf::st_transform(pts_sf, terra::crs(base_raster))

  if (is.null(max_geo_cells)) {
    max_geo_cells <- ceiling(0.7 * terra::ncell(base_raster))
  }

  ## ---- Step 0: Cluster points into groups ----------------------------------
  clusts <- cluster_connectivity(
    x = ibr_matrix,
    pts_sf = pts_sf,
    input = c("features", "distance"),
    fixedClusters = fixedClusters,
    n = n,
    min.nc = min.nc,
    max.nc = max.nc
  )

  ## ---- Step 1: Conservative geographic cores -------------------------------

  ## for inner buffer based re-assignment find the shortest distance between two points in different groups. 
  # then use 40% of this distance for the cap, this ensures nothing can be misclassified. 
  nn <- sf::st_nearest_feature(clusts$clusters)
  nn_x_clusters <- which(clusts$clusters[['ID']] != clusts$clusters[['ID']][nn])
  max_geo_dist <- as.numeric(
    min(
    sf::st_distance(
      clusts$clusters[nn_x_clusters,], 
      clusts$clusters[ nn[nn_x_clusters], ], 
      by_element = TRUE)
      )
  ) * 0.45

  cluster_r <- geographic_core_assignment(
    pop_raster = pop_raster, 
    pts_sf  = clusts$clusters,
    max_dist = max_geo_dist,
    distance_method = distance_method
  )

  ## ---- Step 2: Cheap geographic expansion ----------------------------------
  cluster_r <- expand_geographic_front(
    cluster_r,
    max_cells_per_cluster = max_geo_cells,
    max_iters = geo_iters
  )

  ## ---- Step 3: Fill contested cells from neigbors ---------------------------------
  contested <- find_contested_cells(
    cluster_r,
    pop_rast = pop_raster
  )

  cluster_r <- assign_contested_line(
    cluster_r = contested$cluster_r, 
    contested = contested$contested
  )

  ## randomly fill remaining cells from a neighbor. 
  cluster_r <- finalize_cluster(
    cluster_r = cluster_r,
    pop_raster = pop_raster
  )

  ## format data for return 

  spatialClusters <- terra::as.polygons(cluster_r) |>
    sf::st_as_sf()
  
  # now number the grids in a uniform fashion
  if(!missing(planar_proj)){
    spatialClusters <- sf::st_transform(spatialClusters, planar_proj)
  }
  
  sf::st_agr(spatialClusters) = "constant"
  cents <- sf::st_point_on_surface(spatialClusters)

  cents <- sf::st_transform(cents, sf::st_crs(planar_proj))
  spatialClusters <- sf::st_transform(spatialClusters, sf::st_crs(planar_proj))

    cents <- cents |>
    dplyr::mutate(
      X = sf::st_coordinates(cents)[,1],
      Y = sf::st_coordinates(cents)[,2]
    ) |>
    dplyr::arrange(-Y, X) |>
    dplyr::mutate(ID = seq_len(dplyr::n())) |>
    dplyr::arrange(ID) |>
    dplyr::select(ID, geometry)
  
  sf::st_agr(cents) = "constant"
  ints <- unlist(sf::st_intersects(spatialClusters, cents))
  spatialClusters <- spatialClusters |>
    dplyr::mutate(ID = ints, .before = 1) |>
    dplyr::arrange(ID)

  list(
    points = clusts$clusters,
    geometry = spatialClusters
  )

}


#' @keywords internal
#' @noRd
cluster_connectivity <- function(
  x,
  pts_sf, 
  input = c("features", "distance"),
  fixedClusters = TRUE,
  n = NULL,
  min.nc,
  max.nc
) {

  input <- match.arg(input)

  # drop incomplete cases (rows only)
  x <- x[stats::complete.cases(x), , drop = FALSE]

  # --- CLUSTERING ---
  if (fixedClusters == TRUE) {

    if (is.null(n)) {
      stop("n must be provided when fixedClusters = TRUE")
    }

    d <- stats::as.dist(x)
    hc <- stats::hclust(d, method = 'complete')
    pts_sf$ID <- stats::cutree(hc, n)

  } else {
    coords <- stats::cmdscale(x, k = 2) # have to perform 2d embedding for nbclust method.
    
    NoClusters <- NbClust::NbClust(
      data = coords,
      distance = 'euclidean',
      min.nc = min.nc,
      max.nc = max.nc,
      method = 'complete',
      index = 'silhouette'
    )

    pts_sf$ID <- NoClusters$Best.partition
    hc <- NULL
  }

  list(
    clusters = pts_sf,
    hclust = hc,
    input = input
  )
}

#' @keywords internal
#' @noRd
geographic_core_assignment <- function(
  pop_raster,          # raster of buffered population areas
  pts_sf,              # sf/tibble of points with $ID already assigned
  max_dist,            # maximum distance to assign cells
  distance_method = c("haversine", "cosine")  # geodesic distance
) {
  distance_method <- match.arg(distance_method)

  # --- Setup seed raster ---
  seed_r <- terra::rasterize(
    terra::vect(pts_sf),
    pop_raster,
    field = "ID"
  )

  # --- Compute distances from all pop_raster cells to nearest seed ---
  d <- terra::distance(seed_r, method = distance_method)

  # --- Mask pop_raster to only cells within max_dist ---
  pop_masked <- terra::mask(pop_raster, d <= max_dist)
  pop_masked_cells <- which(!is.na(terra::values(pop_masked)))

  # --- Get centroids of masked raster cells as SpatVector points ---
  cell_vect <- terra::as.points(pop_masked)

  # --- Use terra::nearest() to assign each cell to nearest seed ---
  nearest_seed_idx <- terra::nearest(cell_vect, terra::vect(pts_sf))

  assigned_dists <- terra::distance(
    cell_vect, 
    terra::vect(pts_sf)[nearest_seed_idx$to_id, ],
    pairwise = TRUE
  )

  # Only keep assignments within max_dist
  valid_assignments <- assigned_dists <= max_dist
  cluster_ids <- pts_sf$ID[nearest_seed_idx$to_id]
  cluster_ids[!valid_assignments] <- NA

  # --- Fill masked raster with cluster IDs ---
  cluster_r <- pop_masked
  cluster_r[pop_masked_cells] <- cluster_ids
  names(cluster_r) <- "ID"

  cluster_r
}

#' @keywords internal
#' @noRd
expand_geographic_front <- function(
  cluster_r,
  max_cells_per_cluster,
  max_iters
) {

  for (i in seq_len(max_iters)) {

    adj <- terra::adjacent(
      cluster_r,
      cells = which(!is.na(terra::values(cluster_r))),
      directions = 8,
      pairs = TRUE
    )

    candidates <- adj[is.na(cluster_r[adj[,2]]), , drop = FALSE]
    if (!nrow(candidates)) break

    # Get cluster values for source cells
    cluster_vals <- cluster_r[candidates[,1]]
    
    # Filter out any NA cluster values before creating df
    valid_idx <- !is.na(cluster_vals)
    if (!any(valid_idx)) break
    
    df <- data.frame(
      cell    = candidates[valid_idx, 2],
      cluster = cluster_vals[valid_idx]
    )

    # Now the grouped filter should work
    df <- df |>
      dplyr::group_by(cell) |>
      dplyr::filter(dplyr::n_distinct(.data$cluster) == 1) |>
      dplyr::ungroup()

    if (!nrow(df)) break

    sizes <- table(terra::values(cluster_r))
    df <- df[sizes[as.character(df$cluster)] < max_cells_per_cluster, ]

    if (!nrow(df)) break

    cluster_r[df$cell] <- df$cluster
  }

  cluster_r
}

#' @keywords internal
#' @noRd
find_contested_cells <- function(
  cluster_r,
  pop_rast,
  directions = 8
) {

  cluster_pop <- terra::mask(cluster_r, pop_rast)

  na_cells <- which(
    is.na(terra::values(cluster_pop)) &
    !is.na(terra::values(pop_rast))
  )

  if (!length(na_cells)) {
    return(list(
      cluster_r = cluster_r,
      contested = cluster_r * NA
    ))
  }

  adj <- terra::adjacent(
    cluster_pop,
    cells = na_cells,
    directions = directions,
    pairs = TRUE
  )

  adj_df <- data.frame(
    cell = adj[, 1],
    nbr  = adj[, 2],
    to_cluster = terra::values(cluster_r)[adj[, 2]]
  )

  # drop NA neighbors
  adj_df <- adj_df[!is.na(adj_df$to_cluster), , drop = FALSE]

  # ---- Per-cell neighbor summary ----
  nbr_summary <- adj_df |>
    dplyr::group_by(cell) |>
    dplyr::summarise(
      clusters   = list(unique(to_cluster)),
      n_clusters = dplyr::n_distinct(to_cluster),
      .groups = "drop"
    )

  # ---- Safe assignments ----
  safe <- nbr_summary[nbr_summary$n_clusters == 1, ]
  safe_cells <- safe$cell
  safe_ids   <- vapply(safe$clusters, `[`, numeric(1), 1)

  cluster_out <- cluster_r
  cluster_out[safe_cells] <- safe_ids

  # ---- Contested cells ----
  contested_cells <- nbr_summary$cell[nbr_summary$n_clusters >= 2]

  contested <- cluster_out * NA
  contested[contested_cells] <- 1
  names(contested) <- "contested"

  list(
    cluster_r = cluster_out,
    contested = contested,
    candidates = split(
      nbr_summary$clusters,
      nbr_summary$cell
    )
  )

}

#' @keywords internal
#' @noRd
assign_contested_line <- function(cluster_r, contested) {
  # Identify contested cell indices
  contested_cells <- which(!is.na(terra::values(contested)))
  if (!length(contested_cells)) return(cluster_r)
  
  # Get adjacency graph for contested cells (8 directions)
  adj <- terra::adjacent(cluster_r, cells = contested_cells, directions = 8, pairs = FALSE)
  
  # Build edge list, filtering out NA adjacencies to avoid igraph warnings
  edge_list <- do.call(rbind, lapply(seq_along(adj), function(i) {
    valid_adj <- adj[[i]][!is.na(adj[[i]])]
    if (length(valid_adj) > 0) {
      cbind(contested_cells[i], valid_adj)
    } else {
      NULL
    }
  }))
  
  # If no valid edges, return as-is
  if (is.null(edge_list) || nrow(edge_list) == 0) return(cluster_r)
  
  # Remove any remaining NA rows (belt and suspenders)
  edge_list <- edge_list[complete.cases(edge_list), , drop = FALSE]

  # Check again after NA removal
  if (nrow(edge_list) == 0) return(cluster_r)

  # Find connected components of contested cells
  g <- igraph::graph_from_data_frame(d = edge_list, directed = FALSE)
  comps <- igraph::components(g)$membership
  
  # Split each contiguous component in half
  cluster_out <- cluster_r
  for (comp_id in unique(comps)) {
    cells <- contested_cells[comps == comp_id]
    n <- length(cells)
    half <- ceiling(n/2)
    
    # Get neighboring cluster ids
    nbr_vals <- terra::values(cluster_r)[terra::adjacent(cluster_r, cells = cells, directions = 8)]
    nbr_vals <- nbr_vals[!is.na(nbr_vals)]
    if (length(nbr_vals) < 2) next  # cannot split, just leave assigned later
    
    # Assign first half to one neighbor, second half to the other
    cluster_out[cells[1:half]] <- nbr_vals[1]
    cluster_out[cells[(half+1):n]] <- nbr_vals[2]
  }
  
  cluster_out
}

#' @keywords internal
#' @noRd
finalize_cluster <- function(cluster_r, pop_raster, directions = 8) {
  
  # --- Mask to population raster ---
  cluster_pop <- terra::mask(cluster_r, pop_raster)
  
  # Identify unassigned cells
  na_cells <- which(is.na(terra::values(cluster_pop)) & !is.na(terra::values(pop_raster)))
  
  if (!length(na_cells)) return(cluster_pop)
  
  # --- Assign nearest cluster for each NA cell ---
  for (cell in na_cells) {
    
    nbr_cells <- terra::adjacent(cluster_pop, cells = cell, directions = directions, pairs = FALSE)[[1]]
    
    # Neighbor cluster IDs, drop NA
    nbr_vals <- stats::na.omit(terra::values(cluster_pop)[nbr_cells])
    
    if (length(nbr_vals)) {
      # Randomly pick one if multiple neighbors
      cluster_pop[cell] <- sample(nbr_vals, 1)
    } else {
      # If truly isolated, leave NA (should be rare)
      cluster_pop[cell] <- NA
    }
  }
  
  cluster_pop
}
