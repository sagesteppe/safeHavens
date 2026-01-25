#' Sample a species based on Isolation by Resistance Distance (IBR)
#' 
#' @description Create `n` seed collection areas based on the habitat, and geographic, resistance between points. 
#' @param base_raster a Raster surface to sample the points within, .... 
#' @param resistance_surface a resistance surface raster from  `build_resistance_surface` or `populationResistance`.
#' @param ibr_matrix Distance matrix, from `populationResistance` slot $Rmat
#' @param pts_sf sf/tibble/dataframe of point locations from `populationResistance` $pts_sf 
#' @param n Numeric. the number of clusters desired. 
#' @param fixedClusters Boolean. Defaults to TRUE, which will create n clusters. If False then use NbClust::NbClust to determine the optimal number of clusters.
#' @param min.nc Numeric. Minimum number of clusters to test if fixedClusters=FALSE, defaults to 5. 
#' @param max.nc Numeric. Maximum number of clusters to test if fixedClusters=FALSE, defaults to 20. 
#' @example 
#' 
IBRBasedSample <- function(
  base_raster,
  pop_raster,
  resistance_surface,
  pts_sf, 
  ibr_matrix, 
  fixedClusters = FALSE,
  n = NULL,
  min.nc = 5, 
  max.nc = 20,
  min_cluster_cells = 25,
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

  ## ---- Step 3: Resistance-based assignment ---------------------------------
  contested <- find_contested_cells(
    cluster_r,
    pop_rast = pop_raster
  )

  return(
    list(
      resistance_raster = resistance_surface,
      cluster_r         = contested$cluster_r,
      contested         = contested$contested,
      pts_sf            = clusts$clusters
      )
    )
  
  if (!is.null(resistance_surface)) {
    cluster_r <- assign_by_resistance(
      resistance_raster = resistance_surface,
      cluster_r         = contested$cluster_r,
      contested         = contested$contested,
      pts_sf            = clusts$clusters
    )
  }

  return(cluster_r)

}


#' @keywords internal
#' @tags noRd
cluster_connectivity <- function(
  x,
  pts_sf, 
  input = c("features", "distance"),
  fixedClusters = TRUE,
  n = NULL,
  min.nc = 5,
  max.nc = 10
) {

  input <- match.arg(input)

  # drop incomplete cases (rows only)
  x <- x[stats::complete.cases(x), , drop = FALSE]

  # --- CLUSTERING ---
  if (fixedClusters == TRUE) {

    if (is.null(n)) {
      stop("n must be provided when fixedClusters = TRUE")
    }

    d <- as.dist(x)
    hc <- stats::hclust(d, method = 'complete')
    pts_sf$ID <- stats::cutree(hc, n)

  } else {

    coords <- stats::cmdscale(x, k = 2)  # have to perform 2d embedding for nbclust method. 

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

  nbr_df <- data.frame(
    cell    = adj[, 1],
    nbr_val = terra::values(cluster_pop)[adj[, 2]]
  )

  nbr_df <- nbr_df[!is.na(nbr_df$nbr_val), , drop = FALSE]

  # Per-cell neighbor summary
  summary <- nbr_df |>
    dplyr::group_by(cell) |>
    dplyr::summarise(
      n_clusters = dplyr::n_distinct(nbr_val),
      cluster_id = dplyr::first(nbr_val),
      .groups = "drop"
    )

  # Safe cells: exactly one neighbor cluster
  safe <- summary$cell[summary$n_clusters == 1]
  safe_ids <- summary$cluster_id[summary$n_clusters == 1]

  # Contested cells
  contested_cells <- summary$cell[summary$n_clusters >= 2]

  # ---- Update cluster raster ----
  cluster_out <- cluster_r
  cluster_out[safe] <- safe_ids

  # ---- Contested raster ----
  contested <- cluster_out * NA
  contested[contested_cells] <- 1
  names(contested) <- "contested"

  list(
    cluster_r = cluster_out,
    contested = contested
  )
}

assign_by_resistance <- function(
  resistance_raster,
  cluster_r,
  contested,
  pts_sf
) {

  contested_cells <- which(!is.na(terra::values(contested)))
  if (!length(contested_cells)) return(cluster_r)

  # adjacency: contested cell -> neighboring clusters
  adj <- terra::adjacent(
    cluster_r,
    cells = contested_cells,
    directions = 8,
    pairs = TRUE
  )

  nbr_df <- data.frame(
    cell    = adj[, 1],
    cluster = cluster_r[adj[, 2]]
  )

  nbr_df <- nbr_df[!is.na(nbr_df$cluster), , drop = FALSE]
  if (!nrow(nbr_df)) return(cluster_r)

  candidates <- split(nbr_df$cell, nbr_df$cluster)

  seed_pts <- split(pts_sf, pts_sf$cluster_id)

  # explicit trackers
  best_dist    <- rep(Inf, length(contested_cells))
  best_cluster <- rep(NA_integer_, length(contested_cells))
  names(best_dist) <- names(best_cluster) <- contested_cells

  for (k in names(candidates)) {

    cells_k <- unique(candidates[[k]])
    seeds_k <- seed_pts[[k]]

    if (is.null(seeds_k) || !nrow(seeds_k)) next

    d <- terra::distance(
      resistance_raster,
      terra::vect(seeds_k)
    )

    vals <- terra::values(d)[cells_k]

    idx <- match(cells_k, contested_cells)

    update <- vals < best_dist[idx]
    best_dist[idx[update]]    <- vals[update]
    best_cluster[idx[update]] <- as.integer(k)
  }

  # only assign where we actually found a cluster
  ok <- !is.na(best_cluster)
  if (any(ok)) {
    cluster_r[as.integer(names(best_cluster)[ok])] <- best_cluster[ok]
  }

  cluster_r
}
