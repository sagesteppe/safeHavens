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
  resistance_surface,
  pts_sf, 
  ibr_matrix, 
  fixedClusters = FALSE,
  n = NULL,
  min.nc = 5, 
  max.nc = 20,
  min_cluster_cells = 25,
  max_geo_dist = NULL,
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

  if (is.null(max_geo_dist)) {
    max_geo_dist <- 5 * terra::res(base_raster)[1]
  }

  if (is.null(max_geo_cells)) {
    max_geo_cells <- ceiling(0.7 * terra::ncell(base_raster))
  }

  ## ---- Step 0: Cluster points into groups ----------------------------------
  clusts <- cluster_connectivity(
    input = c("features", "distance"),
    fixedClusters = fixedClusters,
    n = n,
    min.nc = min.nc,
    max.nc = max.nc,
    linkage = NULL
  )

  ## ---- Step 1: Conservative geographic cores -------------------------------
  cluster_r <- geographic_core_assignment(
    seeds_sf        = pts_sf,
    template_raster = base_raster,
    max_dist        = max_geo_dist,
    distance_method = distance_method
  )

  cluster_r <- enforce_min_core_area(
    cluster_r,
    min_cells = min_cluster_cells
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
    width = boundary_width
  )

  if (!is.null(resistance_raster)) {
    cluster_r <- assign_by_resistance(
      resistance_raster = resistance_raster,
      cluster_r         = cluster_r,
      contested         = contested
    )
  }

  ## ---- Step 4: Symmetric trimming loop -------------------------------------

  cluster_r <- trim_interfaces(
    cluster_r,
    min_cells = min_cluster_cells
  )

  cluster_r
}


expand_sparse_to_dense <- function(edges, n) {
  D <- matrix(Inf, n, n)
  diag(D) <- 0

  for (k in seq_len(nrow(edges))) {
    i <- edges$i[k]
    j <- edges$j[k]
    d <- edges$d[k]
    D[i, j] <- d
    D[j, i] <- d
  }

  as.dist(D)
}


#' @keywords internal
#' @tags noRd
cluster_connectivity <- function(
  x,
  input = c("features", "distance"),
  fixedClusters = TRUE,
  n = NULL,
  min.nc = 2,
  max.nc = 10,
  linkage = NULL
) {

  input <- match.arg(input)

  # drop incomplete cases (rows only)
  x <- x[stats::complete.cases(x), , drop = FALSE]

  # --- DISTANCE CONSTRUCTION ---
  if (input == "features") {

    # Euclidean feature space
    w_dist <- stats::dist(x, method = "euclidean")

    # default linkage
    if (is.null(linkage)) {
      linkage <- "ward.D2"
    }

    # sanity check
    if (linkage %in% c("ward.D", "ward.D2") == FALSE) {
      warning("Non-Ward linkage used on feature matrix")
    }

  } else {

    # x is already a distance matrix
    if (!is.matrix(x) || nrow(x) != ncol(x)) {
      stop("For input = 'distance', x must be a square distance matrix")
    }

    if (any(abs(diag(x)) > 1e-8)) {
      stop("Distance matrix must have ~0 diagonal")
    }

    w_dist <- as.dist(x)

    # default linkage
    if (is.null(linkage)) {
      linkage <- "complete"
    }

    # forbid Ward
    if (linkage %in% c("ward.D", "ward.D2")) {
      stop("Ward linkage is not valid for non-Euclidean distance matrices")
    }
  }

  # --- CLUSTERING ---
  if (fixedClusters == TRUE) {

    if (is.null(n)) {
      stop("n must be provided when fixedClusters = TRUE")
    }

    hc <- stats::hclust(w_dist, method = linkage)
    clusterCut <- stats::cutree(hc, n)

  } else {

    NoClusters <- NbClust::NbClust(
      diss = w_dist,
      distance = NULL,
      min.nc = min.nc,
      max.nc = max.nc,
      method = linkage
    )

    clusterCut <- NoClusters$Best.partition
    hc <- NULL
  }

  list(
    clusters = clusterCut,
    hclust = hc,
    distance = w_dist,
    linkage = linkage,
    input = input
  )
}

geographic_core_assignment <- function(
  seeds_sf,
  template_raster,
  max_dist,
  distance_method = "haversine"
) {

  seeds_sf$cluster_id <- seq_len(nrow(seeds_sf))

  seed_r <- terra::rasterize(
    terra::vect(seeds_sf),
    template_raster,
    field = "cluster_id"
  )

  d <- terra::distance(
    template_raster,
    seed_r,
    method = distance_method
  )

  nearest <- terra::nearest(seed_r)

  out <- terra::ifel(d <= max_dist, nearest, NA)
  names(out) <- "cluster_id"
  out
}


enforce_min_core_area <- function(cluster_r, min_cells) {

  tab <- table(terra::values(cluster_r))
  small <- as.integer(names(tab[tab < min_cells]))

  if (length(small)) {
    cluster_r[cluster_r %in% small] <- NA
  }

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

    df <- data.frame(
      cell    = candidates[,2],
      cluster = cluster_r[candidates[,1]]
    )

    df <- df |>
      dplyr::group_by(cell) |>
      dplyr::filter(dplyr::n_distinct(cluster) == 1) |>
      dplyr::ungroup()

    if (!nrow(df)) break

    sizes <- table(terra::values(cluster_r))
    df <- df[sizes[df$cluster] < max_cells_per_cluster, ]

    if (!nrow(df)) break

    cluster_r[df$cell] <- df$cluster
  }

  cluster_r
}


find_contested_cells <- function(cluster_r, width) {

  b <- terra::boundaries(cluster_r, directions = 8)
  terra::buffer(b, width)
}

assign_by_resistance <- function(
  resistance_raster,
  cluster_r,
  contested
) {

  clusters <- na.omit(unique(terra::values(cluster_r)))

  dist_stack <- terra::rast(lapply(clusters, function(k) {
    src <- cluster_r == k
    terra::distance(resistance_raster, src, mask = contested)
  }))

  idx <- terra::which.min(dist_stack)
  cluster_r[contested] <- clusters[idx]

  cluster_r
}

trim_interfaces <- function(cluster_r, min_cells) {

  repeat {

    sizes <- table(terra::values(cluster_r))
    bad <- as.integer(names(sizes[sizes < min_cells]))

    if (!length(bad)) break

    for (k in bad) {
      nbrs <- touching_clusters(cluster_r, k)
      for (j in nbrs) {
        iface <- interface_cells(cluster_r, k, j)
        cluster_r[iface] <- NA
      }
    }
  }

  cluster_r
}
