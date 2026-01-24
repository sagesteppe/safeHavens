#' Identify clusters of populations least separated by landscape resistance
#' 
#' @description 
#' Implements a simplified Isolation by Resistance (IBR) framework. 
#' Given a raster surface parameterizing resistance to movement (e.g., oceans, lakes, rivers, topographic roughness, or habitat suitability), this function calculates pairwise landscape resistance among populations and identifies clusters of populations that are least isolated from each other.
#' 
#' The function can internally build a resistance raster, but for performance across multiple species, it is recommended to provide pre-computed surfaces. 
#' Resistance is treated as isotropic (no slope directionality). 
#' Clustering is based on pairwise IBR distances and can be visualized or used to guide sampling schemes.
#' 
#' @param base_raster SpatRaster. Base raster for the study area. Provides template geometry and resolution.
#' @param populations_sf sf/tibble/data.frame. Coordinates of populations to be clustered. Must contain POINT geometries.
#' @param buffer_dist Numeric. Distance to buffer population points to assign likely occupied cells (optional; ignored if using `sample_population_cells` directly).
#' @param planar_proj Numeric. EPSG code used to project `populations_sf` for buffering operations.
#' @param n_points Numeric (defaults to 50). Number of points to sample for creating the surface. Increasing points drastically increases compute time. 
#' @param resistance_surface SpatRaster. Optional pre-computed resistance raster. If provided, the raster-building arguments are ignored.
#' @param oceans SpatRaster. Binary (0/1) raster for ocean cells. Used to increase movement cost.
#' @param lakes SpatRaster. Binary (0/1) raster for lakes.
#' @param rivers SpatRaster. Binary (0/1) raster for rivers.
#' @param tri SpatRaster. Continuous raster of topographic roughness (TRI). Used to increase cost in mountainous terrain.
#' @param habitat SpatRaster. Continuous raster of habitat suitability. Low values increase cost.
#' @param w_ocean Numeric. Weight applied to oceans (default 2000).
#' @param w_lakes Numeric. Weight applied to lakes (default 200).
#' @param w_rivers Numeric. Weight applied to rivers (default 20).
#' @param w_tri Numeric. Weight applied to TRI (default 1).
#' @param w_habitat Numeric. Weight applied to habitat suitability (default 1).
#' @param graph_method Character. One of `"mst"`, `"gabriel"`, or `"delaunay"`. Defines the spatial graph connecting population points.
#' @param ibr_method Character. One of `"leastcost"` or `"randomwalk"`. Random-walk converts pairwise distances to effective resistance.
#' @param epsilon Numeric. Small value to stabilize conversion to conductance in `randomwalk` method (default 1e-6).
#' 
#' @return A list with the following elements:
#' \describe{
#'   \item{resistance_raster}{The parameterized resistance SpatRaster used for IBR calculations.}
#'   \item{sampled_points}{sf POINT object of population cells used for graph calculations.}
#'   \item{spatial_graph}{igraph object representing the population connectivity graph.}
#'   \item{edge_list}{Data.frame of graph edges used for IBR computation.}
#'   \item{ibr_matrix}{Symmetric numeric matrix of pairwise landscape resistance between populations.}
#' }
#' 
#' @examples
#' \dontrun{
#' # Prepare resistance raster once
#' res <- build_resistance_surface(
#'   base_raster = base_rast,
#'   oceans = ocean_r,
#'   lakes = lakes_r,
#'   rivers = rivers_r,
#'   tri = tri_r
#' )
#' 
#' # Run population resistance clustering
#' out <- population_resistance(
#'   populations_sf = pops_sf,
#'   base_raster = base_rast,
#'   resistance_surface = res,
#'   n_points = 5,
#'   graph_method = "mst",
#'   ibr_method = "leastcost"
#' )
#' 
#' # Access the results
#' ibr_matrix <- out$ibr_matrix
#' graph <- out$spatial_graph
#' }
#' 
#' @export
population_resistance_cluster <- function(
  populations_sf,
  base_raster,
  buffer_dist = 10000,
  planar_proj = NULL,
  n_points = 50,
  resistance_surface = NULL,
  oceans = NULL,
  lakes = NULL,
  rivers = NULL,
  tri = NULL,
  habitat = NULL,
  w_ocean = 2000,
  w_lakes = 200,
  w_rivers = 20,
  w_tri = 1,
  w_habitat = 1,
  graph_method = c("mst", "gabriel", "delaunay"),
  ibr_method = c("leastcost", "randomwalk"),
  epsilon = 1e-6,
  min_resistance = 1L
) {

  graph_method <- match.arg(graph_method)
  ibr_method <- match.arg(ibr_method)

  ## prepare occurrence data for processing. Buffer occurrence points for system
  ## once buffered, each polygon will be extracted and treated as a 'population' to 
  # have representative points sampled from #. 

  # --- buffer occurence data to create population polygons --- #
  if(!is.null(planar_proj)){
      populations_sf <- populations_sf |>
        sf::st_transform(planar_proj) |>
        sf::st_buffer(buffer_dist) |>
        sf::st_cast('POLYGON') |>
        dplyr::mutate(pop_id = seq_len(dplyr::n())) |>
        sf::st_transform(terra::crs(base_raster))
  } else {stop('arg to `planar proj` is required')}

  # --- terra will (sensibly) try to keep the raster on disk and only pull
  # in values as required. we will forcce the raster into memory - which 
  # is safe given how coarse they are. 
  base_raster <- terra::rast(base_raster)
  terra::values(base_raster) <- 0
  names(base_raster) <- "base"

  ## assign population IDs to the base raster. 
  pop_raster <- terra::rasterize(populations_sf, base_raster, field = 'pop_id')
  names(pop_raster) <- 'pop_id'

  # --- build or clone resistance raster ---
  res_rast <- build_resistance_surface(
    base_raster = base_raster,
    resistance_surface = resistance_surface,
    oceans = oceans,
    lakes = lakes,
    rivers = rivers,
    tri = tri,
    habitat = habitat,
    w_ocean = w_ocean,
    w_lakes = w_lakes,
    w_rivers = w_rivers,
    w_tri = w_tri,
    w_habitat = w_habitat,
    min_resistance = min_resistance
  )

  # --- sample population points ---
  pts_sf <- sample_population_cells(
    pop_raster = pop_raster,
    n_total = n_points,
    method = "spread"
  )

  # --- build spatial graph ---
  g_out <- build_spatial_graph(
    pts_sf = pts_sf,
    method = graph_method
  )

  # --- compute pairwise IBR ---
  Rmat <- compute_ibr_edges(
    resistance_raster = res_rast,
    pts_sf = pts_sf,
    edges = g_out$edges,
    method = ibr_method,
    epsilon = epsilon
  )

  list(
    resistance_raster = res_rast,
    sampled_points = pts_sf,
    spatial_graph = g_out$graph,
    edge_list = g_out$edges,
    ibr_matrix = Rmat
  )
}

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


sample_population_cells <- function(
  pop_raster,
  n_total,
  method = "spread"
) {

  # global regular sampling
  pts <- terra::spatSample(
    pop_raster,
    size = n_total,
    method = method,
    as.points = TRUE,
    na.rm = TRUE, 
    xy = TRUE
  ) |> 
    sf::st_as_sf()
}


build_spatial_graph <- function(
  pts_sf,
  method = c("mst", "gabriel", "delaunay"),
  great_circle = TRUE
) {
  method <- match.arg(method)

  if (!inherits(pts_sf, "sf")) {
    stop("pts_sf must be an sf object with POINT geometry")
  }

  if (!all(sf::st_geometry_type(pts_sf) == "POINT")) {
    stop("pts_sf must contain only POINT geometries")
  }

  n <- nrow(pts_sf)
  if (n < 2) {
    stop("At least two points are required to build a graph")
  }

  # --- coordinates ---
  coords <- sf::st_coordinates(pts_sf)

  # --- distance matrix (used for MST weights) ---
  if (great_circle) {
    if (sf::st_is_longlat(pts_sf)) {
      D <- sf::st_distance(pts_sf)
      D <- units::drop_units(D)
    } else {
      warning("great_circle=TRUE but CRS is not lon/lat; using planar distances")
      D <- as.matrix(dist(coords))
    }
  } else {
    D <- as.matrix(dist(coords))
  }

  # --- graph construction ---
  if (method == "mst") {

    # Safety: check finite distances
    if (any(!is.finite(D))) {
      stop("Distance matrix contains non-finite values; cannot compute MST")
    }

    g_full <- igraph::graph_from_adjacency_matrix(
      D,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )

    # Safety: ensure connectivity
    comps <- igraph::components(g_full)
    if (comps$no > 1) {
      stop(
        "Point set is disconnected under the chosen distance metric; ",
        "MST cannot be constructed"
      )
    }

    g <- igraph::mst(g_full)

  } else {

    # --- proximity graphs via spdep ---
    if (method == "delaunay") {
      nb <- spdep::tri2nb(coords)

    } else if (method == "gabriel") {
      nb <- spdep::neigh2nb(spdep::gabrielneigh(coords))
    }

    # Safety: empty neighbor sets
    if (any(lengths(nb) == 0)) {
      warning(
        "Some points have no neighbors under '", method,
        "' graph; graph may be disconnected"
      )
    }

    # Convert nb -> edge list
    edges <- do.call(
      rbind,
      lapply(seq_along(nb), function(i) {
        if (length(nb[[i]]) == 0) return(NULL)
        cbind(i, nb[[i]])
      })
    )

    # Ensure unique undirected edges
    edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]

    g <- igraph::graph_from_edgelist(edges, directed = FALSE)

    # Attach weights
    w <- D[cbind(edges[, 1], edges[, 2])]
    igraph::E(g)$weight <- w
  }

  # --- return both graph + edge table ---
  edge_df <- igraph::as_data_frame(g, what = "edges")

  list(
    graph = g,
    edges = edge_df,
    distances = D
  )
}


compute_ibr_edges <- function(
  resistance_raster,
  pts_sf,
  edges,
  method = c("leastcost", "randomwalk"),
  epsilon = 1e-6
) {
  method <- match.arg(method)

  n <- nrow(pts_sf)
  R <- matrix(NA_real_, n, n)

  pts_v <- terra::vect(pts_sf)

  for (i in unique(edges$from)) {

    src <- pts_v[i]
    d <- terra::distance(resistance_raster, src)

    idx <- which(edges$from == i)
    for (k in idx) {
      j <- edges$to[k]
      R[i, j] <- terra::extract(d, pts_v[j])[, 2]
      R[j, i] <- R[i, j]
    }
  }

  if (method == "randomwalk") {
    W <- 1 / (R + epsilon)
    diag(W) <- 0
    g <- igraph::graph_from_adjacency_matrix(
      W, mode = "undirected", weighted = TRUE
    )
    R <- igraph::resistance_distance(g)
  }

  R
}


build_resistance_surface <- function(
  base_raster,
  resistance_surface = NULL,
  oceans = NULL,
  lakes = NULL,
  rivers = NULL,
  tri = NULL,
  habitat = NULL,
  w_ocean = 1000,
  w_lakes = 200,
  w_rivers = 20,
  w_tri = 1,
  w_habitat = 1,
  min_resistance = 1L
) {

  if (!is.null(resistance_surface)) {
    res <- resistance_surface

    if (!terra::compareGeom(res, base_raster, stopOnError = FALSE)) {
      stop("resistance_surface must match base_raster geometry")
    }

  } else {
    # initialize blank raster
    res <- terra::rast(base_raster)
    terra::values(res) <- 0

    # add weighted features
    if (!is.null(oceans))  res[oceans == 1]  <- res[oceans == 1]  + w_ocean
    if (!is.null(lakes))   res[lakes == 1]   <- res[lakes == 1]   + w_lakes
    if (!is.null(rivers))  res[rivers == 1]  <- res[rivers == 1]  + w_rivers
    if (!is.null(tri))     res <- res + w_tri * terra::setValues(res, scale(terra::values(tri)))
    if (!is.null(habitat)) res <- res + w_habitat * (1 / (habitat + 0.01))
  }

  # --- clamp only if raster has values ---
  vals <- terra::values(res)
  if (!is.null(vals) && length(vals) > 0 && any(!is.na(vals))) {
    res <- terra::clamp(res, lower = min_resistance)
  }

  # convert to int32
  res <- terra::as.int(res)
}
