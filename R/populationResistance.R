#' Identify clusters of populations least separated by landscape resistance
#' 
#' @description 
#' Creates a matrix, and supporting data structures, of a Isolation by Distance distances using a simplified framework. 
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
#' 
#' # this also can run internally in `population resistance`, but for time sakes is best to prep ahead of time
#' # especially if treating multiple species in the same domain. 
#' res <- buildResistanceSurface(
#'   base_raster = base_rast,
#'   oceans = ocean_r,
#'   lakes = lakes_r,
#'   rivers = rivers_r,
#'   tri = tri_r
#' )
#' 
#' # Run population resistance clustering
#' out <- populationResistance(
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
populationResistance <- function(
  populations_sf,
  base_raster,
  buffer_dist = 20000,
  planar_proj = NULL,
  n_points = 150,
  resistance_surface = NULL,
  oceans = NULL,
  lakes = NULL,
  rivers = NULL,
  tri = NULL,
  habitat = NULL,
  w_ocean = 100,
  w_lakes = 50,
  w_rivers = 20,
  w_tri = 1,
  w_habitat = 1,
  graph_method = c("complete", "delaunay"),
  ibr_method = "leastcost",
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
  res_rast <- buildResistanceSurface(
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

  # and now create a dense matrix suitable for use with most stats packages
  # we will impute Inf for values we did not compute pairwise distances with. 
  dense_mat  <- inflate_to_dense(g_out$edges)   # dense

  list(
    resistance_raster = res_rast,
    pop_raster = pop_raster,
    sampled_points = pts_sf,
    spatial_graph = g_out$graph,
    edge_list = g_out$edges,
    ibr_matrix = dense_mat
  )
}

#' @keywords internal
#' @tags noRd
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

#' @keywords internal
#' @tags noRd
#' @keywords internal
#' @tags noRd
build_spatial_graph <- function(
  pts_sf,
  method = "delaunay",
  great_circle = TRUE
) {
  method <- match.arg(method, choices = c("delaunay", "complete"))
  
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
  
  # --- distance matrix ---
  if (great_circle) {
    D <- sf::st_distance(pts_sf)
    D <- units::drop_units(D)
  } else {
    D <- as.matrix(dist(coords))
  }
  
  # --- graph construction ---
  if (method == "delaunay") {
    nb <- spdep::tri2nb(coords)
    
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
    
  } else if (method == "complete") {
    # Complete graph: all pairwise connections
    edges <- t(combn(n, 2))
  }
  
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  
  # Attach weights
  w <- D[cbind(edges[, 1], edges[, 2])]
  igraph::E(g)$weight <- w
  
  # --- return both graph + edge table ---
  edge_df <- igraph::as_data_frame(g, what = "edges")
  
  list(
    graph = g,
    edges = edge_df,
    distances = D
  )
}

#' @keywords internal
#' @tags noRd
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

  R
}

#' @keywords internal
#' @noRd
inflate_to_dense <- function(edges) {
  # edges: data.frame with columns 'from', 'to', 'd' (distance)
  
  # determine number of points from edges
  n <- max(c(edges$from, edges$to))

  finite_max <- max(edges$weight)*2
  D <- matrix(finite_max, n, n)
  diag(D) <- 0

  for (k in seq_len(nrow(edges))) {
    i <- edges$from[k]
    j <- edges$to[k]
    d <- edges$weight[k]
    D[i, j] <- d
    D[j, i] <- d
  }
  
  D
}

