#' Determine which areas in a sample unit should be prioritized
#' 
#' @description This function offers guidance for prioritizing sample locations and orders. 
#' It has two parts, the first essentially creates a very simple 'heatmap' tapering off from the geographic center of each polygon to be sampled for a germplasm collection. 
#' To keep the spatial data 'light' here, we use also hard angled edges (ala vectorization of a raster), and only suggest a few `n_breaks` per polygon. 
#' Each of the levels from this side of the function relate to increasing distances from the geographic center of the polygon. 
#' With cells denoted '1' being the most ideal areas to sample from in each polygon to maintain well spaced distances across the focal taxons range. 
#' 
#' The second part offers a loose order for prioritizing the general order of collections.
#' Using a user specified `metric` it attempts to 'spread' samples across the range to reduce the variance in the distances between samples in case the desired number of samples is not achieved.
#' 
#' within the individual sample units returned by each *Sample function. 
#' The goal of either method is to avoid having collectors teams 'cheat' the system by repeatedly collecting along the border between two or more grid cells. 
#' While we understand that many teams may be collecting closely due to the species biology, land management, or other restrictions, the goal of this function is to try and guide them in dispersing there activity. 
#' 
#' The method used computes the geometric centroid of each region, if the center falls outside of the grid, it is snapped back onto the nearest location by default. 
#' Once the centers of each cell are calculated the remaining area of each grid has distances calculated between the centers and there locations. 
#' In final processing `n_breaks` are applied based on distances from the desired cell center to partition the space into different priority collection units. 
#'
#' Note that if you are submitting data from the `PolygonBasedSample`, the column `n`, must be maintained. 
#' @param x an sf/tibble/dataframe. a set of sample grids from any of the *Sample functions 
#' @param reps Numeric. The number of repetitions used in the sampling design. This is only used for messaging purposes at this time.
#' @param n_breaks Numeric. The number of breaks to return from the function, defaults to 3. Values beyond 5 of are questionable utility.  
#' @param verbose Bool. Whether to print messages to console or not, defaults to TRUE. 
#' @param metric character. The metric to minimize when ordering zones. Options are "var" (variance), "sd" (standard deviation), "energy" (sum of squared distances), and "cv" (coefficient of variation).
#' @return An sf object containing the prioritization zones within each sample unit, with columns: ID, SampleOrder, Level, and geometry. 
#' The ID corresponds to the sample unit, SampleOrder is the order in which to prioritize sampling from each unit,
#' Level is the priority level within each unit (1 being highest), and geometry is the spatial geometry of the prioritization zones.
#' @examples 
#' \dontrun{
#' nc <- sf::st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE) |>
#'   dplyr::select(NAME) |>
#'   sf::st_transform(5070) # should be in planar coordinate system. 
#' 
#' set.seed(1)
#' zones <- EqualAreaSample(nc, n = 20, pts = 1000, planar_proj = 32617, reps = 100)
#' 
#' # the function requires an input sampling strategy to create the prioritization areas
#' ps <- PrioritizeSample(zones$Geometry, n_breaks = 3, metric = 'energy')
#' 
#' ggplot2::ggplot() + 
#'   ggplot2::geom_sf(data = ps[['Geometry']],
#'     ggplot2::aes(fill = factor(Level)), color = 'white', lwd = 1) + 
#'   ggplot2::theme_void() + 
#'   ggplot2::labs(fill = 'Within Zone Priority:', title = 'Focal areas to center sampling within') +
#'   ggplot2::theme(legend.position= 'bottom')
#' 
#' ps[['Geometry']] |> ### to visualize without the priority zones within. 
#'   dplyr::group_by(SampleOrder) |> 
#'   dplyr::summarize(geometry = sf::st_union(geometry)) |>
#' 
#'   ggplot2::ggplot() + 
#'   ggplot2::geom_sf(ggplot2::aes(fill = SampleOrder), color = 'white') +
#'   ggplot2::geom_sf_label(ggplot2::aes(label = SampleOrder), color = 'white', size = 7) + 
#'   ggplot2::labs(fill = 'Sample Order', title = 'Priority guidance for sampling order') + 
#'   ggplot2::theme_void() + 
#'   ggplot2::theme(legend.position= 'bottom')
#' }
#' @export
PrioritizeSample <- function(x, reps, n_breaks = 3, verbose=TRUE, 
  metric = c("var", "sd", "energy", "cv")){

	# Ecoregion based sample is the one method where rows which are not meant to be sampled 
	# can be submitted to the function. We can detect these samples by a unique column name
	# which they feature. We will remove the rows which should not receive samples now. They 
	# will not be reattached to the ouput object as they are by default non-target sample areas. 
  if(any(grepl('^n$', colnames(x))==TRUE)){
    x <- x[x[,grep('^n$', colnames(x))]==1, ]
  }
	
	sf::st_agr(x) <- 'constant'
  POS <- sf::st_point_on_surface(x)
  
  # now we can create points in each row and measure the distance from these
  # points to the desired sample location. We can classify those distances
  # using breaks, creating discrete groups within each feature. 
  dists <- vector(mode = 'list', length = nrow(x))
  ctrlPts <- vector(mode = 'list', length = nrow(x))
  vorons <- vector(mode = 'list', length = nrow(x))
  ints <- vector(mode = 'list', length = nrow(x))
  for (i in seq_len(nrow(x))){
    
    # sample points and place into the polygons
    ctrlPts[[i]] <- sf::st_sample(x[i,], 100, type = 'regular')
    dists[[i]] <- as.numeric(sf::st_distance(ctrlPts[[i]], POS[i,]))
    ctrlPts[[i]] <- sf::st_as_sf(ctrlPts[[i]])
    ctrlPts[[i]] <-  cbind(
      lvl = cut(dists[[i]], breaks = n_breaks, labels = FALSE), 
      geometry = ctrlPts[[i]])
    
    # now chop up the sampling domain
    sf::st_agr(ctrlPts[[i]]) <- 'constant'
    vorons[[i]] <- sf::st_voronoi(sf::st_union(ctrlPts[[i]]))
    vorons[[i]] <- sf::st_intersection(sf::st_cast(vorons[[i]]), sf::st_union(x[i,]))
    
    # and add the pt classifications to each polygon
    vorons[[i]] <- cbind(
      cP = unlist(sf::st_intersects(vorons[[i]], ctrlPts[[i]])), 
      sf::st_as_sf(vorons[[i]]) ) 
    vorons[[i]]$Level <- ctrlPts[[i]]$lvl[vorons[[i]]$cP]
    
  }
  
  melt <- function(x){
    out <- dplyr::group_by(x, Level) |>
      dplyr::summarise(geometry = sf::st_union(x))
  }

  sample_levels <- lapply(vorons, melt)

  ### now perform the prioritization ordering of each seed zone. 
  polygon_orders <- order_by_distance_variance(x, metric = metric)
  SampleOrder = rep( polygon_orders [ x[['ID']] ], each = n_breaks)
  
    ## df, and geometry containing the rough positions to sample in
  sample_levels <- Map(cbind, sample_levels, ID = (seq_along(length(sample_levels)))) |>
    dplyr::bind_rows() |>
    dplyr::mutate(SampleOrder = SampleOrder) |>
    dplyr::select(ID, SampleOrder, Level, geometry) 

  list(Geometry = sample_levels)
}

#' Order zones by minimizing distance variance
#' 
#' @description Order a set of spatial zones based on minimizing distance variance
#' @param x an sf object containing the zones to be ordered
#' @param metric character. The metric to minimize when ordering zones. Options are "var" (variance), "sd" (standard deviation), "energy" (sum of squared distances), and "cv" (coefficient of variation).
#' @return A numeric vector representing the order of zones based on the specified distance variance metric
#' @keywords internal
order_by_distance_variance <- function(x, metric = c("var", "sd", "energy", "cv")) {

  ## create distance matrix
  dist_mat <- as.matrix(sf::st_distance(sf::st_point_on_surface(x)))

  metric <- match.arg(metric)
  n <- nrow(x)
  seed <- which.min(rowSums(dist_mat))

  SampleOrder <- seed
  remaining <- setdiff(seq_len(n), SampleOrder)

  # initialize distance field
  d <- dist_mat[, seed]

  score_fn <- switch(metric,
    var = function(x) stats::var(x),
    sd  = function(x) stats::sd(x),
    cv  = function(x) stats::sd(x) / mean(x),
    energy = function(x) sum(x^2)
  )

  while (length(remaining) > 0) {
    scores <- sapply(remaining, function(cand) {
      d_new <- pmin(d, dist_mat[, cand])
      score_fn(d_new)
    })

    best <- remaining[which.min(scores)]
    SampleOrder <- c(SampleOrder, best)
    d <- pmin(d, dist_mat[, best])
    remaining <- setdiff(remaining, best)
  }

  SampleOrder
}