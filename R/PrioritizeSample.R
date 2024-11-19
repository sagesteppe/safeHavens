#' Determine which areas of a sample unit should be prioritized
#' 
#' @description This function offers two ways to enforce some prioritization within the individual sample units returned by each *Sample function. 
#' The goal of either method is to avoid having collectors teams 'cheat' the system by repeatedly collecting along the border between two or more grid cells. 
#' While we understand that many teams may be collecting closely due to the species biology, land management, or other restrictions, the goal of this function is to try and guide them in dispersing there activity. 
#' 
#' THe method used computes the geometric centroid of each region is computed, if the center falls outside of the grid, it is snapped back onto the nearest location by default. 
#' Once the centers of each cell are calculated the remaining area of each grid has distances calculated between the centers and there locations. 
#' In final processing `n_breaks` are applied based on distances from the desired cell center to partition the space into different priority collection units. 
#'
#' Note that if you are submitting data from the ecoregion based sample, the column `n`, must be maintained. 
#' @param x an sf/tibble/dataframe. a set of sample grids from any of the *Sample functions 
#' @param n_breaks Numeric. The number of breaks to return from the function, defaults to 3. Values much higher than that are untested, and beyond 5 of questionable utility.  
#' @examples /dontrun{
#' nc <- sf::st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE) |>
#' dplyr::select(NAME) |>
#'   sf::st_transform(5070) # should be in planar coordinate system. 
#' 
#' set.seed(1)
#' zones <- EqualAreaSample(nc, n = 20, pts = 1000, planar_projection = 32617, reps = 100)
#' 
#' # the function requires an input sampling strategy to create the prioritization areas
#' ps <- PrioritizeSample(zones$Geometry, method = 'centered', n_breaks = 3)
#' 
#' ggplot2::ggplot() + 
#'  ggplot2::geom_sf(data = ps,  aes(fill = factor(Level))) +
#'  ggplot2::geom_sf(data = zones$Geometry, color = 'red', fill = NA, linewidth = 1) 
#' }
#' @export
PrioritizeSample <- function(x, method, reps, n_breaks, verbose){

	if(missing(n_breaks)){n_breaks <- 3}
	
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
  for (i in seq_along(1:nrow(x))){
    
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
  sample_levels <- Map(cbind, sample_levels, ID = (1:length(sample_levels))) |>
    dplyr::bind_rows() |>
    dplyr::select(ID, Level, geometry)
  
  return(sample_levels)
}