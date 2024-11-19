#' Determine which areas of a sample unit should be prioritized
#' 
#' @description This function offers two ways to enforce some prioritization within the individual sample units returned by each *Sample function. 
#' The goal of either method is to avoid having collectors teams 'cheat' the system by repeatedly collecting along the border between two or more grid cells. 
#' While we understand that many teams may be collecting closely due to the species biology, land management, or other restrictions, the goal of this function is to try and guide them in dispersing there activity. 
#' There are two implemented approaches for achieving these spatial guidelines. 
#' 
#' 'centered', where the geometric centroid of each region is computed, if the center falls outside of the grid, it is snapped back onto the nearest location by default. 
#' Once the centers of each cell are calculated the remaining area of each grid has distances calculated between the centers and there locations.
#' 'distance', where the variance in distance between all grid cells is attempted to be minimized. 
#' This is accomplished by sampling each grid cell `rep` times, and computing the variance between all pairwise distances (e.g. with 1000 replicates and 20 samples desired, place one point in each of the 20 grids 100 times, and calculate the minimum variance in distance between all reps). 
#' The rep with the lowest variance is then chosen for the next step, where the distance between the selected point and the rest of the cell are calculated. 
#' Note that the 'distance' method is many magnitudes of order slower than the centered method. 
#' Both methods end up feeding into the same final processing step. 
#' In final processing `n_breaks` are applied based on distances from the desired cell center to partition the space into different priority collection units. 
#'
#' Note that if you are submitting data from the ecoregion based sample, the column `n`, must be maintained. 
#' @param x an sf/tibble/dataframe. a set of sample grids from any of the *Sample functions 
#' @param method Character String. One of "centered" or "distance" to dispatch the supported methods. If missing defaults to "centered". 
#' @param reps If using method "distance" the number of times to repeat the sampling procedure, defaults to 500. 
#' @param n_breaks Numeric. The number of breaks to return from the function, defaults to 3. Values much higher than that are untested, and beyond 5 of questionable utility.  
#' @param verbose Boolean. Whether to display a progress bar if implementing method distance. 
#' @examples /dontrun{}
#' @export
PrioritizeSample <- function(x, method, reps, n_breaks, verbose){

	if(missing(method)){method <- 'centered'}
	if(missing(reps) & method=='distance'){reps <- 500}
	if(missing(n_breaks)){n_breaks <- 3}
	
	# Ecoregion based sample is the one method where rows which are not meant to be sampled 
	# can be submitted to the function. We can detect these samples by a unique column name
	# which they feature. We will remove the rows which should not receive samples now. They 
	# will not be reattached to the ouput object as they are by default non-target sample areas. 
  if(any(grepl('^n$', colnames(x))==TRUE)){
    x <- x[x[,grep('^n$', colnames(x))]==1, ]
  }
	
	if(method=='centered'){
	
	  sf::st_agr(x) <- 'constant'
	  POS <- sf::st_point_on_surface(x)
	
	} else {
	  
	  # sample each ecoregion reps times. 
	  repetitions <- vector(mode = 'list', length = reps)
	  message('Undergoing sampling repetitions. ')
	  for (i in seq_along(1:reps)){ 
	    repetitions[[i]] <- iter_sample(x)
	    if(verbose ==TRUE){svMisc::progress(i, reps, progress.bar = TRUE)}
	    }
	  return(repetitions)
	}

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
      lvl = cut(dists[[i]], breaks = 4, labels = FALSE), 
      geometry = ctrlPts[[i]])
    
    # now chop up the sampling domain
    sf::st_agr(ctrlPts[[i]]) <- 'constant'
    vorons[[i]] <- sf::st_voronoi(sf::st_union(ctrlPts[[i]]))
    vorons[[i]] <- sf::st_intersection(sf::st_cast(vorons[[i]]), sf::st_union(x[i,]))
    
    # and add the pt classifications to each polygon
    vorons[[i]] <- cbind(
      cP = unlist(sf::st_intersects(vorons[[i]], ctrlPts[[i]])), 
      sf::st_as_sf(vorons[[i]]) ) 
    vorons[[i]]$lvl <- ctrlPts[[i]]$lvl[vorons[[i]]$cP]
    
  }
  
  melt <- function(x){
    out <- dplyr::group_by(x, lvl) |>
      dplyr::summarise(geometry = sf::st_union(x))
  }

  sample_levels <- lapply(vorons, melt)
  sample_levels <- Map(cbind, sample_levels, ID = (1:length(sample_levels))) |>
    dplyr::bind_rows() |>
    dplyr::select(ID, lvl, geometry)
  
  return(sample_levels)
}

# beginning
x <- st_read(system.file("shape/nc.shp", package="sf")) |>
  sf::st_transform(5070)

out <- PrioritizeSample(x, method = 'centered')

# now we can union the polygons together by their levels. 

hist(unlist(get_elements(out, 'vDistance')))
which.max(unlist(get_elements(out, 'vDistance')))
out[[1]]$Points

#' sample points from each row of an sf object, and calculate mean distance between
#' all points, and then take the variance of those measurement. 
iter_sample <- function(x){
  
  obs <- split(x, f = 1:nrow(x))
  pts <- lapply(obs, function(x){
    crds <- sf::st_sample(x, 1) |> 
      sf::st_coordinates() 
    df <- data.frame(X = crds[1], Y = crds[2])
  }) |>
    dplyr::bind_rows() |>
    sf::st_as_sf(coords = c(x = 'X', y = 'Y'), crs = sf::st_crs(x))
  
  d <- sf::st_distance(pts)
  d <- stats::var(apply(d, MARGIN = 2, base::mean, na.rm = TRUE), na.rm = TRUE)
  
  return(
    list(
      Points = pts, 
      vDistance = d
    )
  )
}


#' Recursively grab a named component of a list. 
#' @param x a list of lists
#' @param element the quoted name of the list element to extract. 
#' @keywords internal
get_elements <- function(x, element) { # @ StackOverflow Allan Cameron 
  if(is.list(x))
  {
    if(element %in% names(x)) x[[element]]
    else lapply(x, get_elements, element = element)
  }
}

ggplot() + 
  geom_sf(data = x) + 
  geom_sf(data = out[[49]], aes(fill = lvl))
