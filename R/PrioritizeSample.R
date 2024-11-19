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
#' This is accomplished by sampling each grid cell `rep` times, and computing the variance between all pairwise distances (e.g. with 1000 replicants and 20 samples desired, place one point in each of the 20 grids 1000 times, and calulcate the minimum variance in distance between all reps). 
#' The rep with the lowest variance is then chosen for the next step, where the distance between the selected point and the rest of the cell are calculated. 
#' Both methods end up feeding into the same final processing step. 
#' In final processing `n_breaks` are applied based on distances from the desired cell center to partition the space into different priority collection units. 
#'
#' Note that if you are submitting data from the ecoregion based sample, the column `n`, must be maintained. 
#' @param x an sf/tibble/dataframe. a set of sample grids from any of the *Sample functions 
#' @param method Character String. One of "centered" or "distance" to dispatch the supported methods. If missing defaults to "centered". 
#' @param reps If using method "distance" the number of times to repeat the sampling procedure, defaults to 500. 
#' @param n_breaks Numeric. The number of breaks to return from the function, defaults to 3. Values much higher than that are untested, and beyond 5 of questionable utility.  
#' @examples /dontrun{}
#' @export
PrioritizeSample <- function(x, method, reps, b_breaks){

	if(missing(method)){method <- 'centered'}
	if(missing(reps) & method=='distance'){reps <- 500}
	if(missing(n_breaks)){n_breaks <- 3}
	
	
	# Ecoregion based sample is the one method where rows which are not meant to be sampled 
	# can be sumbitted to the function. We can detect these samples by a unique column name
	# which they feature. We will remove the rows which should not receive samples now. They 
	# will not be reattached to the ouput object as they are by default non-target sample areas. 
	
	
	
	if(method=='centered'){
	
	
	} else {}


}
