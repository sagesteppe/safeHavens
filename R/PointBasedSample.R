#' Design additional collections around already existing collections
nc <- spData::us_states |>
  dplyr::select(NAME) |>
  dplyr::filter(NAME == 'Rhode Island') |>
  sf::st_transform(32617)
existing_collections <- sf::st_sample(nc, size = 5) |>
  sf::st_as_sf() |>
  dplyr::rename(geometry = x)

out <- PointBasedSample(polygon = nc, reps = 150)
head(out$SummaryData)
plot(out$Geometry)

#' Generate a sampling grid using points as the input data - including existing collections
#' 
#' This function utilizes a regular, or nearly so in the case of existing collections, grid of points 
#' to develop a sampling scheme or n polygons. 
#' @param polygon the input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
#' @param n Numeric. The total number of desired collections. Defaults to 20.
#' @param collections an sf point geometry data set of where existing collections have been made.
#' @param reps further arguments passed to np.boot 
#' @BS.reps number of bootstrap replicates for evaluating results. 
#' @examples \donttest{
#' #' Design additional collections around already existing collections
#' nc <- spData::us_states |>
#'   dplyr::select(NAME) |>
#'   dplyr::filter(NAME == 'Rhode Island') |>
#'   sf::st_transform(32617)
#' existing_collections <- sf::st_sample(nc, size = 5) |>
#'   sf::st_as_sf() |>
#'   dplyr::rename(geometry = x)
#'
#' out <- PointBasedSample(polygon = nc, reps = 10, BS.reps = 10) # set very low for example
#' # the function is actually very fast; 150 voronoi reps, with 9999 BS should only take about
#' # 2 seconds per species so not much concern on the speed end of things!
#' head(out$SummaryData)
#' plot(out$Geometry)
#' }
#' @export
PointBasedSample <- function(polygon, n, collections, reps, BS.reps){
  
  if(missing(n)){n = 20}
  if(missing(reps)){reps = 100}
  if(missing(BS.reps)){BS.reps = 9999}
  
  #' Recursively grab a named component of a list. 
  #' @param x a list of lists
  #' @param element the quoted name of the list element to extract. 
  get_elements <- function(x, element) { # @ StackOverflow Allan Cameron 
    if(is.list(x))
    {
      if(element %in% names(x)) x[[element]]
      else lapply(x, get_elements, element = element)
    }
  }
  
  #' Make a voronoi sample of an area n times
  #' 
  #' Split an area up into n polygons of roughly equal area, optionally removing 
  #' some of the default points and replacing them with existing collections to 
  #' build the future collections around. 
  #' 
  #' @param polygon The input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
  #' @param n Numeric. The total number of desired collections. Defaults to 20.
  #' @param collections an sf point geometry data set of where existing collections have been made.
  #' @param reps Numeric. The number of times to rerun the voronoi algorithm, the set of polygons with the most similar sizes, as
  #' measured using their variance of areas will be selected. Defaults to 150, which may accomplish around 100 succesful iterations.  
  VoronoiSampler <- function(polygon, n, collections, reps){
    
    pts <- sf::st_sample(polygon, size = n, type = 'regular') |> 
      sf::st_as_sf() |> 
      dplyr::rename(geometry = x) 
    
    if(!missing(collections)){
      pts <- dplyr::bind_rows(
        collections, 
        pts[-sf::st_nearest_feature(collections, pts),], ) |>
        dplyr::slice_head(n=n)
    } else {pts <- dplyr::slice_head(pts, n=n)}
    
    vorons <- sf::st_voronoi(sf::st_union(pts), sf::st_as_sfc(sf::st_bbox(polygon)))
    vorons <- sf::st_intersection(sf::st_cast(vorons), sf::st_union(polygon))
    variance <- var(as.numeric(sf::st_area(vorons))/10000)
    
    # need to define two slots, one for the variance numeric results, and one for the polygons
    
    return(
      list(
        'Variance' = variance,
        'Polygons' = vorons
      ))
  }
  
  # we apply the voronoi process a number of replicated times, defaults to 100
  if(missing(collections)){
    voronoiPolygons <- replicate(
      reps, 
      VoronoiSampler(polygon = polygon, n = n), 
      simplify = FALSE)} else {
        voronoiPolygons <- replicate(
          reps, 
          VoronoiSampler(polygon = polygon, n = n, collections = collections), 
          simplify = FALSE)
      }
  
  # sf oftentimes gives fewer points than asked for, we will keep only the objects 
  # which have the desired number of points
  voronoiPolygons <- voronoiPolygons[
    lapply(get_elements(voronoiPolygons, 'Polygons'), length) >= n]
  
  # we use variance to determine the configuration of voronoi polygons which have
  # the most equally sized polygons. 
  variance <- unlist(get_elements(voronoiPolygons, 'Variance'))
  
  # but we only select sets of records which actually meet the sample size requirements, 
  # a few weird ones always get through otherwise. 
  SelectedSample <- voronoiPolygons[which.min(variance)][[1]]$Polygons 
  
  SelectedSample <- sf::st_as_sf(SelectedSample) |>
    dplyr::rename(geometry = x)
  
  # Determining the 0.1% quantile for the variance in size of the sampling grids. 
  # Using non-parametric approaches, of bootstrap resampling (replicates = 9999) ,
  # with an 95% confidence level. 
  
  # we can show that the polygon arrangement we have chosen is in the top 1000 of
  # options if npbs[["bca"]][["lower"]] > min(variance) == TRUE . 
  # If the above condition is not meet, we can also say that it is less than the estimate
  # npbs[["t0"]] < min(variance)
  npbs <- nptest::np.boot(
    x = variance, 
    statistic = quantile, 
    R = BS.reps, 
    na.rm = TRUE, 
    probs = c(0.001), 
    level = 0.95) 
  
  # Create an output object containing the bootstrap estimates and the observed variance
  # for the grid, and write out the information on the number of replicates etc. 
  
  output <- list(
    'SummaryData' = data.frame(
      'Metric' = c(
        'variance.observed', 'quantile.0.001', 'lwr.95.CI',
        'upr.95.CI', 'Voronoi.reps.asked', 'Voronoi.reps.received', 'BS.reps'), 
      'Value' = c(
        min(variance), npbs[['t0']], npbs[['bca']][['lower']], 
        npbs[['bca']][['upper']], reps, length(variance),  BS.reps)
    ),
  'Geometry' = SelectedSample)
  
  return(output)
  
} 

