#' Design additional collections around already existing collections
#' 
#' @description This function utilizes a regular, or nearly so in the case of existing collections, grid of points to develop a sampling scheme or n polygons. 
#' @param polygon the input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
#' @param n Numeric. The total number of desired collections. Defaults to 20.
#' @param collections an sf point geometry data set of where existing collections have been made.
#' @param reps further arguments passed to np.boot 
#' @param BS.reps number of bootstrap replicates for evaluating results. 
#' @examples 
#' #' Design additional collections around already existing collections
#' ri <- spData::us_states |>
#'   dplyr::select(NAME) |>
#'   dplyr::filter(NAME == 'Rhode Island') |>
#'   sf::st_transform(32617)
#' existing_collections <- sf::st_sample(ri, size = 5) |>
#'   sf::st_as_sf() |>
#'   dplyr::rename(geometry = x)
#'
#' system.time(
#'   out <- OpportunisticSample(polygon = ri, BS.reps=4999) 
#' ) # set very low for example
#' # the function is actually very fast; 150 voronoi reps, with 9999 BS should only take about
#' # 7 seconds per species so not much concern on the speed end of things.
#' ggplot2::ggplot() + 
#'   ggplot2::geom_sf(data = out$Geometry, ggplot2::aes(fill = ID)) + 
#'   ggplot2::geom_sf(data = existing_collections) 
#'  
#' @return A list containing two sublists, the first of which 'SummaryData' details the number of voronoi polygons generated, and the results of the bootstrap simulations. The second 'Geometry', contains the final spatial data products, which can be written out on your end. See the vignette for questions about saving the two main types of spatial data models (vector - used here, and raster). 
#' @export
OpportunisticSample <- function(polygon, n, collections, reps, BS.reps){
  
  if(missing(n)){n = 20}
  if(missing(reps)){reps = 100}
  if(missing(BS.reps)){BS.reps = 9999}
  
  # we apply the voronoi process a number of replicated times, defaults to 100
  # this allows us to be confident that we get an OK set of results from the function
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
  voronoiPolygons <- voronoiPolygons[lapply(get_elements(voronoiPolygons, 'Polygons'), length) >= 20]
  
  # we use variance to determine the configuration of voronoi polygons which have
  # the most equally sized polygons. 
  variance <- unlist(get_elements(voronoiPolygons, 'Variance'))
  
  # but we only select sets of records which actually meet the sample size requirements, 
  # a few weird ones always get through otherwise. 
  SelectedSample <- voronoiPolygons[which.min(variance)][[1]]$Polygons |>
    sf::st_as_sf()
  
  # try and number all of the polygons from one direction, should help with their
  # usage in the field. 
  ss_cents <- sf::st_point_on_surface(SelectedSample)
  ss_cents <- ss_cents |>
    dplyr::mutate(
      X = sf::st_coordinates(ss_cents)[,1],
      Y = sf::st_coordinates(ss_cents)[,2]
    ) |>
    dplyr::arrange(-Y, X) |>
    dplyr::mutate(ID = 1:dplyr::n()) 
  
  # now assign the arranged notation to the data set overwriting the original 
  # randomly assigned ID's
  
  SelectedSample <- dplyr::mutate(
    SelectedSample, 
    ID = as.numeric(sf::st_intersects(SelectedSample, ss_cents)), .before = 1) |>
    dplyr::rename(dplyr::any_of(c(geometry = 'x', geometry = 'X', geometry = 'geom')))

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
    probs = c(0.001), 
    level = 0.95) 
  
  # Create an output object containing the bootstrap estimates and the observed variance
  # for the grid, and write out the information on the number of replicates etc. 
  list(
    'SummaryData' = data.frame(
      'Metric' = c(
        'variance.observed', 'quantile.0.001', 'lwr.95.CI',
        'upr.95.CI', 'Voronoi.reps.asked', 'Voronoi.reps.received', 'BS.reps'), 
      'Value' = c(
        min(variance), npbs[['t0']], npbs[['bca']][['lower']], 
        npbs[['bca']][['upper']], reps, length(variance),  BS.reps)
    ),
  'Geometry' = SelectedSample) # and of course the spatial data!
  
} 