#' Create equal area polygons over a geographic range
#' 
#' @description This function creates `n` geographic clusters over a geographic area (`x`), typically a species
#' range, using kmeans clustering. 
#' @param x An SF object or terra spatraster. the range over which to generate the clusters.
#' @param n Numeric. the number of clusters desired. Defaults to 20. 
#' @param pts Numeric. the number of points to use for generating the clusters, these will be placed in a grid like fashion across `x`. The exact number of points used may deviate slightly from the user submitted value to allow for equidistant spacing across `x`. Defaults to 5,000.
#' @param planar_projection Numeric, or character vector. An EPSG code, or a proj4 string, for a planar coordinate projection, in meters, for use with the function. For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise a continental scale projection like 5070 See https://projectionwizard.org/ for more information on CRS. The value is simply passed to sf::st_transform if you need to experiment. 
#' @param returnProjected Boolean. Whether to return the data set in the original input CRS (FALSE), or in the new `projection` (TRUE). Defaults to FALSE. 
#' @param reps Numeric. The number of times to rerun the voronoi algorithm, the set of polygons with the most similar sizes, as
#' measured using their variance of areas will be selected. Defaults to 100. 
#' @param BS.reps number of bootstrap replicates for evaluating results. 
#' @examples 
#' nc <- sf::st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE) |>
#' dplyr::select(NAME)
#'
#' set.seed(1)
#' system.time(
#'   zones <- EqualAreaSample(nc, n = 20, pts = 500, planar_projection = 32617, reps = 50)
#' )
#' 
#' plot(nc, main = 'Counties of North Carolina')
#' plot(zones$Geometry, main = 'Clusters')
#' zones$SummaryData
#' 
#' @return A list containing two objects, the first the results of bootstrap simulations.
#' The second an sf dataframe containing the polygons with the smallest amount of variance in size. 
#' @export
EqualAreaSample <- function(x, n = 20, pts = 5000, planar_projection, returnProjected, reps = 100, BS.reps = 9999){
  
  if(missing(planar_projection)){
    message(
      'Argument to `planar_projection` is required. A suitable choice for all of North America is 5070.')
    }
  if(missing(returnProjected)){returnProjected <- FALSE}  
  
  if(returnProjected == TRUE){
    x <- sf::st_transform(x, planar_projection)
    # return the object in the original projection
  } else {
    orig_proj <- sf::st_crs(x)
    x <- sf::st_transform(x, planar_projection)
  }
  
  # now run the request numbed of iterations. 
  voronoiPolygons <- replicate(
    reps, 
    VoronoiSamplerEAS(x = x, kmeans_centers = kmeans_centers, reps = reps, 
                      pts = 100, n = n, planar_projection = planar_projection), 
    simplify = FALSE)
  
  # we use variance to determine the configuration of voronoi polygons which have
  # the most equally sized polygons. 
  variance <- unlist(lapply(voronoiPolygons, 
                             function(x){stats::var(as.numeric(sf::st_area(x))/10000)}))

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
    probs = 0.001, 
    level = 0.95) 
  
  # but we only select sets of records which actually meet the sample size requirements, 
  # a few weird ones always get through otherwise. 
  SelectedSample <- voronoiPolygons[which.min(variance)][[1]]
  
  if(returnProjected == FALSE){
    SelectedSample <- sf::st_transform(SelectedSample, orig_proj)
    # return the object in the original projection
  } 
  
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

#### Voronoi polygons can be formed many times in many ways, We want to run
## a number of iterations, and then determine the configuration which has #
## the smallest amount of variance in the geographic size of each cluster. #
VoronoiSamplerEAS <- function(x, kmeans_centers, reps, pts, n, planar_projection){
  
  # determine which portions of the STZ are likely to be populated by converting
  # the sdm output to vectors and masking the STZ to this area. 
  pts_n <- sf::st_sample(x, size = pts, type = 'regular', by_polygon = F)
  
  kmeans_res <- stats::kmeans(sf::st_coordinates(pts_n), centers = n)
  pts_n$ID <- kmeans_res$cluster
  
  # gather the geographic centers of the polygons. 
  kmeans_centers <- stats::setNames(
    data.frame(kmeans_res['centers'], seq_len(nrow(kmeans_res['centers'][[1]]))), 
    # use the centers as voronoi cores ... 
    c('X', 'Y', 'ID'))
  kmeans_centers <- sf::st_as_sf(kmeans_centers, coords = c('X', 'Y'), crs = planar_projection)
  
  voronoi_poly <- kmeans_centers |> # create polygons surrounding the 
    sf::st_transform(planar_projection) |> # clusters using voronoi polygons
    sf::st_union() |>  
    sf::st_voronoi() |>  
    sf::st_cast() |>  
    sf::st_as_sf() |>
    sf::st_make_valid() 
  
  lkp <- c(geometry = "x") 
  
  voronoi_poly <- sf::st_intersection(
    # reduce the extent of the voronoi polygons to the area of analysis. 
    voronoi_poly, sf::st_union(
      sf::st_transform(x, 
                       sf::st_crs(voronoi_poly)))
  ) |>  
    # we can assign an arbitrary number. 
    dplyr::mutate(ID = seq_len(nrow(voronoi_poly))) |>
    sf::st_make_valid() |>
    sf::st_as_sf() |>
    dplyr::rename(dplyr::any_of(lkp))
  
}
