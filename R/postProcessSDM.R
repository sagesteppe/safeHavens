#' Modify the output rasters from the SDM process to better match target goals
#' 
#' This is the last 'analytical' portion of the SD modelling process. It will
#' produce binary (Yes/ NA) rasters of species suitable habitat based on a three step
#' process. The first step uses dismo::thresholds to determine what feature of the
#' raster we want to maximuize, in this case I want a raster to be *more* than less
#' likely to capture a presence than omit it. Then we can subset the predicted habitat, which
#' has not been dispersed too by comparing nearest neighbor distances of the observed points. 
#' Finally, we can add back in areas to the raster where we know that the species has been 
#' observed, hopefully these are not missed in the original SDM, but there are always
#' some suspicious points which are difficult to fit in light of the inertia 
#' of the rest of the species. 
#' @param x the raw unaltered (except masked) raster predictions. 
#' @param thresh_metric ?dismo::threshold for all options, defaults to 'sensitivity'
#' @param quan_amt the quantile of nearest neighbors distance to use for steps 2 and 3. 
#' defaults to 0.25, using the median nearest neighbor distance of 10 bootstrapping replicates for
#' estimating a buffer to restrict the SDM surface too, and the minimum of the 10 bootstrap reps
#' for adding surface to presence points which were not placed in binary suitable habitat. 

postProcessSDM <- function(rast_cont, thresh_metric, quant_amt){
  
  if(missing(thresh_metric)){thresh_metric <- 'sensitivity'}
  if(missing(quant_amt)){quant_amt <- 0.25}
  
  # determine a threshold for creating a binomial map of the species distribution
  # we want to predict MORE habitat than exists, so we want to maximize sensitivity
  # in our classification. 
  
  test.sf <- sf::st_as_sf(test, coords = c('x', 'y'), crs = 4326) |>
    dplyr::select(occurrence)
  
  test.sf <- terra::extract(rast_cont, test.sf, bind = TRUE) |>
    sf::st_as_sf() |>
    sf::st_drop_geometry() 
  
  eval_ob <- dismo::evaluate(
    p = test.sf[test.sf$occurrence==1,'s0'],
    a = test.sf[test.sf$occurrence==0,'s0']
  )
  thresh <- dismo::threshold(eval_ob)
  cut <- thresh[[thresh_metric]] # ARGUMENT TO FN @PARAM 
  
  m <- matrix( # use this to reclassiy data to a binary raster
    c( # but more simply, turn the absences into NA for easier masking later on? 
      0, cut, NA,
      cut, 1, 1), 
    ncol = 3, byrow= TRUE)
  
  rast_binary <- terra::classify(rast_cont, m) # create a YES/NO raster
  
  rm(eval_ob, cut, m, test.sf)
  # use sf::st_buffer() to only keep habitat within XXX distance from known populations
  # we'll use another set of cv-folds based on all occurrence data 
  # Essentially, we will see how far the nearest neighbor is from each point in each
  # fold
  
  nn_distribution <- function(x, y){
    ob <- unlist(x)
    
    nf <- sf::st_distance(
      y[sf::st_nearest_feature(y[ob, ]), ],
      y[ob, ], by_element = TRUE
    )
  }
  
  pres <- x[ x$occurrence==1, ]
  pres <- sf::st_transform(
    pres, 
    '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs')
  indices_knndm <- CAST::knndm(pres, predictors, k=10)
  
  nn_dist <- lapply(indices_knndm[['indx_train']], nn_distribution, y = pres)
  dists <- unlist(list(lapply(nn_dist, quantile, quant_amt))) 
  
  within_dist <- sf::st_buffer(pres, median(dists)) |>
    dplyr::summarize(geometry = sf::st_union(geometry)) |>
    sf::st_simplify() |>
    sf::st_transform(terra::crs(rast_binary)) |>
    terra::vect()
  
  rast_clipped <- terra::mask(rast_binary, within_dist)
  
  rm(nn_dist, indices_knndm, nn_distribution, within_dist)
  ####### IF WE HAVE POINTS WHICH ARE FLOATING IN SPACE - I.E. POINTS W/O  
  # SUITABLE HABITAT MARKED, THEN LET'S ADD the same amount of suitable habitat 
  # to each of them, that was used as the buffer for clipping suitable habitat to the
  # points above. 
  
  pres <- sf::st_transform(pres, terra::crs(rast_binary))
  outside_binary <- terra::extract(rast_binary, pres, bind = TRUE) |>
    sf::st_as_sf() |>
    dplyr::filter(is.na(s0)) |>
    sf::st_transform(
      '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs') |>
    sf::st_buffer(min(dists)) |>
    dplyr::summarize(geometry = sf::st_union(geometry)) |>
    dplyr::mutate(occurrence = 1) |>
    terra::vect() |>
    terra::project(terra::crs(rast_binary)) |>
    terra::rasterize( rast_binary, field = 'occurrence')
  
  rast_clipped_supplemented <- max(rast_clipped, outside_binary, na.rm = TRUE)
  
  ##########   COMBINE ALL RASTERS TOGETHER FOR A FINAL PRODUCT      #############
  f_rasts <- c(rast_cont, rast_binary, rast_clipped, rast_clipped_supplemented)
  names(f_rasts) <- c('Predictions', 'Threshold', 'Clipped', 'Supplemented')
  
  return(list(f_rasts = f_rasts, thresh = thresh))
}
  
