#' Sample a species based on Isolation by Geographic Distance (IBD)
#' 
#' @description Create `n` seed collection areas based on the distance geographic (great circle) distance between points. 
#' @param x a Raster surface to sample the points within, e.g. the output of SDM$Supplemented. 
#' @param n Numeric. the number of clusters desired. 
#' @param fixedClusters Boolean. Defaults to TRUE, which will create n clusters. If False then use NbClust::NbClust to determine the optimal number of clusters.
#' @param n_pts Numeric. the number of points to use for generating the clusters, these will be randomly sampled within the mask area `mask`. Defaults to 1000, which generally allows for enough points to be split for the KNN training.  
#' @param template Raster. A raster file which can be used as a template for plotting. 
#' @param prop_split Numeric. The proportion of records to be used for training the KNN classifier. Defaults to 
#' 0.8 to use 80% of records for training and 20% for the independent test sample.  
#' @param min.nc Numeric. Minimum number of clusters to test if fixedClusters=FALSE, defaults to 5. 
#' @param max.nc Numeric. Maximum number of clusters to test if fixedClusters=FALSE, defaults to 20. 
#' @param planar_proj Numeric. Optional. A planar projection to use for an sf::st_point_on_surface to ensuare valid spatial operations.
#' @examples 
#' planar_proj <- "+proj=laea +lat_0=-15 +lon_0=-60 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#' 
#' x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
#' x <- x[,c('lon', 'lat')]
#' x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)
#' x_buff <- sf::st_transform(x, planar_proj) |>
#'  sf::st_buffer(125000) |> # we are working in planar metric coordinates, we are
#'  sf::st_as_sfc() |> # buffer by this many / 1000 kilometers. 
#'  sf::st_union()
#'
#' files <- list.files( # note that for this process we need a raster rather than 
#'   path = file.path(system.file(package="dismo"), 'ex'), # vector data to accomplish
#'   pattern = 'grd',  full.names=TRUE ) # this we will 'rasterize' the vector using terra
#' predictors <- terra::rast(files)
#' # this can also be done using 'fasterize'. Whenever
#' # we rasterize a product, we will need to provide a template raster that our vector
#' # will inherit the cell size, coordinate system, etc. from 
#' 
#' x_buff.sf <- sf::st_as_sf(x_buff) |> 
#'   dplyr::mutate(Range = 1) |> 
#'   sf::st_transform(terra::crs(predictors))
#' 
#' # and here we specify the field/column with our variable we want to become 
#' # an attribute of our raster
#' v <- terra::rasterize(x_buff.sf, predictors, field = 'Range') 
#' 
#' # now we run the function demanding 20 areas to make accessions from, 
#' ibdbs <- IBDBasedSample(x = v, n = 20, fixedClusters = TRUE, 
#'    template = predictors, planar_proj = planar_proj)
#' plot(ibdbs[['Geometry']])
#' 
#' @return An simple features (sf) object containing the final grids for saving to computer. See the vignette for questions about saving the two main types of spatial data models (vector - used here, and raster). 
#' @export 
IBDBasedSample <- function(x, n, fixedClusters = TRUE, n_pts = 1000, template, prop_split = 0.8, 
  min.nc = 5, max.nc = 20, planar_proj){
  
  pts <- sf::st_sample(
    sf::st_as_sf(terra::as.polygons(x)), size = n_pts) |>
    sf::st_as_sf() |>
    sf::st_coordinates() |>
    as.matrix()

  pts_sf <- sf::st_as_sf(
    data.frame(pts),
    coords = c("X", "Y"),
    crs = sf::st_crs(x)
  )

  # calculate great circle distances between locations
  geoDist <- matrix(
    as.numeric(
      sf::st_distance(pts_sf, pts_sf)), 
       nrow = nrow(pts_sf), 
       ncol = nrow(pts_sf)
      )

  pts <- as.data.frame(pts)
  if(fixedClusters==TRUE){ # run the clustering processes. 
    
    geoDist_scaled <- stats::dist(base::scale(geoDist), method = 'euclidean') # scale variables
    clusters <- stats::hclust(geoDist_scaled, method = 'ward.D2')
    pts$ID <- stats::cutree(clusters, n)
    
  } else {

    geoDist_scaled <- stats::dist(base::scale(geoDist), method = 'euclidean') # scale variables
    NoClusters <- NbClust::NbClust(
      data = stats::as.dist(geoDist), diss = geoDist_scaled, 
      distance = NULL, min.nc = min.nc, max.nc = max.nc, 
      method = 'complete', index = 'silhouette'
    )
    pts$ID <- NoClusters$Best.partition
  }
  
  pts$ID <- factor(pts$ID)
  pts <- pts[stats::complete.cases(pts),]
  
  KNN <- trainKNN(pts, split_prop = prop_split)
  fit.knn <- KNN$fit.knn
  confusionMatrix <- fit.knn$confusionMatrix
  
  preds <- c(
    terra::init(template, 'x'), 
    terra::init(template, 'y')
  )
  names(preds) <- c('X', 'Y')
  spatialClusters <- terra::mask(preds, x)
  
  spatialClusters <- terra::predict(spatialClusters, model = fit.knn, na.rm = TRUE)
  spatialClusters <- terra::as.polygons(spatialClusters) |>
    sf::st_as_sf()
  
  # now number the grids in a uniform fashion
  if(!is.null(planar_proj)){ ## optional planar proj use for safety
    spatialClusters <- sf::st_transform(spatialClusters, planar_proj)
  }
  
  sf::st_agr(spatialClusters) = "constant"
  cents <- sf::st_point_on_surface(spatialClusters)

  cents <- sf::st_transform(cents, sf::st_crs(x))
  spatialClusters <- sf::st_transform(spatialClusters, sf::st_crs(x))

  cents <- cents |>
    dplyr::mutate(
      X = sf::st_coordinates(cents)[,1],
      Y = sf::st_coordinates(cents)[,2]
    ) |>
    dplyr::arrange(-Y, X) |>
    dplyr::mutate(ID = 1:dplyr::n()) |>
    dplyr::arrange(ID) |>
    dplyr::select(ID, geometry)
  
  sf::st_agr(cents) = "constant"
  ints <- unlist(sf::st_intersects(spatialClusters, cents))
  spatialClusters <- spatialClusters |>
    dplyr::mutate(ID = ints, .before = 1) |>
    dplyr::select(-class) |>
    dplyr::arrange(ID)
  
  list(Geometry = spatialClusters)
}
