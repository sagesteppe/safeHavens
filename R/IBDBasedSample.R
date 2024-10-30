#' Sample a species based on Isolation by Geographic Distance (IBD)
#' 
#' Create n seed collection areas based on the distance geographic (great circle)
#' distance between points. 
#' @param x a Raster surface to sample the points within, e.g. the output of SDM$Supplemented. 
#' @param n Numeric. the number of clusters desired. 
#' @param fixedClusters Boolean. Defaults to TRUE, which will create n clusters. If False then use NbClust::NbClust to determine the optimal number of clusters.
#' @param n_pts Numeric. the number of points to use for generating the clusters, these will be randomly sampled within the mask area `mask`. Defaults to 1000, which generally allows for enough points to be split for the KNN training.  
#' @param prop_split Numeric. The proportion of records to be used for training the KNN classifier. Defaults to 
#' 0.8 to use 80% of records for training and 20% for the independent test sample.  
#' @param min.nc Numeric. Minimum number of clusters to test if fixedClusters=FALSE, defaults to 5. 
#' @param max.nc Numeric. Maximum number of clusters to test if fixedClusters=FALSE, defaults to 20. 
#' @example
#' @return An simple features (sf) object containing the final grids for saving to computer. See the vignette for questions about saving the two main types of spatial data models (vector - used here, and raster). 
#' @export 
IBDBasedSample <- function(x, n, fixedClusters, n_pts, prop_split, min.nc, max.nc){
  
  if(missing(fixedClusters)){fixedClusters <- TRUE}
  if(missing(prop_split)){prop_split <- 0.8}
  if(missing(n_pts)){n_pts <- 1000}
  
  pts <- sf::st_sample(
    sf::st_as_sf(terra::as.polygons(x)), size = n_pts) |>
    sf::st_as_sf() |>
    sf::st_coordinates() |>
    as.matrix()
  
  geoDist <- pointDistance(
    pts, 
    longlat = TRUE) # calculate great circle distances between locations
  
  pts <- as.data.frame(pts)
  if(fixedClusters==TRUE){ # run the clustering processes. 
    
    geoDist_scaled <- stats::dist(scale(geoDist), method = 'euclidean') # scale variables
    clusters <- stats::hclust(geoDist_scaled, method = 'ward.D2')
    pts$ID <- stats::cutree(clusters, n)
    
  } else {
    
    if(missing(min.nc)){min.nc <- 5}
    if(missing(max.nc)){max.nc <- 5}
    
    geoDist_scaled <- stats::dist(scale(geoDist), method = 'euclidean') # scale variables
    NoClusters <- NbClust::NbClust(
      data = stats::as.dist(geoDist), diss = geoDist_scaled, 
      distance = NULL, min.nc = min.nc, max.nc = max.nc, 
      method = 'complete', index = 'silhouette'
    )
    pts$ID <- NoClusters$Best.partition
  }
  
  pts$ID <- factor(pts$ID)
  pts <- pts[complete.cases(pts),]
  
  KNN <- trainKNN(pts, split_prop = prop_split)
  fit.knn <- KNN$fit.knn
  confusionMatrix <- fit.knn$confusionMatrix
  
  preds <- c(
    terra::init(template$Supplemented, 'x'), 
    terra::init(template$Supplemented, 'y')
  )
  names(preds) <- c('X', 'Y')
  spatialClusters <- terra::mask(preds, x)
  
  spatialClusters <- terra::predict(spatialClusters, model = fit.knn, na.rm = TRUE)
  
  # LET'S RETURN AN SF OBJECT INSTEAD !!!!!!!
  return(spatialClusters)
}

out <- IBDBasedSample(template$Supplemented, n = 20, fixedClusters = TRUE)
plot(out)
