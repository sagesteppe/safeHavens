library(tidyverse)
library(sf)

#' place random points in the polygon which will be dissolved with the larger polygons
#' 
#' This function is ran within `GridBasedSample` it will place points throughout polygon geometries
#' which should be merged to larger polygons and assign them to neighboring polygons
#' based on how much area we want to grow these polygons too. 
#' @param neighb_grid 
#' @param focal_grid 
#' @param props 
#' @param nf_pct 
assignGrid_pts <- function(neighb_grid, focal_grid, props, nf_pct){
  
  # Ensure that we are only working on a grid with 100 points
  pts <- sf::st_sample(focal_grid, size = 100, type = 'regular') |>
    sf::st_as_sf() |> 
    dplyr::mutate(ID = 1:dplyr::n())
  
  samp <- 100
  while (nrow(pts) < 100) {
    if (nrow(pts) < 100) {
      pts <- sf::st_sample(focal_grid, size = samp, type = 'regular') |>
        sf::st_as_sf() |> 
        mutate(ID = 1:n())
    }
    samp <- samp + 1
  } 
  pts <- pts[sample(1:nrow(pts), size = 100, replace = F), ]
  rm(samp)
  
  if(nrow(neighb_grid)==1){
    pts$Assigned <- sf::st_drop_geometry(neighb_grid$ID)} else {
      
      # identify the nearest neighbor which the points can be assigned to. 
      # we use these to determine what are the 'neediest' neighbors and assign them 
      # points first 
      pts$nf <- sf::st_nearest_feature(pts, neighb_grid) 
      pts$nf <- sf::st_drop_geometry(neighb_grid)[pts$nf,'ID']
      
      # each grid now receives either the maximum number of nearest neighbors if
      # the desired proportion is lower than the existing nearest neighbors, 
      # or the proportion of nearest points meeting the desired proportion
      dists <- sf::st_distance(pts, neighb_grid)
      dists <- data.frame(apply(dists, 2, as.numeric))
      colnames(dists) <- neighb_grid$ID
      
      # if only one grid remains, assign all remaining points to it. 
      frst_assignments <- assign_pts_frst(dists, props = props, nf_pct)
      if(exists('frst_assignments')){pts$Assigned <- frst_assignments$Assignment}
      rm(nf_pct, dists)
      
      if(any(is.na(frst_assignments$Assignment))){
        
        needAssigned <- pts[is.na(pts$Assigned),]
        dist_final_pts <- data.frame(
          matrix(
            t(sf::st_distance(needAssigned, pts)), ncol = nrow(needAssigned)
          )
        )
        dist_final_pts[dist_final_pts==0] <- NA
        dist_final_pts$ID <- pts$ID
        
        ob <- vector(mode = 'list', length = ncol(dist_final_pts)-1)
        for (i in 1:ncol(dist_final_pts)){
          ob[[i]] <- dist_final_pts[order(dist_final_pts[,i]),  'ID']
        }
      }
      
      rm(frst_assignments)
      # a few points may remain in the center of the object. 
      # We will now try to assign all of these to the groups which do
      # have not yet come close to adequate representation. 
      
      # we will determine what the closest neighbors to the unassigned points are. 
      # to grab them we will buffer the unassigned point and pull out the close neighbors. 
      if(exists('needAssigned')){
        nn <- spdep::knearneigh(pts, k=4)[['nn']][needAssigned$ID,]
        
        for (i in 1:nrow(needAssigned)){
          needAssigned[i, 'Assigned'] <- 
            names(
              which.max(
                table(
                  sf::st_drop_geometry(pts)[
                    unlist(
                      if(length(nn)==4){nn[1:4]} else {nn[i,1:4]}),
                    'Assigned'])
              )
            ) 
        }
        
        pts <- dplyr::filter(pts, ! ID %in% needAssigned$ID) |>
          dplyr::rename(geometry = x) |>
          dplyr::bind_rows(needAssigned) |>
          dplyr::select(Assigned, ID, geometry = x)
      }
      pts <- pts[! sf::st_is_empty(pts), ]
      
      # Determine if there are points which are 'disconnected' from their remaining neighbors
      nn <- spdep::knearneigh(pts, k=4)[['nn']]
      
      indices <- vector(mode = 'list', length = nrow(pts))
      neighs <- vector(mode = 'list', length = nrow(pts))
      focal <- vector(mode = 'list', length = nrow(pts))
      matches <- vector( length = nrow(pts))
      for (i in 1:nrow(pts)){
        
        indices[[i]] <- nn[i,]
        neighs[[i]] <- 
          names(table(sf::st_drop_geometry(pts)[indices[[i]],'Assigned']))
        focal[[i]] <- sf::st_drop_geometry(pts)[i,'Assigned']
        matches[i] <- focal[[i]] %in% neighs[[i]]
        
      }
      
      rm(focal, indices, nn)
      # if a point is disconnected from it's remaining neighbors then assign it to the neighbor which needs more points to reach the ideal sample ratio for it's grid class.
      if(any(matches==F)){
        indices <- neighs[which(matches==F)] 
        
        realized <- table(pts$Assigned)/ nrow(pts) *100
        props [!names(props) %in% names(realized)] <- 0 # if a small grid is missing say so
        diff <- realized - props # the first entry below will gain the point. 
        
        for (i in 1:length(indices)){
          pts[i,'Assigned'] <-
            names(sort(diff[names(diff) %in% unlist(neighs[i])], decreasing = FALSE)[1])
        }
  }
    }
  lkup <- c(geometry = "x")
  
  pts <- pts |>
    dplyr::rename(dplyr::any_of(lkup)) |>
    dplyr::select(Assigned, geometry) |>
    dplyr::mutate(Assigned = as.numeric(Assigned))
  return(pts)
}


# NEED TO FIX SEVERAL ITEMS

# 2) Error in spdep::knearneigh(pts, k = 4)[["nn"]][needAssigned$ID, ] : 
# subscript out of bounds

##################################### SAND BOX ###############################

target <- tigris::states() |> 
  dplyr::filter(NAME == 'Maryland') |>
  sf::st_transform(5070)

TestGridSizes(target)
output <- GridBasedSample(target)

ggplot() + 
  geom_sf(data = target) + 
  geom_sf(data = output, aes(fill = Assigned))

