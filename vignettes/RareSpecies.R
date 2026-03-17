## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

## ----required packages, message=FALSE-----------------------------------------
library(safeHavens)
library(ggplot2) ## plotting 
library(patchwork) ## multiplots
set.seed(99)

## ----import bradypus data for testing, warning = FALSE------------------------
x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)

## ----create map for visualizing results, echo = F, eval = T-------------------
planar_proj <- '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs'

americas <- spData::world

x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)
x_buff <- sf::st_transform(x, planar_proj) |>
  sf::st_buffer(125000) |> 
  sf::st_as_sfc() |> 
  sf::st_union() |>
  sf::st_transform(4326) |>
  sf::st_buffer(100000) |>
  sf::st_bbox()

americas <- spData::world 
americas <- sf::st_crop(americas, x_buff) |>
  dplyr::select(name_long)

bb <- sf::st_bbox(x_buff)

map <- ggplot() + 
  geom_sf(data = americas) + 
  theme(
    legend.position = 'none', 
    panel.background = element_rect(fill = "aliceblue"), 
    panel.grid.minor.x = element_line(colour = "red", linetype = 3, linewidth  = 0.5), 
    axis.ticks=element_blank(),
    axis.text=element_blank(),
    plot.background=element_rect(colour="steelblue"),
    plot.margin=grid::unit(c(0,0,0,0),"cm"),
    axis.ticks.length = unit(0, "pt"))+ 
  coord_sf(xlim = c(bb[1], bb[3]), ylim = c(bb[2], bb[4]), expand = FALSE)

rm(americas, bb, planar_proj, x_buff)

## ----prep bradypus for sampling-----------------------------------------------
n_sites <- nrow(x) 
df <- data.frame(
  site_id = seq_len(n_sites),
  required = FALSE,
  coord_uncertainty = 0, 
  lon = sf::st_coordinates(x)[,1], 
  lat = sf::st_coordinates(x)[,2]
)

knitr::kable(head(df))

## ----calculate distance matrix------------------------------------------------
dist_mat <- sapply(1:nrow(df), function(i) {
   greatCircleDistance(
     df$lat[i], df$lon[i],
     df$lat, df$lon
   )
 })

## ----required points----------------------------------------------------------
dists2c <- greatCircleDistance(
  median(df$lat), 
  median(df$lon), 
  df$lat, 
  df$lon
)
df[order(dists2c)[1],'required'] <- TRUE

## ----simulating coordinate uncertainty----------------------------------------
uncertain_sites <- sample(
  setdiff(seq_len(n_sites), 
  which(df$required)), 
  size = round(n_sites*0.2, 0)
  )
df$coord_uncertainty[uncertain_sites] <- runif(length(uncertain_sites), 1000, 40000) # meters

## ----combine the input data---------------------------------------------------
test_data <- list(
  distances = dist_mat,
  sites = df
  )

str(test_data)

## ----echo = F-----------------------------------------------------------------
rm(x, n_sites, uncertain_sites, dists2c)

## ----run with geographic distances, message=F---------------------------------
st <- system.time( {
    geo_res <- KMedoidsBasedSample( 
       ## reduce params from defaults
       ##  for a quick run. 

      input_data = test_data,
      n = 5,
      n_bootstrap = 10,
      dropout_prob = 0.1,
      n_local_search_iter = 10,
      n_restarts = 2
    )
  }
)

## ----show timing--------------------------------------------------------------
knitr::kable(st)

## ----echo=F-------------------------------------------------------------------
rm(st)

## ----structure of output------------------------------------------------------
str(geo_res)

## -----------------------------------------------------------------------------
knitr::kable(head(geo_res$stability_score))

## -----------------------------------------------------------------------------
knitr::kable(head(geo_res$stability))

## -----------------------------------------------------------------------------
knitr::kable(head(geo_res$input_data))

## -----------------------------------------------------------------------------
knitr::kable(head(geo_res$settings))

## ----first selection of target sites------------------------------------------
map + 
  geom_point(data = geo_res$input_data, 
  aes(
    x = lon, 
    y = lat, 
    shape = required, 
    size = cooccur_strength,
    color = selected
    )
  ) +
 # ggrepel::geom_label_repel(aes(label = site_id), size = 4) + 
  theme_minimal() + 
  labs(title = 'Priority Selection Status of Sites; Geographic Distances')

## ----priority ranking plot----------------------------------------------------
map + 
  geom_point(data = geo_res$input_data, 
    aes(
      x = lon, 
      y = lat, 
      shape = required, 
      size = -sample_rank,
      color = sample_rank
      )
    ) +
 # ggrepel::geom_label_repel(aes(label = sample_rank), size = 4) +
  theme_minimal()   

## ----prep environmental distance matrix---------------------------------------
files <- list.files(
  path = file.path(system.file(package="dismo"), 'ex'), 
  pattern = 'grd',  full.names=TRUE )
predictors <- terra::rast(files) # import the independent variables
rm(files)

## ----run with geographic and environmental distances--------------------------
pts <- terra::spatSample(predictors, 100, na.rm = TRUE)
pts <- pts[, names(pts)!='biome' ] # remove categorical variable for distance calc

pca_results <- stats::prcomp(pts, scale = TRUE)
round(pca_results$sdev^2 / sum(pca_results$sdev^2), 2) # variance explained
pca_raster <- terra::predict(predictors, pca_results)

## ----plot pca layers----------------------------------------------------------
terra::plot(terra::subset(pca_raster, c(1:2))) # prediction of the pca onto a new raster

## ----echo = F-----------------------------------------------------------------
rm(pts, predictors, pca_results)

## ----extract environmental values and calculate distance matrix---------------
env_values <- terra::extract(pca_raster, 
  sf::st_coordinates(
    sf::st_as_sf(
      df, 
      coords = c('lon', 'lat'), 
      crs = 4326
    )
  )
)[,1:2]
plot(env_values, main = 'environmental distance of points from first two PCA axis')

## ----convert pca distances to a matrix----------------------------------------
env_dist_mat <- as.matrix(
    dist(env_values)
  )

## ----echo=F-------------------------------------------------------------------
rm(pca_raster)

## ----run with environmental distance------------------------------------------
test_data <- list(
  distances = env_dist_mat,
  sites = df
  )

st <- system.time( 
  {
    env_res <- KMedoidsBasedSample(  ## reduce some parameters for shorter run time.
      input_data = test_data,
      n = 5,
      n_bootstrap = 10,
      dropout_prob = 0.1,
      n_local_search_iter = 50,
      n_restarts = 2
    )
  }
)

rm(dist_mat, env_dist_mat)

## -----------------------------------------------------------------------------
knitr::kable(st)

## -----------------------------------------------------------------------------
knitr::kable(head(env_res$stability_score))

## -----------------------------------------------------------------------------
map + 
  geom_point(data = env_res$input_data, 
    aes(
      x = lon, 
      y = lat, 
      shape = required, 
      size = cooccur_strength,
      color = selected
      )
    ) +
 # ggrepel::geom_label_repel(aes(label = site_id), size = 4) + 
  theme_minimal() + 
  labs(title = 'Priority Selection Status of Sites; Environmental')

## ----required point method - population centroid - heat map-------------------
dens <- with(df, MASS::kde2d(lon, lat, n = 200))
max_idx <- which(dens$z == max(dens$z), arr.ind = TRUE)[1,]
max_point <- c(dens$x[max_idx[1]], dens$y[max_idx[2]])

pops_centre <- sweep(df[c('lon', 'lat')], 2, max_point, "-")
pop_centered_id <- which.min(rowSums(abs(pops_centre^2)))

rm(dens, max_idx, max_point, pops_centre)

## ----required point method - environmental centroid---------------------------
env_centered <- sweep(env_values, 2, sapply(env_values, median), "-")
env_centered_id <- which.min(rowSums(abs(env_centered^2)))

rm(env_values)

## ----plot alternative centroids-----------------------------------------------
# geographic centroid was pt 47
centers <- df[ c(env_centered_id, pop_centered_id, 47), ] 
centers$type <- c('Environmental', 'Population', 'Geographic')

map +
  geom_point(
    data = df, 
    aes(x = lon, y = lat)
    ) + 
  geom_point(
    data = centers,  
    aes(x = lon, y = lat),
    col = '#FF1493', size = 4
    ) + 
 # ggrepel::geom_label_repel(
 #   data = centers, 
 #   aes(label = type, x = lon, y = lat)
 #   ) + 
  theme_minimal() + 
  labs(title = 'Possbilities for centers')

## ----echo=F-------------------------------------------------------------------
rm(env_centered_id, env_centered, pop_centered_id)

