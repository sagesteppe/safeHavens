## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

## ----attach packages, message = F, warning = F--------------------------------
library(safeHavens)
library(sf) ## vector operations
library(terra) ## raster operations
library(geodata) ## environmental variables
library(dplyr) ## general data handling
library(tidyr) ## general data handling
library(ggplot2) ## plotting 
library(patchwork) ## multiplots
set.seed(22)

## ----Download species occurrence data-----------------------------------------
cols = c('decimalLatitude', 'decimalLongitude', 'dateIdentified', 'species', 'acceptedScientificName', 'datasetName', 
  'coordinateUncertaintyInMeters', 'basisOfRecord', 'institutionCode', 'catalogNumber')

## download species data using scientificName, can use keys and lookup tables for automating many taxa. 
hemi <- rgbif::occ_search(scientificName = "Helianthella microcephala")
hemi <- hemi[['data']][,cols]  |>
  drop_na(decimalLatitude, decimalLongitude) |> # any missing coords need dropped. 
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) |> # no dupes can be present
  st_as_sf(coords = c( 'decimalLongitude', 'decimalLatitude'), crs = 4326, remove = F) 

western_states <- spData::us_states |> ## for making a quick basemap. 
  dplyr::filter(NAME %in% 
    c('Utah', 'Arizona', 'Colorado', 'New Mexico', 'Wyoming', 'Nevada', 'Idaho', 'California')) |>
  dplyr::select(NAME, geometry) |>
  st_transform(4326)

bb <- st_bbox(
  c(
    xmin = -116, 
    xmax = -105, 
    ymax = 44, 
    ymin = 33.5),
    crs = st_crs(4326)
    )

western_states <- st_crop(western_states, bb)

bmap <- ggplot() + 
    geom_sf(data = western_states) + 
    geom_sf(data = hemi) +
    theme_minimal() +
    coord_sf(
        xlim = c(bb[['xmin']], bb[['xmax']]), 
        ylim = c(bb[['ymin']], bb[['ymax']])
        )  +
    theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
      )

bmap +
    labs(title = 'Helianthella microcephala\noccurrence records')

## ----remove cols and western states, echo = F---------------------------------
rm(cols, western_states)

## ----download data------------------------------------------------------------
# Download WorldClim bioclim at ~10 km
bio_current <- worldclim_global(var="bioc", res=2.5)
bio_future <- cmip6_world(
  model = "CNRM-CM6-1", ## modelling method
  ssp   = "245", ## "Middle of the Road" scenario
  time  = "2041-2060", # time period
  var   = "bioc", # just use the bioclim variables
  res   = 2.5
)

# Crop to domain - use a large BB to accomodate range shift
# under future time points. 
# but going too large will dilute the absence records
bbox <- ext(bb)

bio_current <- crop(bio_current, bbox)
bio_future <- crop(bio_future, bbox)

## ----standardize raster layer names-------------------------------------------
simplify_names <- function(x){
    paste0('bio_', sprintf("%02d", as.numeric(gsub('\\D+','', names(x)))))
}

names(bio_current) <- gsub('^.*5m_', '', names(bio_current))
names(bio_future) <- gsub('^.*2060_', '', names(bio_future))

names(bio_current) <- simplify_names(bio_current)
names(bio_future) <- simplify_names(bio_future)

# TRUE means all names match, and are in the same position. 
all(names(bio_current) == names(bio_future))

## we will also drop some variables that are essentially identical in the study area
drops <- c(
  'bio_05' , 'bio_02',
  'bio_11', 'bio_12'
  ) 

bio_current <- subset(bio_current, negate = TRUE, drops)
bio_future <- subset(bio_future, negate = TRUE, drops)

## ----diff the raster layers to show areas of difference-----------------------
pct_diff <- function(x, y){((x - y)/((x + y)/2)) * 100}
difference <- pct_diff(bio_current, bio_future)
plot(difference)

## ----remove diff and simplify_names, echo = F---------------------------------
rm(difference, simplify_names, drops)

## ----fit elasticSDM for species,message = FALSE-------------------------------
## note we subset to just the geometries for the prediction records. 
# this is because we feed in all columns to the elastic net model internally. 
hemi <- select(hemi, geometry) 

# as always verify our coordinate reference systems (crs) match. 
st_crs(bio_current) == st_crs(hemi) 

# and fit the elasticnet model 
eSDM_model <- elasticSDM(
    x = hemi, 
    predictors = bio_current, 
    planar_projection = 5070, 
    PCNM = FALSE ## set to FALSE for this workstream!!!!!
    )

# Test with just 3 runs to see variability quickly
run1 <- elasticSDM(hemi, bio_current, 5070)
run2 <- elasticSDM(hemi, bio_current, 5070)
run3 <- elasticSDM(hemi, bio_current, 5070)

# Compare coefficients
coef1 <- as.matrix(coef(run1$Model))
coef2 <- as.matrix(coef(run2$Model))
coef3 <- as.matrix(coef(run3$Model))

# Quick comparison
all_vars <- unique(c(rownames(coef1), rownames(coef2), rownames(coef3)))
compare <- data.frame(
  var = all_vars,
  run1 = coef1[match(all_vars, rownames(coef1)), 1],
  run2 = coef2[match(all_vars, rownames(coef2)), 1],
  run3 = coef3[match(all_vars, rownames(coef3)), 1]
)
compare[is.na(compare)] <- 0
compare$cv <- apply(compare[,-1], 1, function(x) sd(x)/mean(x))
print(compare)

## ----current sdm threshold results--------------------------------------------
knitr::kable(eSDM_model$ConfusionMatrix$byClass[
  c('Sensitivity', 'Specificity', 'Recall', 'Balanced Accuracy')])
plot(eSDM_model$RasterPredictions)

## ----rescale rasters----------------------------------------------------------
bio_current_rs <- RescaleRasters(
    model = eSDM_model$Model, 
    predictors = eSDM_model$Predictors,
    training_data = eSDM_model$Train,
    pred_mat = eSDM_model$PredictMatrix
    )

plot(bio_current_rs$RescaledPredictors)

bio_future_rs <- rescaleFuture(
  eSDM_model$Model, 
  bio_future, 
  eSDM_model$Predictors,
  training_data = eSDM_model$Train,
  pred_mat = eSDM_model$PredictMatrix
)

plot(bio_future_rs)

## ----rescaled raster difference-----------------------------------------------
difference_rs <- pct_diff(bio_current_rs$RescaledPredictors, bio_future_rs)
plot(difference_rs)

## ----remove rs differnce surfaces, echo = F-----------------------------------
rm(difference_rs, bio_future_rs)

## ----identify current clusters, message = FALSE, warning = FALSE--------------
threshold_rasts <- PostProcessSDM(
  rast_cont = eSDM_model$RasterPredictions, 
  test = eSDM_model$TestData,
  train = eSDM_model$TrainData,
  planar_proj = 5070,
  thresh_metric =
     'equal_sens_spec', 
  quant_amt = 0.25
  )

plot(threshold_rasts$FinalRasters)
knitr::kable(threshold_rasts$Threshold$equal_sens_spec)

bmap + 
  geom_sf(data = 
  sf::st_as_sf(
    terra::as.polygons(
      threshold_rasts$FinalRasters['Threshold'])
      ), fill = 'cornsilk'
    ) + 
    geom_sf(data = hemi) 


## ----EnvironmentalBasedSample, message=F, warning = F, results='hide'---------
ENVIbs <- EnvironmentalBasedSample(
  pred_rescale = bio_current_rs$RescaledPredictors, 
  write2disk = FALSE, # we are not writing, but showing how to provide some arguments
  f_rasts = threshold_rasts$FinalRasters, 
  coord_wt = 0.001, 
  fixedClusters = FALSE,
  lyr = 'Threshold',
  n_pts = 500, 
  planar_proj = "epsg:5070",
  buffer_d = 3,
  prop_split = 0.8,
  min.nc = 5, 
  max.nc = 15
  )

## ----message=F,warning=F------------------------------------------------------
bmap + 
  geom_sf(data = ENVIbs$Geometry, aes(fill = factor(ID))) + 
  geom_sf(data = hemi) + 
  theme(legend.position = 'bottom') + 
  labs(fill = 'Cluster', title = 'Current')

## ----classify future climate surface with clusters, warning = FALSE, message = FALSE----
future_clusts <- projectClusters(
  eSDM_object = eSDM_model, 
  current_clusters = ENVIbs,
  future_predictors = bio_future,
  current_predictors = bio_current,
  thresholds = threshold_rasts$Threshold, 
  planar_proj = "epsg:5070", 
  thresh_metric = 'equal_sens_spec',
  n_sample_per_cluster = 20
)

## ----plot the mess surfaces,message=F,warning=F-------------------------------
plot(future_clusts$mess)
bmap +
  geom_sf(data = 
    st_as_sf(terra::as.polygons(future_clusts$novel_mask)), 
    fill = 'red') + 
  labs(title = 'MESS regions')

## ----identify additional clusters under future climate, message=F,warning=F----
current = bmap + 
  geom_sf(data = ENVIbs$Geometry, aes(fill = factor(ID))) +
  labs(title = 'Current', fill = 'Cluster') +
  theme(legend.position = "none")

future = bmap + 
  geom_sf(data = future_clusts$clusters_sf, aes(fill = factor(ID))) + 
  labs(title = '2041-2060', fill = 'Cluster') 

current + future

## ----determine analog current climate clusters for the novel groups-----------
knitr::kable(future_clusts$novel_similarity)

## ----cluster changes----------------------------------------------------------
future_clusts$changes

