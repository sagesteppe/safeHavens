## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

## ----load packaeges, message=FALSE, warning=FALSE-----------------------------
library(safeHavens)

## ----Prepare data for a Species Distribution Model----------------------------
x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)

planar_proj <- 3857 # Web Mercator for planar distance calcs

files <- list.files(
  path = file.path(system.file(package="dismo"), 'ex'), 
  pattern = 'grd',  full.names=TRUE )
predictors <- terra::rast(files) # import the independent variables

## ----Create Species Distribution Model, message=F, warning=F------------------
sdModel <- elasticSDM(
  x = x,
  predictors = predictors,
  quantile_v = 0.025,
  planar_proj = planar_proj
  )

## ----save sdmodel for getting started vignette, echo =F, eval=F---------------
# #sdModel$RasterPredictions <- terra::wrap(sdModel$RasterPredictions)
# #sdModel$Predictors <- terra::wrap(sdModel$Predictors)
# #sdModel$PCNM <- terra::wrap(sdModel$PCNM)
# ### this is only for vignette setup!!!!!!!!!!!!!!!!!!!
# #setwd('~/Documents/assoRted/safeHavens/vignettes')
# #saveRDS(sdModel, '../inst/extdata/sdModel.rds')
# 
# sdModel$RasterPredictions <- terra::unwrap(sdModel$RasterPredictions)
# sdModel$Predictors <- terra::unwrap(sdModel$Predictors)
# sdModel$PCNM <- terra::unwrap(sdModel$PCNM)

## ----Explore SDM output - Confusion Matrix------------------------------------
sdModel$ConfusionMatrix

## ----Explore SDM output - Map of Predictions----------------------------------
terra::plot(sdModel$RasterPredictions)

## ----Threshold the SDM output-------------------------------------------------
threshold_rasts <- PostProcessSDM(
  rast_cont = sdModel$RasterPredictions, 
  test = sdModel$TestData,
  train = sdModel$TrainData,
  planar_proj = planar_proj,
  thresh_metric = 'sensitivity', 
  quant_amt = 0.5
  )

## ----Compare Threshold Results------------------------------------------------
terra::plot(threshold_rasts$FinalRasters)

## ----save the ouputs for vignettes, echo = F, eval = F------------------------
# #terra::writeRaster(threshold_rasts$FinalRasters, filename = file.path('..', 'inst', 'extdata', 'SDM_thresholds.tif'))

## ----Rescale Predictor Variables----------------------------------------------
rr <- RescaleRasters(
  model = sdModel$Model,
  predictors = sdModel$Predictors, 
  training_data = sdModel$TrainData, 
  pred_mat = sdModel$PredictMatrix)

terra::plot(rr$RescaledPredictors)

## ----Save SDM results, eval = F-----------------------------------------------
# bp <- '~/Documents/assoRted/StrategizingGermplasmCollections'
# 
# writeSDMresults(
#   cv_model = sdModel$CVStructure,
#   pcnm = sdModel$PCNM,
#   model = sdModel$Model,
#   cm = sdModel$ConfusionMatrix,
#   coef_tab = rr$BetaCoefficients,
#   f_rasts = threshold_rasts$FinalRasters,
#   thresh = threshold_rasts$Threshold,
#   file.path(bp, 'results', 'SDM'), 'Bradypus_test')
# 
# # we can see that the files were placed here using this.
# list.files( file.path(bp, 'results', 'SDM'), recursive = TRUE )

## ----Clean up SDM variables, warning=FALSE, echo = FALSE----------------------
rm(rr, predictors, files, sdModel, bp)

