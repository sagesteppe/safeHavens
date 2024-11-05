#' Create a quick SDM using elastic net regression
#' 
#' @description This function quickly creates a SDM using elastic net regression, and it will 
#' properly format all data for downstream use with the `safeHavens` workflow. 
#' Note that elastic net models are used for a couple very important reasons: 
#' they rescale all input independent variables before modelling, allowing us 
#' to combine the raw data with the beta coefficients to use for clustering 
#' algorithms downstream. They also allowing for 'shrinking' of terms from models
#' by shrinking terms from models we are able to get levels of ecological inference
#' prohibited by older model selection frameworks. 
#' @param x A (simple feature) sf data set of occurrence data for the species. 
#' @param predictors a terra 'rasterstack' of variables to serve as indepedent predictors. 
#' @param planar_projection Numeric, or character vector. An EPSG code, or a proj4 string, for a planar coordinate projection, in meters, for use with the function. For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise a continental scale projection like 5070 See https://projectionwizard.org/ for more information on CRS. The value is simply passed to sf::st_transform if you need to experiment. 
#' @param domain Numeric, how many times larger to make the entire domain of analysis than a simple bounding box around the occurrence data in `x`. 
#' @param quantile_v Numeric, this variable is used in thinning the input data, e.g. quantile = 0.05 will remove records within the lowest 5% of distance to each other iteratively, until all remaining records are further apart than this distance from each other. If you want essentially no thinning to happen just supply 0.01. Defaults to 0.025. 
#' @examples \dontrun{
#' 
#'  x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
#'  x <- x[,c('lon', 'lat')]
#'  x <- sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326)
#' 
#'  files <- list.files(
#'    path = file.path(system.file(package="dismo"), 'ex'), 
#'    pattern = 'grd',  full.names=TRUE )
#'  predictors <- terra::rast(files)
#' 
#' sdModel <- elasticSDM(
#'    x = x, predictors = predictors, quantile_v = 0.025,
#'    planar_projection =
#'      '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs')
#'      
#'  terra::plot(sdModel$RasterPredictions)
#' }
#' @returns A list of 12 objects, each of these subsequently used in the downstream  SDM Post processing sequence, or which we think are best written to disk. 
#' The actual model prediction on a raster surface are present in the first list 'RasterPredictions', the indepedent variables used in the final model are present in 'Predictors, and just the global PCNM/MEM raster surfaces are in 'PCNM'. 
#' The fit model is in 'Model', while the cross validation folds are stored in 'CVStructure', results from a single test/train partition in 'ConfusionMatrix', and the two data split in 'TrainData' and 'TestData' finally the 'PredictMatrix' which was used for classifying the test data for the confusion matrix. 
#' @export
elasticSDM <- function(x, predictors, planar_projection, domain, quantile_v){
  
  if(missing(quantile_v)){quantile_v <- 0.025}
  
  pts_plan <-   sf::st_transform(x, planar_projection)
  
  bb <- sf::st_bbox(pts_plan)
  buff_dist <- as.numeric( # here we get the mean distance of the XY distances of the bb
    ((bb[3] - bb[1]) + (bb[4] - bb[2])) / 2 
  ) / 2 # the mean distance * 0.25 is how much we will enlarge the area of analysis. 
  
  bb1 <- sf::st_union(pts_plan) |>
    sf::st_buffer(buff_dist) |>
    terra::vect() |>
    terra::project(terra::crs(predictors)) |>
    terra::ext()
  
  # THIS EXTENDING THE BB IS NOT DONE YET ##
  
  #######
  #######
  #######
  
  # Step 1 Select Background points - let's use SDM package `envidist` for this
  pa <- sdm::background(x = predictors, n = nrow(x), sp = x, method = 'eDist') |>
    dplyr::select(lon = x,  lat = y) |>
    sf::st_as_sf(coords = c('lon', 'lat'), crs = terra::crs(predictors)) |>
    dplyr::mutate(occurrence = 0)
  
  x$occurrence <- 1
  x <- dplyr::bind_rows(x, pa) |> # combine the presence and pseudoabsence points
    dplyr::mutate(occurrence = factor(occurrence))
  
  sp.coords <- data.frame(Species = 'Species', data.frame(sf::st_coordinates(x)))
  dists <- sf::st_distance(x[sf::st_nearest_feature(x), ], x, by_element = TRUE)
  thinD <- as.numeric(quantile(dists, c(quantile_v), na.rm = TRUE) / 1000) # ARGUMENT TO FN @PARAM 
  
  # Step 2 thin points to ensure there are not too many too close to each other. 
  # at worst these points are total duplicates 
  thinned <- spThin::thin(
    loc.data = sp.coords, thin.par = thinD,
    spec.col = 'Species',
    lat.col = 'Y', long.col = 'X', reps = 100, 
    locs.thinned.list.return = TRUE, 
    verbose = FALSE, 
    write.files = FALSE, 
    write.log.file = FALSE)
  
  thinned <- data.frame(thinned[ which.max(unlist(lapply(thinned, nrow)))]) |>
    sf::st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326)
  x <- x[lengths(sf::st_intersects(x, thinned))>0,]
  
  # Step 3 - Extract data to points for modelling
  x <- terra::extract(predictors, x, bind = TRUE) |>
    sf::st_as_sf() 
  
  # Step 4 - create a data split for testing the residuals of the glmnet model
  # It's not ideal to do a simple split of these data, because spatial autocorrelation
  # will mean that our results could be overly optimistic. 
  index <- unlist(caret::createDataPartition(x$occurrence, p=0.8)) # @ ARGUMENT TO FN @PARAM
  train <- x[index,]
  test <- x[-index,]

  # Develop CV folds for modelling
  indices_knndm <- CAST::knndm(train, predictors, k=5)
  
  # Recursive feature elimination using CAST developed folds
  train_dat <- sf::st_drop_geometry(train[, -which(names(train) %in% c("occurrence"))])
  ctrl <- caret::rfeControl(
    method = "LG0CV",
    repeats = 5,
    number = 10,
    functions = caret::lrFuncs,
    index = indices_knndm$indx_train,
    verbose = FALSE)
  
  lmProfile <- caret::rfe(
    method = 'glmnet',
    sizes = c(3:ncol(train_dat)), 
    x = train_dat,
    y = sf::st_drop_geometry(train)$occurrence,
    rfeControl = ctrl)
  
  # Step 4. fit model 
  
  # GLMNET converts all of it's input features in some fashion, where using 1/0 
  # for a boolean doesn't work so we need to make a quick intermediate here with the
  # Y/N. 
  train1 <- dplyr::mutate(
    train, # requires a YES OR NO or T/F just not anything numeric alike. 
    occurrence = dplyr::if_else(occurrence==1, 'YES', 'NO'))
  
  cv_model <- train( # we will do a quick pass through glmnet and drop
    # predictors which are totally shrunk out. 
    x = sf::st_drop_geometry(train1[,predictors(lmProfile)]), 
    sf::st_drop_geometry(train)$occurrence, 
    method = "glmnet", 
    family = 'binomial', 
    index = indices_knndm$indx_train) 
  
  sub <- train_dat[,predictors(lmProfile)]
  
  # now fit the model just using glmnet::glmnet in order that we can get the 
  # type of response for type='prob' rather than log odds or labelled classes
  # which we need to work with terra::predict. 
  mod <- glmnet::glmnet(
    x = sub, 
    sf::st_drop_geometry(train)$occurrence, 
    family = 'binomial', 
    keep = TRUE, 
    lambda = cv_model$bestTune$lambda, alpha = cv_model$bestTune$alpha
  )
  
  obs <- createPCNM_fitModel(
    x = train, 
    ctrl = ctrl, 
    indices_knndm = indices_knndm, 
    planar_proj = planar_projection, 
    sub = sub
    )
  
  mod <- obs$mod; pcnm <- obs$pcnm
  
  predictors <- c(predictors, pcnm) 
  # get the variables to extract from the rasters for creating a matrix for 
  # predictions, glmnet predict is kind of wonky and needs exact matrix dimensions. 
  vars <- rownames(coef(mod)); vars <- vars[2:length(vars)] 
  
  # now we need just the COORDINATES FOR TEST and will extract the data from
  # this set of predictors to them... 
  predict_mat <- predictors[[vars]]
  predict_mat <- as.matrix(
    terra::extract(predict_mat, test, ID = FALSE) 
  )
  
  cm <- caret::confusionMatrix(
    data = as.factor(predict(mod, newx = predict_mat, type = 'class')), 
    reference = test$occurrence,
    positive = "1")
  
  ## Predict our model onto a gridded surface (raster) ## This will allow for downstream
  # use with the rest of the safeHavens workflow. 
  preds <- predictors[[vars]]
  predfun <- function(model, data, ...){
    predict(model, newx=as.matrix(data), type = 'response')
  }
  
  rast_cont <- terra::predict(preds, model = mod, fun=predfun, na.rm=TRUE)
  
  return( # many objects were made in this function! Given that a sampling
    # schema may have a very long lifetime, it is likely best to save all of them. 
    list(
      RasterPredictions = rast_cont, 
      Predictors = predictors, 
      PCNM = pcnm,
      Model = mod,
      CVStructure = obs$cv_model, 
      ConfusionMatrix = cm, 
      TrainData = train,
      TestData = test,
      PredictMatrix = obs$pred_mat
    )
  )
  
}

