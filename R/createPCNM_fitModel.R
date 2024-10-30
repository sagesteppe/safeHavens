#' Create global/regional PCNM surfaces and fit elastic net regression to all covariates
#' 
#' This function does most of the lifting in the SDM workflow. It will create PCNM/MEM
#' surfaces and subset them to the local/global maps. It can then use Thin plate regression
#' to predict those onto an actual raster surface which can be used for prediction
#' downstream. 
#' 
#' It will then use cross validation to determine a suitable glmnet model alpha and
#' lambda, and fit them using glmnet. It returns three objects which are spit 
#' out into the environment, 1) pcnm, surfaces for only those eigenvectors used in 
#' the glmnet model (including if shrunk out), 2) the glmnet model 3) all fitting
#' information from carets process. 
#' @param x should be the training data as an sf/tibble/dataframe
#' @param planar_proj Numeric, or character vector. An EPSG code, or a proj4 string, for a planar coordinate projection, in meters, for use with the function. For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise a continental scale projection like 5070 See https://projectionwizard.org/ for more information on CRS. The value is simply passed to sf::st_transform if you need to experiment. 
createPCNM_fitModel <- function(x, planar_proj){
  
  train_planar <- sf::st_transform(x, planar_proj) 
  
  dis <- sf::st_distance(train_planar)
  dis <- apply(dis, 2, as.numeric)
  xypcnm <- vegan::pcnm(dis)
  xypcnm.df <- data.frame(xypcnm$vectors)[,1:20]
  
  pcnmProfile <- caret::rfe( 
    method = 'glmnet', # https://doi.org/10.1111/gean.12054
    x = xypcnm.df, 
    sizes = 1:5,
    sf::st_drop_geometry(x)$occurrence,
    rfeControl = ctrl, 
    index = indices_knndm$indx_train
  )
  
  xypcnm.df <- xypcnm.df[,caret::predictors(pcnmProfile)]
  preds <- cbind(sub, xypcnm.df)
  
  if(is.numeric(xypcnm.df)){
    colnames(preds)[length(preds)] <- caret::predictors(pcnmProfile)}
  
  rm(xypcnm, dis, train_planar)
  # trying to refit the glmnet
  
  cv_model <- caret::train( # let's extract the model performance for the top alpha/lambda info here. 
    x = preds, 
    sf::st_drop_geometry(x)$occurrence, 
    method = "glmnet", 
    metric = 'Accuracy',
    family = 'binomial', 
    index = indices_knndm$indx_train) 
  
  mod <- glmnet::glmnet(
    x = preds, 
    sf::st_drop_geometry(x)$occurrence, 
    family = 'binomial', 
    keep = TRUE,
    lambda = cv_model$bestTune$lambda, alpha = cv_model$bestTune$alpha
  )
  
  # to predict onto the confusion matrix, we now need to add the PCNM/MEM values
  # for the relevant layers onto our independent test data, this will require us
  # to create PCNM raster surfaces
  
  xypcnm.sf <- cbind(xypcnm.df, dplyr::select(x, geometry)) |> 
    sf::st_as_sf()
  
  if(is.data.frame(xypcnm.df)){
    pcnm2raster <- function(x){
      
      fit <- fields::Tps(sf::st_coordinates(xypcnm.sf), x)
      p <- terra::rast(predictors[[1]])
      pcnm <- terra::interpolate(p, fit)
      pcnm <- terra::mask(pcnm, predictors[[1]])
      
      return(pcnm)
    }
    
    pcnm <- lapply(xypcnm.df, pcnm2raster)
    pcnm <- terra::rast(pcnm)
    names(pcnm) <- caret::predictors(pcnmProfile)
    
  } else if(is.numeric(xypcnm.df)){
    
    fit <- fields::Tps(sf::st_coordinates(xypcnm.sf), xypcnm.df)
    p <- terra::rast(predictors[[1]])
    pcnm <- terra::interpolate(p, fit)
    pcnm <- terra::mask(pcnm, predictors[[1]])
    names(pcnm) <- caret::predictors(pcnmProfile)
    
    rm(fit, p)
  }  
  
  rm(pcnm2raster, pcnmProfile)
  
  return(list(
    mod = mod,
    pred_mat = preds,
    cv_model = cv_model, 
    pcnm = pcnm))
}
