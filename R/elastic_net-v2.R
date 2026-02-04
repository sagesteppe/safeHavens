#' Create a quick SDM using elastic net regression
#'
#' @description This function quickly creates a SDM using elastic net regression, and it will
#' properly format all data for downstream use with the `safeHavens` workflow.
#' Note that elastic net models are used for a couple very important reasons:
#' they rescale all input independent variables before modelling, allowing us
#' to combine the raw data with the beta coefficients to use for clustering
#' algorithms downstream. They also allow for 'shrinking' of terms from models
#' by shrinking terms from models we are able to get levels of ecological inference
#' prohibited by older model selection frameworks.
#' @param x A (simple feature) sf data set of occurrence data for the species.
#' @param predictors A terra 'rasterstack' of variables to serve as independent predictors.
#' @param planar_projection Numeric, or character vector. An EPSG code, or a proj4 string, for a planar coordinate projection, in meters, for use with the function. For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). 
#' Otherwise a continental scale projection like 5070 See https://projectionwizard.org/ for more information on CRS. 
#' The value is simply passed to sf::st_transform if you need to experiment.
#' @param domain Numeric, how many times larger to make the entire domain of analysis than a simple bounding box around the occurrence data in `x`.
#' @param quantile_v Numeric, this variable is used in thinning the input data, e.g. quantile = 0.05 will remove records within the lowest 5% of distance to each other iteratively, until all remaining records are further apart than this distance from each other. 
#' If you want essentially no thinning to happen just supply 0.01. Defaults to 0.025.
#' @param PCNM Boolean, defaults to TRUE. Whether to use PCNM surfaces for fitting model.
#' For current scenarios (e.g. EnvironmentalBasedsample) use TRUE, for predictive provenance use FALSE. 
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
#' @returns A list of 12 objects, each of these subsequently used in the downstream SDM post processing sequence, or which we think are best written to disk.
#' The actual model prediction on a raster surface are present in the first list 'RasterPredictions', the independent variables used in the final model are present in 'Predictors', and just the global PCNM/MEM raster surfaces are in 'PCNM'.
#' The fit model is in 'Model', while the cross validation folds are stored in 'CVStructure', results from a single test/train partition in 'ConfusionMatrix', and the two data split in 'TrainData' and 'TestData' finally the 'PredictMatrix' which was used for classifying the test data for the confusion matrix.
#' @export
elasticSDM <- function(
  x,
  predictors,
  planar_projection,
  domain = NULL,
  quantile_v = 0.025,
  PCNM = TRUE
  ){
  
  # Calculate study extent
  study_extent <- calculate_study_extent(
    x,
    planar_projection,
    domain,
    predictors
  )

  # Generate background points
  pa <- generate_background_points(predictors, x)

  # Combine presence and pseudo-absence
  x$occurrence <- 1
  x <- dplyr::bind_rows(x, pa) |>
    dplyr::mutate(occurrence = factor(occurrence))

  # Spatially thin points
  x <- thin_occurrence_points(x, quantile_v)

  # Extract predictor values
  x <- extract_predictors_to_points(x, predictors)

  # Create train/test split
  index <- unlist(caret::createDataPartition(x$occurrence, p = 0.8))
  train <- x[index, ]
  test <- x[-index, ]

  # Create spatial CV folds
  indices_knndm <- create_spatial_cv_folds(train, predictors, k = 5)

  # Perform recursive feature elimination
  lmProfile <- perform_feature_selection(train, indices_knndm)

  # Fit elastic net model
  model_results <- fit_elastic_net_model(
    train,
    predictors(lmProfile),
    indices_knndm
  )

  if(PCNM){
    # Create PCNM variables and refit model
    obs <- createPCNM_fitModel(
      x = train,
      ctrl = caret::rfeControl(
        method = "LGOCV",
        repeats = 5,
        number = 10,
        functions = caret::lrFuncs,
        index = indices_knndm$indx_train,
        verbose = FALSE
      ),
      indices_knndm = indices_knndm,
      planar_proj = planar_projection,
      sub = model_results$selected_data,
      predictors = predictors
    )

    mod <- obs$mod
    pcnm <- obs$pcnm
    cv_structure <- obs$cv_model
    pred_matrix <- obs$pred_mat

    
    # Combine predictors with PCNM
    predictors <- c(predictors, pcnm)

  } else { 
    # this is the branch for predictive provenance because we can not
    # else use the already fit model. nothing else to add.  
    # you cannot have PCNM surfaces of theoretical scenarios.  
    pcnm = NULL
    mod <- model_results$glmnet_model
    cv_structure <- model_results$cv_model
    pred_matrix <- model_results$selected_data

  }

  # Get variables from final model
  vars <- rownames(stats::coef(mod))
  vars <- vars[2:length(vars)]

  # Evaluate model on test data
  cm <- evaluate_model_performance(mod, test, predictors, vars)

  # Create spatial predictions
  rast_cont <- create_spatial_predictions(mod, predictors, vars)

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
}

#' Calculate study domain extent from occurrence points
#'
#' @param x SF object with occurrence points
#' @param planar_projection Projection string or EPSG code
#' @param domain Multiplier for domain size (currently unused but reserved for future)
#' @param predictors Raster stack of predictors
#' @return Terra extent object in predictors' CRS
#' @keywords internal
#' @noRd
calculate_study_extent <- function(
  x,
  planar_projection,
  domain = NULL,
  predictors
) {
  pts_plan <- sf::st_transform(x, planar_projection)

  bb <- sf::st_bbox(pts_plan)
  buff_dist <- as.numeric(
    ((bb[3] - bb[1]) + (bb[4] - bb[2])) / 2
  ) /
    2

  sf::st_union(pts_plan) |>
    sf::st_buffer(buff_dist) |>
    terra::vect() |>
    terra::project(terra::crs(predictors)) |>
    terra::ext()
}

#' Generate background (pseudo-absence) points using environmental distance
#'
#' @param predictors Raster stack of environmental predictors
#' @param occurrences SF object with occurrence points
#' @return SF object with pseudo-absence points (occurrence = 0)
#' @keywords internal
#' @noRd
generate_background_points <- function(predictors, occurrences) {
  sdm::background(
    x = predictors,
    n = nrow(occurrences),
    sp = occurrences,
    method = 'eDist'
  ) |>
    dplyr::select(lon = x, lat = y) |>
    sf::st_as_sf(coords = c('lon', 'lat'), crs = terra::crs(predictors)) |>
    dplyr::mutate(occurrence = 0)
}

#' Spatially thin occurrence points to reduce spatial autocorrelation
#'
#' @param x SF object with presence and pseudo-absence points
#' @param quantile_v Quantile threshold for thinning distance (default 0.025)
#' @return Thinned SF object
#' @keywords internal
#' @noRd
thin_occurrence_points <- function(x, quantile_v = 0.025) {
  sp.coords <- data.frame(
    Species = 'Species',
    data.frame(sf::st_coordinates(x))
  )

  dists <- sf::st_distance(
    x[sf::st_nearest_feature(x), ],
    x,
    by_element = TRUE
  )

  thinD <- as.numeric(stats::quantile(dists, quantile_v, na.rm = TRUE) / 1000)

  thinned <- suppressMessages(
    spThin::thin(
      loc.data = sp.coords,
      thin.par = thinD,
      spec.col = 'Species',
      lat.col = 'Y',
      long.col = 'X',
      reps = 100,
      locs.thinned.list.return = TRUE,
      verbose = FALSE,
      write.files = FALSE,
      write.log.file = FALSE
    )
  )

  # Select the thinning result with the most points retained
  thinned <- data.frame(thinned[which.max(unlist(lapply(thinned, nrow)))]) |>
    sf::st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326)

  x[lengths(sf::st_intersects(x, thinned)) > 0, ]
}

#' Extract predictor values to point locations
#'
#' @param occurrences SF object with occurrence points
#' @param predictors Raster stack of predictors
#' @return SF object with predictor values extracted
#' @keywords internal
#' @noRd
extract_predictors_to_points <- function(occurrences, predictors) {
  terra::extract(predictors, occurrences, bind = TRUE) |>
    sf::st_as_sf()
}

#' Create spatial cross-validation folds using k-fold nearest neighbor distance matching
#'
#' @param train_data SF object with training data
#' @param predictors Raster stack of predictors
#' @param k Number of folds (default 5)
#' @return CAST knndm object with CV fold indices
#' @keywords internal
#' @noRd
create_spatial_cv_folds <- function(train_data, predictors, k = 5) {
  CAST::knndm(train_data, predictors, k = k)
}

#' Perform recursive feature elimination to select important predictors
#'
#' @param train_data SF object with training data
#' @param cv_indices Cross-validation fold indices from knndm
#' @return rfe object with selected features
#' @keywords internal
#' @noRd
perform_feature_selection <- function(train_data, cv_indices) {
  train_dat <- sf::st_drop_geometry(
    train_data[, -which(names(train_data) %in% "occurrence")]
  )

  ctrl <- caret::rfeControl(
    method = "LGOCV",
    repeats = 5,
    number = 10,
    functions = caret::lrFuncs,
    index = cv_indices$indx_train,
    verbose = FALSE
  )

  suppressMessages(
    caret::rfe(
      method = 'glmnet',
      sizes = 3:ncol(train_dat),
      x = train_dat,
      y = sf::st_drop_geometry(train_data)$occurrence,
      rfeControl = ctrl
    )
  )
}

#' Fit elastic net model with selected features
#'
#' @param train_data SF object with training data
#' @param selected_predictors Character vector of selected predictor names
#' @param cv_indices Cross-validation fold indices
#' @return List containing caret model and glmnet model
#' @keywords internal
#' @noRd
fit_elastic_net_model <- function(train_data, selected_predictors, cv_indices) {
  # Prepare data with YES/NO for caret's glmnet
  train1 <- dplyr::mutate(
    train_data,
    occurrence = dplyr::if_else(occurrence == 1, 'YES', 'NO')
  )

  # Fit with caret to find optimal hyperparameters
  cv_model <- suppressMessages(
    caret::train(
      x = sf::st_drop_geometry(train1[, selected_predictors]),
      y = sf::st_drop_geometry(train_data)$occurrence,
      method = "glmnet",
      family = 'binomial',
      index = cv_indices$indx_train
    )
  )

  # Refit with glmnet directly for better predict compatibility
  sub <- sf::st_drop_geometry(train_data[, selected_predictors])

  mod <- glmnet::glmnet(
    x = as.matrix(sub),
    y = sf::st_drop_geometry(train_data)$occurrence,
    family = 'binomial',
    keep = TRUE,
    lambda = cv_model$bestTune$lambda,
    alpha = cv_model$bestTune$alpha
  )

  list(
    cv_model = cv_model,
    glmnet_model = mod,
    selected_data = sub
  )
}

#' Evaluate model performance on test data
#'
#' @param model Fitted glmnet model
#' @param test_data SF object with test data
#' @param predictors Raster stack with selected predictors
#' @param selected_vars Character vector of variable names in model
#' @return Confusion matrix from caret
#' @keywords internal
#' @noRd
evaluate_model_performance <- function(
  model,
  test_data,
  predictors,
  selected_vars
) {
  predict_mat <- predictors[[selected_vars]]
  predict_mat <- as.matrix(
    terra::extract(predict_mat, test_data, ID = FALSE)
  )

  predictions <- stats::predict(model, newx = predict_mat, type = 'class')

  caret::confusionMatrix(
    data = as.factor(predictions),
    reference = test_data$occurrence,
    positive = "1"
  )
}

#' Create spatial predictions from fitted model
#'
#' @param model Fitted glmnet model
#' @param predictors Raster stack with selected predictors
#' @param selected_vars Character vector of variable names in model
#' @return SpatRaster with continuous predictions (0-1 probability)
#' @keywords internal
#' @noRd
create_spatial_predictions <- function(model, predictors, selected_vars) {
  preds <- predictors[[selected_vars]]

  predfun <- function(model, data, ...) {
    stats::predict(model, newx = as.matrix(data), type = 'response')
  }

  terra::predict(preds, model = model, fun = predfun, na.rm = TRUE)
}