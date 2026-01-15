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
#' @param planar_proj Numeric, or character vector. An EPSG code, or a proj4 string, for a planar coordinate projection, in meters, for use with the function. 
#' For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). 
#' Otherwise a continental scale projection like 5070 See https://projectionwizard.org/ for more information on CRS.
#'  The value is simply passed to sf::st_transform if you need to experiment. 
#' @param ctrl the control object created by character in the SDM function. 
#' @param indices_knndm from sdm function
#' @param sub the subset predictors from elasticSDM
#' @param predictors Raster stack of environmental predictors (used as template)
#' @keywords internal
#' @noRd
createPCNM_fitModel <- function(x, planar_proj, ctrl, indices_knndm, sub, predictors) {
  
  # Calculate distance matrix
  dis <- calculate_distance_matrix(x, planar_proj)
  
  # Create PCNM eigenvectors
  pcnm_df <- create_pcnm_vectors(dis, n_vectors = 20)
  
  # Select important PCNM features
  occurrence_data <- sf::st_drop_geometry(x)$occurrence
  selected_pcnm <- select_pcnm_features(pcnm_df, occurrence_data, ctrl, indices_knndm)
  
  # Combine environmental and spatial predictors
  preds <- combine_predictors(sub, pcnm_df, selected_pcnm)
  
  # Fit combined model
  model_results <- fit_combined_model(preds, occurrence_data, indices_knndm)
  
  # Create PCNM raster surfaces
  pcnm_rast <- create_pcnm_rasters(
    pcnm_df, 
    selected_pcnm, 
    x, 
    template_raster = terra::subset(predictors, 1)
  )
  
  list(
    mod = model_results$glmnet_model,
    pred_mat = preds,
    cv_model = model_results$cv_model,
    pcnm = pcnm_rast
  )
}


#' Calculate distance matrix for spatial points
#'
#' @param x SF object with occurrence data
#' @param planar_proj Planar projection for accurate distance calculation
#' @return Numeric distance matrix
#' @keywords internal
#' @noRd
calculate_distance_matrix <- function(x, planar_proj) {
  train_planar <- sf::st_transform(x, planar_proj)
  dis <- sf::st_distance(train_planar)
  apply(dis, 2, as.numeric)
}

#' Create PCNM (Principal Coordinates of Neighbor Matrices) eigenvectors
#'
#' @param distance_matrix Numeric distance matrix from calculate_distance_matrix
#' @param n_vectors Number of PCNM vectors to extract (default 10)
#' @return Data frame of PCNM eigenvectors
#' @keywords internal
#' @noRd
create_pcnm_vectors <- function(distance_matrix, n_vectors = 10) {

  if (is.null(distance_matrix) || length(distance_matrix) == 0) {
    stop("Cannot compute PCNM: distance matrix is empty")
  }

  xypcnm <- vegan::pcnm(distance_matrix)
  total_vectors <- ncol(xypcnm$vectors)

  if(total_vectors < n_vectors){
    message('Requested eigenvectors is greater than available, returning all eigenvectors.')
    data.frame(xypcnm$vectors)[,1:total_vectors, drop = FALSE]
  } else {
    data.frame(xypcnm$vectors)[,1:n_vectors, drop = FALSE]
  }
}



#' Select important PCNM vectors using forward feature elimination
#'
#' @param pcnm_df Data frame of PCNM eigenvectors
#' @param occurrence_data Factor vector of occurrence data (0/1)
#' @param ctrl Control object from caret::rfeControl
#' @param cv_indices Cross-validation indices from knndm
#' @return Character vector of selected PCNM variable names
#' @keywords internal
#' @noRd
select_pcnm_features <- function(pcnm_df, occurrence_data, ctrl, cv_indices) {

  ctrl <- caret::trainControl(method="cv", index = spatial_cv$index)
  ffsmodel <- suppressMessages(
    ffs(
      predictors = pcnm_df,
      response = occurrence_data,
      method = "rf",
      trControl = ctrl,
      ntree = 20,
      seed = 1,
      cores = 16
    )
  )
  
  ffsmodel$selectedvars
}

#' Combine environmental and spatial predictors
#'
#' @param env_predictors Data frame of environmental predictors
#' @param pcnm_df Data frame of PCNM vectors
#' @param selected_pcnm Character vector of selected PCNM names
#' @return Combined data frame with proper column names
#' @keywords internal
#' @noRd
combine_predictors <- function(env_predictors, pcnm_df, selected_pcnm) {

  if (is.numeric(pcnm_df)) { # check in case it has collapsed to numeric vector.
    pcnm_df <- data.frame(PCNM1 = pcnm_df)
  }

  pcnm_subset <- pcnm_df[, selected_pcnm, drop = FALSE]
  preds <- cbind(env_predictors, pcnm_subset)
  
  # Handle case where only one PCNM variable selected
  if (is.numeric(pcnm_subset) && !is.data.frame(pcnm_subset)) {
    colnames(preds)[ncol(preds)] <- selected_pcnm
  }
  
  preds
}

#' Fit elastic net model with combined environmental and spatial predictors
#'
#' @param combined_predictors Data frame with env + PCNM predictors
#' @param occurrence_data Factor vector of occurrence data
#' @param cv_indices Cross-validation indices
#' @return List with cv_model and glmnet_model
#' @keywords internal
#' @noRd
fit_combined_model <- function(combined_predictors, occurrence_data, cv_indices) {

  if (ncol(combined_predictors) < 2) {
    warning("glmnet requires ≥2 predictors. Returning NULL.")
    return(list(cv_model = NULL, glmnet_model = NULL))
  }
  
  cv_model <- suppressMessages(
    caret::train(
      x = combined_predictors,
      y = occurrence_data,
      method = "glmnet",
      metric = 'Accuracy',
      family = 'binomial',
      index = cv_indices$indx_train
    )
  )
  
  mod <- suppressMessages(
    glmnet::glmnet(
      x = combined_predictors,
      y = occurrence_data,
      family = 'binomial',
      keep = TRUE,
      lambda = cv_model$bestTune$lambda,
      alpha = cv_model$bestTune$alpha
    )
  )
  
  list(cv_model = cv_model, glmnet_model = mod)
}

#' Interpolate PCNM vector to raster surface using thin plate spline
#'
#' @param pcnm_vector Numeric vector of PCNM values
#' @param coordinates Matrix of spatial coordinates
#' @param template_raster Template raster for interpolation extent/resolution
#' @return SpatRaster of interpolated PCNM values
#' @keywords internal
#' @noRd
interpolate_pcnm_to_raster <- function(pcnm_vector, coordinates, template_raster) {
  # Ensure enough variation
  if (length(unique(pcnm_vector)) < 2) {
    warning("PCNM vector is constant: returning raster of constant values")
    p <- terra::rast(template_raster)
    terra::values(p) <- unique(pcnm_vector)
    return(p)
  }

  # Ensure coordinates are not degenerate
  if (any(apply(coordinates, 2, function(z) length(unique(z)) < 2))) {
    warning("Coordinates are degenerate: returning raster of NA")
    p <- terra::rast(template_raster)
    terra::values(p) <- NA
    return(p)
  }

  fit <- suppressMessages(fields::Tps(coordinates, pcnm_vector))
  p <- terra::rast(template_raster)
  pcnm <- terra::interpolate(p, fit)
  terra::mask(pcnm, template_raster)
}


#' Create raster surfaces for selected PCNM vectors
#'
#' @param pcnm_df Data frame of selected PCNM vectors
#' @param selected_pcnm Character vector of selected PCNM names
#' @param spatial_data SF object with spatial locations
#' @param template_raster Template raster for extent/resolution
#' @return SpatRaster stack of PCNM surfaces
#' @keywords internal
#' @noRd
create_pcnm_rasters <- function(pcnm_df, selected_pcnm, spatial_data, template_raster) {
  # Combine PCNM data with geometry

  if (is.numeric(pcnm_df)) {
    pcnm_df <- data.frame(PCNM1 = pcnm_df)
    selected_pcnm <- names(pcnm_df)  # override selected_pcnm if numeric
  }

  pcnm_subset <- pcnm_df[, selected_pcnm, drop = FALSE]
  xypcnm.sf <- cbind(pcnm_subset, dplyr::select(spatial_data, geometry)) |>
    sf::st_as_sf()
  
  coords <- sf::st_coordinates(xypcnm.sf)
  
  # Handle single vs multiple PCNM vectors
  if (is.data.frame(pcnm_subset)) {
    # Multiple PCNM vectors
    pcnm_list <- lapply(pcnm_subset, function(x) {
      suppressWarnings(interpolate_pcnm_to_raster(x, coords, template_raster))
    })
    pcnm_rast <- terra::rast(pcnm_list)
    names(pcnm_rast) <- selected_pcnm
  } else {
    # Single PCNM vector (numeric)
    pcnm_rast <- suppressWarnings(interpolate_pcnm_to_raster(pcnm_subset, coords, template_raster))
    names(pcnm_rast) <- selected_pcnm
  }
  
  pcnm_rast
}

